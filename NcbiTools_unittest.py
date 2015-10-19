import unittest

from NcbiTools import *



#def test(mapObj, data, cmpFunc=cmp, testFunc="__getitem__"):
#    for key, value in data:
#        sys.stdout.write("%s\t" % key)
#        try:
#            returned = getattr(mapObj, testFunc)(key)
#        except Exception as e:
#            import traceback
#            sys.stderr.write(traceback.format_exc())
#            sys.stdout.write("\tFailed. Exception raised: \"%s\". "
#                             "See log file for more information\n" % str(e))
#            continue
#        sys.stdout.write(str(returned))
#        if cmpFunc(returned, value) == 0:
#            sys.stdout.write("\t\033[92mOK\033[0m\n")
#        else:
#            sys.stdout.write("\t\033[91mFailed. Unexpected value: %s\033[0m\n" % returned)

email = "fheeger@mi.fu-berlin.de"

class noCachUsageTest(unittest.TestCase):
    def setUp(self):
        self.baseMap = NcbiMap(email, useCache=False)
    
    def runTest(self):
        with self.assertRaises(CacheNotUsedError):
            self.baseMap.save()

class mapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value)

class SetMapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertSetEqual(self.map[key], value)

class SciNam2IdTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("Homo sapiens", "9606"),
                        ("Mus musculus", "10090"),
                        ("Human immunodeficiency virus 1", "11676"),
                        ("Escherichia coli", "562"),
                        ("Clavariopsis aquatica", "253306"),
                       ]
        self.map = SpeciesName2TaxId(email, retry=3)

class TaxIdToSpeciesNameTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("9606", "Homo sapiens"),
                        ("10090", "Mus musculus"),
                        ("11676", "Human immunodeficiency virus 1"),
                        ("562", "Escherichia coli"),
                        ("253306", "Clavariopsis aquatica"),
                       ]
        self.map = TaxIdToSpeciesNameMap(email)


class Nucl2taxNameTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [(297822720, "Arabidopsis lyrata subsp. lyrata"),
                        (237825467, "Clavariopsis aquatica"),
                        (21694053, "Homo sapiens"),
                        (34559759, "Rupicapra rupicapra"),
                        (666359714, "Ursus maritimus"),
                       ]
        self.map = NuclId2SpeciesNameMap(email)

class TaxNodeNametoIdTest(SetMapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("Clavariopsis aquatica", set(["253306"])),
                        ("leotiomyceta", set(["716546"])),
                        ("Escherichia coli", set(["562"])),
                        ("Bacteria", set(["2", "629395"])),
                        ("cellular organisms", set(["131567"]))
                       ]
        self.map = TaxonomyNodeName2TaxId(email)

class Nucl2taxIdTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [(297822720, "81972"),
                        (237825467, "253306"),
                        (21694053, "9606"),
                        (34559759, "34869"),
                        (666359714, "29073"),
                       ]
        self.map = NuclId2TaxIdMap(email)

class Prot2nameTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [(684157356, "hypothetical protein HMPREF1120_00314 [Exophiala dermatitidis NIH/UT8656]"),
                        (918400562, "hemoglobin [Homo sapiens]"),
                        (190663721, "laccase, partial [Clavariopsis aquatica]"),
                        (284794136, "Chain F, Crystal Structure Of Zaire Ebola Vp35 Interferon Inhibitory Domain K339a Mutant")
                       ]
        self.map = ProtId2ProtNameMap(email)

class TaxParentTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [(666681, "1055487"),
                        (4890, "451864"),
                        (63221, "9606"),
                        (9606, "9605"),
                        (153057, "327045"),
                       ]
        self.map = TaxonomyParentMap(email)

class LcnTest(unittest.TestCase):
    def setUp(self):
        self.testSet = [([1491027, 341667, 577457, 42310], "42310"),
                        ([109873, 4751, 451459, 4805], "109873"),
                        ([2, 414713,85026, 1760], "1760"),
                        ]
        self.tree = NcbiTaxonomyTree(email,
                                     "/home/heeger/data/ncbi_tax/tax2parent.db")

    def runTest(self):
        for inData, outData in self.testSet:
            self.assertEqual(self.tree.lcn(inData), outData)

class LineageMapTest(unittest.TestCase):
    def setUp(self):
        self.map = LineageMap(email)

    def test_succeed(self):
        #polar bear <-> brown bear
        lineages = zip(self.map[9644], self.map[29073])
        for tupA, tupB in lineages[:-1]:
            self.assertTupleEqual(tupA, tupB)
        self.assertNotEqual(*lineages[-1])
        #E. coli <-> E. alberti
        lineages = zip(self.map[562], self.map[208962])
        for tupA, tupB in lineages[:-1]:
            self.assertTupleEqual(tupA, tupB)
        self.assertNotEqual(*lineages[-1])
    
    @unittest.expectedFailure
    def test_fail(self):
        #mouse <-> rat
        for tupA, tupB in zip(self.map[10090], self.map[10116]):
            self.assertTupleEqual(tupA, tupB)



if __name__ == '__main__':
    unittest.main()
