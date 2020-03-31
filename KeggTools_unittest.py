import unittest
from warnings import catch_warnings

from KeggTools import *

class mapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value)
            
# 2016-02-08 kegg service seems to be broken
#class NcbiGiToKeggMapTest(mapTest, unittest.TestCase):
#    def setUp(self):
#        self.testSet = [("16130957", None),
#                        ("256810887", "mfe:Mefer_0938"),
#                        ("83590826", "mta:Moth_1996"),
#                        ("317410866", None)
#                       ]
#        self.map = NcbiGiToKeggMap(retry=3)
        
class KeggGeneTopathwayMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("eco:b3061", set(["eco00630"])),
                        ("mfe:Mefer_0938", set(["mfe00680", "mfe01100",
                                                "mfe01120", "mfe01200"])),
                        ("mta:Moth_1996", set(["mta00220", "mta01100", "mta01230"])),
                       ]
        self.map = KeggGeneToPathwayMap(retry=3, useCache=False)

class KeggPathwayIdToNameMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("eco00630", "Glyoxylate and dicarboxylate metabolism"),
                        ("mta00220", "Arginine biosynthesis"),
                        ("mfe00680", "Methane metabolism")
                        ]
        self.map = KeggPathwayIdToNameMap(retry=3, useCache=False)

class KeggKoIdToDefMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("K00262", "glutamate dehydrogenase (NADP+) [EC:1.4.1.4]"),
                        ("K14805" , "ATP-dependent RNA helicase DDX24/MAK5 [EC:3.6.4.13]"),
                        ("K13953", "alcohol dehydrogenase, propanol-preferring [EC:1.1.1.1]")
                       ]
        self.map = KeggKoIdToDefMap()

class KeggReactionIdToEcMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("R00577", set(["2.3.1.68"])),
                        ("R04558", set(["2.4.2.-","4.1.3.-"])),
                        ("R04066", set(["1.14.13.5"])),
                       ]
        self.map = KeggReactionIdToEcMap()

class KeggProteinToKoMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("spo:SPCC4B3.01", set(["ko:K01011"])),
                        ("mcc:695842", set(["ko:K01011"])),
                        ("mpur:MARPU_10910", set(["ko:K11181"])),
                        ("cal:CAALFM_C701250WA", set(["ko:K14753"])),
                       ]
        self.map = KeggProteinToKoMap(retry=3, useCache=False)

class KeggKoToPathwayMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("ko:K00939", set(["ko00230", "ko01100", "ko01110",
                                           "ko01130"])),
                        ("ko:K07679" , set(["ko02020", "ko05133"])),
                        ("ko:K07678", set(["ko02020"]))
                       ]
        self.map = KeggKoToPathwayMap()

class KeggKoToEnzymeMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("ko:K01961", set(["6.4.1.2", "6.3.4.14"])),
                        ("ko:K07679", set(["2.7.13.3"])),
                        ("ko:K07678", set(["2.7.13.3"]))
                       ]
        self.map = KeggKoToEnzymeMap()

class KeggEcToPathwayMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("ec:1.1.1.388", set(["map00030"])),
                        ("ec:2.7.4.3", set(["map00230", "map00730", "map01100", "map01110", "map01130"])),
                        ("ko:1.10.3.2", set([]))
                       ]
        self.map = KeggEcToPathwayMap()

class KeggEcToKoMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("ec:1.1.1.100", set(["K00059", "K11610"])),
                        ("ec:1.1.1.101", set(["K06123"])),
                        ("ec:1.1.1.201", set([]))
                       ]
        self.map = KeggEcToKoMap()

if __name__ == '__main__':
    unittest.main()
