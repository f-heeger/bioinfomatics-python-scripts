import unittest

from OttTools import *

#class mapTest(object):
#    def runTest(self):
#        for key, value in self.testSet:
#            self.assertEqual(self.map[key], value)

class MapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value)

class Nam2IdTest(MapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("Candida", set([4085890])),
                        ("Clavariopsis", set([686240])),
                        ("Bacteria", set([844192, 5021155])),
                        ("Ursus maritimus", set([10732]))
                       ]
        self.map = OttName2IdMap()

class Nam2IdContextTest(unittest.TestCase):
    
    def test_valid(self):
        self.map = OttName2IdMap(context="Insects")
        self.assertSetEqual(self.map["Bacteria"],  set([5021155]))

    def test_invalid(self):
        with self.assertRaises(ValueError):
            self.map = OttName2IdMap(context="INVALID")

#TODO test CachedOttName2IdMap

class Id2otherTaxonomyTest(MapTest, unittest.TestCase):
    def setUp(self):
        self.map = OttId2otherTaxonomyMap()
        self.testSet = [(686240, {u"ncbi": u"253305", u"if": u"7692", 
                                 u"gbif": u"2614823", u"irmng": u"1373703"}),
                        (10732, {u"ncbi": u"29073", u"gbif": u"2433451", 
                                 u"irmng": u"11061491"})
                       ]

class OttId2NcbiTaxIdTest(MapTest, unittest.TestCase):
    def setUp(self):
        self.map = OttId2NcbiTaxIdMap()
        self.testSet = [(686240, u"253305"),
                        (10732, u"29073"),
                       ]

#TODO test CachedOttId2NcbiTaxIdMap

if __name__ == '__main__':
    unittest.main()
