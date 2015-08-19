import unittest

from OttTools import *

#class mapTest(object):
#    def runTest(self):
#        for key, value in self.testSet:
#            self.assertEqual(self.map[key], value)

class MapTest(object):
    def runTest(self):
        print(self.__class__)
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value, "Problem with key: "+key)

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
        print(self.__class__)
        self.map = OttName2IdMap(context="Insects")
        self.assertSetEqual(self.map["Bacteria"],  set([5021155]))

    def test_invalid(self):
        print(self.__class__)
        with self.assertRaises(ValueError):
            self.map = OttName2IdMap(context="INVALID")

#TODO test CachedOttName2IdMap

class Id2otherTaxonomyTest(MapTest, unittest.TestCase):
    def setUp(self):
        self.map = OttId2otherTaxonomyMap()
        self.testSet = [(686240, {"ncbi": 253305, "if": 7692, 
                                 "gbif": 2614823, "irmng": 1373703}),
                        (10732, {"ncbi": 29073, "gbif": 2433451, 
                                 "irmng": 11061491})
                       ]

if __name__ == '__main__':
    unittest.main()
