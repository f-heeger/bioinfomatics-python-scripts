import unittest

from KeggTools import *

class mapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value)
            
class NcbiGiToKeggMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("16130957", "eco:b3061"),
                        ("256810887", "mfe:Mefer_0938"),
                        ("83590826", "mta:Moth_1996"),
                        ("317410866", None)
                       ]
        self.map = NcbiGiToKeggMap(retry=3)
        
class KeggGeneTopathwayMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("eco:b3061", set(["eco00630"])),
                        ("mfe:Mefer_0938", set(["mfe00680", "mfe01100",
                                                "mfe01120", "mfe01200"])),
                        ("mta:Moth_1996", set(["mta00220", "mta01100"])),
                       ]
        self.map = KeggGeneToPathwayMap(retry=3)

class KeggPathwayIdToNameMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("eco00630", "Glyoxylate and dicarboxylate metabolism"),
                        ("mta00220", "Arginine biosynthesis"),
                        ("mfe00680", "Methane metabolism")
                        ]
        self.map = KeggPathwayIdToNameMap(retry=3)

class KeggReactionIdToEcMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("R00577", set(["2.3.1.68"])),
                        ("R04558", set(["2.4.2.-","4.1.3.-"])),
                        ("R04066", set(["1.14.13.5"])),
                       ]
        self.map = KeggReactionIdToEcMap()

if __name__ == '__main__':
    unittest.main()
