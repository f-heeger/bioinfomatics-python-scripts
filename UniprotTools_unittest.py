import unittest

from UniprotTools import *


class GoTest(unittest.TestCase):
    def setUp(self):
        self.goTestData = [("N4UXJ0",[{'function': 'ATP binding', 
                                       'evidenceCode': 'IEA', 
                                       'aspect': 'F', 
                                       'goId': 'GO:0005524', 
                                       'evidenceSource': 'InterPro'}, 
                                      {'function': 'protein kinase activity', 
                                       'evidenceCode': 'IEA', 
                                       'aspect': 'F', 
                                       'goId': 'GO:0004672', 
                                       'evidenceSource': 'InterPro'}
                                     ]
                          ),
                          ]
    def runTest(self):
        u=UniprotInterface(contact="fheeger@mi.fu-berlin.de")
        for key, value in self.goTestData:
            u.getData(key)
            self.assertEqual(u.readGoInfo(), value)

class CazyTest(unittest.TestCase):
    def setUp(self):
        self.cazyTestData = [
            ("A2QYE5", {"id": "GH28", 
                        "name": "Glycoside Hydrolase Family 28"}),
            ("Q39Z37", {"id": "GT2", 
                        "name": "Glycosyltransferase Family 2"}),
        ]
        
    def runTest(self):
        u=UniprotInterface(contact="fheeger@mi.fu-berlin.de")
        for key, value in self.cazyTestData:
            u.getData(key)
            self.assertEqual(u.readCazyInfo(), value)

class KeggTest(unittest.TestCase):
    def setUp(self):
        self.testData = [("Q55F68", ["ddi:DDB_G0267376"]),
                         ("P54202", ["ani:AN3741.2"]),
                         ("Q9UUE1", ["spo:SPBC17G9.11c"]),
                        ]
                        
    def runTest(self):
        u=UniprotInterface(contact="fheeger@mi.fu-berlin.de")
        for key, value in self.testData:
            u.getData(key)
            self.assertEqual(u.readKeggInfo(), value)

class mapTest(object):
    def runTest(self):
        for key, value in self.testSet:
            self.assertEqual(self.map[key], value)
            
class UniprotIdMapTest(mapTest, unittest.TestCase):
    def setUp(self):
        self.testSet = [("684157356", "H6BMR7"),
                        ("915106749", "A0A072PTG6"),
                        ("358058139", None),
                       ]
                       
        self.map=UniprotIdMap("P_GI","ACC", retry=3, 
                              contact="fheeger@mi.fu-berlin.de",
                              returnNone=True)

if __name__ == '__main__':
    unittest.main()
