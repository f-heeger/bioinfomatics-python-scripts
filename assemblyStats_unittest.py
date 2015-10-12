import unittest
from StringIO import StringIO

from assemblyStats import *


class Test1(unittest.TestCase):
    def setUp(self):
        self.fasta = StringIO(""">c1 length=40
TACGGGGGCGAGAGACTTGCCGACAAGACAATAGCCGCGT
>c2 length=40
CGAACCAGGGGAGCAAGAGGGCTACCGTTT
>c3 length=10
TCGGCCTCAG
>c4 length=10
GTTCAAAGCC
>c5 length=10
ATAGCCGCGT
""")

    def runTest(self):
        n50, l50 = assemblyStats(readFasta(self.fasta))
        self.assertEqual(n50, 2)
        self.assertEqual(l50, 30)
        
        
class Test2(unittest.TestCase):
    def setUp(self):
        self.fasta = StringIO(""">c1 length=40
TACGGGGGCGAGAGACTTGCCGACAAGACAATAGCCGCGT
>c2 length=40
CGAACCAGGGGAGCAAGAGGGCTACCGTTTATAGCCGCGT
>c3 length=10
TCGGCCTCAG
>c4 length=10
GTTCAAAGCC
>c5 length=10
ATAGCCGCGT
""")

    def runTest(self):
        n50, l50 = assemblyStats(readFasta(self.fasta))
        self.assertEqual(n50, 2)
        self.assertEqual(l50, 40)
        
class Test3(unittest.TestCase):
    def setUp(self):
        self.fasta = StringIO(""">c1 length=40
TACGGGGGCGAGAGACTTGCCGACAAGACAATAGCCGCGT
>c2 length=10
TCGGCCTCAG
>c3 length=10
GTTCAAAGCC
>c4 length=10
ATAGCCGCGT
>c5 length=10
CGAACCAGGG
>c6 length=10
TACGGGGGCG
""")

    def runTest(self):
        n50, l50 = assemblyStats(readFasta(self.fasta))
        self.assertEqual(n50, 2)
        self.assertEqual(l50, 10)
        
        
if __name__ == '__main__':
    unittest.main()
