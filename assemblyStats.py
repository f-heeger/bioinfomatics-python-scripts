from optparse import OptionParser
try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False

from Bio.SeqIO import parse
from Bio.SeqUtils import GC

def assemblyStats(contLenTable):
    """Compute assembly stats from a table of contig IDs and length"""
    sTable = sorted(contLenTable, cmp=lambda x,y: cmp(x[1],y[1]), reverse=True)
#    import pdb; pdb.set_trace()
    lenSum = sum([c[1] for c in contLenTable])
    i=0
    iSum=0
    while i<len(contLenTable) and iSum < lenSum/2:
        iSum += sTable[i][1]
        i += 1
    return (i, sTable[i-1][1])

def readFasta(inStream):
    """Read fasta file and save a table with sequence IDs and sequence length"""
    out = []
    for record in parse(inStream, "fasta"):
        out.append((record.id, len(record), GC(record.seq), record.seq.count("N")))
    return out
    
if __name__ == "__main__":

    usage = "usage: %prog [options] input.fasta[.gz]"
    epi = "Computes basic assembly statistics given that input.fasta is a" \
          " multi fasta file with one contig per record."
    parser = OptionParser(usage, epilog=epi)
    
    if gzImported:
        #if gzip library could be imported offer this option
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
        parser.add_option("-t", "--tabular",
                           action="store_true", dest="tabular",
                           default=False, help="write tabular output",)
    (options, args) = parser.parse_args()
    

    
    if gzImported and options.gzip:
        inStream = gzip.open(args[0], "r")
    else:
        inStream = open(args[0], "r")
    
    tab = readFasta(inStream)
    n50, l50 = assemblyStats(tab)
    lenSum = sum([c[1] for c in tab])
    maxLen = max([c[1] for c in tab])
    minLen = min([c[1] for c in tab])
    avgGC = sum([c[2] for c in tab])/float(len(tab))
    nSum = sum([c[3] for c in tab])
    if options.tabular:
        print("%s\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%i" 
              % (args[0], len(tab), l50, n50, lenSum, maxLen, minLen, gc, nSum))
    else:
        print("Number of contigs:           %15i" % len(tab))
        print("N50 length (normally used):  %15i bp" % l50)
        print("N50 rank (rarely used):      %15i" % n50)
        print("Total length of all contigs: %15i bp" % lenSum)
        print("Longest contig:              %15i bp" % maxLen)
        print("Shortest contig:             %15i bp" % minLen)
        print("GC content:                  %15.2f %%" % avgGC)
        print("Total number of Ns:          %15i" % nSum)
