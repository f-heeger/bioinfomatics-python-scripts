#!/usr/bin/env python

import sys
from optparse import OptionParser

try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False
try:
    from Bio.SeqIO import parse
except ImportError:
    print("Could not import Biopython. Please make sure it is installed.")
    exit(1)
    
def removSeqsWithN(threshold, inStream1, outStream1, inStream2=None, 
                   outStream2=None, fileForm="fastq"):
    """Write all sequences with equal or less than threshold Ns to the 
outStream.
    
Optinal sequences with more Ns can be written to a different stream.
    """
    t=0
    n=0
    strIt1 = parse(inStream1, fileForm)
    if not inStream2 is None:
        strIt2 = parse(inStream2, fileForm)
    while True:
        t+=1
        try:
            r1 = strIt1.next()
        except StopIteration:
            break
        remove = r1.seq.count("N") > threshold
        if not inStream2 is None:
            r2 = strIt2.next()
            remove |= r2.seq.count("N") > threshold
        if remove:
            n+=1
        else:
            outStream1.write(r1.format(fileForm))
            if not inStream2 is None:
                outStream2.write(r2.format(fileForm))
    return t, n
    
if __name__ == "__main__":

    usage = "usage: %prog [options] inputfile1 [inputfile2]"
    epi = "Remove sequences with more then a certain number of Ns"
    parser = OptionParser(usage, epilog=epi)
    
    if gzImported:
        #if gzip library could be imported offer this option
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
    parser.add_option("-a", "--fasta",
                      action="store_true", dest="fasta",
                      default=False, help="set input and output format to"
                      " fasta [default: fastq]",)
    parser.add_option("-n", "--max-n",
                      action="store", type="int", dest="maxNs", default=5, 
                      help="remove all sequences with more than X Ns "
                      "[default: 5]", metavar="X")
#    parser.add_option("-r", "--print-removed",
#                       action="store_true", dest="removed", default=False, 
#                       help="write files with removed sequences",)
    (options, args) = parser.parse_args()

    if gzImported and options.gzip:
        tOpen = gzip.open
        maxsplit = 2
        sys.stderr.write("Writing and reading with gzip\n")
    else:
        maxsplit = 1
        tOpen = open

    #determine file streams and open if neccessary
    if len(args) < 1:
        sys.stderr.write("Please give at least one input file.\n")
        exit(1)
    sys.stderr.write("Reading from file: %s\n" % args[0])
    inStream1 = tOpen(args[0])
    splitName = args[0].rsplit(".", maxsplit)
    outStream1 = tOpen("%s_Ns_removed.%s" % (splitName[0], 
                                             ".".join(splitName[1:])),
                       "w")
    paired = False
    if len(args) > 1:
        paired = True
        sys.stderr.write("Running in paired mode\n")
        sys.stderr.write("Reading from second file: %s\n" % args[1])
        inStream2 = tOpen(args[1])
        splitName2 = args[1].rsplit(".", maxsplit)
        outStream2 = tOpen("%s_Ns_removed.%s" % (splitName2[0], 
                                                 ".".join(splitName2[1:])),
                           "w")
    
#    remStream = None
#    if not options.removed is None:
#        sys.stderr.write("Writing removed sequences to file: %s\n"
#                          % options.removed)
#        remStream = open(options.removed, "w")
    #setting fiel format
    
    fileFormat = "fastq"
    if options.fasta:
        fileFormat = "fasta"
    sys.stderr.write("Reading and writing in %s format\n" % fileFormat)
    #running actual function
    sys.stderr.write("Start filtering...\n")
    if paired:
        total, removed = removSeqsWithN(options.maxNs, inStream1, outStream1,
                                        inStream2, outStream2, fileFormat)
    else:
        total, removed = removSeqsWithN(options.maxNs, inStream1, outStream1,
                                        fileForm=fileFormat)
                                    
    #close all open streams
    inStream1.close()
    outStream1.close()
    if paired:
        inStream2.close()
        outStream2.close()

    sys.stderr.write("Done. Removed %i of %i sequneces\n" % (removed, total))

