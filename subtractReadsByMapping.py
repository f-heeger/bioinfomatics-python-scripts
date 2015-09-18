import sys
import gzip
from optparse import OptionParser

from Bio.SeqIO import parse

def readSamFile(samStream, exchange, threshold=None):
    """Read a sam file and put the IDs of all mapping reads into a dict"""
    if not threshold is None:
        raise NotImplemented("Threshold parameter is not used "
                             "in this implementation.")
    mapped = {}
    for line in samStream:
        if line.startswith("@"):
            #skip header lines
            continue
        qname, flag, sname, _ = line.split("\t", 3)
        #check bit flag for mapping (see sam specification)
        if int(flag) & 4:
            #unmapped
            continue
        if exchange:
            mapped[sname] = None
        else:
            mapped[qname] = None
    return mapped
    
def readTabularBlastOutput(blastStream, exchange, th=0.000001):
    """Read a blast tabular file and put the sequence IDs of any sequences
    that have a good e-value hit into a dict. "Good" is defined by the e-value
    threshold."""
    mapped = {}
    for line in blastStream:
        #standart blast tabular file is tab sepearted and has the evlaue at 
        # poition 10
        lineArr = line.split("\t")
        if float(lineArr[10])<=th:
            if exchange:
                mapped[lineArr[1]] = None
            else:
                mapped[lineArr[0]] = None
    return mapped
    
def subtractReads(readStream, format, mapped, outStream, outFormat="fastq", 
                  mappedStream=None):
    """Iterate over a stream of reads from a fasta or fastq file and only 
    output the ones that are NOT in the 'mapped' dictionary"""
    for rec in parse(readStream, format):
        if rec.id in mapped or (rec.id[-2]=="/" and rec.id[-2] in mapped):
            if not mappedStream is None:
                mappedStream.write(rec.format(outFormat))
        else:
            outStream.write(rec.format(outFormat))

def subtractReadsBySamFile(samPath, readPath1, outPath1, readPath2=None, 
                           outPaht2=None, inForm="fastq", outForm="fastq", 
                           mappedPath1=None, mappedPath2=None, gziped=False, 
                           mgziped=False, log=None, exchange=False):
    """Convenience function to call sbtractReadsByFile with a sam file"""
    subrtractReadsByFile(samPath, readSamFile, readPath1, outPath1, readPath2, 
                         outPaht2, inForm, outForm, mappedPath1, mappedPath2, 
                         gziped, mgziped, log, exchange)

def subtractReadsByBlastFile(blastPath, readPath1, outPath1, readPath2=None,
                             outPaht2=None, inForm="fastq", outForm="fastq", 
                             mappedPath1=None, mappedPath2=None, gziped=False, 
                             mgziped=False, log=None, exchange=False, threshold=None):
    """Convenience function to call sbtractReadsByFile with a tabular blast file
    """
    subrtractReadsByFile(blastPath, readTabularBlastOutput, readPath1, outPath1,
                         readPath2, outPaht2, inForm, outForm, mappedPath1, 
                         mappedPath2, gziped, mgziped, log, exchange, threshold)

def subrtractReadsByFile(filePath, readFunc, readPath1, outPath1, 
                         readPath2=None, outPath2=None, inForm="fastq", 
                         outForm="fastq", mappedPath1=None, mappedPath2=None, 
                         gziped=False, mgziped=False, log=None, exchange=False, threshold=None):
    """Function to subtract reads that were "mapped" from a fasta/fastq file
    
    A "mapping" file is read with a function that has to be supplied. Either one
    input file with single end reads or two input files with paired end reads 
    can be provided. Read and mapping files can be gzip if the according switch 
    is set. A threshold for "mapping" quality can be given if applicable.
    """
    if mgziped:
        tOpen = gzip.open
    else:
        tOpen = open
        if log:
            log.write("Reading mapping file...\n")
            if mgziped:
                log.write("Using gzip for reading mapping file.\n")
    with tOpen(filePath) as inFile:
        if threshold is None:
            mappedReads = readFunc(inFile, exchange)
        else:
            mappedReads = readFunc(inFile, exchange, threshold)
    if log:
        log.write("%i reads mapped\n" % len(mappedReads))
    rPath = [readPath1, readPath2]
    oPath = [outPath1, outPath2]
    mPath = [mappedPath1, mappedPath2]
    
    if gziped:
        tOpen = gzip.open
    else:
        tOpen = open
    if log:
        log.write("Subtracting Reads...\n")
        if gziped:
            log.write("Using gzip compression for in and out files.\n")
    for i in range(2):
        if rPath[i] != None:
            read =  tOpen(rPath[i]) 
            out = tOpen(oPath[i], "w")
            
            if not mPath[i] is None:
                mappedOut = open(mPath[i], "w")
            else:
                mappedOut = None
            try:
                subtractReads(read, inForm, mappedReads, out, outForm, 
                              mappedOut)
            finally:
                if mappedOut:
                    mappedOut.close()

if __name__ == "__main__":

    usage = "usage: %prog [options] mappingfile read1.fastx out1.fastx " \
            "read2.fastx out2.fastx"

    parser = OptionParser(usage)
    
    parser.add_option("-q", "--quiet",
                      action="store_true", dest="quiet", default=False, 
                      help="do not print status messages to the screen",)
    parser.add_option("-a", "--fasta",
                      action="store_true", dest="fasta",
                      default=False, help="input file(s) is/are fasta",)
    parser.add_option("-m", "--mapped1",
                      action="store", type="string", dest="mapped1", 
                      default=None, 
                      help="write mapped reads from read 1 to this files",
                      metavar="X")
    parser.add_option("-n", "--mapped2",
                      action="store", type="string", dest="mapped2", 
                      default=None, 
                      help="write mapped reads from read 2 to this files",
                      metavar="X")
    parser.add_option("-z", "--gzip",
                      action="store_true", dest="gzip",
                      default=False, help="input file(s) is/are gzipped",)
    parser.add_option("-y", "--mapping-gzip",
                      action="store_true", dest="m_gzip",
                      default=False, help="mapping file is gzipped",)
    parser.add_option("-b", "--blast",
                      action="store_true", dest="blast", default=False, 
                      help="mapping file is tabular blast output instead "
                           "of sam file",)
    parser.add_option("-x", "--exchange",
                      action="store_true", dest="exchange", default=False, 
                      help="exchange identity of query ad source in mapping file",)
    parser.add_option("-t", "--threshold",
                      action="store", type="float", dest="threshold", 
                      default=None, 
                      help="consider reads with an e-value lower than X as"
                           " \"mapped\". (Can only be used in blast mode) "
                           "[default: 0.000001]",
                      metavar="X")
                       
    (options, args) = parser.parse_args()
    
    if (options.threshold!=None and not options.blast):
        parser.error("Options -t can only be used in blast mode (-b).")
    
    if options.quiet:
        log = None
    else:
        log = sys.stderr
    
    inForm = "fastq"
    if options.fasta:
        inForm = "fasta"
    
    if len(args) > 3:
        myArgs=args[0:5] + [inForm, inForm, options.mapped1, options.mapped2, 
                            options.gzip, options.m_gzip, log, options.exchange]
    else:
        myArgs=args[0:3] + [None, None, inForm, inForm, options.mapped1, None, 
                            options.gzip, options.m_gzip, log, options.exchange]
    if options.blast:
        if log:
            log.write("Using BLAST file as mapping\n")
        subtractReadsByBlastFile(*myArgs, threshold=options.threshold)
    else:
        if log:
            log.write("Using SAM file as mapping\n")
        subtractReadsBySamFile(*myArgs)
        
