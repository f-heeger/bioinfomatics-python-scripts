import sys
from optparse import OptionParser

from Bio.SeqIO import parse

def _trimNStreams(trimFunc ,readFiles, outFiles, outFormat, stats, 
                  *additionalPars):
    discarded = 0
    overall = 0
    while 1:
        try:
            #make a list with the next read from each read file
            reads = []
            for readFile in readFiles:
                reads.append(readFile.next())
        except StopIteration:
            break
        overall += 1
        outReads = []
        for r in reads:
            trimmedRead = trimFunc(r, stats, *additionalPars)
            outReads.append(trimmedRead)
        if all(outReads):
            for i, read in enumerate(outReads):
                outFiles[i].write(read.format(outFormat))
        else:
            discarded += 1
    return (discarded, overall)

def trimByQuality(path1, outPath1, path2=None, outPath2=None, minQual=20, 
    minLen=75, bothEnds=True, outFormat="fasta", stats=None):
    """Trim fastq reads according to a quality threshold.
    
    If path two is given path1 and path2 should contain reads from read1 and 
    read2 of a paired end library.
    Reads with length below minLen after trimming will be discarded (together 
    with there mate if path2 is given).
    If bothEnds is true, nucleotides with quality below minQual will also be 
    removed from the begining of the read.
    """
    
    files = [parse(open(path1), "fastq")]
    if outPath1 is None:
        outFiles = [sys.stdout]
    else:
        outFiles = [open(outPath1, "w")]
    if path2:
        files.append(parse(open(path2), "fastq"))
        outFiles.append(open(outPath2, "w")) 
        
    return _trimNStreams(_trimReadByQuality, files, outFiles, outFormat, stats, 
           minQual, minLen, bothEnds)

def _trimReadByQuality(r, stats, minQual, minLen, bothEnds):
    
    if minQual == 0:
        return r
    
    start = len(r)/3 # <-- I assume that if I want to trim at the start 
                     #     I can ignore the last two thirds of the read
    end = len(r)/3 # <---- I assume that if I want to trim at the end I can 
                   #       ignore the first third of the read 
    qual = r._per_letter_annotations["phred_quality"]
    #find first nucleotide so that the rest of the read as an average 
    # quality below the threshold
    while end < len(r) and float(sum(qual[end:]))/(len(r)-end) > minQual:
        end += 1
    #push the end further out as long as the nucleotide in question has a good
    # enough quality
    while end < len(r) and qual[end] > minQual:
        end += 1
    #if trimmingfrom the beginning is active
    if bothEnds:
        #find the first nucleotide so that the average quality up to that
        # nucleotid is below the threshold
        while start > 0 and float(sum(qual[:start]))/start > minQual:
            start -= 1
        #push the start further to the begining if the as long as the quality 
        # of the nuclotide is good enough
        while start > 0 and qual[start] > minQual:
            start -= 1
    else:
        #if no trimming at the start is done the first nucleotide to use
        # is the 0th
        start = 0    
    if end-start < minLen:
        return None
    else:
        if not stats is None:
            stats.append((start, len(r)-end))
        return r[start:end]

def trimToErrorProb(inPath1, outPath1, inPath2=None, outPath2=None, prob=0.95,
                    minLen=75, outFormat="fasta", stats=None):  
    files = [parse(open(inPath1), "fastq")]
    if outPath1 is None:
        outFiles = [sys.stdout]
    else:
        outFiles = [open(outPath1, "w")]
    if inPath2:
        files.append(parse(open(inPath2), "fastq"))
        outFiles.append(open(outPath2, "w")) 
        
    return _trimNStreams(_trimReadToErrorProb, files, outFiles, outFormat, 
                         stats, prob, minLen)

def _trimReadToErrorProb(r, stats, maxErrorProb, minLen):
    i = 0
    qual = r._per_letter_annotations["phred_quality"]
    p = []
    while i<len(r) and sum(p)/float(i+1) < maxErrorProb:
        p.append(10**(float(qual[i])/-10))
        i+=1
    if len(r)-i < minLen:
        return None
    else:
        if not stats is None:
            stats.append((0, len(r)-i))
        return r[:i]

def trimByConstant(inPath1, outPath1, inPath2=None, outPath2=None, cEnd=20, 
                   minLen=75, cBegin=None, outFormat="fasta", inFormat="fastq",
                   stats=None):
    """Trim constant number of nucleotides from the end of each read in a 
    fastq file.
    
    If cBegin is given it will be used as the number of reads to cut from the 
    beginning of a read.
    """
    files = [parse(open(inPath1), "fastq")]
    if outPath1 is None:
        outFiles = [sys.stdout]
    else:
        outFiles = [open(outPath1, "w")]
    if inPath2:
        files.append(parse(open(path2), "fastq"))
        outFiles.append(open(outPath2, "w")) 
    return _trimNStreams(_trimReadByConstant, files, outFiles, outFormat, 
                         stats, cEnd, cBegin, minLen)
    
def _trimReadByConstant(r, stats, cEnd, cBegin, minLen):
    if len(r) - cEnd - cBegin < minLen:
        return None
    else:
        if not stats is None:
            stats.append((cBegin, cEnd))
        return r[cBegin:cEnd]
    
def trimToConstant(inPath1, outPath1, inPath2=None, outPath2=None, length=75, 
                   minLen=75, outFormat="fasta", inFormat="fastq", stats=None):
    """Trim each read in a fastq file to a certain length.
    
    Reads that are shorter than the desired length but at least min. length will
    be papped by Ns. Reads shorter than the min. length will be discareded.
    """
    files = [parse(open(inPath1), "fastq")]
    if outPath1 is None:
        outFiles = [sys.stdout]
    else:
        outFiles = [open(outPath1, "w")]
    if inPath2:
        files.append(parse(open(path2), "fastq"))
        outFiles.append(open(outPath2, "w")) 
    return _trimNStreams(_trimReadToConstant, files, outFiles, outFormat, 
                         stats, length, minLen)

def _trimReadToConstant(r, stats, length, minLen):
    if len(r) < minLen:
        return None
    else:
        if len(r) < length:
            if not stats is None:
                stats.append((0, len(r)-length))
            return r[:length] + "N"*(length - len(r))
        else:
            if not stats is None:
                stats.append((0, len(r)-length))
            return r[:length]

if __name__ == "__main__":

    usage = "usage: %prog [options] input1.fastq [output1.fasta input2.fastq "\
            "output2.fasta]"

    parser = OptionParser(usage)
    
    parser.add_option("-q", "--quite",
                       action="store_true", dest="quiet", default=False, 
                       help="do not print status messages to the screen",)
    parser.add_option("-u", "--fastq",
                       action="store_true", dest="fastq", default=False, 
                       help="set output format to fastq [default: fasta]",)
    parser.add_option("-t", "--min-quality",
                       action="store", type="int", dest="minQuality",
                       default=0, help="quality threshold",
                       metavar="X")
    parser.add_option("-l", "--min-length",
                       action="store", type="int", dest="minLength", default=0, 
                       help="minimal length for reads after trimming"
                            " [default: 0]",
                       metavar="X")
    parser.add_option("-p", "--max-error-prob",
                       action="store", type="float", dest="maxProb", 
                       default=None, 
                       help="maximal over all error probability in one read",
                       metavar="X")
    parser.add_option("-c", "--const",
                       action="store", type="int", dest="const",
                       default=0, 
                       help="remove X bases from the end of the read",
                       metavar="X")
    parser.add_option("-b", "--begin-const",
                       action="store", type="int", dest="beginConst",
                       default=0, 
                       help="remove X bases from the begining of the read",
                       metavar="X")
    parser.add_option("-r", "--crop",
                       action="store", type="int", dest="crop",
                       default=0, 
                       help="cut all reads to length X",
                       metavar="X")
                       
    (options, args) = parser.parse_args()
    
    if options.quiet:
        log = None
    else:
        log = sys.stderr
    
    inPath1 = None
    inPath2 = None
    outPath1 = None
    outPath2 = None
    
    if len(args) < 1:
        parser.error("Please give at least an input file.")
    if len(args) > 0:
        inPath1 = args[0]
    if len(args) > 1:
        outPath1 = args[1]
    if len(args) > 2:
        inPath2 = args[2]
    if len(args) > 3:
        outPath2 = args[3]
    if len(args) > 4 and log:
        log.write("Additional arguments will be ignored!\n")
            
    format = "fasta"
    
    if options.fastq:
        format = "fastq"
    
    stat = []
    d = o = None
    if options.const > 0 or options.beginConst > 0:
        if log:
            log.write("Removing %i bases form the beginning and %i from the "
                      "end of each read\n" % (options.beginConst, options.const)
                      )
        d, o = trimByConstant(inPath1, outPath1, inPath2, outPath2, 
                              cEnd=options.const, cBegin=options.beginConst, 
                              outFormat="fastq", inFormat="fastq")
        
    elif not options.crop is None:
        if log:
            log.write("Cropping all reads to %i.\n" % (options.crop))
        d, o = trimToConstant(inPath1, outPath1, inPath2, outPath2, 
                              length=options.crop, outFormat=format, stats=stat)
    elif not options.maxProb is None:
        if log:
            log.write("Quality trimming to overall error probability of %f. "
                      "Minimum length after timming %i.\n" \
                       % (options.maxProb, options.minLength))
        d, o = trimToErrorProb(inPath1, outPath1, inPath2, outPath2, 
                               prob=options.maxProb,  minLen=options.minLength,
                               outFormat=format, stats=stat)
    elif options.minQuality != 0:
        if log:
            log.write("Quality trimming removing from the end and beginning"
                      " what has average quality below %f. Minimum length after"
                      " timming %i.\n" % (options.minQuality, options.minLength)
                     )
        d, o = trimByQuality(inPath1, outPath1, inPath2, outPath2, 
                             minQual=options.minQuality, 
                             minLen=options.minLength, 
                             bothEnds=True, outFormat=format, stats=stat)
    if log:
        log.write("Removed %.2f%% (%i) of %i pairs\n" % (float(d)/o*100, d,o))
        
    if stat:
        out = open("stat.tsv","w")
        for line in stat:
            out.write("%s\t%s\n" % line)
        out.close()
