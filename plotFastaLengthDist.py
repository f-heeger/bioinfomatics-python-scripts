import os, sys
try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False
from optparse import OptionParser
from subprocess import Popen, PIPE

from Bio.SeqIO import parse


def getStats(stream, inFormat, log=None):
    """Extract sequence ID and length of the sequecene from fast(a|q)"""
    stats = {"length": []}
    for rec in parse(stream, inFormat):
        stats["length"].append((rec.id, len(rec)))
    return stats


def plotLengthDist(stat, outFile, marks=None, logY=False):
    """Generates a temporary R script and calls R, pipes in the data generate 
    plots this way"""
    # template for R script with solots for outfile, marks and logarithmic y axis
    r_tmpl = """library(ggplot2)
library(scales)

data = read.table(file("stdin"), header=F)
colnames(data) = c("seqId", "len")

pdf("%(out)s")
ggplot(data, aes(x=len))%(mark)s + geom_histogram(binwidth=1, color="black")%(log)s
ggplot(data, aes(x=len))%(mark)s + geom_density()
dev.off()
"""

    #buil dictonary with additional R code for plotting options
    opt = {"out": outFile}
    opt["mark"] = ""
    if not marks is None:
        for i, m in enumerate(marks):
            opt["mark"] += "+ geom_vline(xintercept=%(m)i, color=palette()[2+%(i)i], alpha=1/2)  + annotate(\"text\", label=\"%(m)i\", x=%(m)i, y=0, size=4, color=palette()[2+%(i)i])" % {"m": m, "i": i}
    if logY:
        opt["log"] = "+ scale_y_log10(breaks = trans_breaks(\"log10\", function(x) 10^x), labels = trans_format(\"log10\", math_format(10^.x)))"
    else:
        opt["log"] = ""
    with open("r_plot.tmp", "w") as rScript:
        rScript.write(r_tmpl % opt)
    #start R process with the script
    p = Popen(["Rscript", "r_plot.tmp"], stdin=PIPE)
    #pipe in data in tabular format
    for seqId, length in stats["length"]:
        p.stdin.write("%s\t%i\n" % (seqId, length))
    p.stdin.close()
    p.wait()
    os.remove("r_plot.tmp")

if __name__ == "__main__":

    usage = "usage: %prog [options] input.fasta"

    parser = OptionParser(usage)
    
    parser.add_option("-u", "--fastq",
                       action="store_true", dest="fastq",
                       default=False, help="input file is fastq",)
    if gzImported:
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
    parser.add_option("-o", "--out-folder",
                      action="store", type="string", dest="outFolder",
                      default=".", help="write output files to folder OUT/PATH",
                      metavar="OUT/PATH")
    parser.add_option("-m", "--mark-value",
                      action="append", type="int", dest="marks",
                      default=None, help="mark postion X on the x-Axis",
                      metavar="X")
    parser.add_option("-t", "--text-output",
                      action="store_true", dest="textout",
                      default=False, help="also write text output",)
    parser.add_option("-l", "--log-yaxis",
                      action="store_true", dest="ylog",
                      default=False, help="create plot with logarithmic y-axis",)

    (options, args) = parser.parse_args()
    
    fileType = "fasta"
    if options.fastq:
        fileType = "fastq"
    
    if gzImported and options.gzip:
        inStream = gzip.open(args[0], "r")
    else:
        inStream = open(args[0], "r")
    
    stats = getStats(inStream, fileType, sys.stderr)
    pdfPath = os.path.join(options.outFolder, "%s_lengthDist.pdf" \
                           % args[0].rsplit(".fast", 1)[0].rsplit("/", 1)[-1])
    try:
        plotLengthDist(stats, pdfPath, options.marks, options.ylog)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.stderr.write("Problem starting R to do the plotting\n")
            exit(-1)
        else:
            raise
    #additionally write table of sequence ID and length to a file
    if options.textout:
        textPath = os.path.join(options.outFolder, 
                                "%s_lengthDist.tsv" \
                                % args[0].rsplit(".fast", 1)[0].rsplit("/", 1)[-1])
        with open(textPath, "w") as outFile:
            for seqId, length in stats["length"]:
               outFile.write("%s\t%i\n" % (seqId, length))
    
