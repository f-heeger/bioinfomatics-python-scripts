import os, sys
try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False
from optparse import OptionParser
from subprocess import Popen, PIPE

from Bio.SeqIO import parse

from NcbiTools import SingleLevelLineageMap, CachedNuclId2TaxIdMap


def getData(stream, inFormat, ranks, log=None):
    """Extract sequence ID and length of the sequecene from fast(a|q)"""
    GI2TAX_DBPATH = "/home/heeger/data/dbs/gi_taxid_nucl.db"
    mapping = {}
    gi2spec = CachedNuclId2TaxIdMap(GI2TAX_DBPATH)
    spec2rank = {}
    for rank in ranks:
        spec2rank[rank] = SingleLevelLineageMap(rank, cachePath="spec2%s.csv" % rank, retry=1)
        mapping[rank] = []
    for rec in parse(stream, inFormat):
        gi = rec.id.split("|", 2)[1]
        try:
            spec = gi2spec[gi]
        except KeyError:
            for rank in ranks:
                mapping[rank].append((rec.id, len(rec), "NA"))
        else:
            for rank in ranks:
                d = spec2rank[rank][spec]
                if len(d) > 1:
                    raise ValueError("More than one taxonomy group of "
                                     "rank '%s' for gi '%s'" % (rank, gi))
                if len(d) == 0:
                    tax = "NA"
                else:
                    tax = d[0][2]
                mapping[rank].append((rec.id, len(rec), tax))
    for lineageMap in spec2rank.values():
        lineageMap.save()
    return mapping


def plotLengthDist(stat, outFile, marks=None, logY=False, title=None):
    """Generates a temporary R script and calls R, pipes in the data, generate 
    plots this way"""
    # template for R script with solots for outfile, marks and logarithmic y axis
    r_tmpl = """library(ggplot2)
library(scales)

data = read.table(file("stdin"), header=F)
colnames(data) = c("seqId", "len", "taxon")

pdf("%(out)s")
ggplot(data, aes(x=len, color=taxon, fill=taxon))%(mark)s + geom_histogram(binwidth=10)%(log)s
ggplot(data, aes(x=len, color=taxon))%(mark)s + geom_density(aes(len, ..count..), alpha=.5)
dev.off()
"""

    #buil dictonary with additional R code for plotting options
    opt = {"out": outFile}
    opt["title"] = ""
    if not title is None:
        opt["title"] = "+ ggtitle(%s)" % title
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
    for seqId, length, tax in stat:
        p.stdin.write("%s\t%i\t%s\n" % (seqId, length, tax))
    p.stdin.close()
    p.wait()
#    os.remove("r_plot.tmp")

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
    parser.add_option("-r", "--rank",
                      action="append", type="string", dest="ranks",
                      default=None, 
                      help="color length plots by this taxonomic level [default: kingdom]",
                      metavar="RANK")

    (options, args) = parser.parse_args()
    
    log = sys.stderr
    
    fileType = "fasta"
    if options.fastq:
        fileType = "fastq"
    
    if not log is None:
        log.write("File type is %s\n" % fileType)
    
    if gzImported and options.gzip:
        inStream = gzip.open(args[0], "r")
        if not log is None:
            log.write("Reading compressed file.\n")
    else:
        inStream = open(args[0], "r")
    
    if len(options.ranks) == 0:
        if not log is None:
            log.write("No rank given defaulting to 'kingdom'\n")
        ranks = ["kingdom"]
    else:
        ranks = options.ranks
    if not log is None:
        log.write("Reading data...\n")
    mapping = getData(inStream, fileType, ranks, sys.stderr)
    for rank, d in mapping.items():
        if not log is None:
            log.write("Plotting %i data points for rank '%s'...\n" % (len(d),rank))
        pdfPath = os.path.join(options.outFolder, "%s_%s_lengthDist.pdf" \
                               % (args[0].rsplit(".fast", 1)[0].rsplit("/", 1)[-1],
                                  rank))
        try:
            plotLengthDist(d, pdfPath, options.marks, options.ylog, rank)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                sys.stderr.write("Problem starting R to do the plotting\n")
                exit(-1)
            else:
                raise
    #additionally write table of sequence ID and length to a file
    if options.textout:
        if not log is None:
            log.write("Writing text output...\n")
        textPath = os.path.join(options.outFolder, 
                                "%s_lengthDist.tsv" \
                                % args[0].rsplit(".fast", 1)[0].rsplit("/", 1)[-1])
                                  
        with open(textPath, "w") as outFile:
            for rank, d in mapping.items():
                for seqId, length, tax in d:
                   outFile.write("%s\t%i\t%s\t%s\n" % (seqId, length, rank, tax))
    
