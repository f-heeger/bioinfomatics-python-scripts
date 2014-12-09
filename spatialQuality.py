import sys
import os
from optparse import OptionParser
try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False

from Bio.SeqIO import parse
try:
    from rpy2 import robjects
    import rpy2.robjects.lib.ggplot2 as ggplot
    rpyImported = True
except ImportError:
    rpyImported = False
    sys.stderr.write("Could not import rpy2. No plotting will be available")

class Log:
    def __init__(self, logStream=None):
        self.stream = logStream
        
    def write(self, string):
        if not self.stream is None:
            self.stream.write(string)

def getReadStats(inStream):
    readStats = {"rId": [], "lane": [], "tile": [], "x": [], "y": [], 
                 "qual": [], "n_count": [], "length": []}
    for rec in parse(inStream, "fastq"):
        #position
        idArr = rec.id.split(" ")[0].split(":")
        lane = int(idArr[3])
        tile = int(idArr[4])
        x = int(idArr[5])
        y = int(idArr[6])
        readStats["rId"].append(rec.id)
        readStats["lane"].append(lane)
        readStats["tile"].append(tile)
        readStats["x"].append(x)
        readStats["y"].append(y)
        #mean quality
        qual = float(sum(rec._per_letter_annotations["phred_quality"]))/len(rec)
        readStats["qual"].append(qual)
        #number of Ns
        nCount = rec.seq.count("N")
        readStats["n_count"].append(nCount)
        #length
        length = len(rec)
        readStats["length"].append(length)
    return readStats

def makeDataFrame(stat):
    #prepare data
    rId = robjects.StrVector(stats["rId"])
    tile = robjects.IntVector(stats["tile"])
    x = robjects.IntVector(stats["x"])
    y = robjects.IntVector(stats["y"])
    q = robjects.FloatVector(stats["qual"])
    n = robjects.IntVector(stats["n_count"])
    data = robjects.DataFrame({"rid": rId, "tile": tile, "x": x, "y": y,
                               "qual": q, "n_count": n})
    return data

def plotStats(data, outFolder, prop="qual", prefix="",
              high="blue", low="red"):
    #overview plot
    p = ggplot.ggplot(data) 
    p = p + ggplot.aes_string(x="x", y="y", col=prop) \
        + ggplot.geom_point(size=0.5) \
        + ggplot.facet_wrap(robjects.Formula("~ tile")) \
        + ggplot.scale_colour_gradient(high=high, low=low) \
        + ggplot.ggtitle("Overview %s" % (prop))
    fileName = "%soverview_%s.png" % (prefix, prop)
    p.save(os.path.join(outFolder, fileName), scale=2)
    
    #detail plots
    for t in set(stats["tile"]):
        p = ggplot.ggplot(data.rx(data.rx2("tile").ro == t, True))
        p = p + ggplot.aes_string(x="x", y="y", col=prop) \
            + ggplot.geom_point(size=1) \
            + ggplot.facet_wrap(robjects.Formula("~ tile")) \
            + ggplot.scale_colour_gradient(high=high, low=low) \
            + ggplot.ggtitle("%i %s" % (t, prop))
        fileName = "%s%i_%s.pdf" % (prefix, t, prop)
        p.save(os.path.join(outFolder, fileName), scale=2)
        fileName = "%s%i_%s.png" % (prefix, t, prop)
        p.save(os.path.join(outFolder, fileName), scale=2)

def writeStats(s, outFile):
    with open(outFile, "w") as out:
        for i in range(len(s["rId"])):
            out.write("%s\t%s\t%i\t%i\t%f\t%i\n" 
                      % (s["rId"][i], s["tile"][i], s["x"][i], s["y"][i],
                         s["qual"][i], s["n_count"][i]))

if __name__ == "__main__":

    usage = "usage: %prog [options] input.fastq[.gz] [input2.fastq[.gz] ...]"
    epi = "Plot properties of reads according to position in the flow cell"
    parser = OptionParser(usage, epilog=epi)
    
    parser.add_option("-q", "--quite",
                      action="store_true", dest="quiet", default=False, 
                      help="do not print status messages to the screen",)
    parser.add_option("-o", "--output-folder",
                      action="store", type="string", dest="outFolder",
                      default="spatialQual", 
                      help="write results to folder X [default: ./saptialQual]",
                      metavar="X")
    if gzImported:
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
    if rpyImported:
        parser.add_option("-n", "--n-count",
                          action="store_true", dest="ncount", default=False, 
                          help="also make plots for N count",)
    
    (options, args) = parser.parse_args()
    
    if options.quiet:
        log = Log(None)
    else:
        log = Log(sys.stderr)
    
    if gzImported and options.gzip:
        log.write("Using gzip compression for input files...\n")
        tOpen = gzip.open
    else:
        tOpen = open
    
    if not os.path.exists(options.outFolder):
        log.write("Creating output folder...\n")
        os.makedirs(options.outFolder)
    
    for inFile in args:
        log.write("Working on file %s...\n" % inFile)
        log.write("Gathering statistics...\n")
        stats = getReadStats(tOpen(inFile))
        log.write("Writing statistics...\n")
        writeStats(stats, os.path.join(options.outFolder, "spatialStats.tsv"))
        if rpyImported:
            log.write("Building data frame...\n")
            dataFrame = makeDataFrame(stats)
            log.write("Plotting quality...\n")
            plotStats(dataFrame, options.outFolder, prop="qual")
            if options.ncount:
                log.write("Plotting n counts...\n")
                plotStats(dataFrame, options.outFolder, prop="n_count", 
                          high="red", low="lightblue")
    
