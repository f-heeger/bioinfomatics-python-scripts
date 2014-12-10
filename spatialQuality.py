import sys
import os
import re
from optparse import OptionParser
try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False

try:
    from Bio.SeqIO import parse
except ImportError:
    sys.stderr.write("ERROR: Could not import biopython.")
    exit(1)
try:
    from rpy2 import robjects
    import rpy2.robjects.lib.ggplot2 as ggplot
    rpyImported = True
except ImportError:
    rpyImported = False
    sys.stderr.write("WARNING: Could not import rpy2. "
                     "No plotting will be available.")

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

def makeDataFrame(stats):
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

#def plotStatsAlpha(data, outFolder, prop="qual", prefix="",
#              high=0.9, low=0.1, pdf=False):
#    #overview plot
#    p = ggplot.ggplot(data) 
#    p = p + ggplot.aes_string(x="x", y="y", alpha=prop) \
#        + ggplot.geom_point(size=0.5) \
#        + ggplot.facet_wrap(robjects.Formula("~ tile")) \
#        + ggplot.scale_alpha(range=robjects.FloatVector([high, low])) \
#        + ggplot.ggtitle("Overview %s" % (prop))
#    fileName = "%soverview_%s.png" % (prefix, prop)
#    p.save(os.path.join(outFolder, fileName), scale=2)
#    
#    #detail plots
#    for t in set(stats["tile"]):
#        p = ggplot.ggplot(data.rx(data.rx2("tile").ro == t, True))
#        p = p + ggplot.aes_string(x="x", y="y", alpha=prop) \
#            + ggplot.geom_point(size=1) \
#            + ggplot.scale_alpha(range=robjects.FloatVector([high, low])) \
#            + ggplot.ggtitle("%i %s" % (t, prop))
#        fileName = "%s%i_%s.png" % (prefix, t, prop)
#        p.save(os.path.join(outFolder, fileName), scale=2)
#        if pdf:
#            fileName = "%s%i_%s.pdf" % (prefix, t, prop)
#            p.save(os.path.join(outFolder, fileName), scale=2)

def plotStats(data, outFolder, tiles, prop="qual", prefix="",
              high="yellow", low="blue", pdf=False, detail=True):
    #overview plot
    p = ggplot.ggplot(data) 
    p = p + ggplot.aes_string(x="x", y="y", col=prop) \
        + ggplot.geom_point(size=0.1) \
        + ggplot.facet_wrap(robjects.Formula("~ tile")) \
        + ggplot.scale_colour_gradient(high=high, low=low) \
        + ggplot.ggtitle("Overview %s" % (prop))
    if prefix:
        fileName = "%s_overview_%s.png" % (prefix, prop)
    else:
        fileName = "overview_%s.png" % (prop)
    p.save(os.path.join(outFolder, fileName), scale=2)
    
    #detail plots
    if detail:
        for t in tiles:
            p = ggplot.ggplot(data.rx(data.rx2("tile").ro == t, True))
            p = p + ggplot.aes_string(x="x", y="y", col=prop) \
                + ggplot.geom_point(size=1) \
                + ggplot.facet_wrap(robjects.Formula("~ tile")) \
                + ggplot.scale_colour_gradient(high=high, low=low) \
                + ggplot.ggtitle("%i %s" % (t, prop))
            if prefix:
                fileName = "%s_%i_%s.png" % (prefix, t, prop)
            else:
                fileName = "%i_%s.png" % (t, prop)
            p.save(os.path.join(outFolder, fileName), scale=2)
            if pdf:
                fileName = "%s%i_%s.pdf" % (prefix, t, prop)
                p.save(os.path.join(outFolder, fileName), scale=2)

def writeStats(s, outFile):
    with open(outFile, "w") as out:
        for i in range(len(s["rId"])):
            out.write("%s\t%s\t%i\t%i\t%f\t%i\n" 
                      % (s["rId"][i], s["tile"][i], s["x"][i], s["y"][i],
                         s["qual"][i], s["n_count"][i]))

def main(inStream, outFolder, name, pdf, detail, ncount, log):
    log.write("Gathering statistics...\n")
    stats = getReadStats(inStream)
    log.write("Writing statistics...\n")
    writeStats(stats, os.path.join(outFolder, "%s_spatialStats.tsv" % name))
    if rpyImported:
        log.write("Building data frame...\n")
        dataFrame = makeDataFrame(stats)
        log.write("Plotting quality...\n")
        tiles = set(stats["tile"])
        plotStats(dataFrame, outFolder, tiles, prefix=name, prop="qual",pdf=pdf)
        if ncount:
            log.write("Plotting n counts...\n")
            plotStats(dataFrame, outFolder, tiles, prefix=name, prop="n_count", 
                      high="blue", low="yellow", pdf=pdf, detail=detail)

if __name__ == "__main__":

    usage = "usage: %prog [options] input.fastq[.gz] [input2.fastq[.gz] ...]"
    epi = "Plot properties of reads according to position on the flow cell"
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
        parser.add_option("-d", "--detail-plot",
                          action="store_true", dest="detail", default=False, 
                          help="make detail plot for each tile",)
        parser.add_option("-p", "--pdf",
                      action="store_true", dest="pdf", default=False, 
                      help="output pdf files additional to png files",)
    
    (options, args) = parser.parse_args()
    
    if options.quiet:
        log = Log(None)
    else:
        log = Log(sys.stderr)
    
    if gzImported and options.gzip:
        log.write("Using gzip compression for input...\n")
        tOpen = gzip.open
    else:
        tOpen = open
    
    if not os.path.exists(options.outFolder):
        log.write("Creating output folder...\n")
        os.makedirs(options.outFolder)
    
    if len(args) == 0:
        log.write("Reading from stdin...\n")
        inStream = sys.stdin
        main(inStream, options.outFolder, "stdin", options.pdf, 
             options.ncount, log)
    else:
        for i, inFile in enumerate(args):
            log.write("Reading from file %s...\n" % inFile)
            inStream = tOpen(inFile)
            m = re.match("(.*)_S\d*_L\d{3}_(R[12])_\d{3}.fastq(.gz)?", 
                         os.path.basename(inFile))
            try:
                name = "%s_%s" % (m.group(1), m.group(2))
            except AttributeError:
                log.write("WARNING. Could not parse file name '%s'. "
                          "Falling back to position (%i) in command line for "
                          "output naming." % (inFile, i))
                name = str(i)
            main(inStream, options.outFolder, name, options.pdf, options.detail,
                 options.ncount, log)
            inStream.close()
    log.write("Everything done...\n")
