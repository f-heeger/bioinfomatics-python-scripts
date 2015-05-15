import json
import subprocess
import sys

from ruffus import *
from ruffus.combinatorics import *

from plotFastaLengthDist import getStats

def loadConfig(path="config.json"):
    #TODO: handle parsing errors gracefully
    return json.load(open(path))

primerFiles = ["lsu_ssu.txt"]

@transform(primerFiles, suffix(".txt"), ".out", loadConfig())
def runTntBlast(primerFile, tntOut, conf):
    for dbPath in conf["path"]["dbs"]:
        subprocess.call([conf["path"]["tntblast"], 
                         "-l", str(conf["param"]["max_len"]),
                         "-e", str(conf["param"]["min_temp"]),
                         "-E", str(conf["param"]["max_temp"]),
                         "-i", primerFile,
                         "-d", dbPath,
                         "-o", tntOut
                        ])

#extract amplicons
@transform(runTntBlast, suffix(".out"), ".amplicon.fasta")
def extractAmplicons(tntOut, ampliconFasta):
    with open(tntOut, "r") as inFile, open(ampliconFasta, "w") as outFile:
        while True:
            try:
                inFile.next()
                if line[0] == ">":
                    outFile.write(line)
                    try:
                        seqLine = inFile.next()
                        outFile.write(seqLine)
                    except StopIteration:
                        raise ValueError("TntBlast output file contains "
                                         "incomplete sequence data. No sequence"
                                         " found after header line: %s" % line)
            except StopIteration:
                break

#get length statistics
@transform(extractAmplicons, suffix(".amplicon.fasta"), "_lengthDist.tsv")
def getLengthStats(inFasta, outTsv):
    stats = getStats(open(inFasta, "r"), "fasta", sys.stderr)
    with open(outTsv, "w") as out:
        for seqId, length in stats["length"]:
               outFile.write("%s\t%i\n" % (seqId, length))

#generate R script for plotting length distibution
@transform(getLengthStats, suffix(".tsv"), ".R")
def generateRScript(inTsv, outR): 
    r_tmpl = """library(ggplot2)
library(scales)

args = commandArgs(trailingOnly = TRUE)

data = read.table(file("stdin"), header=F)
colnames(data) = c("seqId", "len")

pdf(%s.pdf)
ggplot(data, aes(x=len)) + geom_histogram(binwidth=1, color="black")
ggplot(data, aes(x=len)) + geom_density()
dev.off()
"""
    with open(outR, "w") as rScript:
        rScript.write(r_tmpl % outR.rsplit(".", 1)[0])
        
#plot length distibution with R
@transform(generateRScript, suffix(".R"), ".pdf")
def plotWithR(inScript, outPlot):
    subprocess.call(["RScript", inScript, outPlot])

    #get taxonomy distribution
    #plot taxonomy distribution
    
#conf = loadConfig()
pipeline_run(multiprocess=4)
