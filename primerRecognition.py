import sys
if sys.version_info[0] < 3:
    raise RuntimeError("Please use Python 3 to run this script")

import argparse
import gzip
import glob

from Bio import SeqIO

iupac = {
"A": "A",
"C": "C",
"G": "G",
"T": "T",
"AC": "M",
"AG": "R",
"AT": "W",
"CG": "S",
"CT": "Y",
"GT": "K",
"ACG": "V",
"ACT": "H",
"AGT": "D",
"CGT": "B",
"ACGT": "N"
}

def readData(inFileList, length, gz=True, maxReads=0):
    data = {}
    if maxReads == 0:
        maxReads = float("Inf")
    for inFile in inFileList:
        data[inFile] = []
        for i in range(length):
            data[inFile].append({"A": 0, "C": 0, "G": 0, "T": 0, "N": 0})
        if gz:
            oFunc = gzip.open
        else:
            oFunc = open
        with oFunc(inFile, "rt") as inStream:
            for r, rec in enumerate(SeqIO.parse(inStream, "fastq")):
                if r > maxReads:
                    break
                for i in range(length):
                    data[inFile][i][rec.seq[i]] += 1
    
    return data

def primerPerFile(data, r=0.9):
    primers = {}
    for inFile, posData in data.items():
        primers[inFile] = ""
        for counts in posData:
            maxChar = set()
            countSum = sum(counts.values())
            acum = 0.0
            for char, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
                if char == "N":
                    #Ns are save to ignore as they will not change the outcome
                    continue
                acum += count/countSum
                maxChar.add(char)
                if acum >= r:
                    break
            primers[inFile] += iupac["".join(sorted(list(maxChar)))]
    return primers

def saveData(data, opt, path):
    with open(path, "w") as out:
        if opt[1] > 0:
            out.write("# nucleotide distribution in the first %i bases, based on %i reads per sample. The IUPAC code in each position represents %.2f%% of the bases.\n" % opt)
        else:
            out.write("# nucleotide distribution in the first %i bases, based on all available reads per sample. The IUPAC code in each position represents %.2f%% of the bases.\n" % (opt[0], opt[2]))
        for inFile, posData in data.items():
            for pos, counts in enumerate(posData):
                for char, count in counts.items():
                    out.write("%s\t%i\t%s\t%i\n" % (inFile, pos, char, count))

def loadData(path):
    data = {}
    for line in open(path):
        if line[0] == "#":
            continue
        inFile, posStr, char, countStr = line.strip("\n").split("\t")
        pos = int(posStr)
        count = int(countStr)
        if inFile not in data:
            data[inFile] = []
        if len(data[inFile]) <= pos:
            data[inFile].append(dict(zip("ACGTN",[None]*5)))
        data[inFile][pos][char] = count
    return data


if __name__ == "__main__":
    def proportion(n):
        try:
            f=float(n)
        except ValueError:
            raise argparse.ArgumentError("Parameter -r: '%s' can not be convereted into a floating point number" % n)
        if not 0<f<=1:
            raise argparse.ArgumentError("Parameter -r is '%f', but should be between 0 and 1.")
        return f

    aParser = argparse.ArgumentParser(description="Try to recognize a possible primer sequence from the starting bases of a read file in fastq[.gz] format")
    aParser.add_argument("length", type=int,
                         help="up to which position should be analysed for the primer")
    aParser.add_argument("path", help="where to look for fastq files")
    aParser.add_argument("-o", "--out", help="write results to this file",
                         default=None)
    aParser.add_argument("-g", "--gz", help="input files are gzipped", 
                         action="store_true", default=True)
    aParser.add_argument("-p", "--pre-computed", dest="pre",
                         help="give a file of precomputed data here",
                         default=None)
    aParser.add_argument("-m", "--max-reads", dest="max", type=int,
                         help="maximum number of reads to read per file (0 for all)",
                         default=0)
    aParser.add_argument("-r", "--represented", type=proportion, default=0.9,
                         help="Which proportion of the nucleotides in a postiotion should be represented by the primer sequence (0-1, deafault 0.9)")
    
    args = aParser.parse_args()
    
    if args.out is None:
        out = sys.stdout
    else:
        out = open(args.out, "w")
    try:
        #load/read data
        if args.pre:
            data = loadData(args.pre)
        else:
            if args.gz:
                fList = [f for f in glob.glob(args.path) if f.endswith(".fastq.gz")]
            else:
                fList = [f for f in glob.glob(args.path) if f.endswith(".fastq")]
            data = readData(fList, args.length, args.gz, args.max)
        
        #save data
        saveData(data, (args.length, args.max, args.represented*100), "nucDist.tsv")
        
        #recognize primers
        primers = primerPerFile(data, args.represented)
        fileByPrimer = {}
        for fName, pSeq in primers.items():
            try:
                fileByPrimer[pSeq].append(fName)
            except KeyError:
                fileByPrimer[pSeq] = [fName]
        
        if args.max > 0:
            out.write("# Primer estimation based on the first %i bases of %i reads per sample. The IUPAC code in each position represents %.2f%% of the bases.\n" % (args.length, args.max, args.represented*100))
        else:
            out.write("# Primer estimation based on the first %i bases of all available reads per sample. The IUPAC code in each position represents %.2f%% of the bases.\n" % (args.length, args.represented*100))
        for pSeq, fList in fileByPrimer.items():
            out.write("%s\t%i\t%s\n" % (pSeq, len(fList), ",".join(fList)))
        
    finally:
        if not args.out is None:
            out.close()
