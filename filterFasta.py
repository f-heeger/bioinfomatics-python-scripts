#! /usr/bin/env python

import sys, re
from warnings import warn

try:
    import gzip
    gzImported = True
except ImportError:
    gzImported = False
from random import sample
from optparse import OptionParser

from Bio import SeqIO

def filterFasta(inStream, outPath, minLength=None, idList=None, 
                random=None, fastq=False, regex=False, neg=False, 
                noWarn=False, ignoreAt=False, log=sys.stderr):
    if fastq:
        format = "fastq"
    else:
        format = "fasta"
    if random:
        sampleRandom(inStream, outPath, format, random, log)
    else:
        filterLengthIdList(inStream, outPath, format, minLength, idList, 
                           regex, neg, noWarn, ignoreAt, log)
    

def filterLengthIdList(inStream, outPath, format, minLength=None, 
                       idList=None, regex=False, neg=False, noWarn=False, 
                       ignoreAt=False, log=sys.stderr):
    if not idList is None:
        if ignoreAt:
            log.write("Ignoring all \"@\" signs at the start of IDs")
            idList = [rId.lstrip("@") for rId in idList]
        if regex:
            idRes = [re.compile(x) for x in idList]
        else:
            #use a dict to have random acces in O(1)
            idDict = dict(zip(idList, [0]*len(idList)))
    try:
        out = open(outPath, "w")
    except TypeError: 
        #if the outPath paramerter is not a path but already a stream
        out = outPath
    if log:
        log.write("Filtering and writing...\n")
    try:
        l = 0
        n = 0
        for rec in SeqIO.parse(inStream, format):
            l+=1
            if log and l%1000==0:
                log.write("\r%i records done" % l)
            write = True
            if not minLength is None and len(rec)<minLength: 
                #if a min length was set skip the record if is to short
                write = False
            elif not idList is None:
                if regex:
                    if all([r.match(rec.id) is None for r in idRes]):
                        #if the record id does not match any of the given REs
                        # i.e. all results of matching are None
                        write = False
                else:
                    if rec.id not in idDict:
                        #if a ID list was given skip the record if it is not in it
                        write = False
                    else:
                        idDict[rec.id] += 1
            if (not neg and write) or (neg and not write) :
                out.write(rec.format(format))
                n+=1
        notFound = 0
        multiFound = 0
        if not idList is None and not regex:
            for recId, found in idDict.items():
                if found == 0 and not noWarn:
                    warn("ID '%s' was not found in the file." % recId)
                    notFound +=1
                if found > 1 and not noWarn:
                    warn("ID '%s' encountered %i times in the file" \
                         % (recId, found))
                    multiFound +=1
        if log:
            if notFound > 0:
                log.write("Of the %i given IDs %i were NOT found in the file\n"
                          % (len(idDict), notFound))
            if multiFound > 0:
                log.write("Of the %i given IDs %i were found MULTIPLE TIMES"
                          " in the file\n" % (len(idDict), multiFound))
            log.write("\n%i out of %i sequences remained after filtering.\n"
                      % (n,l))
    finally:
        out.close()   

def sampleRandom(inStream, outPath, format, number, log):
    try:
        try:
            out = open(outPath, "w")
        except TypeError:
            out = outPath
        if log:
            log.write("Sampling %i sequences at random...\n" % number)
            log.write("Counting sequences in file...\n")
        l = 0
        for rec in SeqIO.parse(inStream, format):
            l+=1
        inStream.seek(0, 0)
        if number>l:
            raise ValueError("Sample size(%i) is bigger than number of "
                             "sequences in input file(%i)." % (number, l))
        if log:
            log.write("Sampling %i from %i sequences...\n" % (number, l))
        n = 0 #record number
        i = 0 #index of next sampled record
        ls = sample(range(0,l), number)
        ls.sort()
        for rec in SeqIO.parse(inStream, format):
            if n == ls[i]: #if this record is the nextsampled
                out.write(rec.format(format))
                i+=1
                if i>=len(ls):
                    #stop iteration over the input file if all
                    # sampled records were written
                    break
            n+=1    
    finally:
        out.close()   

if __name__ == "__main__":

    usage = "usage: %prog [options] input.fasta [output.fasta]"

    parser = OptionParser(usage)
    
    parser.add_option("-q", "--quiet",
                       action="store_true", dest="quiet", default=False, 
                       help="do not print status messages to the screen",)
    parser.add_option("-w", "--supress-warnings",
                       action="store_true", dest="noWarn", default=False, 
                       help="do not print warnings if a ID was found more or "
                            "less than onece",)
    parser.add_option("-u", "--fastq",
                       action="store_true", dest="fastq",
                       default=False, help="input file is fastq",)
    if gzImported:
        parser.add_option("-z", "--gzip",
                           action="store_true", dest="gzip",
                           default=False, help="input file is gzipped",)
    parser.add_option("-l", "--min-length",
                       action="store", type="int", dest="minLength",
                       default=None, 
                       help="write only sequence with lengths at least X",
                       metavar="X")
    parser.add_option("-i", "--id-list",
                       action="store", type="string", dest="idList",
                       default=None, 
                       help="write only sequence with an ID from this list. "
                            "List can be comma separated string of IDs or a "
                            "path to a file with a line separated list of IDs",
                       metavar="X")
    parser.add_option("-r", "--random",
                      action="store", type="int", dest="random",
                      default=None, 
                      help="randomly sample X sequence from input file",
                      metavar="X")
    parser.add_option("-e", "--regexp",
                       action="store_true", dest="regexp",
                       default=False, 
                       help="use regular expression instead of exact "
                            "matching for IDs",)
    parser.add_option("-a", "--ignore-at",
                       action="store_true", dest="ignore_at",
                       default=False, 
                       help="ignore the first letter of the query IDs if it is"
                            "an @ (this is for more convinent filter list "
                            "creation from fastq files)",)
    parser.add_option("-n", "--negative",
                       action="store_true", dest="neg",
                       default=False, 
                       help="do exactly the opposite of what would normally "
                            "be done",)
    (options, args) = parser.parse_args()
    
    if (options.idList and options.random):
        parser.error("Options -i and -r are mutually exclusive.")
    if (options.minLength and options.random):
        parser.error("Options -l and -r are mutually exclusive.")
    if (options.regexp and not options.idList):
        parser.error("Options -e can only be used with -i.")
    if (options.random and options.neg):
        parser.error("Negative mode does not work with random mode.")
    
    if options.quiet:
        log = None
    else:
        log = sys.stderr
        
    if len(args) < 1:
        if gzImported and options.gzip:
            parser.error("Pipe mode (no input file argument) does not work together with -z (gzipped input).")
        if log:
            log.write("NOTE: Running in pipe mode. Waiting for input from stdin.\n")
            log.write("Will be writing to stdout.\n")
        out = sys.stdout
    elif len(args) == 1:
        #if no output file was given write to std out
        if log:
            log.write("Will be writing to stdout.\n")
        out = sys.stdout
    else:
        out = args[1]
    
    if len(args)>2:
        if log:
            log.write("Additional arguments will be ignored!\n")
    if options.neg:
        if log:
            log.write("NOTE: Running in negative mode.\n")
    idList = None
    if options.idList:
        idList = []
        try:
            iListFile = open(options.idList)
            for line in iListFile:
                idList.append(line.strip())
            iListFile.close()
        except IOError:
            pass
            if log:
                log.write("ID list parameter is not a valid path. Assume it to "
                          "be comma separated string.\n")
            idList = options.idList.strip().split(",")
    if log:
        if options.regexp and log:
            log.write("Using list of Regular Expression on IDs to filter.\n")
        elif idList and log:
            log.write("Using list of IDs to filter.\n")
                    
    if len(args) == 0:
        inStream = sys.stdin
    elif gzImported and options.gzip:
        inStream = gzip.open(args[0], "r")
    else:
        inStream = open(args[0], "r")
    
    try:                
        filterFasta(inStream, out, options.minLength, idList, options.random, 
                    options.fastq, options.regexp, options.neg, options.noWarn, 
                    options.ignore_at, log=log)
    finally:
        if len(args) > 0:
            inStream.close()
    
