import re
import gzip as gz
import argparse

from Bio import SeqIO

def splitFasta(fastaPath, pieces=2, number=None, sequential=False, form="fasta", gzip=False):
    
    if gzip:
        tOpen = gz.open
    else:
        tOpen = open
        
    if not number is None:
        if not pieces is None:
            raise ValueError("Either piece or number parameter has to be None")
        step = number[0]
    else:
        l=0
        for r in SeqIO.parse(tOpen(fastaPath), form):
            l+=1
        step = (l/pieces[0])+1
    start = 0
    end = step
    p = 0
    
    r=0
    rec = SeqIO.parse(tOpen(fastaPath), form)
    while True:
        if gzip:
            path = fastaPath.rsplit(".", 2)[0]
            ext = form + ".gz"
        else:
            path = fastaPath.rsplit(".", 1)[0]
            ext = form
        if sequential:
            path += "_%i.%s" % (p, ext)
        else:
            path += "_%i.%s" % (start, ext)
        out = open(path ,"w")
        p+=1
        try:
            while r<end:
                out.write(next(rec).format(form))
                r+=1
        except StopIteration:
            break
        finally:
            out.close()
        start=end
        end+=step

def splitMultiFasta(fastaPath, form="fasta", gzip=False):
    
    if gzip:
        tOpen = gz.open
        ext = "%s.gz" % form
    else:
        tOpen = open
        ext = form
    
    for rec in SeqIO.parse(tOpen(fastaPath), form):
        valid_name = "".join(i if i not in "\/:*?<>|" else "_" for i in rec.id)
        out = tOpen("%s.%s" % (valid_name, form), "w")
        try:
            out.write(rec.format(form))
        finally:
            out.close()

if __name__ == "__main__":

#    usage = "usage: python %prog [options] INFILE"

    parser = argparse.ArgumentParser(description="Split fast[a|q] file into multiple files.")

    parser.add_argument("infile", metavar="INFILE", nargs=1,
                        help="fast[a|q] file to split")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-p", "--pieces",
                      action="store", nargs=1, default=None, type=int,
                      help="split sequence to X pieces", metavar="X")
    mode.add_argument("-n", "--number",
                      action="store", nargs=1, type=int, dest="number",
                      default=None, help="split sequence to pieces with X reads",
                      metavar="X")
    mode.add_argument("-s", "--single-record",
                      action="store_true", dest="each",
                      default=False, help="split a every record",)
    parser.add_argument("-u", "--fastq",
                       action="store_true",default=False, dest="fastq",
                       help="input file is in fastq format",)
    parser.add_argument("-z", "--gzip",
                      action="store_true", dest="gzip",
                      default=False, help="input file is gzipped",)

    args = parser.parse_args()
    
#    if len(args)!=1:
#        parser.error("Please give an input file.")
        
#    if sum([not options.pieces is None, options.each, not options.number is None])>1:
#        parser.error("options -p, -n and -s are mutally exclusive")
    
    form = "fasta"
    if args.fastq:
        form = "fastq"
    
    if args.each:
        splitMultiFasta(args.infile[0], form, args.gzip)
    else:
        splitFasta(args.infile[0], args.pieces, args.number, args.sequential, form, args.gzip)
