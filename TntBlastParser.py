def tntBlastParser(inStream):
    """Smiple parser for standart output of tntblast
    
    returns per hit a tuple of the form: (Ncbi GI, aplification start, 
    amplification end, amplified sequence)"""
    l=0 #line number INSIDE THE CURRENT ENTRY
    #line number 0: the first line in the file is filled with #
    # the first line of every other entry will be empty
    for line in inStream:
        if l == 1:
            #primer name (ignored for now)
            name = line.split("=")[1].strip()
        elif l == 23:
            #amplification start and end
            alnRange = line.split("=")[1].strip()
            alnStart, alnEnd = alnRange.split(" .. ")
        elif l == 24:
            #amplificon length
            length = line.split("=")[1].strip()
        elif l == 33:
            #gi
            seqId = line[1:].split(" ")[0]
            gi = seqId.split("|")[1]
        elif l == 34:
            #sequence
            seq = line.strip()
            #trigger yield since this is the end of one entry
            yield (gi, alnStart, alnEnd, seq)
            l = -1 #set line to -1, because it will be incremented and is than 0
        l += 1
