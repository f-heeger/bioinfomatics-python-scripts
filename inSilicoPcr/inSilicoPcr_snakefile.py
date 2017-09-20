import json
import subprocess
import sys

from snakemake.utils import R

from NcbiTools import CachedNuclId2TaxIdMap
from TntBlastParser import tntBlastParser
from plotFastaLengthDist import getStats

def getTaxIds(ampl, dbpath, email):
    gi2tax = CachedNuclId2TaxIdMap(dbpath, email)
    perTax = {}
    for rec in ampl:
        try:
            tax = gi2tax[rec[0]]
        except KeyError:
            tax = -1
        try:
            perTax[tax].append(rec)
        except KeyError:
            perTax[tax] = [rec]
    return perTax

configfile: "config.json"

#print(str(config))

primer_pairs = list(config["primer_pairs"].keys())
databases = list(config["databases"].keys())

#print(config)

rule all:
    input: expand(["{pair_name}/{pair_name}_vs_{dbname}_lengthDist.pdf", "{pair_name}/{pair_name}_vs_{dbname}.krona.html"], pair_name=primer_pairs, dbname=databases)

rule creatPrimerFile:
    output: "{pair_name}/{pair_name}.txt"
    run:
        pName = wildcards.pair_name
        cmd = 'echo "%s\t%s\t%s\n" > {output}' % (pName, config["primer_pairs"][pName]["forward_seq"], config["primer_pairs"][pName]["reverse_seq"])
        shell(cmd)

rule tntBlast:
    input: primerFile="{pair_name}/{pair_name}.txt"
    output: "{pair_name}/{pair_name}_vs_{dbname}.out"
    params: database="{dbname}"
    threads: 6
    message: "Running tnt blast for primer pair {wildcards.pair_name} on db {wildcards.dbname} with {threads} threads"
    run:
        name = wildcards.pair_name
        props = config["primer_pairs"][name]
        if threads>1:
            cmd = "mpirun -np {threads} "
        else:
            cmd = ""
        cmd += config["tntblast"] + " -l %(max_len)i -e %(min_temp)i -x %(max_temp)i -i {input.primerFile}" % props + " -d %s -o {output}" % config["databases"][params.database]
        shell(cmd)
    
rule extractAmplicons:
    input: "{pair_name}/{pair_name}_vs_{dbname}.out"
    output: "{pair_name}/{pair_name}_vs_{dbname}_amplicon.fasta"
    run:
        with open(input[0], "r") as inFile, open(output[0], "w") as outFile:
            while True:
                try:
                    line = inFile.readline()
                    if not line:
                        break #end of file reached
                    if line[0] == ">":
                        outFile.write(line)
                        seqLine = inFile.readline()
                        if len(seqLine.strip()) == 0:
                            raise ValueError("TntBlast output file contains "
                                             "incomplete sequence data. No sequence"
                                             " found after header line: %s" % line)
                        outFile.write(seqLine)
                except StopIteration:
                    break


rule getLengthStats:
    input: "{pair_name}/{pair_name}_vs_{dbname}_amplicon.fasta"
    output: "{pair_name}/{pair_name}_vs_{dbname}_lengthDist.tsv"
    run: 
        stats = getStats(open(input[0], "r"), "fasta", sys.stderr)
        with open(output[0], "w") as out:
            for seqId, length in stats["length"]:
                   out.write("%s\t%i\n" % (seqId, length))


rule generateLengthPlot:
    input: "{pair_name}/{pair_name}_vs_{dbname}_lengthDist.tsv"
    output: "{pair_name}/{pair_name}_vs_{dbname}_lengthDist.pdf"
    run:
        R("""
        library(ggplot2)
        library(scales)
        library(gridExtra)

        data =  tryCatch({{return(read.table("%(infile)s",  header=F))}}, error=function(e) {{return(data.frame())}})
        
        if (dim(data)[1] == 0) {{
            pdf("%(outfile)s")
            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 8, "No data", cex=2)
            dev.off()
        }} else {{
        
            colnames(data) = c("seqId", "len")

            plotList = list()
            plotList[[1]] = ggplot(data, aes(x=len)) + geom_histogram(binwidth=1, color="black") + ggtitle(expression(atop("Length distribution of predicted amplicons", atop("with primer pair: %(primers)s on database: %(db)s", atop("min. Temp.: %(minT).2f, max. Temp.: %(maxT).2f, max. Len.: %(len)i"))))) + xlab("length")
            plotList[[2]] = ggplot(data, aes(x=len)) + geom_density() + ggtitle(expression(atop("Length distribution of predicted amplicons", atop("with primer pair: %(primers)s on database: %(db)s", atop("min. Temp.: %(minT).2f, max. Temp.: %(maxT).2f, max. Len.: %(len)i"))))) + xlab("length")
            pdf("%(outfile)s")
            print(plotList)
            dev.off()
        }}
        """ % {"infile": input[0], "outfile": output[0], 
               "primers": wildcards.pair_name, 
               "db": wildcards.dbname, 
               "len": config["primer_pairs"][wildcards.pair_name]["max_len"], 
               "maxT":config["primer_pairs"][wildcards.pair_name]["max_temp"],
               "minT": config["primer_pairs"][wildcards.pair_name]["min_temp"],
               }
         )
        
rule getTaxonomyInfo:
    input: "{pair_name}/{pair_name}_vs_{dbname}.out"
    output: "{pair_name}/{pair_name}_vs_{dbname}.taxDist.tsv"
    resources: dbaccess=1
    run:
        ampl = list(tntBlastParser(open(input[0])))
        taxDist = getTaxIds(ampl, config["gi2tax_db"], config["email"])
        with open(output[0], "w") as out:
            for tax, amplList in taxDist.items():
                for gi, start, end, seq in amplList:
                    out.write("%s_%s-%s\t%s\n" % (gi, start, end, tax))
                
rule makeKronaPlots:
    input: "{pair_name}/{pair_name}_vs_{dbname}.taxDist.tsv"
    output: "{pair_name}/{pair_name}_vs_{dbname}.krona.html"
    shell:
        "%(krona)s -o {output} {input}" % config
        
