from subprocess import Popen, PIPE, check_output

from Bio import SeqIO
from snakemake.utils import R

configfile: "config.json"

rule all:
    input: expand(["{sample}.full.good.unique.abund.otus.krona.html", "{sample}.full.good.unique.abund.otus.rarefaction.pdf"], sample=config["samples"])

rule extractFull:
    input: "{sample}.full_and_partial.fasta"
    output: "{sample}.full.fasta"
    run:
        with open(output[0], "w") as out:
            for rec in SeqIO.parse(open(input[0]), "fasta"):
                if rec.description.split(" ",1)[1].startswith("Extracted Full"):
                    out.write(rec.format("fasta"))

rule goodReads:
    input: "{sample}.full.fasta"
    output: "{sample}.full.good.fasta"
    shell:
        "%(mothur)s \"#screen.seqs(fasta={input}, maxambig=0, maxn=0);\"" % config
    
rule uniqueReads:
    input: "{sample}.full.good.fasta"
    output: "{sample}.full.good.unique.fasta", "{sample}.full.good.names"
    shell:
        "%(mothur)s \"#unique.seqs(fasta={input});\"" % config 
    
rule nonSigletonReads:
    input: fasta = "{sample}.full.good.unique.fasta", names = "{sample}.full.good.names"
    output: "{sample}.full.good.unique.abund.fasta", "{sample}.full.good.abund.names"
    shell:
        "%(mothur)s \"#split.abund(fasta={input.fasta}, name={input.names}, cutoff=1);\"" % config
    
rule countReads:
    input: "{sample}.full.good.abund.names"
    output: "{sample}.full.good.abund.count_table"
    shell:
        "%(mothur)s \"#count.seqs(name={input})\"" % config
        
rule prepareForUsearch:
    input: fasta = "{sample}.full.good.unique.abund.fasta", counts = "{sample}.full.good.abund.count_table"
    output: "{sample}.full.good.unique.abund.size.fasta"
    run:
        seqCount = {}

        with open(input.counts) as counts:
            #skip header
            next(counts)
            for line in counts:
                seqId, count = line.strip().split("\t")
                seqCount[seqId] = int(count)
            
        reads = list(SeqIO.parse(open(input.fasta), "fasta"))
        sortedReads = sorted(reads, key=lambda x: seqCount[x.id], reverse=True)

        with open(output[0], "w") as outFasta:
            for rec in sortedReads:
                rec.id = "%s;size=%i" % (rec.id, seqCount[rec.id])
                outFasta.write(rec.format("fasta"))

rule usearch:
    input: "{sample}.full.good.unique.abund.size.fasta"
    output: fasta= "{sample}.full.good.unique.abund.otus.size.fasta", details = "{sample}.full.good.unique.abund.uparseout.up"
    shell:
        "%(usearch)s -cluster_otus {input} -otus {output.fasta} -sizein -sizeout -uparseout {output.details}" % config
        
rule cleanUsearchFasta:
    input: fasta = "{sample}.full.good.unique.abund.otus.size.fasta"
    output: fasta = "{sample}.full.good.unique.abund.otus.fasta"
    run:
        with open(output.fasta,"w") as out:
            for rec in SeqIO.parse(open(input.fasta), "fasta"):
                rec.id = rec.id.split(";")[0]
                out.write(rec.format("fasta"))
        
rule usearchStats:
    input: names = "{sample}.full.good.abund.names", uout = "{sample}.full.good.unique.abund.uparseout.up"
    output: oList = "{sample}.full.good.unique.abund.otus.list", counts = "{sample}.full.good.unique.abund.otus.counts.tsv", readOtuMap = "{sample}.full.good.unique.abund.read_otu.mapping.tsv"
    run:
        listTab = {}
        namesTab = {}

        with open(input.names) as names:
            for line in names:
                name, reads = line.strip().split("\t")
                namesTab[name] = reads.split(",")

        with open(input.uout) as up:
            for line in up:
                seqIdsize, cls, ident, chimScr, otu = line.strip().split("\t")[:5]
                if cls == "otu":
                    seqId = seqIdsize.split(";")[0]
                    listTab[seqId] = namesTab[seqId]
                elif cls == "match":
                    seqId = seqIdsize.split(";")[0]
                    listTab[otu].extend(namesTab[seqId])
                else:
                    print("chimera")
                    
        with open(output.oList, "w") as out:
            out.write("0.03\t%i\t" % len(listTab))
            out.write("\t".join([",".join(seqList) for seqList in listTab.values()]))
            
        with open(output.counts, "w") as out:
            for name, seqList in listTab.items():
                out.write("%s\t%i\n" % (name, len(seqList)))
                
        with open(output.readOtuMap, "w") as out:
            for otu, seqList in listTab.items():
                for seq in seqList:
                    out.write("%s\t%s\n" % (seq, otu))

rule computRarefaction:
    input: "{sample}.full.good.unique.abund.otus.list"
    output: "{sample}.full.good.unique.abund.otus.rarefaction"
    shell:
        "%(mothur)s \"#rarefaction.single(list={input})\"" % config
        
rule plotRarefaction:
    input: "{sample}.full.good.unique.abund.otus.rarefaction"
    output: "{sample}.full.good.unique.abund.otus.rarefaction.pdf"
    run:
        R("""library(ggplot2)
d=read.table("%s", header=T)

p = ggplot(d, aes(numsampled)) + geom_point(aes(y=X0.03), colour="blue", size=1) + geom_ribbon(aes(ymin=lci, ymax=hci), alpha=0.2) + xlab("number of reads sampled") + ylab("number of OTUs observed") + ggtitle("Rarefaction curve - %s")
ggsave("%s", p)""" % (input[0], "{wildcards.sample}", output[0]))

rule blast:
    input: "{sample}.full.good.unique.abund.otus.fasta"
    output: "{sample}.full.good.unique.abund.otus.blast.unite.out"
    threads: 6
    shell:
        "export BLASTDB=%(blast_db_folder)s; %(blast_bin_folder)s/blastn -query {input} -task blastn -db unite_db -out {output} -outfmt 6 -evalue %(blast_threshold)f -num_threads {threads}" % config

rule megan:
    input: blast = "{sample}.full.good.unique.abund.otus.blast.unite.out", fasta = "{sample}.full.good.unique.abund.otus.fasta"
    output: megan = "{sample}.full.good.unique.abund.otus.assignment.rma", assign = "{sample}.full.good.unique.abund.otus.assignment"
    log: "logs/{sample}_megan.log"
    run:
        read_number = int(check_output('grep -c ">" %s' % input.fasta, shell=True))
        #TODO check how the identity filter influences the result
        megan_cmd = "load synonymsFileTaxonomy='%s'; import blastFile='%s' fastaFile='%s' meganFile='%s' maxMatches=100 minScore=50.0 maxExpected=0.01 topPercent=10.0 minSupport=1 minComplexity=0.44 useMinimalCoverageHeuristic=false useSeed=false useCOG=false useKegg=false paired=false useIdentityFilter=true textStoragePolicy=Embed blastFormat=BlastTAB mapping='Taxonomy:SYNONYMS_MAP=true';\nset totalReads=%i;\ncollapse rank='Species';\nselect nodes=all;\nexport what=CSV format=readname_taxonid separator=tab file='%s';\nquit;\n" % (config["synonym_file"], input.blast, input.fasta, output.megan, read_number, output.assign)
        #print(megan_cmd)
        with open(log, "w") as out:
            proc = Popen(["xvfb-run", "--auto-servernum", config["megan"], "-g", "-E"], stdin=PIPE, stderr=out.buffer, stdout=out.buffer)
            proc.communicate(megan_cmd.encode("ascii"))
rule krona:
    input: assign = "{sample}.full.good.unique.abund.otus.assignment", counts = "{sample}.full.good.unique.abund.otus.counts.tsv"
    output: "{sample}.full.good.unique.abund.otus.krona.html"
    shell:
        "%(krona)s -o {output} {input.assign}:{input.counts}" % config
