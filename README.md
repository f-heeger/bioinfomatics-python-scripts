# Collection of Bioinformatic Python Scripts

## Scripts

### assemblyStats
> Usage: assemblyStats.py [options] input.fasta[.gz]  
>
>Options:  
>  -h, --help     show this help message and exit  
>  -z, --gzip     input file is gzipped  
>  -t, --tabular  write tabular output  
>
>Computes basic assembly statistics given that input.fasta is a multi fasta  
>file with one contig per record.

Normally write human readable output. With the -t option a tabular format of the same information will be written, that can be used in statistic evalutation of multiple asssemblies.

### filterFasta

>Usage: filterFasta.py [options] input.fasta [output.fasta]  

>Options:  
>  -h, --help            show this help message and exit  
>  -q, --quite           do not print status messages to the screen  
>  -u, --fastq           input file is fastq  
>  -z, --gzip            input file is gzipped  
>  -l X, --min-length=X  write only sequence with lengths at least X  
>  -i X, --id-list=X     write only sequence with an ID from this list. List  
>                        can be comma separated string of IDs or a path to a  
>                        file with a line separated list of IDs  
>  -r X, --random=X      randomly sample X sequence from input file  
>  -e, --regexp          use regular expression instead of exact matching for  
>                        IDs  
>  -n, --negative        do exactly the opposite of what would normally be done  

Filter fasta/fastq files in different way.

1. Filter by minimum length (`-l`)  
Only write sequences of a certain length to the output.
2. Filter by list of IDs (`-i`)  
Only write sequences with an ID from the list given. The list can either be given as a comma separated list of IDs or as in a file with one ID per line. With the `-e` option the given "IDs" will be used as regular expression isntead of exact matching (using python regex).
3. "Filter" randomly
Write a random subset of the sequences of the given size.

The `-n` option will switch to negative mode. Meaning the script will do exactly the oposite it normally does.

Input data can be provided as a file (first argument) or be piped in.

### plotFastaLengthDist

>Usage: plotFastaLengthDist.py [options] input.fasta  
>
>Options:  
>  -h, --help            show this help message and exit  
>  -u, --fastq           input file is fastq  
>  -z, --gzip            input file is gzipped  
>  -o OUT/PATH, --out-folder=OUT/PATH  
>                        write output files to folder OUT/PATH  
>  -m X, --mark-value=X  mark position X on the x-Axis  
>  -t, --text-output     also write text output  
>  -l, --log-yaxis       create plot with logarithmic y-axis  

Plot the length distribution of the sequences in the input file. 
Two plots in one pdf file weill be drawn.
A barplot of the length counts (NOT a histogram, there will be no binning) and a smoth density plot.
With the `-t` option the collected data (sequence ID and sequence length) will also be written to a text file.
This scripts needs R to do the plotting.
It will write a temporary R script and pipe the data into a R process executing this script.

### removeSeqsWithN

>Usage: removeSeqsWithN.py [options] inputfile1 [inputfile2]  
>
>Options:  
>  -h, --help       show this help message and exit  
>  -z, --gzip       input file is gzipped  
>  -a, --fasta      set input and output format to fasta [default: fastq]  
>  -n X, --max-n=X  remove all sequences with more than X Ns [default: 5]  
>
>Remove sequences with more then a certain number of Ns  

Write all sequences that have less than a certain number of Ns to a different file.

### spatialQuality

>Usage: spatialQuality.py [options] input.fastq[.gz] [input2.fastq[.gz] ...]  
>
>Options:  
>  -h, --help            show this help message and exit  
>  -q, --quite           do not print status messages to the screen  
>  -o X, --output-folder=X  
>                        write results to folder X [default: ./spatialQual]  
>  -z, --gzip            input file is gzipped  
>  -n, --n-count         also make plots for N count  
>  -d, --detail-plot     make detail plot for each tile  
>  -p, --pdf             output pdf files additional to png files  
>
>Plot properties of reads according to position on the flow cell  

Quality of Illumina reads can be linked to the physical position on the flow cell.
This script plots read properties accorind to their position on the flow cell.
Default only one quality overview per file will be plotted as png file.
Additional options allow to also plot N-count per read (`-n`), produce detail plots per tile (`-d`) and to produce additiional pdf files (`-p`).

## Modules

### NcbiTools
Collection of tools to query [Ncbi](http://www.ncbi.nlm.nih.gov/). Most classes work as python dictonaries.

 | NcbiSoapMap | Base class for querying NCBI via soap | 
 | SpeciesName2TaxId | Map a scientific species name to a NCBI taxonomy ID | 
 | LineageMap | Map a NCBI taxonomy ID to the complete taxonomic lineage defined by the NCBI taxonomy tree. Return value will be a list of 3-tuples representing NCBI taxonomy nodes. The tuples will contain: rank, taxonomy ID, scientific name | 
 | NuclId2TaxIdMap | Map the GI number of a NCBI nucleotide record to the according NCBI taxonomy ID | 
 | NuclId2SpeciesNameMap | Map the GI number of a NCBI nucleotide record to sietific name of the according species | 
 | CachedNuclId2TaxIdMap | Map the GI number of a NCBI nucleotide record to the according NCBI taxonomy ID. Use a sqlite3 database as persistent cache to reduce requests to NCBI | 

### UniprotTools
Collection of tools to query [Uniprot](http://www.uniprot.org/)

 | CachedUniprotIdMap | Map protein IDs via the Uniprot mapping service. Uses a multi-layer cache (RAM and sqlite3 database). Can be configured to map different ID types to each other as long as they are supported by Uniprot. The sqlite3 cache can be initialized with a Uniprot flat file to reduce web requests. | 
 | UniprotIdMap | Map proteins IDs via the Uniprot mapping service. Will send one request per mapping to Uniprot. Can be configured to map different ID types to each other as long as they are supported by Uniprot. | 
 | UniprotInterface | Class to query the Uniprot web service. Comes with functions to read CAZy and Gene Ontology informations. | 
 
### MultiCachedDict
Module for multi-layer cached dictionaries. Multi-layer cache means that the dictionary will beside the normal python dictionary also have other cache layers, that will be queried one after another until a value for the key can be found. Cache layers should be ordered by speed.

 | MultiCachedDict | Base class that implements the multi-layer cache. | 
 | SqliteCache | Cache class to use as a cache layer in a `MultiCachedDict`. It will use a sqlite3 database to implement a local, persistent mapping. | 
 | PersistantDict | Simple cached dictionary. It uses only a `SqliteCache` (beside the normal RAM dictionary) to save mappings between program calls. | 
 
### TntBlastParser
Simple parser for Thermonuclear Blast standart output files.
 | tntBlastParser | Generator that will yield one amplicon at a time in the form of a tuple. consisting of: (NCBI GI number, start of amplicon in the data base sequence, end of amplicon in the data base sequence, amplicon sequece) | 
