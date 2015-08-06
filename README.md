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

Normally write human readable output. With the -t option a tabular format of the same information will be written, that can be used in statistic evaluation of multiple assemblies.

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
Only write sequences with an ID from the list given. The list can either be given as a comma separated list of IDs or as in a file with one ID per line. With the `-e` option the given "IDs" will be used as regular expression instead of exact matching (using python regex).
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
>  -f FORMAT, --img-format=FORMAT
>                        set plot format to FORMAT [default: pdf]  


Plot the length distribution of the sequences in the input file. 
Two plots in one pdf file will be drawn.
A bar plot of the length counts (NOT a histogram, there will be no binning) and a smoothed density plot.
With the `-t` option the collected data (sequence ID and sequence length) will also be written to a text file.
Specify the image format with `-f`. Possible formats are: 'pdf', 'png', 'jpeg', 'bmp', 'postscript'.
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
This script plots read properties according to their position on the flow cell.
Default only one quality overview per file will be plotted as png file.
Additional options allow to also plot N-count per read (`-n`), produce detail plots per tile (`-d`) and to produce additional pdf files (`-p`).

### subtractReadsByMapping

>Usage: subtractReadsByMapping.py [options] mappingfile read1.fastx out1.fastx read2.fastx >out2.fastx  
>  
>Options:  
>  -h, --help           show this help message and exit  
>  -q, --quite          do not print status messages to the screen  
>  -a, --fasta          input file(s) is/are fasta  
>  -m X, --mapped1=X    write mapped reads from read 1 to this files  
>  -n X, --mapped2=X    write mapped reads from read 2 to this files  
>  -z, --gzip           input file(s) is/are gzipped  
>  -y, --mapping-gzip   mapping file is gzipped  
>  -b, --blast          mapping file is tabular blast output instead of sam  
>                       file  
>  -t X, --threshold=X  consider reads with an e-value lower than X as  
>                       "mapped". (Can only be used in blast mode) [default:  
>                       0.000001]  

From a set of single end or paired end reads in a fasta or fastq file (or two for paired end) remove all reads that were mapped in a mapping result given as a sam file or a blast tabular file (`-b`). Reads files can be gziped (`-z`) as well as the mapping file (`-y`). For blast results a minimal e-value can be given for a match to be considered as a "mapping" (`-t`).

## Modules

### NcbiTools
Collection of tools to query [Ncbi](http://www.ncbi.nlm.nih.gov/). Most classes work as python dictionaries.

#### NcbiMap 
Base class for querying NCBI via BioPyhton and the NCBI web interface. Do not use directly!
#### SpeciesName2TaxId 
Map a scientific species name to a NCBI taxonomy ID
#### TaxonomyNodeName2TaxId
Map a Ncbi Taxonomy Node name to its ID. This is slitely different from the `SpeciesName2TaxId` map: It always returns a list, which contains multiple IDs for ambiguous node names and it also works for higher taxonomic levels and nodes that do not have a standard rank (like sub-phylum or no rank).
#### LineageMap
Map a NCBI taxonomy ID to the complete taxonomic lineage defined by the NCBI taxonomy tree. Return value will be a list of 3-tuples representing NCBI taxonomy nodes. The tuples will contain: rank, taxonomy ID, scientific name
#### SingleLevelLineageMap
Map a NCBI taxonomy ID to specific taxonomic level from the NCBI taxonomy tree. Return value will be a list of 3-tuples representing NCBI taxonomy nodes. The tuples will contain: rank, taxonomy ID, scientific name
#### TaxonomyParentMap
Map a NCBI taxonomy ID to its parent node (in the NCBI taxonomy tree) NCBI taxonomy ID.
#### NuclId2TaxIdMap
Map the GI number of a NCBI nucleotide record to the according NCBI taxonomy ID 
#### NuclId2SpeciesNameMap
Map the GI number of a NCBI nucleotide record to scientific name of the according species
#### CachedNuclId2TaxIdMap
Map the GI number of a NCBI nucleotide record to the according NCBI taxonomy ID. Use a sqlite3 database as persistent cache to reduce requests to NCBI 
#### NcbiTaxonomyTree
A class representing the tree given by the NCBI taxonomy database. If a cache path is given to the constructor a database at this path will be used as persistent cache. The object can be initialized by the function with the same name. It takes a node.dmp file of the NCBI taxonomy file dump as input. This option is only available if a cache is used (otherwise there is no place to store the initialized data). Missing data will be loaded directly from the NCBI data base via SOAP request, but only once it is needed.  
The object can be queried for information on the tree with NCBI taxonomy IDs representing a node. The function include: parent of a node, full path to the root, lowest common ancestor of two or more nodes and a variant of the lowest common ancestor that for a set of nodes returns the lowest node for which all input nodes are either ancestors or descendants of the output node (called lowest common node (LCN) here).

### UniprotTools
Collection of tools to query [Uniprot](http://www.uniprot.org/)

#### CachedUniprotIdMap
Map protein IDs via the Uniprot mapping service. Uses a multi-layer cache (RAM and sqlite3 database). Can be configured to map different ID types to each other as long as they are supported by Uniprot. The sqlite3 cache can be initialized with a Uniprot flat file to reduce web requests.
#### UniprotIdMap 
Map proteins IDs via the Uniprot mapping service. Will send one request per mapping to Uniprot. Can be configured to map different ID types to each other as long as they are supported by Uniprot. 
#### UniprotInterface
Class to query the Uniprot web service. Comes with functions to read CAZy and Gene Ontology informations.
 
### MultiCachedDict
Module for multi-layer cached dictionaries. Multi-layer cache means that the dictionary will beside the normal python dictionary also have other cache layers, that will be queried one after another until a value for the key can be found. Cache layers should be ordered by speed.

#### MultiCachedDict
Base class that implements the multi-layer cache.
#### SqliteCache
Cache class to use as a cache layer in a `MultiCachedDict`. It will use a sqlite3 database to implement a local, persistent mapping.
#### SqliteListCache
Same as `SqliteCache`, but can deal with values that are lists.
#### PersistantDict
Simple cached dictionary. It uses only a `SqliteCache` (beside the normal RAM dictionary) to save mappings between program calls.
 
### TntBlastParser
Simple parser for Thermonuclear Blast standard output files.
#### tntBlastParser 
Generator that will yield one amplicon at a time in the form of a tuple. consisting of: (NCBI GI number, start of amplicon in the data base sequence, end of amplicon in the data base sequence, amplicon sequece)

### EolTools
Colelction of tools to query the [Encyclopedia of Life](eol.org)

#### EolName2IdMap
Dictionary that maps names to EOL IDs using the `search` EOL web service.
Will return a list of IDs.
Search parameters can be configured in the `config` dictionary.
Note: exact search is switched on by default.

#### EolInterface
Loads data for an EOL entry once with the `page` EOL web service.
Data is kepped in memory and can be querried directly in the `data` member of this class or by specialized functions.
