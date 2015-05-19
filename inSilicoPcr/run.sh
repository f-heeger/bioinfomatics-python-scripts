if [ -n "$1" ]; then 
    THREADS=$1
else
    THREADS=6
fi

snakemake -s inSilicoPcr_snakefile.py --forceall --dag | dot -Tpdf > pipeline.pdf
snakemake -s inSilicoPcr_snakefile.py -j ${THREADS} --resources dbaccess=1
