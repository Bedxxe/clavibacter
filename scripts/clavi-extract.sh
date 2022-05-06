#!/bin/sh
# This program is going to extract the Clavibacter reads from a set of samples


# EXTRACTING THE CLAVIBACTER READS

# With the next piece of code, the reads clasiffied as from the genus "Clavibacter", will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Clavibacter by kraken2: 1573

mkdir -p reads/clavi/

cat metadata/run-labels.txt | while read line; do echo "\nExtracting Clavibacter reads from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/clavi/$line-clav-1.fq -o2 reads/clavi/$line-clav-2.fq -t 1573; done

