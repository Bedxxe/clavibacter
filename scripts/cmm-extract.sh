#!/bin/sh
# This program is going to extract the Clavibacter michiganensis michiganensis reads from a set of samples

mkdir -p reads/cmm/

# EXTRACTING THE CLAVIBACTER MICHIGANESIS-MICHIGANENSIS READS

# With the next piece of code, the reads clasiffied as from the Clavibacter michiganensis michiganensis, will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Cmm by kraken2: 33013

cat metadata/run-labels.txt | while read line; do echo "\nExtracting Cmm reads from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/cmm/clavi-$line-1.fq -o2 reads/cmm/clavi-$line-2.fq -t 33013; done


