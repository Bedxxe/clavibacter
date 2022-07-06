#!/bin/sh
# This program is going to extract the Clavibacter michiganensis michiganensis reads from a set of samples

# The program ask you to give two pieces of information:
# 1) A first prefix that is the name of the author were the data has been extracted
# 2) A second prefix that is the name of the host plant (in the case of the project were this script was created)

aut=$1 #A prefix to name some of the files. In this case, the author name.
pref=$2 #The first prefix. In this case the name of the plant's host

# CREATING NEEDED FOLDERS
mkdir -p reads/cmm/

# EXTRACTING THE CLAVIBACTER MICHIGANESIS-MICHIGANENSIS READS
# With the next piece of code, the reads clasiffied as from the Clavibacter michiganensis michiganensis, will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Cmm by kraken2: 33013

cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Cmm reads in fasta format from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/cmm/cmm-$pref-$aut-$line-1.fasta -o2 reads/cmm/cmm-$pref-$aut-$line-2.fasta -t 33013 --include-children; done
