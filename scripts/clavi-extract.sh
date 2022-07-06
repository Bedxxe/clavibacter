#!/bin/sh
# This program is going to extract the Clavibacter reads from a set of samples.
# You need to specify the program 1)the name of the author or a first prefix to name 
# the output files, and 2) a second prefix that is advised to be the name of the host plant.

aut=$1 #A prefix to name some of the files. In this case, the author name.
pref=$2 #The first prefix. In this case the name of the plant's host

# EXTRACTING THE CLAVIBACTER READS

# With the next piece of code, the reads clasiffied as from the genus "Clavibacter", will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Clavibacter by kraken2: 1573

# MAKING THE NECESSARY FOLDERS TO ALLOCATE THE OUTPUTS
mkdir -p reads/clavi/
mkdir -p reads/fasta-clavi

#EXTRACT THE READS IN FASTQ FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Clavibacter reads in fastq from sample:" $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/clavi/$pref-$aut-$line-clav-1.fq -o2 reads/clavi/$pref-$aut-$line-clav-2.fq -t 1573 --fastq-output --include-children; done

#EXTRACT THE READS IN FASTA FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Clavibacter reads in fasta from sample:" $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/fasta-clavi/$pref-$aut-$line-clav-1.fasta -o2 reads/fasta-clavi/$pref-$aut-$line-clav-2.fasta -t 1573 --include-children; done
