#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to run the next programs.

# This program requires that you give 3 input data. 1) where this 
#SraRunTable is located, 2) where the kraken database has been saved, and 
# 3) a sufix that you want for the files to have (from the biom file)

metd=$1 #Location to the SraRunTable.txt
kdat=$2 #Location of the kraken2 database
sufx=$3 #The choosen suffix for some files
root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
# Now we will define were the reads are:
runs='reads'

# CREATING NECCESARY FOLDERS
mkdir reads
mkdir -p reads/clavi/
mkdir -p reads/cmm/
mkdir -p taxonomy/kraken
mkdir -p taxonomy/taxonomy-logs/scripts
mkdir -p taxonomy/kraken/reports
mkdir -p taxonomy/kraken/krakens
mkdir -p taxonomy/biom-files


# DOWNLOADING THE DATA

#Let's use the next piece of code to download the data
cat $metd  |  sed -n '1!p' | while read line;  do read=$(echo $line | cut -d',' -f1); fasterq-dump -S $read -p -e 8 -o $read ; done
mv *.fastq reads/

# MANAGING THE DATA

# We will change the names of the reads files. They have a sufix that makes impossible
#to be read in a loop
ls $runs | while read line ; do new=$(echo $line | sed 's/_/-/g'); mv $runs/$line $runs/$new; done

# Now, we will create a file where the information of the run labes can be located
cat $metd  | while read line; do read=$(echo $line | cut -d',' -f1); echo $read ; done > run-labels.txt
mv run-labels.txt metadata/

# TAXONOMIC ASSIGNATION WITH KRAKEN2

cat metadata/run-labels.txt | while read line; do file1=$(echo $runs/$line-1.fastq); file2=$(echo $runs/$line-2.fastq) ; echo '\n''working in run' "$line"\ 
#kraken2 --db $kdat --threads 12 --paired $file1 $file2 --output taxonomy/kraken/krakens/$line.kraken --report taxonomy/kraken/reports/$line.report \ 
echo '#!/bin/sh''\n''\n'"kraken2 --db $kdat --threads 12 --paired" "$runs/$line"'-1.fastq' "$runs/$line"'-2.fastq' "--output taxonomy/kraken/krakens/$line.kraken --report taxonomy/kraken/reports/$line.report" > taxonomy/taxonomy-logs/scripts/$line-kraken.sh; sh taxonomy/taxonomy-logs/scripts/$line-kraken.sh; done

#CREATING THE BIOM FILE

# Now we will create the biom file using kraken-biom 
kraken-biom taxonomy/kraken/reports/* --fmt json -o taxonomy/biom-files/$sufx.biom

# EXTRACTING THE CLAVIBACTER READS

# With the next piece of code, the reads clasiffied as from the genus "Clavibacter", will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Clavibacter by kraken2: 1573

cat metadata/run-labels.txt | while read line; do echo "\nExtracting Clavibacter reads from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/clavi/clavi-$line-1.fq -o2 reads/clavi/clavi-$line-2.fq -t 1573; done

# EXTRACTING THE CLAVIBACTER MICHIGANESIS-MICHIGANENSIS READS

# With the next piece of code, the reads clasiffied as from the Clavibacter michiganensis michiganensis, will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Cmm by kraken2: 33013

cat metadata/run-labels.txt | while read line; do echo "\nExtracting Cmm reads from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/cmm/clavi-$line-1.fq -o2 reads/cmm/clavi-$line-2.fq -t 33013; done




