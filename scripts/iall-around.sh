#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to download, trim and move the libraries information. 
# It gives the user the option to extract the reads from a specific genus 
# of bacterial lineages. 

# REQUIREMENTS:
#a) A kraken2 database in a known location
#b) The "names.dmp" file from NCBI: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/
#c) The SraRunTable of the experiments of insterest

# This program requires that you give 6 input data: 
#1) where the SraRunTable is located 
#2) where the kraken database has been assembled
#3) the name of the author of the library
#4) A prefix for some of the output files, in this case the name of the plant's host.
#5) The localization of the "names.dmp" file
#6) The genus name of interest. This is for the extraction of this genus reads from the different libraries to analize.

#ASSIGNATIONS
metd=$1 #Location to the SraRunTable.txt
kdat=$2 #Location of the kraken2 database
aut=$3 #A first prefix to name some of the files. In this case, the author name.
pref=$4 #The second prefix. In this case the name of the plant's host
ndmp=$5 #Location of the "names.dmp" file
igen=$6 #Genus of interest

root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
# Now we will define were the reads are:
runs='reads'

# CREATING NECCESARY FOLDERS
mkdir reads
mkdir -p taxonomy/kraken
mkdir -p taxonomy/taxonomy-logs/scripts
mkdir -p taxonomy/kraken/reports
mkdir -p taxonomy/kraken/krakens
mkdir -p taxonomy/biom-files

mkdir -p reads/$igen
mkdir -p reads/fasta-$igen

# DOWNLOADING THE DATA

#Let's use the next piece of code to download the data
cat $metd  |  sed -n '1!p' | while read line;  do read=$(echo $line | cut -d',' -f1); fasterq-dump -S $read -p -e 8 -o $read ; done
mv *.fastq reads/
# The -e flag can be customized. This indicates the number of threads used to do this task.

# MANAGING THE DOWNLADED DATA

# We will change the names of the reads files. They have a sufix that makes impossible
#to be read in a loop
ls $runs | while read line ; do new=$(echo $line | sed 's/_/-/g'); mv $runs/$line $runs/$new; done

# Now, we will create a file where the information of the run labes can be located
cat $metd  | sed -n '1!p' | while read line; do read=$(echo $line | cut -d',' -f1); echo $read ; done > run-labels.txt
mv run-labels.txt metadata/

# TAXONOMIC ASSIGNATION WITH KRAKEN2

cat metadata/run-labels.txt | while read line; do file1=$(echo $runs/$line-1.fastq); file2=$(echo $runs/$line-2.fastq) ; echo '\n''working in run' "$line"\ 
#kraken2 --db $kdat --threads 6 --paired $file1 $file2 --output taxonomy/kraken/krakens/$line.kraken --report taxonomy/kraken/reports/$line.report \ 
echo '#!/bin/sh''\n''\n'"kraken2 --db $kdat --threads 6 --paired" "$runs/$line"'-1.fastq' "$runs/$line"'-2.fastq' "--output taxonomy/kraken/krakens/$line.kraken --report taxonomy/kraken/reports/$line.report" > taxonomy/taxonomy-logs/scripts/$line-kraken.sh; sh taxonomy/taxonomy-logs/scripts/$line-kraken.sh; done

#CREATING THE BIOM FILE

# Now we will create the biom file using kraken-biom 
kraken-biom taxonomy/kraken/reports/* --fmt json -o taxonomy/biom-files/$aut.biom

# OBTAINING THE TAX-ID NUMBER FROM THE "names.dmp" FILE
taxid=$(grep "$(printf '\t')$igen$(printf '\t')" $ndmp | grep -o '[0-9]*')

# Putting this information in a file:
grep "$(printf '\t')$igen$(printf '\t')" $ndmp | while read line; 
do tid=$(echo $line | cut -d' ' -f1); oname=$(echo $line | cut -d' ' -f3);
echo "Genus""\t""tax-id" > taxonomy/genus-interest.txt;
echo $oname"\t"$tid >> taxonomy/genus-interest.txt; 
done

# EXTRACTING THE READS FROM THE GENUS OF INTEREST
# With the next piece of code, the reads clasiffied as the Genus of interest, will be separated from the main reads. 

#EXTRACT THE READS IN FASTQ FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\n"Extracting $igen reads in fastq from sample: $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/$igen/$pref-$aut-$igen-$line-1.fq -o2 reads/$igen/$pref-$aut-$igen-$line-2.fq -t $taxid --fastq-output --include-children; done

#EXTRACT THE READS IN FASTA FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\n"Extracting $igen reads in fasta from sample: $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/fasta-$igen/$pref-$aut-$igen-$line-1.fasta -o2 reads/fasta-$igen/$pref-$aut-$igen-$line-2.fasta -t $taxid --include-children; done

