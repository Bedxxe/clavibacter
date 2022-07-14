---
source: md
title: "Automatization of the pipeline"
---

# Automatization of the pipeline

<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >


<img src="/clavibacter/figures/fuji-mountains-in-clear-weather-1831.jpg" >

## Extracting the tax-id number

Around the chapters, I have been creating a pipeline that can take the SRA 
data to download the data, do the taxonomic identification with kraken2, 
create the biom files, and extract _Clavibacter_ reads. But what about 
doing the same set of steps for a different genus? 

As I explained before, what we need to extract the specific reads of a specific 
lineage I need the **tax-id** number from NCBI. For the genus _CLavibacter_ it 
was the **1573**. 

First, I will search for a file where the information for all the NCBI taxonomic 
information is allocated inside the [NCBI-taxonomy page](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/). There are a lot of options, but I will use the files that are 
inside [`new_taxdump`](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/). 
Specifically I will use the file `names.dmp` because all the OTUs are sorted 
in tax-id ascendent order:

~~~
$ head -n20 names.dmp
~~~
{: .language-bash}

~~~
1       |       all     |               |       synonym |
1       |       root    |               |       scientific name |
2       |       Bacteria        |       Bacteria <bacteria>     |       scientific name |
2       |       bacteria        |               |       blast name      |
2       |       eubacteria      |               |       genbank common name     |
2       |       Monera  |       Monera <bacteria>       |       in-part |
2       |       Procaryotae     |       Procaryotae <bacteria>  |       in-part |
2       |       Prokaryotae     |       Prokaryotae <bacteria>  |       in-part |
2       |       Prokaryota      |       Prokaryota <bacteria>   |       in-part |
2       |       prokaryote      |       prokaryote <bacteria>   |       in-part |
2       |       prokaryotes     |       prokaryotes <bacteria>  |       in-part |
6       |       Azorhizobium Dreyfus et al. 1988 emend. Lang et al. 2013        |               |    authority        |
6       |       Azorhizobium    |               |       scientific name |
7       |       Azorhizobium caulinodans Dreyfus et al. 1988    |               |       authority    |
7       |       Azorhizobium caulinodans        |               |       scientific name |
7       |       Azotirhizobium caulinodans      |               |       equivalent name |
9       |       Acyrthosiphon pisum symbiont P  |               |       includes        |
9       |       Buchnera aphidicola Munson et al. 1991  |               |       authority       |
9       |       Buchnera aphidicola     |               |       scientific name |
9       |       primary endosymbiont of Schizaphis graminum     |               |       includes     |
~~~
{: .output}

Next, I need a way to extract the tax-id number at the genus level. Since each 
character string is separated by a tab separators, I will use this feature to 
extract the genus specific line of _Nostoc_:

~~~
$ grep "$(printf '\t')Nostoc$(printf '\t')" names.dmp
~~~
{: .language-bash}

~~~
1177    |       Nostoc  |               |       scientific name |
~~~
{: .output}

That is a great output, now I want to obtain just the first numbers:

~~~
$ grep "$(printf '\t')Nostoc$(printf '\t')" names.dmp |  grep -o '[0-9]*'
~~~
{: .language-bash}

~~~
1177
~~~
{: .output}

## The all-around pipeline

Until now, the `all-around.sh` program have the next lines of code:

~~~
#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to download, trim and move the libraries information.

# This program requires that you give 4 input data: 1) where this 
#SraRunTable is located, 2) where the kraken database has been saved,
# 3) the name of the author of the library, and 4) A prefix for some of the
# output files, in this case the name of the plant's host.

#ASSIGNATIONS
metd=$1 #Location to the SraRunTable.txt
kdat=$2 #Location of the kraken2 database
aut=$3 #A first prefix to name some of the files. In this case, the author name.
pref=$4 #The second prefix. In this case the name of the plant's host

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

mkdir -p reads/clavi/
mkdir -p reads/cmm/
mkdir -p reads/fasta-clavi

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

# EXTRACTING THE CLAVIBACTER READS
# With the next piece of code, the reads clasiffied as from the genus "Clavibacter", will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Clavibacter by kraken2: 1573

#EXTRACT THE READS IN FASTQ FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Clavibacter reads in fastq from sample:" $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/clavi/$pref-$aut-$line-clav-1.fq -o2 reads/clavi/$pref-$aut-$line-clav-2.fq -t 1573 --fastq-output --include-children; done

#EXTRACT THE READS IN FASTA FORMAT
cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Clavibacter reads in fasta from sample:" $line;
extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/fasta-clavi/$pref-$aut-$line-clav-1.fasta -o2 reads/fasta-clavi/$pref-$aut-$line-clav-2.fasta -t 1573 --include-children; done

# EXTRACTING THE CLAVIBACTER MICHIGANESIS-MICHIGANENSIS READS
# With the next piece of code, the reads clasiffied as from the Clavibacter michiganensis michiganensis, will be separated from the main reads. 
# The number needed for the extraction is the numeric-ID given to Cmm by kraken2: 33013

cat metadata/run-labels.txt | while read line;
do echo "\nExtracting Cmm reads in fasta format from sample:" $line; extract_kraken_reads.py -k taxonomy/kraken/krakens/$line.kraken -r taxonomy/kraken/reports/$line.report -s1 reads/$line-1.fastq -s2 reads/$line-2.fastq -o reads/cmm/cmm-$pref-$aut-$line-1.fasta -o2 reads/cmm/cmm-$pref-$aut-$line-2.fasta -t 33013 --include-children; done
~~~
{: .language-bash}

I am going to modify it to ask for a new input: The genus of interest. I also 
will modify this to use this tax-id to extract the reads from this genus.

The new input lines will be as follows:

~~~
ndmp=$5 #Location of the "names.dmp" file
igen=$6 #Genus of interest
~~~
{: .language-bash}

And the original last 15 lines of `all-around.sh`, will change for the next piece of code:

~~~
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

~~~
{: .language-bash}

This program will generate an output file `genus-interest.txt` inside the `taxonomy/` folder.


In the end, we will have a new scrit that we will call `iall-around.sh` because of the little interactivity 
of enabling the user to choose the genus of interest:

~~~
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
~~~
{: .language-bash}

As all the other scripts, this will be located in the [scripts folder](https://github.com/Bedxxe/clavibacter/tree/main/scripts)
of the repository.



<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
