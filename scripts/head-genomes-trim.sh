#!/bin/sh
# This program is to change the headers of a fasta file in order for all the 
#headers to be the same for the entire file
# This is useful when the user want to concatenate all the genomes from a 
#database and want to know from which of the initial genomes is the 
#best match

#This program ask from the user to specify 1 input: 1) the location 
#of the fasta files without the "/" cahracter at the end

gem=$1 #Location to the genomes in fasta format without the "/" character at the end
# CREATE A DIRECTORY TO STORE THE TRIMMED GENOME
mkdir -p $gem/trim-header-genomes

ls $gem | grep -v 'trim-header-genomes' | while read line;
do fas=$(echo $line);
#CREATE A .TXT FILE TO STORE THE ORIGINAL HEADERS
grep '>' $gem/$fas > heads.txt;
cp $gem/$fas $gem/trim-header-genomes/htrim-$fas;
#REPLACE ALL THE HEADERS FOR THE SUFFIX NAME 
cat heads.txt | while read line;
do sed "s/$line/>$fas/" $gem/trim-header-genomes/htrim-$fas > $gem/trim-header-genomes/temp.txt;
mv $gem/trim-header-genomes/temp.txt $gem/trim-header-genomes/htrim-$fas;
done;done

#REMOVING TEMP FILES
rm heads.txt


