#!/bin/sh
# This program is to change the headers of a fasta file in order for all the 
#headers to be the same for the entire file
# This is useful when the user want to concatenate all the genomes from a 
#database and want to know from which of the initial genomes is the 
#best match

#This program ask from the user to specify 2 things: 1) the location 
#of the fasta file and 2) the suffix that all the identifiers are going 
#to be changed for.

gem=$1 #Location to the genome in fasta format
sufx=$2 #The new header without the '>' symbol at the begginng

#CREATE A .TXT FILE TO STORE THE ORIGINAL HEADERS
grep '>' $gem > heads.txt

cp $gem $sufx.fna

#REPLACE ALL THE HEADERS FOR THE SUFFIX NAME 
cat heads.txt | while read line; do sed "s/$line/>$sufx/" $sufx.fna > temp.txt; mv temp.txt $sufx.fna ;done

#REMOVING TEMP FILES
rm heads.txt


