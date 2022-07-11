# !/bin/bash

#This program is ment to do a Blast search of a set of files on a database.
# The user needs to specify 1) where is located that database and 2) the 
#location of the sequences that are going to be submitted to blast

db=$1 #Location of the database to do the blast
seq=$2 #Location of the sequences to blast

mkdir output-blast

ls $seq | while read line; do name=$(echo $line | cut -d'.' -f1);
for i in 0.000001 ; do mkdir -p output-blast/$i;
blastn -db $db -query $seq/$line -outfmt "6 sseqid slen qstart qend bitscore evalue sseq " -evalue $i -num_threads 12 -out output-blast/$i/$name;
done;done