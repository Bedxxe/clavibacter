# !/bin/bash

#This program is ment to separate the reads of a fasta file. In order to 
#work, the read header must begin with a '>' as the next example:
# >SRR13319511.872821 872821 length=101

#This little program asks for the folder where the fasta files are located

fas=$1 # folder where the fasta files are located

#CREATE A FOLDER TO CONTAIN THE READS
mkdir -p $fas/ind-reads

#EXTRACTING THE INDIVIDUAL READS
ls $fas | grep -v 'ind-reads'| while read line; do name=$(echo $line | cut -d'.' -f1);
file=$(echo $line);
grep '>' $fas/$file | while read line; do ref=$(echo $line| cut -d' ' -f1| cut -d'.' -f2);
grep "$line" $fas/$file -A 3 > ind-temp.fasta;
echo ">"$name-$ref > $name-$ref.fasta;
grep 'G' ind-temp.fasta >> $name-$ref.fasta;
rm ind-temp.fasta;
mv $name-$ref.fasta $fas/ind-reads;
done;done


