---
source: md
title: "Searching for Cmm"
---


# Searching for _Clavibacter michiganensis_ michiganensis (Cmm) specific sequences

<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >


<img src="/clavibacter/figures/the-fields-of-sekiya-by-the-sumida-river.jpg" >

After extracting the specific _Clavibacter michiganensis_ michiganensis (Cmm) reads 
from all the libraries, I need to evaluate if these sequences are truly from 
Cmm. The reads separation was made by the `kraken2` algorithm using its database. 
The main issue is if this sequences are really from Cmm and not from one of the 
other subspecies fo _Clavibacter michiganensis_. 

I downloaded a set of genomes of the michiganensis subspecies from the 
NCBI:

|-------------------+-------------------------------------------------------------------------------------------------------------------|
| _Cm_ subspecies | NCBI-link | Number of genomes |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| capsici | [link](https://ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=64934) | 4 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| tessellarius | [link](https://ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=64935) | 3 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| insidiosus | [link](https://ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=64932) | 7 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| nebraskensis | [link](https://ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=64933) | 12 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| sepedonicus | [link](https://ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=64931) | 3 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|

The reference genomes are located inside the `data/reference/clavi-genomes/` 
[folder of the repository](https://github.com/Bedxxe/clavibacter/tree/main/data/reference/clavi-genomes).


I also used a reference genome of [Cmm](https://github.com/Bedxxe/clavibacter/blob/main/data/reference/cmm/cmm-contigs.fa)

With this amassed database, I faced a problem. Each one of them has its own 
header/label nomenclature for each read/contig/scaffold. For example, the 
reference Cmm genome has this headers:


~~~
$ grep '>' cmm-contigs.fa
~~~
{: .language-bash}

~~~
>NC_009480.1
>NC_009478.1
>NC_009479.1
~~~
{: .output}

To change this, I created a little program called [changing-headers.sh](https://github.com/Bedxxe/clavibacter/blob/main/scripts/changing-headers.sh) 
which takes the headers of a fasta file and change them for a desired one. 
It is important to highlight that all the headers are going to be the same.

~~~
$ cat changing-headers.sh
~~~
{: .language-bash}

~~~
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
~~~
{: .language-bash}

I will run it in a new directory called `trim-c-genomes`:

~~~
$ mkdir trim-c-genomes
$ sh changing-headers.sh ../cmm/cmm-contigs.fa Cmm
$ grep '>' Cmm.fna
~~~
{: .language-bash}

~~~
>Cmm
>Cmm
>Cmm
~~~
{: .output}

With this result, I will apply this to all the _Clavibacter_ genomes that I 
downloaded previously:

~~~
$ ls clavi-genomes/
~~~
{: .language-bash}

~~~
Clavibacter_michiganensis_subsp_capsici_1101.fna
Clavibacter_michiganensis_subsp_capsici_1106.fna
Clavibacter_michiganensis_subsp_capsici_1207.fna
Clavibacter_michiganensis_subsp_capsici_PF008.fna
Clavibacter_michiganensis_subsp_insidiosus_ATCC_10253.fna
Clavibacter_michiganensis_subsp_insidiosus_CFBP_1195.fna
Clavibacter_michiganensis_subsp_insidiosus_CFBP_2404.fna
Clavibacter_michiganensis_subsp_insidiosus_CFBP_6488.fna
Clavibacter_michiganensis_subsp_insidiosus_LMG_3663.fna
Clavibacter_michiganensis_subsp_insidiosus_R1-1.fna
Clavibacter_michiganensis_subsp_insidiosus_R1-3.fna
Clavibacter_michiganensis_subsp_nebraskensis_419B.fna
Clavibacter_michiganensis_subsp_nebraskensis_44.fna
Clavibacter_michiganensis_subsp_nebraskensis_61-1.fna
Clavibacter_michiganensis_subsp_nebraskensis_7580.fna
Clavibacter_michiganensis_subsp_nebraskensis_A6096.fna
Clavibacter_michiganensis_subsp_nebraskensis_CFBP_7577.fna
Clavibacter_michiganensis_subsp_nebraskensis_CIBA.fna
Clavibacter_michiganensis_subsp_nebraskensis_DOAB_395.fna
Clavibacter_michiganensis_subsp_nebraskensis_DOAB_397.fna
Clavibacter_michiganensis_subsp_nebraskensis_HF4.fna
Clavibacter_michiganensis_subsp_nebraskensis_NCPPB_2581.fna
Clavibacter_michiganensis_subsp_nebraskensis_SL1.fna
Clavibacter_michiganensis_subsp_sepedonicus_ATCC33113.fna
Clavibacter_michiganensis_subsp_sepedonicus_CFIA-Cs3N.fna
Clavibacter_michiganensis_subsp_sepedonicus_CFIA-CsR14.fna
Clavibacter_michiganensis_subsp_tessellarius_ATCC_33566.fna
Clavibacter_michiganensis_subsp_tessellarius_DOAB_609.fna
~~~
{: .output}

The line of code that I use to change all the headers from all the files, 
is the next one:
~~~
$ cd trim-c-genomes
$ ls clavi-genomes/*.fna | while read line; do name=$(echo $line | cut -d'/' -f2 |cut -d'_' -f4,5,6| cut -d'.' -f1); sh changing-headers.sh $line $name; done
$ ls
~~~
{: .language-bash}

~~~
capsici_1101.fna           insidiosus_R1-1.fna         nebraskensis_DOAB_397.fna
capsici_1106.fna           insidiosus_R1-3.fna         nebraskensis_HF4.fna
capsici_1207.fna           nebraskensis_419B.fna       nebraskensis_NCPPB_2581.fna
capsici_PF008.fna          nebraskensis_44.fna         nebraskensis_SL1.fna
Cmm.fna			           nebraskensis_61-1.fna       sepedonicus_ATCC33113.fna
insidiosus_ATCC_10253.fna  nebraskensis_7580.fna       sepedonicus_CFIA-Cs3N.fna
insidiosus_CFBP_1195.fna   nebraskensis_A6096.fna      sepedonicus_CFIA-CsR14.fna
insidiosus_CFBP_2404.fna   nebraskensis_CFBP_7577.fna  tessellarius_ATCC_33566.fna
insidiosus_CFBP_6488.fna   nebraskensis_CIBA.fna       tessellarius_DOAB_609.fna
insidiosus_LMG_3663.fna    nebraskensis_DOAB_395.fna
~~~
{: .output}


As I explored the data, I became aware that in the _Capsicum_ and _Tuberosum_ 
samples, one of the dominant OTUs was 
[_Agrobacterium tumefaciens_](https://github.com/Bedxxe/clavibacter/blob/main/data/outgrups/Agrobacterium_tumefaciens_LN-1.fna) :

<img src="/clavibacter/figures/05-01-dominantGenus.png" >

In order to see if this obtained reads do not belong to this dominant population, 
I add the _Agrobacterium tumefaciens_ to the database.

I will also add an outgroup-like genome from [_Nostoc desertorum_](https://github.com/Bedxxe/clavibacter/blob/main/data/outgrups/Nostoc_desertorum_CM1_VF14.fna).

With these genomes I will create a Blast database. First I will concatenate all 
this genomes in one big file and then I will run the next lines of code to do 
the BLAST database:

~~~
$  cat trim-c-genomes/* outgroups/* > all-genomes.fasta
$  makeblastdb -in ../reference/all-genomes.fasta -dbtype nucl -out all-clavi-g/all-genomes
~~~
{: .language-bash}

Now that we have all we need to try a Blast of the "Cmm" extracted-reads, I will 
take the `Capsicum/choi-2020` data to try this. I will create a directory at the same level 
as `all-clavi-g/` directory, and I will name it `capsi-choi/`. Here I will copy 
the extracted Cmm from `choi-2020/` inside a folder named `extracted-cmm-choi/`. If we explore one of these files, we will see that is composed for a set of reads:

~~~
$ head -n15 extracted-cmm-choi/cmm-capsi-choi-2020-SRR12778013-1.fasta
~~~
{: .language-bash}

~~~
>SRR12778013.32167716 32167716 length=151
TCGTCATCGACTGCCCCAACCGGCAGGGCGGCCCCCTCACCCTCTCCGCCCTCAACGCCG
CCGACACCGTCGTCTACGCCGCCACCCCCAGCGGCGACGGCGTCGACGGCGTCGCCGGCG
CCCGCCGCACCGTCGACCAGTTCCGGATCAA
>SRR12778013.32168284 32168284 length=151
TCGTCATCGACTGCCCCAACCGGCAGGGCGGCCCCCTCACCCTCTCCGCCCTCAACGCCG
CCGACACCGTCGTCTACGCCGCCACCCCCAGCGGCGACGGCGTCGACGGCGTCGCCGGCG
CCCGCCGCACCGTCGACCAGTTCCGGATCAA
>SRR12778013.32987997 32987997 length=151
TTGTGATCGGAGACCTCGTCGCGCAGCCGGGCGGCCAACTCGAACTTGAGTTCGCCGGCG
GCGATGAGCATCTGGGCGTTGAGATCCTCGATGATCGCCTCGAGCTCGGCACCGCCCGAG
GCGCCCTTCGCGCCCGAGCCCAGCGACGGGG
>SRR12778013.34496501 34496501 length=151
ACGTCGTCGTGTACTGCTGGGGCCCCGGCTGCAACGGCAGCACCCGCGCCGGCCTCGCCC
TCGCGCTCGCGGGCTACGGCCGCGTGAAGGAGCTGGTCGGCGGCTACGAGTACTGGGTCC
~~~
{: .output}

In order to separe each sequence in a file, I have prepared a little program 
`dividing-ind-reads.sh`

~~~
$ cat dividing-ind-reads.sh
~~~
{: .language-bash}

~~~
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
~~~
{: .language-bash}

I will use this program to divide the sequences of each of the `fasta` files inside 
the `extracted-cmm-choi/` direcotry:

~~~
$ sh dividing-ind-reads.sh extracted-cmm-choi
$ ls extracted-cmm-choi/ind-reads/ | wc -l
~~~
{: .language-bash}

~~~
7392
~~~
{: .output}

Now we are ready to do the BLAST with the mock microbiome that we amassed. I 
prepared another little program to do the blast, `sec-blast.sh`:

~~~
$ cat sec-blast.sh
~~~
{: .language-bash}

~~~
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
~~~
{: .language-bash}

It is important to note that one of the advantages of the `-outfmt 6` of BLAST is 
that one can [customize the output](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) . We are asking for 6 parameters:

* sseqid:	Subject Seq-id
* slen:		Subject sequence length		
* qstart:	Start of alignment in query
* qend:		End of alignment in query
* bitscore:	Bit score
* evalue:	Expect value
* sseq: 	Aligned part of subject sequence

~~~
$ sh sec-blast.sh ../all-clavi-g/all-genomes extracted-cmm-choi/ind-reads/
~~~
{: .language-bash}

It will take a long time since we got a great amount of sequences

~~~

~~~
{: .language-bash}

~~~

~~~
{: .language-bash}

~~~

~~~
{: .language-bash}

~~~
$ makeblastdb -in ../reference/cmm/cmm.fna -dbtype nucl -out cmm/ndatabase/cmm
~~~
{: .language-bash}


~~~


Building a new DB, current time: 05/23/2022 16:15:51
New DB name:   /home/betterlab/diego/clavibacter/blast/cmm/ndatabase/cmm
New DB title:  ../reference/cmm/cmm.fna
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 3053 sequences in 0.0727699 seconds.
~~~
{: .output}

~~~
cat
cat headers.txt | while read line; do sed "s/$line/>Cmm/" cmm-contigs.fa > temp.txt; mv temp.txt cmm-contigs.fa ;done
~~~
{: .language-bash}

~~~
$ bwa mem index/cmm-contigs.fa ../../../blast/cmm/shotg-seq/clavi-SRR13319511-1.fq ../../../blast/cmm/shotg-seq/clavi-SRR13319511-2.fq -o clavi-SRR13319511.sam
~~~
{: .language-bash}

To delete zero-blast results
~~~
ls -l output-blast/0.000001/ |  grep " 0 May 27" | while read line; do file=$(echo $line | cut -d' ' -f9); rm output-blast/0.000001/$file; done
~~~
{: .language-bash}

~~~

~~~
{: .language-bash}


~~~

~~~
{: .language-bash}

~~~

~~~
{: .output}

~~~

~~~
{: .language-bash}

~~~

~~~
{: .output}

~~~

~~~
{: .language-bash}

~~~

~~~
{: .output}



<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
