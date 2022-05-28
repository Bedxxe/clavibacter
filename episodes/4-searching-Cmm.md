# Clavibacter

<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >

## Searching for _Clavibacter michiganensis_ michiganensis (Cmm) specific genes

<img src="/clavibacter/figures/shimomeguro.jpg" >

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

I also used a reference genome of [Cmm](https://github.com/Bedxxe/clavibacter/blob/main/data/reference/cmm/cmm-contigs.fa)

With this amassed database, I faced a problem. Each one of them has its own 
header/label nomenclature for each read/contig/scaffold. For example, the 
reference Cmm genome has this headers:


~~~
$ grep '>' cmm-contigs.fa
~~~
{: .bash}

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
{: .bash}

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
{: .output}

I will run it in a new directory called `trim-c-genomes`:

~~~
$ mkdir trim-c-genomes
$ sh changing-headers.sh ../cmm/cmm-contigs.fa Cmm
$ grep '>' Cmm.fna
~~~
{: .bash}

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
{: .bash}

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
{: .bash}

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

I will take the _Capsicum_ data to do a first try. As I explored the data, I 
became aware that in the _Capsicum_ samples, one of the dominant OTUs was 
[_Agrobacterium tumefaciens_](https://github.com/Bedxxe/clavibacter/blob/main/data/outgrups/Agrobacterium_tumefaciens_LN-1.fna) :

<img src="/clavibacter/figures/shimomeguro.jpg" >
<em> Figure 1. Barplot of dominant OTUs from _Capsicum_ and _Tuberosum_ libraries <em/>


I will add also an outgroup-like genome from [_Nostoc desertorum_](https://github.com/Bedxxe/clavibacter/blob/main/data/outgrups/Nostoc_desertorum_CM1_VF14.fna).

With these genomes I will create a Blast database. First I will concatenate all 
this genomes in one big file and then I will run the code to do the database:

~~~
$  cat trim-c-genomes/* outgroups/* > all-genomes.fasta
$  makeblastdb -in ../reference/all-genomes.fasta -dbtype nucl -out all-clavi-g/all-genomes
~~~
{: .bash}


~~~

~~~
{: .output}

~~~

~~~
{: .bash}

~~~

~~~
{: .output}

~~~
$ makeblastdb -in ../reference/cmm/cmm.fna -dbtype nucl -out cmm/ndatabase/cmm
~~~
{: .bash}


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
{: .bash}

~~~
$ bwa mem index/cmm-contigs.fa ../../../blast/cmm/shotg-seq/clavi-SRR13319511-1.fq ../../../blast/cmm/shotg-seq/clavi-SRR13319511-2.fq -o clavi-SRR13319511.sam
~~~
{: .bash}

To delete zero-blast results
~~~
ls -l output-blast/0.000001/ |  grep " 0 May 27" | while read line; do file=$(echo $line | cut -d' ' -f9); rm output-blast/0.000001/$file; done
~~~
{: .bash}

~~~

~~~
{: .bash}
