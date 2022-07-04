---
source: md
title: "Taxonomic assignation"
---

# Taxonomic assignation with kraken2
<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >


<img src="/clavibacter/figures/fuji-from-the-platform-of-sasayedo.jpg" >

## Preparing the database

For the taxonomic assigantion I will going to use [kraken2](https://github.com/DerrickWood/kraken2). The algorith that this software use has demostrated to be 
accurate and reliable. 

According to the kraken2 manual, it is imperative that first we assemble a databese. I downloaded the next set of databases because I want to have the capability to 
discern if the the source of the reads  that we are obtaining:
* Archea
* Bacteria
* Fungi
* Viral

The database at the end occupies close to 209 Gb of disk space and needs at least 58 
Gb of RAM memory to be able to do the taxonomic assignation. I put the information 
of the database in a location I can easily remember because I will need to write 
the directory's location in during the next steps:

~~~
$ pwd
~~~
{: .language-bash}

~~~
/mnt/d/programs/kraken/database
~~~
{: .output}

I run the next set lines of code to assemble the databese. First, I downloaded each 
set of information of the four groups I mentioned before:

~~~
$ kraken2-build --download-library archaea --db datab03242022
$ kraken2-build --download-library bacteria --db datab03242022
$ kraken2-build --download-library fungi --db datab03242022
$ kraken2-build --download-library viral --db datab03242022
~~~
{: .language-bash}

I named the database with the date it was assembled to have a record of this and 
be aware when it is convenient to update it. Secondly, I downloaded the taxonomy 
information:

~~~
$ kraken2-build --download-taxonomy --db datab03242022
~~~
{: .language-bash}

Finally, I build the database:

~~~
$ kraken2-build --build --db datab03242022
~~~
{: .language-bash}

All this process took my hardware (SSD, 64 RAM memory, 16 cores) close to 36 hours.

## Running kraken2

With the database already created, we can return to the directory where we have the 
information from [Choi-2020](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA667562). 

First, we will create a set of folders where we can save all the 
outputs that we will obtain from kraken:

~~~
$ mkdir -p taxonomy/kraken/reports
$ mkdir -p taxonomy/kraken/krakens
~~~
{: .language-bash}

We will use the next line to run the `kraken2` algorith with, using the "recently" 
created database:
~~~
$ kraken2 --db /mnt/d/programs/kraken/database/database/datab03242022 --threads 8 --paired reads/SRR12778013-1.fastq reads/SRR12778013-2.fastq --output taxonomy/kraken/krakens/SRR12778013.kraken --report taxonomy/kraken/krakens/SRR12778013.report
~~~
{: .language-bash}

~~~
reads/SRR12778013-1.fastq reads/SRR12778013-2.fastq
Loading database information... done.
60417 sequences (29.85 Mbp) processed in 3.700s (979.6 Kseq/m, 483.93 Mbp/m).
  60403 sequences classified (99.98%)
  14 sequences unclassified (0.02%)
~~~
{: .output}

Now, we have our first result from the first library:

~~~
$ head -n 15 taxonomy/kraken/reports/SRR12778013.report
~~~
{: .language-bash}

~~~
 95.30  65820224        65820224        U       0       unclassified
  4.70  3245750 10996   R       1       root
  4.50  3104957 55642   R1      131567    cellular organisms
  2.37  1638696 0       D       2759        Eukaryota
  2.37  1638693 0       D1      33154         Opisthokonta
  2.37  1638693 0       K       33208           Metazoa
  2.37  1638693 0       K1      6072              Eumetazoa
  2.37  1638693 0       K2      33213               Bilateria
  2.37  1638693 0       K3      33511                 Deuterostomia
  2.37  1638693 0       P       7711                    Chordata
  2.37  1638693 0       P1      89593                     Craniata
  2.37  1638693 0       P2      7742                        Vertebrata
  2.37  1638693 0       P3      7776                          Gnathostomata
  2.37  1638693 0       P4      117570                          Teleostomi
  2.37  1638693 0       P5      117571                            Euteleostomi
~~~
{: .output}


## Continue with the program construction

As I began on the last episode, I will update the `all-around` program with this 
new step of running the `kraken2` in all the files that we have. The program has the name `kraken-reads.sh` and can be located inside the [scripts folder](https://github.com/Bedxxe/clavibacter/blob/main/scripts/kraken-reads.sh) of this repository.
Let's see what is new inside:

~~~
$ cat kraken-reads.sh
~~~
{: .language-bash}


~~~
#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to download, trim and move the libraries information.

# This program requires that you give 2 input data: 1) where this 
#SraRunTable is located, and 2) where the kraken database has been saved

#ASSIGNATIONS
metd=$1 #Location to the SraRunTable.txt
kdat=$2 #Location of the kraken2 database

root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
# Now we will define were the reads are:
runs='reads'

# CREATING NECCESARY FOLDERS
mkdir reads
mkdir -p taxonomy/kraken
mkdir -p taxonomy/taxonomy-logs/scripts
mkdir -p taxonomy/kraken/reports
mkdir -p taxonomy/kraken/krakens

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

~~~
{: .output}

As it says in the first lines, now you will need to provide it with 1) the location of your `SraRunTable.txt` file, and 2) the location of your `kraken2` database. This program will create in your actual directory the folders where the outputs are going to be allocated and little pieces of code in the `taxonomy-logs/scripts` folder that were used to run `kraken2` on each 
pair of reads.

The 26th line is commented because I could not made `kraken2` to run inside the `while loop`. If you can find any solutions or any improvements, please let me know.
(diego.garfias@cinvestav.mx)

You can always change the number of threads to be used.

So, in my case I am going to run it with the following command:
~~~
$ sh kraken-reads.sh metadata/SraRunTable.txt /mnt/d/programs/kraken/database/database/datab03242022
~~~
{: .language-bash}

~~~
working in run SRR12778013
reads/SRR12778013-1.fastq reads/SRR12778013-2.fastq
Loading database information... done.
60417 sequences (29.85 Mbp) processed in 3.700s (979.6 Kseq/m, 483.93 Mbp/m).
  60403 sequences classified (99.98%)
  14 sequences unclassified (0.02%)

working in run SRR12778014
reads/SRR12778014-1.fastq reads/SRR12778014-2.fastq
Loading database information... done.
68119 sequences (33.72 Mbp) processed in 3.412s (1197.9 Kseq/m, 592.97 Mbp/m).
  68109 sequences classified (99.99%)
  10 sequences unclassified (0.01%)
~~~
{: .output}

The resulting output is telling us that `kraken2` is running on the files that 
we gave to the program. It will take a several minutes to a couple of hours to 
process this set of 24 samples.

By the end of the process, we will have a new set of content inside the `taxonomy/` folder. Let's see what it is inside.

~~~
$ tree -L 2
~~~
{: .language-bash}

~~~
.
├── kraken
│   ├── krakens
│   │   ├── SRR12778013.kraken
│   │   ├── SRR12778014.kraken
│   │   ├── SRR12778015.kraken
│   │   ├── SRR12778016.kraken
│   │   ├── SRR12778017.kraken
│   │   ├── SRR12778018.kraken
│   │   ├── SRR12778019.kraken
│   │   ├── SRR12778020.kraken
│   │   ├── SRR12778021.kraken
│   │   ├── SRR12778022.kraken
│   │   ├── SRR12778023.kraken
│   │   └── SRR12778024.kraken
│   └── reports
│       ├── SRR12778013.report
│       ├── SRR12778014.report
│       ├── SRR12778015.report
│       ├── SRR12778016.report
│       ├── SRR12778017.report
│       ├── SRR12778018.report
│       ├── SRR12778019.report
│       ├── SRR12778020.report
│       ├── SRR12778021.report
│       ├── SRR12778022.report
│       ├── SRR12778023.report
│       └── SRR12778024.report
└── taxonomy-logs
    └── scripts
        ├── SRR12778013-kraken.sh
        ├── SRR12778014-kraken.sh
        ├── SRR12778015-kraken.sh
        ├── SRR12778016-kraken.sh
        ├── SRR12778017-kraken.sh
        ├── SRR12778018-kraken.sh
        ├── SRR12778019-kraken.sh
        ├── SRR12778020-kraken.sh
        ├── SRR12778021-kraken.sh
        ├── SRR12778022-kraken.sh
        ├── SRR12778023-kraken.sh
        └── SRR12778024-kraken.sh

6 directories, 36 files
~~~
{: .output}

Inside the `reports/` folder are all the reports from all the samples that the 
script processed. The same goes for the `krakens/` folder

Now, we can do this same process with each of the `capsicum` folders. Inside `miscelaneous-capsicum` and `newberry-2020` folders, I will use the next line of code to obtain the desired results:

~~~
$ sh kraken-reads.sh metadata/SraRunTable.txt /mnt/d/programs/kraken/database/database/datab03242022
~~~
{: .language-bash}

If we expore what we obtained, we will see the results are 
there:

~~~
$ for i in miscelaneous-capsicum newberry-2020; do echo -e "\n"; tree $i -L 3; done
~~~
{: .language-bash}

~~~
miscelaneous-capsicum
├── kraken-reads.sh
├── metadata
│   ├── run-labels.txt
│   ├── SraRunTable.txt
│   └── SRR_Acc_List.txt
├── reads
│   ├── ERR5639101-1.fastq
│   ├── ERR5639101-2.fastq
│   ├── SRR13319509-1.fastq
│   ├── SRR13319509-2.fastq
│   ├── SRR13319510-1.fastq
│   ├── SRR13319510-2.fastq
│   ├── SRR13319511-1.fastq
│   └── SRR13319511-2.fastq
└── taxonomy
    ├── kraken
    │   ├── krakens
    │   └── reports
    └── taxonomy-logs
        └── scripts

8 directories, 12 files


newberry-2020
├── kraken-reads.sh
├── metadata
│   ├── run-labels.txt
│   ├── SraRunTable.txt
│   └── SRR_Acc_List.txt
├── reads
│   ├── SRR10527387-1.fastq
│   ├── SRR10527387-2.fastq
│   ├── SRR10527388-1.fastq
│   ├── SRR10527388-2.fastq
│   ├── SRR10527389-1.fastq
│   └── SRR10527389-2.fastq
└── taxonomy
    ├── kraken
    │   ├── krakens
    │   └── reports
    └── taxonomy-logs
        └── scripts

8 directories, 10 files
~~~
{: .output}

With this, we can progress to do the same for all the set of libraries that we want 
to process with `kraken2`




<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
