---
source: md
title: "First exploration of the data "
---

# First exploration of the data
<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >


<img src="/clavibacter/figures/bw-hokusai.jpg" alt="Black and white representation of a wave made by Hokusai." >

## Introductions

I was invited to participate in a project to uncover the presence of a group of bacterial lineajes from the genus 
_Clavibacter_ in different cultivars. As planned by my Sensei (Professor Nelly Selem-Mojica), we will use 
both approximation of metagenomics (_i.e._ Metabarcoding and Shotgun) to try to obtain some answers 
from the questions that were formulated because of the widely presence of these lineages.

## The data

I want to start this project with the exploration of some of the shotgun data. The 
next two lists of data was given by the team:

### The first one with a literature exploration:

|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| Plant                                               |  Number of samples     | Type of data  | Link to the data | Paper| 
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| _Solanum lycopersicum_| 13 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=603603)| Barajas et al., 2020|
|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| _Zea mays_| 5| Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA645385) |Xiong et al., 2021|
|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| _Zea mays_| 5| Shotgun | [Link](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA645371) |Akinola et al., 2021 |
|-------------------+-------------------------------------------------------------------------------------------------------------------| 
| _Solanum tuberosum_ | 20 | Shotgun| [Link](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA477767) |Shi et al., 2019 |
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| _Solanum tuberosum_ | 39 | 16S |[Link](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA477767) | Shi et al., 2019 |
|-------------------+-----------------------------------------------------------------------------------------------------------------| 
| _Capsicum annuum_ | 12 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA667562) | Choi et al., 2020|
|-------------------+-----------------------------------------------------------------------------------------------------------------| 
| _Solanum lycopersicum_ in _Capsicum annuum_ soils | 15| Shotgun| [Link](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA590717) | Newberry et al., 2020 | 
|-------------------+-----------------------------------------------------------------------------------------------------------------| 
| Soil of _Zea mays_ and _Triticum_ rotation |27 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA640885) | X. Wu et al., 2021| 
|-------------------+-----------------------------------------------------------------------------------------------------------------| 
| _Triticum_ | - | - | No data provided by authors of the paper Quiza et al., 2021| Quiza et al., 2021|
|-------------------+-----------------------------------------------------------------------------------------------------------------| 
| _Zea mays_ grown in _Triticum_ soil | 42 | 16S | [Link](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA745034) | M. Wu et al., 2021|
|-------------------+-----------------------------------------------------------------------------------------------------------------| 

### A second set of data obtained by a general exploration on the NCBI database:

|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| Plant                                               |  Number of samples     | Type of data  | Link to the data | 
|-------------------+-------------------------------------------------------------------------------------------------------------------|
| _Medicago sativa_| 18 | 16S  | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=Medicago+sativa+metagenomic+NOT+(amplicon)+NOT+(RNA)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| _Zea mays_ | 41 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=metagenomics+zea+mays+NOT+(amplicon)+NOT+(RNA)+NOT+(PCR)+NOT+(454)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 
| _Capsicum annuum_ | 4 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=capsicum+annuum+metagenome+NOT+(amplicon)+NOT+(RNA)+NOT+(454)+NOT+(PCR)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 
| _Titricum_ | 54 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=triticum+metagenome+NOT+(amplicon)+NOT+(RNA)+NOT+(454)+NOT+(PCR)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 

## First exploration of *Capsicum* dataset from Choi-2020

I will use the data from [Choi-2020](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA667562) to begin with the analysis. 
I have created a folder structure inside `clavibacter` main folder as follows:

~~~
$ tree
~~~
{: .bash}

~~~
.
├── 16S
└── shotgun
    ├── capsicum
    ├── lycopersicum
    ├── triticum
    ├── tuberosum
    └── zea-mays
~~~
{: .output}

Inside the `capsicum` folder, I have created the next folder organization:

~~~
$ tree -L 1
~~~
{: .bash}

~~~
.
├── choi-2020
├── metadata
├── miscelaneous-capsicum
└── newberry-2020

4 directories, 0 files
~~~
{: .output}

I will enter to the `choi-2020` folder and create a `metadata/` folder and 
allocated all the metadata information right there to give structure to out 
working directory.

Let's first downloaded them from the [SRA](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_627021272c1b2e5f9507d1a7&o=acc_s%3Aa) repository in NCBI and move it into the `metadata/` folder. 
I will download the metadata table and the Accesion table. With the Accseion table I will use SRA-toolkit to download the shotgun reads.

This files are located inside the its [portion](https://github.com/Bedxxe/clavibacter/tree/main/data/shotgun/capsicum/choi-2020/metadata) of the `data` folder.

I will use the next command to obtain the reads:

~~~
$ cat metadata/SRR_Acc_List.txt | while read line; do fasterq-dump $line -S -p -e 12 -o $line; done
$ mkdir reads/
$ mv *.fastq reads/
$ ls reads/*.fastq | wc -l
~~~
{: .bash}

~~~
24
~~~
{: .output}

Now, we have both the forward and reverse reads to begin to work with them.
Next, I want to explore the diversity inside each of these metagenomes. I will use [kraken2](https://github.com/DerrickWood/kraken2) to this purpose. I want to use the information 
inside the metadata table `SraRunTable.txt` to run the kraken commands, and also 
to use it to run all the other programs that I will be using along this analysis.
Let's see the structure of the file:
~~~
$ head -n 5 metadata/SraRunTable.txt
~~~
{: .bash}

~~~
Run,Assay Type,AvgSpotLen,Bases,BioProject,BioSample,BioSampleModel,Bytes,Center Name,Collection_date,Consent,DATASTORE filetype,DATASTORE provider,DATASTORE region,Experiment,geo_loc_name_country,geo_loc_name_country_continent,geo_loc_name,HOST,Instrument,Isolation_Source,Lat_Lon,Library Name,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,ReleaseDate,replicate,Sample Name,SRA Study
SRR12778013,OTHER,302,20857924148,PRJNA667562,SAMN16378122,Metagenome or environmental,6662690206,CHINESE ACADEMY OF SCIENCES,2018-08-23,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",SRX9247618,China,Asia,China: Huishui,pepper,Illumina NovaSeq 6000,SU_epidermis_D,not applicable,10,PAIRED,other,METAGENOMIC,plant metagenome,ILLUMINA,2020-11-05T00:00:00Z,biological replicate 1,Sample_106,SRP286471
SRR12778014,OTHER,302,24003640104,PRJNA667562,SAMN16378121,Metagenome or environmental,7547690247,CHINESE ACADEMY OF SCIENCES,2018-08-23,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",SRX9247617,China,Asia,China: Huishui,pepper,Illumina NovaSeq 6000,SU_epidermis_H,not applicable,9,PAIRED,other,METAGENOMIC,plant metagenome,ILLUMINA,2020-11-05T00:00:00Z,biological replicate 3,Sample_99,SRP286471
SRR12778015,OTHER,302,20960803468,PRJNA667562,SAMN16378120,Metagenome or environmental,6757603248,CHINESE ACADEMY OF SCIENCES,2018-08-23,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",SRX9247616,China,Asia,China: Huishui,pepper,Illumina NovaSeq 6000,SU_epidermis_H,not applicable,8,PAIRED,other,METAGENOMIC,plant metagenome,ILLUMINA,2020-11-05T00:00:00Z,biological replicate 2,Sample_98,SRP286471
SRR12778016,OTHER,302,22817685198,PRJNA667562,SAMN16378119,Metagenome or environmental,7251222351,CHINESE ACADEMY OF SCIENCES,2018-08-23,public,"fastq,sra","gs,ncbi,s3","gs.US,ncbi.public,s3.us-east-1",SRX9247615,China,Asia,China: Huishui,pepper,Illumina NovaSeq 6000,SU_epidermis_H,not applicable,7,PAIRED,other,METAGENOMIC,plant metagenome,ILLUMINA,2020-11-05T00:00:00Z,biological replicate 1,Sample_97,SRP286471
~~~
{: .output}

If I use the next piece of code, I can obtain the first column of all the rows, 
which is the information inside the `Run` column, the same name that each forward and reverse reads files has.

~~~
$ cat metadata/SraRunTable.txt| sed -n '1!p' | while read line; do read=$(echo $line | cut -d',' -f1); echo $read;done
~~~
{: .bash}

~~~
SRR12778013
SRR12778014
SRR12778015
SRR12778016
SRR12778017
SRR12778018
SRR12778019
SRR12778020
SRR12778021
SRR12778022
SRR12778023
SRR12778024
~~~
{: .output}

## Creation of the all-around program

One of the great goals of this project is to create a program that can take the 
information that I obtained in this little chapter and autimatically do all the 
needed steps to process the data. I will do the first lines of that program here. 
This first little program will take the information found inside the 
`SraRunTable.txt` and will download the reads from all the libraries listed 
inside. I will name it `down-reads.sh`. This and all the scripts will be locaten 
inside the [scripts-folder](https://github.com/Bedxxe/clavibacter/tree/main/scripts) for this repository:

~~~
$ cat down-reads.sh
~~~
{: .bash}

~~~
#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to download, trim and move the libraries information.

# This program requires that you give 1 input data. 1) where this 
#SraRunTable is located.

metd=$1 #Location to the SraRunTable.txt

root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
# Now we will define were the reads are:
runs='reads'

# CREATING NECCESARY FOLDERS
mkdir reads

# DOWNLOADING THE DATA

#Let's use the next piece of code to download the data
cat $metd  |  sed -n '1!p' | while read line;  do read=$(echo $line | cut -d',' -f1); fasterq-dump -S $read -p -e 8 -o $read ; done
mv *.fastq reads/
# The -e flag can be customized. This indicates the number of threads used to do this task.

# MANAGING THE DATA

# We will change the names of the reads files. They have a sufix that makes impossible
#to be read in a loop
ls $runs | while read line ; do new=$(echo $line | sed 's/_/-/g'); mv $runs/$line $runs/$new; done

# Now, we will create a file where the information of the run labes can be located
cat $metd  | sed -n '1!p' | while read line; do read=$(echo $line | cut -d',' -f1); echo $read ; done > run-labels.txt
mv run-labels.txt metadata/
~~~
{: .output}


I will like to write about some of the steps of this program. According to my 
experience, the presence of a `_` in the names can cause issues with some programs. This is why the second step of this script changes the names fo all the downloaded 
files. The last line is to create a file where we can locate the name of each 
library. This information will be located inside the `run-labels.txt` file.


<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
