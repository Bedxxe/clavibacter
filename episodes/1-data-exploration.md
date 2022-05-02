# Clavibacter project

<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >


## First steps on exploring the data

<a href="{{ page.root }}/figures/bw-hokusai.jpg">
  <img src="{{ page.root }}/figures/bw-hokusai.jpg" alt="Black and white representation of a wave made by Hokusai" />
</a>

### Introductions

I was invited to participate in a project to uncover the presence of a group of bacterial lineajes from the genus 
_Clavibacter_ in different cultivars. As planned by my Sensei (Professor Nelly Selem-Mojica), we will use 
both approximation of metagenomics (_i.e._ Metabarcoding and Shotgun) to try to obtain some answers 
from the questions that were formulated because of the widely presence of these lineages.

### The data

I want to start this project with the exploration of some of the shotgun data. The 
next two lists of data was given by the team:

#### The first one with a literature exploration:

|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| Plant                                               |  Number of samples     | Type of data  | Link to the data | Paper| 
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

#### A second set of data obtained by a general exploration on the NCBI database:

|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| Plant                                               |  Number of samples     | Type of data  | Link to the data | 
| _Medicago sativa_| 18 | 16S  | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=Medicago+sativa+metagenomic+NOT+(amplicon)+NOT+(RNA)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------|   
| _Zea mays_ | 41 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=metagenomics+zea+mays+NOT+(amplicon)+NOT+(RNA)+NOT+(PCR)+NOT+(454)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 
| _Capsicum annuum_ | 4 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=capsicum+annuum+metagenome+NOT+(amplicon)+NOT+(RNA)+NOT+(454)+NOT+(PCR)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 
| _Titricum_ | 54 | Shotgun | [Link](https://www.ncbi.nlm.nih.gov/sra/?term=triticum+metagenome+NOT+(amplicon)+NOT+(RNA)+NOT+(454)+NOT+(PCR)) | 
|-------------------+-------------------------------------------------------------------------------------------------------------------| 

### First exploration of *Capsicum* dataset from Choi-2020

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

#### Running kraken2

I prepared a little program that is capable to extract this information, create an organization for the folders where the results are going to be located, and 
run the `kraken2` in all the files that we have. The program has the name `kraken-reads.sh` and can be located inside the [scripts folder](https://github.com/Bedxxe/clavibacter/blob/main/scripts/kraken-reads.sh) of this repository.
Let's see what is inside:

~~~
$ cat kraken-reads.sh
~~~
{: .bash}

~~~
#!/bin/sh
# This is a program that is going to pick a SraRunTable of metadata and 
#extract the run label to run the next programs.

# This program requires that you give 2 input data. 1) where this 
#SraRunTable is located, and 2) where the kraken database has been saved.

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
mkdri -p taxonomy/biom-files

# DOWNLOADING THE DATA

#Let's use the next piece of code to download the data
cat $metd  |  sed -n '1!p' | while read line; do fasterq-dump -S $line -p -e 8 -o $line ; done
mv *.fastq reads/

# MANAGING THE DATA

# We will change the names of the reads files. They have a sufix that makes impossible
#to be read in a loop
ls $runs | while read line ; do new=$(echo $line | sed 's/_/-/g'); mv $runs/$line $runs/$new; done

# Now, we will create a file where the information of the run labes can be located
cat $metd  | while read line; do read=$(echo $line | cut -d',' -f1); echo $read ; done > run-labels.txt
mv run-labels.txt metadata/

# TAXONOMIC ASSIGNATION WITH KRAKEN2

cat metadata/run-labels.txt | while read line; do mkdir taxonomy/kraken/$line; file1=$(echo $runs/$line-1.fastq); file2=$(echo $runs/$line-2.fastq) ; echo '\n''working in run' "$line"\ 
#kraken2 --db $kdat --threads 12 --paired $file1 $file2 --output taxonomy/kraken/$line/$line.kraken --report taxonomy/kraken/$line/$line.report \ 
echo '#!/bin/sh''\n''\n'"kraken2 --db $kdat --threads 12 --paired" "$runs/$line"'-1.fastq' "$runs/$line"'-2.fastq' "--output taxonomy/kraken/$line/$line.kraken --report taxonomy/kraken/$line/$line.report" > taxonomy/taxonomy-logs/scripts/$line-kraken.sh; sh taxonomy/taxonomy-logs/scripts/$line-kraken.sh; cp taxonomy/kraken/$line/$line.report taxonomy/kraken/reports;done
~~~
{: .output}

There are some things that I would like to highlight on this little program. As it says in the first lines, you will need to provide it with 1) the location of your 
`SraRunTable.txt` file, 2) the location of your `kraken2` database, and 3) location where you located your reads. This program will create in your actual directory 
the folders where the outputs are going to be allocated and little pieces of code in the `taxonomy-logs/scripts` folder that were used to run `kraken2` on each 
pair of reads.

The 26th line is commented because I could not made `kraken2` to run inside the `while loop`. If you can find any solutions or any improvements, please let me know.

You can always change the number of threads to be used.

So, in my case I am going to run it with the following command:
~~~
$ sh kraken-reads.sh metadata/SraRunTable.txt /home/betterlab/kraken2/database/db_kraken2
~~~
{: .bash}

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

By the end of the process, I have a new folder called `taxonomy/`. Let's see what 
it is inside.

~~~
$ tree -L 2
~~~
{: .bash}

~~~
.
├── kraken
│   ├── reports
│   ├── SRR12778013
│   ├── SRR12778014
│   ├── SRR12778015
│   ├── SRR12778016
│   ├── SRR12778017
│   ├── SRR12778018
│   ├── SRR12778019
│   ├── SRR12778020
│   ├── SRR12778021
│   ├── SRR12778022
│   ├── SRR12778023
│   └── SRR12778024
└── taxonomy-logs
    └── scripts

16 directories, 0 files
~~~
{: .output}

Each folder inside `kraken/` contains it's respective `.kraken` and `.report` file. 
Inside the `reports/` folder are all the reports from all the samples that the 
script processed.

Now, we can do this same process with each of the `capsicum` folders. Inside `miscelaneous-capsicum` and `newberry-2020` folders, I will use the next line of code to obtain the desired results:

~~~
$ sh kraken-reads.sh metadata/SraRunTable.txt /home/betterlab/kraken2/database/db_kraken2
~~~
{: .bash}

And for the `newberry-2020`, the next one.
If we expore what it is inside these two folders, we will see the results are 
there:

~~~
$ for i in miscelaneous-capsicum newberry-2020; do echo -e "\n"; tree $i -L 3; done
~~~
{: .bash}

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
    │   ├── ERR5639101
    │   ├── reports
    │   ├── SRR13319509
    │   ├── SRR13319510
    │   └── SRR13319511
    └── taxonomy-logs
        └── scripts

11 directories, 12 files


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
    │   ├── reports
    │   ├── SRR10527387
    │   ├── SRR10527388
    │   └── SRR10527389
    └── taxonomy-logs
        └── scripts

10 directories, 10 files
~~~
{: .output}



