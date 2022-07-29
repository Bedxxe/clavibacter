# Clavibacter Project

![grecas-mitla1](https://user-images.githubusercontent.com/67386612/177206178-e3f3a9ff-608c-4fce-bd58-aa4b29ef6990.png)

I was invited to participate in a project to uncover the presence of a group of bacterial lineajes from the genus 
_Clavibacter_ in different cultivars. As planned by my Sensei (Professor Nelly Selem-Mojica), we will use 
both approximation of metagenomics (_i.e._ Metabarcoding and Shotgun) to try to obtain some answers 
from the questions that were formulated because of the widely presence of these lineages.

This repository is to manage, create, and allocate the data from the Clavibacter project (Selem-Mojica, et al. 2022)

The main page:
https://bedxxe.github.io/clavibacter/

Episodes:

https://bedxxe.github.io/clavibacter/episodes/1-data-exploration.html

https://bedxxe.github.io/clavibacter/episodes/2-Taxonomic-assignation.html

https://bedxxe.github.io/clavibacter/episodes/3-Taxonomic-exploration-R.html

https://bedxxe.github.io/clavibacter/episodes/4-Clavi-extraction.html

https://bedxxe.github.io/clavibacter/episodes/5-Searching-for-Cmm.html

https://bedxxe.github.io/clavibacter/episodes/6-16S-database-selection.html

https://bedxxe.github.io/clavibacter/episodes/7-16S-Taxonomic-exploration-R.html

https://bedxxe.github.io/clavibacter/episodes/8-Automatization.html


![changuito-tavehua](https://user-images.githubusercontent.com/67386612/166823738-87a8da81-11d4-4dcb-88f9-06206c5bd824.png)


After the completion of all the episodes, I have amassed all the data into the 
`betterlab` server from the Selem-Lab. In the next lines I will try to describe 
where these files are. All the files are located in the next directory:

~~~
$ /home/betterlab/diego/clavibacter
~~~
{: .language-bash}

Inside this folder, there are a set of 7 directories that contain all the 
information:

~~~
$ tree -L 1
~~~
{: .language-bash}

~~~
.
├── 16S
├── blast
├── ex-reads-clavi
├── ex-reads-cmm
├── reference
├── scripts
└── shotgun

7 directories, 0 files
~~~
{: .output}

## 16S/

Here is where all the information from the 16S sets of data:

* Episode 6: [16-S database selection](https://bedxxe.github.io/clavibacter/episodes/6-16S-database-selection.html)
* Episode 7: [16-S taxonomic exploration with R](https://bedxxe.github.io/clavibacter/episodes/7-16S-Taxonomic-exploration-R.html
)

Inside the `16S/` folder, there are three subdirectories. Each of them amasses the 
data from a different plant host and its metadata files:

~~~
$ tree -L 2
~~~
{: .language-bash}

~~~
.
├── medicago-sativa
│   ├── miscellaneous
│   ├── miscellaneous-SraRunTable.txt
│   └── miscellaneous-SRR_Acc_List.txt
├── tuberosum
│   ├── miscellaneous-tuberosum
│   ├── miscellaneous-tuberosum-SraRunTable.txt
│   ├── miscellaneous-tuberosum-SRR_Acc_List.txt
│   ├── shi-2019
│   ├── shi-2019-SraRunTable.txt
│   └── shi-2019-SRR_Acc_List.txt
└── zea-mayz
    ├── m-wu-2021
    ├── m-wu-2021-SraRunTable.txt
    └── m-wu-2021-SRR_Acc_List.txt

7 directories, 8 files
~~~
{: .output}

Inside the two _tuberosum_ directories (_i.e._ miscellaneous-tuberosum and shi-2019), and the _zea-mays_ directory m-wu-2021, the user can find the next 
folder structure:

~~~
├── metadata
├── reads
└── taxonomy

3 directories, 0 files
~~~
{: .language-bash}

Inside each the `taxonomy/` folder, the user can find all the information from 
the `kraken2` taxonomic assignation:

~~~
tuberosum/miscellaneous-tuberosum/taxonomy/
├── biom-files
│   └── tuber.biom
├── kraken
│   ├── krakens
│   └── reports
└── taxonomy-logs
    └── scripts

6 directories, 1 file
~~~
{: .language-bash}

The `medicago-sativa/miscellaneous/` is the only one which has a different folder 
structure. Beyond the `metadata/` `reads/`, and `taxonomy/` folders, the user 
will find a set of files that are the output of the [Episode 6](https://bedxxe.github.io/clavibacter/episodes/6-16S-database-selection.html) and 
that were useful for the 16S-database selection. 

## blast/

Inside the `blast/` folder is the information used in the next episode:

* Episode 5: [Searching for Cmm](https://bedxxe.github.io/clavibacter/episodes/5-Searching-for-Cmm.html)


Inside the data `databeses/` subfolder is the files created by the command 
`makeblastdb` used in the mentioned episode.


If we explore the contents of the `capsi-choi/` subfolder, we will found 
three scripts that were used in the Episode 5, and two output folders: 
~~~
├── blast.sh
├── dividing-ind-reads.sh
├── enhanced-database-output
├── extracted-cmm-choi
├── first-database-output
└── sec-blast.sh

3 directories, 3 files
~~~
{: .language-bash}

The `extracted-cmm-choi/` contains the files that we used as an example. Indide 
this folder there is a subdirectory that contains the individual reads that 
were used to run the blast:

~~~
$ ls capsi-choi/extracted-cmm-choi/
~~~
{: .language-bash}

~~~
cmm-capsi-choi-2020-ERR5639101-1.fasta   cmm-capsi-choi-2020-SRR12778019-2.fasta
cmm-capsi-choi-2020-ERR5639101-2.fasta   cmm-capsi-choi-2020-SRR12778020-1.fasta
cmm-capsi-choi-2020-SRR10527387-1.fasta  cmm-capsi-choi-2020-SRR12778020-2.fasta
cmm-capsi-choi-2020-SRR10527387-2.fasta  cmm-capsi-choi-2020-SRR12778021-1.fasta
cmm-capsi-choi-2020-SRR10527389-1.fasta  cmm-capsi-choi-2020-SRR12778021-2.fasta
cmm-capsi-choi-2020-SRR10527389-2.fasta  cmm-capsi-choi-2020-SRR12778022-1.fasta
cmm-capsi-choi-2020-SRR12778013-1.fasta  cmm-capsi-choi-2020-SRR12778022-2.fasta
cmm-capsi-choi-2020-SRR12778013-2.fasta  cmm-capsi-choi-2020-SRR12778023-1.fasta
cmm-capsi-choi-2020-SRR12778014-1.fasta  cmm-capsi-choi-2020-SRR12778023-2.fasta
cmm-capsi-choi-2020-SRR12778014-2.fasta  cmm-capsi-choi-2020-SRR12778024-1.fasta
cmm-capsi-choi-2020-SRR12778015-1.fasta  cmm-capsi-choi-2020-SRR12778024-2.fasta
cmm-capsi-choi-2020-SRR12778015-2.fasta  cmm-capsi-choi-2020-SRR13319509-1.fasta
cmm-capsi-choi-2020-SRR12778016-1.fasta  cmm-capsi-choi-2020-SRR13319509-2.fasta
cmm-capsi-choi-2020-SRR12778016-2.fasta  cmm-capsi-choi-2020-SRR13319510-1.fasta
cmm-capsi-choi-2020-SRR12778017-1.fasta  cmm-capsi-choi-2020-SRR13319510-2.fasta
cmm-capsi-choi-2020-SRR12778017-2.fasta  cmm-capsi-choi-2020-SRR13319511-1.fasta
cmm-capsi-choi-2020-SRR12778018-1.fasta  cmm-capsi-choi-2020-SRR13319511-2.fasta
cmm-capsi-choi-2020-SRR12778018-2.fasta  ind-reads
cmm-capsi-choi-2020-SRR12778019-1.fasta
~~~
{: .output}


## ex-reads-clavi/

Here is all the _Clavibacter_ sequences extracted from the entire libraries. 
This was made in the next episode:

* Episode 4: [Clavi extraction](https://bedxxe.github.io/clavibacter/episodes/4-Clavi-extraction.html)

If we explore inside the contents of this folder, we will see that the outputs 
are in `fasta` and `fastq` format:

~~~
$ tree -L 2 ex-reads-clavi/
~~~
{: .language-bash}

~~~
ex-reads-clavi/
├── fasta
│   ├── capsicum
│   ├── lycopersicum
│   ├── triticum
│   ├── tuberosum
│   └── zea-mays
└── fastq
    ├── capsicum
    ├── lycopersicum
    ├── triticum
    ├── tuberosum
    └── zea-mays

12 directories, 0 files
~~~
{: .output}


## ex-reads-cmm/

Here is all the _Clavibacter michiganensis_ sequences extracted from the entire libraries. 
This was made in the next episode:

* Episode 4: [Clavi extraction](https://bedxxe.github.io/clavibacter/episodes/4-Clavi-extraction.html)

If we explore inside the contents of this folder, we will see that the outputs 
are in `fasta` format:

~~~
$ tree -L 1 ex-reads-cmm/
~~~
{: .language-bash}

~~~
├── capsicum
├── lycopersicum
├── triticum
├── tuberosum
└── zea-mays

5 directories, 0 files
~~~
{: .output}

## reference/

Inside this directory, we can found the data that was used to create the 
database for the blast run. It was covered in the next episode:

* Episode 5: [Searching for Cmm](https://bedxxe.github.io/clavibacter/episodes/5-Searching-for-Cmm.html)


~~~
├── all-genomes.fasta
├── change-head.sh
├── clavi-genomes
│   ├── Clavibacter_michiganensis_subsp_capsici_1101.fna
│   ├── Clavibacter_michiganensis_subsp_capsici_1106.fna
│   ├── Clavibacter_michiganensis_subsp_capsici_1207.fna
│   ├── Clavibacter_michiganensis_subsp_capsici_PF008.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_ATCC_10253.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_CFBP_1195.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_CFBP_2404.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_CFBP_6488.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_LMG_3663.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_R1-1.fna
│   ├── Clavibacter_michiganensis_subsp_insidiosus_R1-3.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_419B.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_44.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_61-1.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_7580.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_A6096.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_CFBP_7577.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_CIBA.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_DOAB_395.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_DOAB_397.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_HF4.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_NCPPB_2581.fna
│   ├── Clavibacter_michiganensis_subsp_nebraskensis_SL1.fna
│   ├── Clavibacter_michiganensis_subsp_sepedonicus_ATCC33113.fna
│   ├── Clavibacter_michiganensis_subsp_sepedonicus_CFIA-Cs3N.fna
│   ├── Clavibacter_michiganensis_subsp_sepedonicus_CFIA-CsR14.fna
│   ├── Clavibacter_michiganensis_subsp_tessellarius_ATCC_33566.fna
│   ├── Clavibacter_michiganensis_subsp_tessellarius_DOAB_609.fna
│   ├── cmm-contigs.fa
│   └── trim-c-genomes
├── cmm
│   ├── cmm-selem
│   └── other-cmm
├── enhanced-datab.fasta
├── head-genomes-trim.sh
├── outgroups
│   ├── Agrobacterium.fna
│   ├── Nostoc.fna
│   └── raw
├── plasmids
│   ├── Cryobacterium-sp.-LW097-plasmid-unnamed1.fasta
│   ├── Microbacterium-hominis-strain-PDNC016-plasmid-unnamed.fasta
│   ├── Rathayibacter-sp.-VKM-Ac-2760-plasmid-unnamed1.fasta
│   ├── Rathayibacter-sp.-VKM-Ac-2760-plasmid-unnamed2.fasta
│   ├── Rhodococcus-hoagii-JID03-56-plasmid-pVAPN03-56.fasta
│   └── trim-header-genomes
└── plasmids.fasta

9 directories, 41 files
~~~
{: .output}

Inside the `clavi-genomes/` are located all the _Clavibacter_ lineajes that are 
not _michiganensis_. 

The `cmm/` directory, allocates the _Clavibacter michiganensis_ that 
was used as reference (`cmm-selem/`) and genomes from other Cmm (`other-cmm/`).

In the `outgroups/` folder, one can find the sequences of the two bacterial 
lineages that were used as outgroups.

Inside the `plasmids/` folder are the sequences corresponding to plasmids found to 
gave a good match for the "Cmm" sequences, when they were compared against the 
NCBI database (as explained in Episode 5)

The `fasta` files are the resulting data from concatenating all the sequences 
that are part of the final database.

## scripts/

This folder is made of the scripts used along all the episodes. It is analogous to 
the `scripts/` folder found in the [GitHub repository](https://github.com/Bedxxe/clavibacter)

## shotgun/

Inside this folder, we can find a structure that is very similar to what we 
found inside `16S/`. The data that is located here is the result of the next 
episodes:

* Episode 1: [Data exploration](https://bedxxe.github.io/clavibacter/episodes/1-data-exploration.html)
* Episode 2: [Taxonomic assignation](https://bedxxe.github.io/clavibacter/episodes/2-Taxonomic-assignation.html)
* Episode 4: [Clavi extraction](https://bedxxe.github.io/clavibacter/episodes/4-Clavi-extraction.html)

As well as in the `16S/` folder, there is a subdirectory for each host plant 
with the same structure described previously:

~~~
$  tree -L 3
~~~
{: .language-bash}

~~~
.
├── capsicum
│   ├── choi-2020
│   │   ├── metadata
│   │   ├── reads
│   │   └── taxonomy
│   ├── metadata
│   │   ├── choi-2020-SraRunTable.txt
│   │   ├── choi-2020-SRR_Acc_List.txt
│   │   ├── miscelaneous-capsicum-SraRunTable.txt
│   │   ├── miscelaneous-capsicum-SRR_Acc_List.txt
│   │   ├── newberry-2020-SraRunTable.txt
│   │   └── newberry-2020-SRR_Acc_List.txt
│   ├── miscelaneous-capsicum
│   │   ├── metadata
│   │   ├── reads
│   │   └── taxonomy
│   └── newberry-2020
│       ├── metadata
│       ├── reads
│       └── taxonomy
├── lycopersicum
│   ├── barajas-2020
│   │   ├── metadata
│   │   ├── reads
│   │   └── taxonomy
│   ├── metadata
│   │   ├── barajas-2020-SraRunTable.txt
│   │   ├── barajas-2020-SRR_Acc_List.txt
│   │   ├── newberry-2020-SraRunTable.txt
│   │   └── newberry-2020-SRR_Acc_List.txt
│   └── newberry-2020
│       ├── kraken-reads.sh
│       ├── metadata
│       ├── reads
│       └── taxonomy
├── tuberosum
│   ├── metadata
│   │   ├── shi-2019-SraRunTable.txt
│   │   └── shi-2019-SRR_Acc_List.txt
│   └── shi-2019
│       ├── metadata
│       ├── reads
│       └── taxonomy
└── zea-mays
    ├── akinola-2021
    │   ├── metadata
    │   ├── reads
    │   └── taxonomy
    ├── metadata
    │   ├── akinola-2021-SraRunTable.txt
    │   ├── akinola-2021-SRR_Acc_List.txt
    │   ├── miscellaneous-mays-SraRunTable.txt
    │   ├── miscellaneous-mays-SRR_Acc_List.txt
    │   ├── xiong-2021-SraRunTable.txt
    │   ├── xiong-2021-SRR_Acc_List.txt
    │   ├── x-wu-2021-SraRunTable.txt
    │   └── x-wu-2021-SRR_Acc_List.txt
    ├── miscellaneous-mays
    │   ├── metadata
    │   ├── reads
    │   └── taxonomy
    ├── xiong-2021
    │   ├── metadata
    │   ├── reads
    │   └── taxonomy
    └── x-wu-2021
        ├── metadata
        ├── reads
        └── taxonomy

48 directories, 21 files
~~~
{: .output}


![grecas-mitla1](https://user-images.githubusercontent.com/67386612/177206200-34922b72-228a-43cc-9ecf-4a0ffab66036.png)

