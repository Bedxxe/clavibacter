---
source: md
title: "Taxonomic exploration with R"
---

# Taxonomic exploration with R
<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
## First steps on exploring the data

<img src="/clavibacter/figures/the-fuji-seen-from-the-mishima-pass.jpg" >

## Correction of kraken2 report output

In the last chapter, we processed all the data from _Capsicum_. I processed all the 
_Lycopersicum_ libraries as well in order to have two sets of data from two 
different host-plants to compare. First, I want to repeat that we are interested 
in the different populations of _Clavibacter_ species that are present in the 
holobiont (plant). If we take a look in how `kraken2` assignated the names of the 
_Clavibacter_ OTUs, we will see that we need to do some data trimming before we 
can continue:

~~~
$  grep 'Clavibacter' capsicum/choi-2020/taxonomy/kraken/reports/SRR12778013.report
~~~
{: .language-bash}

~~~
  0.00  544     84      G       1573                    Clavibacter
  0.00  398     151     S       28447                     Clavibacter michiganensis
  0.00  62      62      S1      1874630                     Clavibacter michiganensis subsp. capsici
  0.00  49      49      S1      31965                       Clavibacter michiganensis subsp. tessellarius
  0.00  33      33      S1      31964                       Clavibacter michiganensis subsp. sepedonicus
  0.00  32      32      S1      31963                       Clavibacter michiganensis subsp. nebraskensis
  0.00  29      29      S1      1401995                     Clavibacter michiganensis subsp. californiensis
  0.00  25      25      S1      33014                       Clavibacter michiganensis subsp. insidiosus
  0.00  17      15      S1      33013                       Clavibacter michiganensis subsp. michiganensis
  0.00  2       2       S2      443906                        Clavibacter michiganensis subsp. michiganensis NCPPB 382
  0.00  34      34      S       2768071                   Clavibacter zhangzhiyongii
  0.00  28      0       G1      2626594                   unclassified Clavibacter
  0.00  28      28      S       2860285                     Clavibacter sp. A6099
~~~
{: .output}

If we take a look to the information regarding this [output](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown), we can see that most of the 
_Clavibacter_ species that we are interested on, are clasiffied as subspecies of 
_Clavibacter michiganensis_ (_e.g_ _capsici_, and _tessellarius_). I have created 
a program that will correct this issue. The `trim-clavi-reports.sh` program is located 
inside the [scripts-folder](https://github.com/Bedxxe/clavibacter/tree/main/scripts)


~~~
$ cat trim-clavi-reports.sh
~~~
{: .language-bash}

~~~
#!/bin/bash

#This program is to trim the Clavibacter michiganensis identifiers from a kraken.report
#file

#The program will ask you 1 thing. a) The name of the report file 

repo=$1 #Name of the report file
sufx=$(echo $repo |cut -d'.' -f1)

#Creating output directory
mkdir -p trim-reports


#Obtaining the values before Cmm
sed -n '/Clavibacter michiganensis/q;p' $repo > before-$sufx.txt

#Obtaining the C. michiganensis values
cat $repo | grep  'Clavibacter michiganensis' > cm-$sufx.txt

#Obtaining the values after Cmm
sed -n '/Clavibacter michiganensis/,$p' $repo | grep -v 'Clavibacter michiganensis'> after-$sufx.txt

# Taking the value of C. michiganensis
val1=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f2)

#Obtaining the value that need to be substracted to the C. michiganensis field
i=0
grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1!p' | while read line;
do cou=$(echo $line | cut -d' ' -f2); i=$(($i + $cou)); echo $i ; done > temp
val2=$(tail -n1 temp)

##The new value
val3=$(($val1 - $val2))

#Making the new file of Clavi

##The first line with the unclassified Cm
a=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f1)
b=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f3)
c=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f4)
d=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f5)
e=$(grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1p' | cut  -f6)
echo " "$a"\t"$val3"\t"$b"\t"$c"\t"$d"\t""                  ""unclassified"$e > clavi.txt

##The rest of the species
grep 'michiganensis' cm-$sufx.txt | grep -v 'S2' | sed -n '1!p' | while read line;
do ta=$(echo $line | cut -d' ' -f1); tb=$(echo $line | cut -d' ' -f2);
tc=$(echo $line | cut -d' ' -f3); td=$(echo $line | cut -d' ' -f4);
te=$(echo $line | cut -d' ' -f5); tf=$(echo $line | cut -d' ' -f6);
tg=$(echo $line | cut -d' ' -f9);
echo "  "$ta"\t"$tb"\t"$tc"\t""S""\t"$te"\t""                  "$tf" "$tg >> clavi.txt ;
done

##The last line of the sub-sub species of Cmm
ka=$(grep 'michiganensis' cm-$sufx.txt | tail -n1| cut  -f1)
kb=$(grep 'michiganensis' cm-$sufx.txt | tail -n1| cut  -f2)
kc=$(grep 'michiganensis' cm-$sufx.txt | tail -n1| cut  -f3)
#kd
ke=$(grep 'michiganensis' cm-$sufx.txt | tail -n1| cut  -f5)
kf=$(grep 'michiganensis' cm-$sufx.txt | tail -n1| cut  -f6)
echo " "$ka"\t"$kb"\t"$kc"\t""S1""\t"$ke"\t""                   "$kf >> clavi.txt

#Creating the trimmed report
cat before-$sufx.txt > t-$sufx.report
cat clavi.txt >> t-$sufx.report
cat after-$sufx.txt >> t-$sufx.report

#Moving the new report to trim-reports
mv t-$sufx.report trim-reports/

#Removing temporary files
rm temp
rm cm-$sufx.txt
rm after-$sufx.txt
rm before-$sufx.txt
rm clavi.txt
~~~
{: .language-bash}

Inside each of the `reports` folder, I will use the next line to run the program on all the outputs:
~~~
$ ls *.report | while read line; do sh parche.sh $line; done
~~~
{: .language-bash}

## Using kraken-biom to process the reports files

I will use [kraken-biom](https://github.com/smdabdoub/kraken-biom) to process all the reports and put them into a `biom` file. I will do an example inside the 
`capsicum/choi-2020/` folder:

~~~
$ mkdir biom-files/
$ kraken-biom reports/* --fmt json -o biom-files/choi-2020.biom
~~~
{: .language-bash}

We will repear the same process for all the rest of the author's folders.

## Adjusting the all-aroun program

We will add this last step to the script that we have been constructing. We will 
the user for a new input: 3) the name of a prefix, in this case the name of the 
author. And we will get a new output inside the `biom-files/` folder. The 
`kraken-biom.sh` script inside the [scripts folder](https://github.com/Bedxxe/clavibacter/blob/main/scripts/kraken-reads.sh)

~~~
$ cat kraken-biom.sh
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
aut=$3 #A prefix to name some of the files. In this case, the author name.


root=$(pwd) #Gets the path to the directory of this file, on which the outputs ought to be created 
# Now we will define were the reads are:
runs='reads'

# CREATING NECCESARY FOLDERS
mkdir reads
mkdir -p taxonomy/kraken
mkdir -p taxonomy/taxonomy-logs/scripts
mkdir -p taxonomy/kraken/reports
mkdir -p taxonomy/kraken/krakens
mkdir -p taxonomy/biom-files

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

#CREATING THE BIOM FILE

# Now we will create the biom file using kraken-biom 
kraken-biom taxonomy/kraken/reports/* --fmt json -o taxonomy/biom-files/$aut.biom
~~~
{: .language-bash}

Now, we will obtain the `biom-file` for all the author's libraries in the next 
time we run the script for the other host plants.

## Using R for the analysis

### Loading the packages 


~~~

~~~
{: .language-r}





<img src="/clavibacter/figures/grecas-mitla1.png" alt="Picture of the fretwork on the ruins in Mitla, Oaxaca." >
