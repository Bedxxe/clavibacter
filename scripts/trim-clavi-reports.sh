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