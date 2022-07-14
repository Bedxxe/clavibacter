# Analysis of the data from Solena suministred to Professor Nelly Selem Mojica

## Kraken-biom

Inside the folder that was suministred: `nelly`, I created three new folders. 
I moved all the `*.report.txt` files inside the folder `reports`, and the two 
metadata files to a `metadata` folder:

~~~
$ mkdir reports/
$ mkdir biom-files/
$ mkdir metadata/
$ mv 
$ mv kraken_braken/*.report.txt reports/
~~~
{: .language-bash}

Then, I trimmed all the report files to make a shorter set of names:

~~~
$ ls | while read line; do name=$(echo $line| cut -d'.' -f1,4) ; mv $line $name; done
~~~
{: .language-bash}

There is a problems with the allocation of subspecies in the kraken and bracken 
reports: The subspecies of _Clavibacter_ lineajes are not detected by 
the programs that we usually use to precess them. Inside our _Clavibacter_ 
repository, we developed a [chapter]() explaining how we solve this problem. We developed the 
next program to correct the reports:

~~~
$ cat trim-clavi.sh
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
{: .output}

Inside the `reports` folder, I will use the next line to run the program on all the reads:

~~~
$ ls *.report | while read line; do sh parche.sh $line; done
~~~
{: .language-bash}


I used [kraken-biom](https://github.com/smdabdoub/kraken-biom) to process all the 
reports and put them into a `biom` file: `solena.biom`.

~~~
$ kraken-biom reports/* --fmt json -o biom-files/solena.biom
~~~
{: .language-bash}

This now can be read by RStudio and the package [`phyloseq`](https://joey711.github.io/phyloseq/)

## Analyzing the data with R

### Loading and trimming the data

I will load the needed packages to work with this data:

~~~
> library("phyloseq")
> library("ggplot2")
> library("edgeR")
> library("DESeq2")
> library("pheatmap")
> library("RColorBrewer")
> library("stringr")
> library("tidyverse")
> library("vegan")
~~~
{: .language-r}

The data seems quite large:

~~~
> solena <- import_biom("../../solena/nelly/biom-files/solena.biom")
> solena
~~~
{: .language-r}
~~~
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 12233 taxa and 37 samples ]
tax_table()   Taxonomy Table:    [ 12233 taxa by 7 taxonomic ranks ]
~~~
{: .output}

We will trim the data:

~~~
> solena@tax_table@.Data <- substring(solena@tax_table@.Data, 4)
> colnames(solena@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
> colnames(solena@otu_table@.Data) <- substring(sample_names(solena), 3)
~~~
{: .language-r}

And trim the names of the samples that are too long to be practical:

~~~
> colnames(solena@otu_table@.Data) <- substring(sample_names(solena), 15)
~~~
{: .language-r}

### Metadata

We will load the metadata that was provided by Solena:
~~~
> met.sol <- read.csv2("../../solena/nelly/metadata/fastp_metadat.csv",
           header =  TRUE, row.names = 1, sep = ",")
~~~
{: .language-r}

We will rename the rows in order to them to match the ones that we generated in the 
`solena` object.
~~~
> rownames(met.sol) <- substring(rownames(met.sol),first = 13)
~~~
{: .language-r}

~~~
> met.sol <- sample_data(met.sol)
~~~
{: .language-r}

And finally we will merge the metadata with the `solena` object that has already 
all the taxonomic assignation data:

~~~
> solena<- merge_phyloseq(solena, met.sol)
> solena
~~~
{: .language-r}

~~~
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 12238 taxa and 37 samples ]
sample_data() Sample Data:       [ 37 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 12238 taxa by 7 taxonomic ranks ]
~~~
{: .output}

### Emploring the depth of the data

We will create an object to see the depth of the sequenciation:

~~~
> dep.sol <- data.frame(Samples = as.factor(colnames(solena@otu_table@.Data)),
                    Reads = sample_sums(solena),
                    Cultivo = as.factor(solena@sam_data@.Data[[3]]))
~~~
{: .language-r}

In order to have all the samples from the same plant together, we will change the 
order of factors in the `Cultivo` column:

~~~
> dep.sol <- arrange(.data = dep.sol, dep.sol$Cultivo)
> dep.sol$Samples <- factor(dep.sol$Samples, levels = dep.sol$Samples)
~~~
{: .language-r}


And we will plot it with the next line of code:
~~~
> ggplot(data = dep.sol, mapping = aes(x = Samples, y = Reads, fill = Cultivo))+
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#d95f02","#1b9e77","#7570b3")) +
    theme(text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

~~~
{: .language-r}

<img src="/clavibacter/figures/sol-01.png" >

It seems that the depth of reads varies a lot. We will need to normalize the data 
for futere analysis.

### Relative abundance and dominant Phyla

We will transform the data to relative abundance in order to compare the libraries 
against each other:
~~~
> sol.rel <- transform_sample_counts(solena, function(x) x*100/ sum(x))
~~~

We need to cancatenate all the OTUs that posses the belon to the same phylum, thus 
helping the visualization of the data:

{: .language-r}

~~~
> top.sol <- tax_glom(sol.rel,taxrank = rank_names(sol.rel)[2])
~~~
{: .language-r}

And finally maintain the dominant phyla (more than 0.3% of abundance in all the samples)
~~~
> top.sol <- filter_taxa(top.sol, function(x) mean(x) > 0.1, TRUE)
~~~
{: .language-r}

With this, we will prepare tohe data.frame to plot the results:
~~~
> sol.ftot <- psmelt(top.sol)
# Trimming the factors to cluster the libraries from the same plant
> sol.ftot$Sample <- as.factor(sol.ftot$Sample)
> sol.ftot <- arrange(.data = sol.ftot, sol.ftot$Cultivo)
> sol.ftot$Sample <- factor(sol.ftot$Sample, levels = unique(sol.ftot$Sample))
~~~
{: .language-r}

And we can plot the result with the next line of code:

~~~
> ggplot(data = sol.ftot, mapping = aes(x= Sample, y = Abundance, fill = Phylum, color = Cultivo))+
    geom_bar(stat = "identity", position = "stack", size=1) +
    scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3"))+
    theme(text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

~~~
{: .language-r}

<img src="/clavibacter/figures/sol-02.png" >

### _Clavibacter_ 

To observe the abundance of _Clavibacter_ inside the data, we will first extract 
the information belonging to this group of OTUs to a new object from the `solena` 
dataset:

~~~
> sol.cla <- subset_taxa(solena, Genus == "Clavibacter")
~~~
{: .language-r}

We want to know the number of reads that belongs to this:

~~~
> sol.c.fra <- psmelt(sol.cla)
# reordering the factors
> sol.c.fra $Sample <- as.factor(sol.c.fra $Sample)
> sol.c.fra  <- arrange(.data = sol.c.fra , sol.c.fra $Cultivo)
> sol.c.fra $Sample <- factor(sol.c.fra $Sample, levels = unique(sol.c.fra $Sample))
~~~
{: .language-r}

~~~
> ggplot(data = sol.c.fra, mapping = aes(y= Abundance, x = Sample, fill = Species, color = Cultivo))+
    geom_bar(position = "stack", stat = "identity", size=1)+
    scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3"))+
    #facet_wrap(~Cultivo)
    theme(text = element_text(size= 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

~~~
{: .language-r}

<img src="/clavibacter/figures/sol-03.png" >

If we want to see how much is this of the total number of reads:
~~~
> sol.c.fra$ClaRelative <- sample_sums(sol.cla)/sample_sums(solena)
~~~
{: .language-r}

Plotting the obtained data:

~~~
> ggplot(data = sol.c.fra, mapping = aes(y= ClaRelative, x = Sample, fill = Species, color = Cultivo))+
    geom_bar(position = "stack", stat = "identity", size=1.5)+
    scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3"))+
    #facet_wrap(~Cultivo)
    theme(text = element_text(size= 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
~~~
{: .language-r}

<img src="/clavibacter/figures/sol-04.png" >

We can also see how much they differ according to the zone from where the samples 
were taken:
~~~
> ggplot(data = sol.c.fra, mapping = aes(y= ClaRelative, x = Sample, fill = Species, color = Cultivo))+
    geom_bar(position = "stack", stat = "identity", size=1.5)+
    scale_color_manual(values = c("#d95f02","#1b9e77","#7570b3"))+
    facet_wrap(~Origen) +
    theme(text = element_text(size= 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

~~~
{: .language-r}

<img src="/clavibacter/figures/sol-05.png" >

### Multivatiate analysis

First, we will re-assign the metadata to a nwe object to use it in the analysis:
~~~
> meta.sol <- read.csv2("../../solena/nelly/metadata/fastp_metadat.csv",
                      header =  TRUE, row.names = 1, sep = ",")
# Trimming the sample names:
> rownames(meta.sol) <- substring(rownames(meta.sol),first = 13)
~~~
{: .language-r}

Also, the NAs in the metadata will be removed to unknown data so as to preclude an 
error from the `adonis` function:

~~~
meta.sol[is.na(meta.sol)] = "Unknown"
~~~
{: .language-r}

We will evaluate the data with a PERMANOVA test around all the taxonomica 
levels. We first begin with Phylum. We will use the `tax_glom()` function to 
cut up to the phylum level. From this new `phyloseq` object, we will extract the 
OTU information from the `otu_table` and we will transpose it:

~~~
> p.sol <- tax_glom(solena,taxrank = rank_names(solena)[2])
> p.solena <- t(p.sol@otu_table@.Data)
~~~
{: .language-r}

We will use the `adonis()` function that comes with the `vegan` package to do the 
multivariate analysis:

~~~
adonis(p.solena ~ Cultivo , data = meta.sol, permutations = 9999, method = "bray")
~~~
{: .language-r}


~~~
Call:
adonis(formula = p.solena ~ Cultivo, data = meta.sol, permutations = 9999,      method = "bray") 

Permutation: free
Number of permutations: 9999

Terms added sequentially (first to last)

          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
Cultivo    2   0.27408 0.137041  1.7609 0.09386 0.1261
Residuals 34   2.64596 0.077822         0.90614       
Total     36   2.92004                  1.00000   
~~~
{: .output}

It seems that the diversity of OTUs is not well correlated with the plant that 
was extracted at the phylum level.

Since we have the information of the phenological stages when the plant was sampled 
we can try to see how much the OTU diversity change with a nested design where the 
phenological stages are nested inside the plant host (Cultivo). We will use the 
`strata` option from `adonis` to accomplish this.

~~~
adonis(p.solena ~ Cultivo/Etapa_fenologica , data = meta.sol, permutations = 9999, 
       method = "bray", strata = meta.sol$Cultivo)
~~~
{: .language-r}

~~~
Call:
adonis(formula = p.solena ~ Cultivo/Etapa_fenologica, data = meta.sol,      permutations = 9999, method = "bray", strata = meta.sol$Cultivo) 

Blocks:  strata 
Permutation: free
Number of permutations: 9999

Terms added sequentially (first to last)

                         Df SumsOfSqs  MeanSqs F.Model
Cultivo                   2   0.27408 0.137041  1.8186
Cultivo:Etapa_fenologica  7   0.61133 0.087332  1.1589
Residuals                27   2.03463 0.075357        
Total                    36   2.92004                 
                              R2 Pr(>F)
Cultivo                  0.09386 0.1931
Cultivo:Etapa_fenologica 0.20936 0.1931
Residuals                0.69678       
Total                    1.00000   
~~~
{: .output}

Next, there are two tables where the results for each taxonomic level is shown:


<img src="/clavibacter/figures/sol-06.png" >


<img src="/clavibacter/figures/sol-07.png" >


### Beta diversity

That the data is is not significantly separated by the host-plant origin can be 
seen in a NMDS ordination. First, we need to obtain the ordination with some 
of the `phyloseq` tools:

~~~
> sol.ord <- ordinate(physeq = solena, method = "NMDS", distance = "bray")
~~~
{: .language-r}

With this result, we can plot the NMDS space:

~~~
> plot_ordination(physeq = solena, ordination = sol.ord,
                  color = "Cultivo") +
    stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = Cultivo))+
    geom_point(size=5 , alpha = 0.5)+
    theme(text = element_text(size = 18))
~~~
{: .language-r}


<img src="/clavibacter/figures/sol-08.png" >
