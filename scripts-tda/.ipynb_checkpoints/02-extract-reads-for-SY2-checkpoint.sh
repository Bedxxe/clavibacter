fol=$1
file=$2
#This script will run over all SY to extraxct all reads from *sam file's
base2=$(basename ${file} .txt)
for infile in output-tda/$fol/sorted/*.sam
do 
	base=$(basename ${infile} .sam)
	head -4 ${infile} > output-tda/$fol/sam_extract/${base2}-${base}.bam  #get the head the .sam file's
	samtools view ${infile} | grep -P "${file}\t" >> output-tda/$fol/sam_extract/${base2}-${base}.bam # append the read selected 
	done
