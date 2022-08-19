file=$1

#This script will run in a while over all reference genomes

#file is a reference genome
#reference=subsp_capsici_1101  ; usar para extraer las lineas de archivo sed '1q;d' index-genomes.txt
for infile in ex-reads-clavi/fastq/fastq/capsicum/*-clav-1.fq
do 
	base=$(basename ${infile} -clav-1.fq)
	dir=$(dirname ${infile} )
	bowtie2 -x genomas-clavi/clavi-genomes/index-bowtie/${file} -1 ${dir}/${base}-clav-1.fq -2 ${dir}/${base}-clav-2.fq  --no-unal -p 12 -S genomas-clavi/clavi-genomes/out-bowtie/${file}-${base}.sam
done
