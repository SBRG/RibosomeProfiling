#!/bin/bash
#
#if [ "$#" -ne 1 ] || ! [ -d "$1" ]; then
#  echo "Usage: $0 bowtie_index_dir <untrimmed_fastq_1> <untrimmed_fastq_2>...<untrimmed_fastq_N>" &> log.txt
#  exit 1
#fi
#
################
#make sure that bowtie index is one folder up from the sample fastq files
#Also, might want to push stderr and stdout to a logfile by using &> log.txt


for input in "$@"
do
	fastq="$input"
	echo 'Now processing '$fastq
	cutadapt -a AGATCGGAAGAGCACACGTCT -m 15 $fastq > ${fastq%.*}_trimmed.fastq

	for index in ec_all_btw ec_trna_btw ec_rRNA_ext_btw ec_RNA_other_btw 
	do
		echo 'Now mapping against the following index: '$index    
		bowtie -n 1 -q -p 4 -S ../$index ${fastq%.*}_trimmed.fastq ${fastq%.*}_${index}_aligned.sam
	done
	
	echo 'Now mapping against all RNA species index'
	bowtie -n 1 -q -p 4 ../ec_RNA_all_btw ${fastq%.*}_trimmed.fastq --un ${fastq%.*}_${index}_unaligned.fastq > reads.fastq
	rm reads.fastq
	
	echo 'Now mapping unaligned reads to reference genome'
	bowtie -n 1 -q -p 4 -S ../ec_all_btw ${fastq%.*}_${index}_unaligned.fastq ${fastq%.*}_non_RNA_aligned.sam

    python ribo_profiling_read_len_dist.py ${fastq%.*}*trna*.sam ${fastq%.*}*rRNA*.sam ${fastq%.*}*RNA_other*.sam ${fastq%.*}*non_RNA*.sam ${fastq%.*}*all*.sam
    rm ${fastq%.*}*non_RNA*.gff
    python makegff.py --flip ${fastq%.*}*non_RNA*.sam
    python makegff.py --flip ${fastq%.*}*all_btw*.sam
done


