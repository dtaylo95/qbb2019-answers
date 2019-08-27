#!/bin/bash

GENOME=../genomes/BDGP6.fa
ANNOTATION=../genomes/BDGP6.Ensembl.81.gtf
THREADS=4

for SAMPLE in SRR072893 SRR072903 SRR072905 SRR072915
do
	echo "*** Processing $SAMPLE"
	echo "*** Running fastqc"
	fastqc -t $THREADS ../../data/rawdata/$SAMPLE.fastq
	echo "*** Running hisat2"
	hisat2 -x ../genomes/BDGP6 -U ../../data/rawdata/$SAMPLE.fastq -p 4 -S $SAMPLE.sam
	echo "*** Running samtools sort"
	samtools sort -@ 4 -O BAM $SAMPLE.sam > $SAMPLE.bam
	echo "*** Running samtools index"
	samtools index -@ 4 -b $SAMPLE.bam
	echo "*** Running StringTie"
	stringtie $SAMPLE.bam -G $ANNOTATION -o $SAMPLE.gtf -p 4 -e -B
done