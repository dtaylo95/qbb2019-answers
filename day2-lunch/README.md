## EXERCISE 1

# To create a fastq file that only contains the first 10,000 reads from SRR072893.fastq
head -40000 ../rawdata/SRR072893.fastq > SRR072893.10k.fastq

# To verify that the resulting file has 10,000 reads
grep "^@SRR072893" SRR072893.10k.fastq | wc -l

# To create a quality control report for the truncated .fq file
fastqc -t 4 SRR072893.10k.fastq

# To align the reads against the generated indices
hisat2 -x ../genomes/BDGP6 -U SRR072893.10k.fastq -p 4 -S SRR072893.10k.sam

# To convert .sam file to a sorted .bam file
samtools sort -@ 4 -O BAM SRR072893.10k.sam > SRR072893.10k.bam

# To create an index .bai file from the .bam file
samtools index -b -@ 4 SRR072893.10k.bam

# To quantitate the sorted .bam file using StringTie
stringtie SRR072893.10k.bam -G ../genomes/BDGP6.Ensembl.81.gtf -o SRR072893.10k.gtf -p 4 -e -B


## EXERCISE 3

# The really slow way to find the number of alignments to each chromosome is:
grep "^SRR072893" SRR072893.sam | cut -f 3 | sort | uniq -c | sort -r -k 1 -g > chromosome_alignment_counts.txt

