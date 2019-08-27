# EXERCISE 1

### To create a fastq file that only contains the first 10,000 reads from SRR072893.fastq
head -40000 ../rawdata/SRR072893.fastq > SRR072893.10k.fastq

### To verify that the resulting file has 10,000 reads
grep "^@SRR072893" SRR072893.10k.fastq | wc -l

### To create a quality control report for the truncated .fq file
fastqc -t 4 SRR072893.10k.fastq

### To align the reads against the generated indices
hisat2 -x ../genomes/BDGP6 -U SRR072893.10k.fastq -p 4 -S SRR072893.10k.sam

### To convert .sam file to a sorted .bam file
samtools sort -@ 4 -O BAM SRR072893.10k.sam > SRR072893.10k.bam

### To create an index .bai file from the .bam file
samtools index -b -@ 4 SRR072893.10k.bam

### To quantitate the sorted .bam file using StringTie
stringtie SRR072893.10k.bam -G ../genomes/BDGP6.Ensembl.81.gtf -o SRR072893.10k.gtf -p 4 -e -B


# EXERCISE 3

### The really slow way to find the number of alignments to each chromosome is:
grep "^SRR072893" SRR072893.sam | cut -f 3 | sort | uniq -c | sort -r -k 1 -g > chromosome_alignment_counts.txt


# EXERCISE 4

### The differences between lines of different column-length are as follows:
All lines have the required 11 fields:
1. Name of query sequence
2. Bitwise flags
3. Name of aligned reference (or * if not available)
4. Start position of aligned read in reference (0 if not aligned)
5. MAPQ mapping quality (0 if not aligned)
6. CIGAR string( * if not available)
7. RNEXT ( * if not available)
8. PNEXT(0 if not available)
9. TLEN (0 for single-segment template)
10. Sequence of aligned read
11. ASCII quality (base-by-base) of aligned read

They each also have the YT optional tag

Lines with 13 columns have an additional YF tag

Lines with 20 or more columns have the following tags:
1. AS:i:<alignment_score>
