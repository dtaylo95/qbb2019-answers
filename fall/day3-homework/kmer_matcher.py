#!/usr/bin/env Python3

"""
Finds all exact matches of a certain size of a query sequence against a set of target sequences.

Usage: ./kmer_matcher.py <target.fa> <query.fa> <k>

Commands:
	k		The length of the seeds used to find exact matches
"""

import sys
from fasta import FASTAReader

k = int(sys.argv[3])

target_kmers = {}


for ident, sequence in FASTAReader(open(sys.argv[1])):
	target_kmers[ident] = {}
	sequence = sequence.upper()
	for i in range(0, len(sequence) - k + 1):
		kmer = sequence[i:i+k]
		target_kmers[ident].setdefault(kmer,set())
		target_kmers[ident][kmer].add(str(i))

query_file = open(sys.argv[2])
query_seq = ''
for line in query_file:
	if line.startswith('>'):
		continue
	query_seq += line.strip()
 
query_kmers = {}
for i in range(0, len(query_seq) - k + 1):
	kmer = query_seq[i:i+k]
	query_kmers.setdefault(kmer,set())
	query_kmers[kmer].add(str(i))

for query_kmer in query_kmers:
	for seq_id in target_kmers:
		if query_kmer in target_kmers[seq_id]:
			prtln = seq_id + '\t' + ','.join(list(target_kmers[seq_id][query_kmer])) + '\t' + ','.join(list(query_kmers[query_kmer])) + '\t' + query_kmer
			print(prtln)
