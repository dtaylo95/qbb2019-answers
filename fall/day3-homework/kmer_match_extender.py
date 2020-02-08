#!/usr/bin/env Python3

"""
Finds longest exact matches of a query sequence against a set of target sequences.

Usage: ./kmer_matcher.py <target.fa> <query.fa> <k>

Commands:
	k		The length of the seeds used to find exact matches
"""

# ______________________________________________________________________________
# This section is copied from the kmer_matcher. Finds all matches

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

matches = []

for query_kmer in query_kmers:
	for seq_id in target_kmers:
		if query_kmer in target_kmers[seq_id]:
			for query_pos in query_kmers[query_kmer]:
				for target_pos in target_kmers[seq_id][query_kmer]:
					matches.append((seq_id, int(target_pos), int(query_pos)))

# ______________________________________________________________________________


# This dictionary will store tuples containing query_start, query_end of each match for each target_id
longest_matches = {}

while matches != []:
	match_info = matches[0]
	matches.remove(match_info)
	seq_id = match_info[0]
	target_pos = match_info[1]
	query_pos = match_info[2]
	min_target_pos = target_pos
	min_query_pos = query_pos
	max_query_pos = query_pos
	i = 1

	# This while loop will incrementally search for tuples from the same target seq whose 
	# target position and query position are one higher than previously
	while True:
		if (seq_id, target_pos + i, query_pos + i) in matches:
			matches.remove((seq_id, target_pos + i, query_pos + i))
			max_query_pos += 1
			i += 1
		else:
			break
	
	# This while loop will do the same, except searching backwards
	j = 1
	while True:
		if (seq_id, target_pos - j, query_pos - j) in matches:
			matches.remove((seq_id, target_pos - j, query_pos - j))
			min_query_pos -= 1
			min_target_pos -= 1
			j += 1
		else:
			break

	# Once the longest possible match is found at these locations in the target and query,
	# the start and end positions of the match in the query sequence, as well as the start
	# position of the match in the target sequence will be added (as a tuple) to a list
	# of matches in a dictionary organized by target_id
	longest_matches.setdefault(seq_id,[])
	longest_matches[seq_id].append((min_query_pos, max_query_pos+k, min_target_pos))


# This section parses through each target_seq and incrementally finds the match with the
# longest length, then prints the target_id, target start position, query start position,
# and the sequence to the console
print('target_id\ttarget_start_pos\tquery_start_pos\tmatch_seq')
for seq_id in longest_matches:
	while longest_matches[seq_id] != []:
		query_start = 0
		query_end = 0
		target_start = 0
		for match_len in longest_matches[seq_id]:
			if match_len[1] - match_len[0] > query_end - query_start:
				query_start = match_len[0]
				query_end = match_len[1]
				target_start = match_len[2]
		print(seq_id + '\t' + str(target_start) + '\t' + str(query_start) + '\t' + query_seq[query_start:query_end])
		longest_matches[seq_id].remove((query_start, query_end, target_start))










