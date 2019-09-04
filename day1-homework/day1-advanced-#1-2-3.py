#!/usr/bin/env Python3

import sys
import re

def parse_CIGAR(c_string):
	numerals = '1234567890'
	parsed_CIGAR = []
	count = ''
	identity = ''
	for char in c_string:
		if char in numerals:
			count += char
		else:
			identity = char
			parsed_CIGAR.append(count)
			parsed_CIGAR.append(identity)
			count = ''
			identity = ''
	return parsed_CIGAR


if len(sys.argv) > 1:
	f = open(sys.argv[1])
else:
	f = sys.stdin

reverse_comp_count = 0
bit_mask = 0b10000
over30score_count = 0
indel1 = 0
indel2 = 0
indel3 = 0
indel4 = 0
indelbig = 0

for line in f:
	fields = line.rstrip('\n').split('\t')
	flag = int(fields[1])
	if flag & bit_mask:
		reverse_comp_count += 1
	
	qual_score = fields[10]
	chars = 0
	total_scores = 0
	for char in qual_score:
		chars += 1
		total_scores += ord(char) - 33
	avg_score = float(total_scores)/chars
	if avg_score > 30:
		over30score_count += 1

	CIGAR_string = fields[5]
	parsed_CIGAR = parse_CIGAR(CIGAR_string)
	if 'I' in parsed_CIGAR:
		indel_score = int(parsed_CIGAR[parsed_CIGAR.index('I') - 1])
		if indel_score == 1:
			indel1 += 1
		elif indel_score == 2:
			indel2 += 1
		elif indel_score == 3:
			indel3 += 1
		elif indel_score == 4:
			indel4 += 1
		elif indel_score > 4:
			indelbig += 1
	if 'D' in parsed_CIGAR:
		indel_score = int(parsed_CIGAR[parsed_CIGAR.index('D') - 1])
		if indel_score == 1:
			indel1 += 1
		elif indel_score == 2:
			indel2 += 1
		elif indel_score == 3:
			indel3 += 1
		elif indel_score == 4:
			indel4 += 1
		elif indel_score > 4:
			indelbig += 1		


print('Advanced Exercise #1\n', reverse_comp_count)
print('Advanced Exercise #2\n', over30score_count)
print('Advanced Exercise #3\n', indel1, indel2, indel3, indel4, indelbig)