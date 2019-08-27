#!/usr/bin/env Python3

"""
From the given fly.txt file containing the various accession numbers, 
parses out only the FlyBase AC's and the UniProt AC's into a tab-delimited
file with FlyBase AC's in the first column and their associated UniProt AC's
in the second column

Usage:	./parse_mapping.py <fly.txt> <out_file>
"""

import sys

map_file = open(sys.argv[1],'r')
out_file = open(sys.argv[2],'w')

start_map_test = False

for line in map_file:
	if line.startswith('_'):
		start_map_test = True
	else:
		if start_map_test:
			if line.startswith('-'):
				start_map_test = False
				continue
			fields = line.strip().split()
			if len(fields) != 4 or 'DROME' not in fields[1]:
				continue
			out_file.write('%s\t%s\n' %(fields[3], fields[2]))

map_file.close()
out_file.close()