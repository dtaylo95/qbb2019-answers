#!/usr/bin/env Python3

"""
Converts gene_id's in a t_data.ctab file from FlyBase to UniProt

Usage:	./id_mapping.py <t_data_file> <map_file> <out_file> [i|d]

Commands:
[i|d]		Determines protocol for genes without a valid FlyBase AC.
			"i" will ignore genes without a valid FlyBase AC in the output file.
			"d" will replace the gene_id with a default '*' value.
			(Defaults to 'd')
"""

import sys

t_data = open(sys.argv[1])
mapping = open(sys.argv[2])
out_file = open(sys.argv[3],'w')
no_map_test = 'd'
if len(sys.argv) == 5:
	no_map_test = sys.argv[4]

gene_id_mapping = {}

for line in mapping:
	line = line.strip().split()
	flybase = line[0]
	uniprot = line[1]
	gene_id_mapping[flybase] = uniprot

for i, line in enumerate(t_data):
	if i == 0:
		out_file.write(line)
		continue
	fields = line.strip().split()
	fb_id = fields[8]
	if fb_id in gene_id_mapping:
		fields[8] = gene_id_mapping[fb_id]
		new_line = '\t'.join(fields) + '\n'
		out_file.write(new_line)
	else:
		if no_map_test.lower() == 'i':
			continue
		elif no_map_test.lower() == 'd':
			fields[8] = '*'
			new_line = '\t'.join(fields) + '\n'
			out_file.write(new_line)

t_data.close()
mapping.close()
out_file.close()
