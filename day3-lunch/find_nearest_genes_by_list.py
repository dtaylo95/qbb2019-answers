#!/usr/bin/env Python3


"""
Finds the closest genes to a set of positions (as an input file)

Usage: ./find_nearest_genes_by_list.py <gtf_file> <pos_file>
"""


import sys


def binary_search(sorted_list,value,iteration_count=1):
	if len(sorted_list) == 2:
		if abs(sorted_list[0]-value) <= abs(sorted_list[1] - value):
			return sorted_list[0], iteration_count 
		else:
			return sorted_list[1], iteration_count
	if len(sorted_list) == 1:
		return sorted_list[0], iteration_count

	median_i = int(len(sorted_list)/2)
	new_bin = 0
	if sorted_list[median_i] > value:
		new_bin = sorted_list[:median_i+1]
	elif sorted_list[median_i] < value:
		new_bin = sorted_list[median_i:]
	elif sorted_list[median_i] == value:
		return value, iteration_count
	return binary_search(new_bin,value,iteration_count+1)


target_locs = []
with open(sys.argv[2]) as pos_file:
	for line in pos_file:
		fields = line.strip().split()
		target_locs.append((fields[0], int(fields[1])))

gene_locs_dict = {}
gene_id_w_locs = {}
all_locs = {}

for line in open(sys.argv[1]):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	if fields[2] == 'gene' and 'gene_biotype "protein_coding"' in line:
		gene_id = fields[13][1:-2]
		chromosome = fields[0]
		gene_id_w_locs.setdefault(chromosome, {})
		gene_id_w_locs[chromosome][gene_id] = (int(fields[3]),int(fields[4]))
		gene_locs_dict.setdefault(chromosome, {})
		gene_locs_dict[chromosome][fields[3]] = gene_id
		gene_locs_dict[chromosome][fields[4]] = gene_id
		all_locs.setdefault(chromosome, [])
		all_locs[chromosome].append(int(fields[3]))
		all_locs[chromosome].append(int(fields[4]))

for locs in all_locs:
	all_locs[locs].sort()

for target_loc in target_locs:
	target_chromosome = target_loc[0]
	pos = target_loc[1]

	closest_loc, necessary_iterations = binary_search(all_locs[target_chromosome], int(pos))
	closest_loc = str(closest_loc)
	closest_gene = gene_locs_dict[target_chromosome][closest_loc]

	print('For the position: %s %s, the closest gene is %s' %(target_chromosome, pos, closest_gene))
	if int(pos) in range(gene_id_w_locs[target_chromosome][closest_gene][0],gene_id_w_locs[target_chromosome][closest_gene][1]+1):
		print('The target region is within the gene.')
	elif int(pos) < gene_id_w_locs[target_chromosome][closest_gene][0]:
		print('The target region is %d bases before the start of the gene' %(gene_id_w_locs[target_chromosome][closest_gene][0]-int(pos)))	
	elif int(pos) > gene_id_w_locs[target_chromosome][closest_gene][1]:
		print('The target region is %d bases after the end of the gene' %(int(pos)-gene_id_w_locs[target_chromosome][closest_gene][1]))	

	print('It took %d iterations to find this gene\n' %(necessary_iterations))
