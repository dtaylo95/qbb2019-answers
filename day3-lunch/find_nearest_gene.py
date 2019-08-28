#!/usr/bin/env Python3

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


##########################################################################################################

gtf = open(sys.argv[1])
pos = sys.argv[2]
iterations = 1

if len(sys.argv) > 3:
	iterations = int(sys.argv[3])

gene_locs_dict = {}
gene_id_w_locs = {}
all_locs = []

for line in gtf:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	if (fields[2] == 'gene' and fields[0] == '3R') and 'gene_biotype "protein_coding"' in line:
		gene_id = fields[13][1:-2]
		gene_id_w_locs[gene_id] = (int(fields[3]),int(fields[4]))
		gene_locs_dict[fields[3]] = gene_id
		gene_locs_dict[fields[4]] = gene_id
		all_locs.append(int(fields[3]))
		all_locs.append(int(fields[4]))

all_locs.sort()

for i in range(1,iterations+1):

	closest_loc, necessary_iterations = binary_search(all_locs, int(pos))
	closest_loc = str(closest_loc)
	closest_gene = gene_locs_dict[closest_loc]

	print('\nThe #%s closest protein-coding gene is %s.' %(str(i),closest_gene))

	if int(pos) in range(gene_id_w_locs[closest_gene][0],gene_id_w_locs[closest_gene][1]+1):
		print('The target region is within the gene.')
	elif int(pos) < gene_id_w_locs[closest_gene][0]:
		print('The target region is %d bases before the start of the gene' %(gene_id_w_locs[closest_gene][0]-int(pos)))	
	elif int(pos) > gene_id_w_locs[closest_gene][1]:
		print('The target region is %d bases after the end of the gene' %(int(pos)-gene_id_w_locs[closest_gene][1]))	

	print('It took %d iterations to find this gene\n' %(necessary_iterations))
	all_locs.remove(gene_id_w_locs[closest_gene][0])
	all_locs.remove(gene_id_w_locs[closest_gene][1])