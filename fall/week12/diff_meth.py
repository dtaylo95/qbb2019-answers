"""
Usage: python diff_meth.py <sample1_CpG_context.txt> <sample2_CpG_context.txt>
"""

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

sample1_methylation  = {}
sample2_methylation = {}

sample1_up_meth = []
sample2_up_meth = []

for i, line in enumerate(open(file1)):
	if i == 0:
		continue
	fields = line.strip().split()
	site = int(fields[3])
	methylation = fields[4]
	if methylation == 'Z':
		methylation = 1
	else:
		methylation = 0
	sample1_methylation[site] = methylation

for i, line in enumerate(open(file2)):
	if i == 0:
		continue
	fields = line.strip().split()
	site = int(fields[3])
	methylation = fields[4]
	if methylation == 'Z':
		methylation = 1
	else:
		methylation = 0
	sample2_methylation[site] = methylation

for site in sample1_methylation:
	if site in sample2_methylation:
		if sample1_methylation[site] == 1 and sample2_methylation[site] == 0:
			sample1_up_meth.append(site)
		elif sample1_methylation[site] == 0 and sample2_methylation[site] == 1:
			sample2_up_meth.append(site)

sample1_name = file1[12:-21]
sample2_name = file2[12:-21]

all_sites = sorted(sample1_up_meth + sample2_up_meth)

# sample1_out = open(sample1_name + '_up-methylated_sites.txt','w')
# for site in sample1_up_meth:
# 	sample1_out.write(str(site) + '\n')
# sample1_out.close()

# sample2_out = open(sample2_name + '_up-methylated_sites.txt','w')
# for site in sample2_up_meth:
# 	sample2_out.write(str(site)+ '\n')
# sample2_out.close()

out_file = open('differential_methylation.tsv','w')
out_file.write('Chromosome\tPosition\t{}_meth\t{}_meth\n'.format(sample1_name, sample2_name))
for start in all_sites:
	if start in sample1_up_meth:
		out_file.write('chr19\t{}\ty\tn\n'.format(str(start)))
	else:
		out_file.write('chr19\t{}\tn\ty\n'.format(str(start)))

out_file.close()

