"""
Usage: python gen_ManhattanPlot.py <Phenotype_ID_file>
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

###################################################################################################

colors = ['#42d4f4', '#911eb4']
sig_color = '#bfef45'

desired_pval = 0.01

###################################################################################################

phenotype_file = open(sys.argv[1])
phenotypes = []
for i, line in enumerate(phenotype_file):
	if i == 0:
		phenotypes = line.strip().split()
		break

snp_count = 0
for file_name in os.listdir(os.getcwd()):
	if file_name.endswith('.qassoc'):
		for i, line in enumerate(open(file_name)):
			if i != 0:
				snp_count +=1
	break

num_tests = len(phenotypes)*snp_count
corrected_pval = float(desired_pval)/snp_count
minuslog10threshold = -1*np.log10(corrected_pval)

for file_name in os.listdir(os.getcwd()):
	if file_name.endswith('.qassoc'):
		phenotype_id = int(file_name.split('.')[1][1:])
		pvals_by_chr = {}
		qassoc_file = open(file_name)
		
		for i, line in enumerate(qassoc_file):
			if i == 0:
				continue
			fields = line.strip().split()
			chromosome = fields[0]
			if chromosome == '26':
				chromosome = 'chrM'
			if chromosome == '23':
				chromosome = 'chrX'
			p_val = fields[-1]
			if p_val == 'NA':
				continue
			pvals_by_chr.setdefault(chromosome,[])
			pvals_by_chr[chromosome].append(-1*np.log10(float(p_val)))

		xticks = []
		xtick_labels = []
		fig, ax = plt.subplots()
		prev_points = 0
		# ax.axhline(y=minuslog10threshold,color='red', zorder=0)
		for i, chromosome in enumerate(sorted(pvals_by_chr.keys())):
			x_vals = [x+prev_points for x in range(len(pvals_by_chr[chromosome]))]
			y_vals = pvals_by_chr[chromosome]
			base_color = colors[i%2]
			marker_colors = [base_color if y_val < minuslog10threshold else sig_color for y_val in y_vals]
			ax.scatter(x_vals, y_vals, color=marker_colors, s=2, zorder=10)
			xticks.append((x_vals[-1] + x_vals[0])/2)
			prev_points += len(pvals_by_chr[chromosome])
			xtick_labels.append(chromosome)
		ax.set_xticks(xticks)
		ax.set_xticklabels(xtick_labels,rotation=90,size=7)
		ax.set_xlabel('Chromosomes')
		ax.set_ylabel('-log10(pval)')
		plt.title(phenotypes[phenotype_id - 1])
		plt.tight_layout()
		fig.savefig(file_name + '.png')
		plt.close()
