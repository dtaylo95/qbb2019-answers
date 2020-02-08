"""
python parse_vcf.pg <variants.vcf>
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

vcf_file = open(sys.argv[1])

qualities = []
depths = []
variant_types = {}
allele_freqs = []
allele_nums = []

for line in vcf_file:
	if line.startswith('#'):
		continue
	fields = line.strip().split()

	quality = float(fields[5])
	qualities.append(quality)

	info = fields[7].split(';')

	allele_num = info[4][3:].split(',')
	for num in allele_num:
		allele_nums.append(int(num))

	read_depth = info[7][3:].split(',')
	for depth in read_depth:
		depths.append(int(depth))

	allele_freq = info[3][3:].split(',')
	for freq in allele_freq:
		allele_freqs.append(float(freq))

	annotations = info[-1].split('|')
	for i in range(len(annotations)):
		if annotations[i].startswith('ANN') or annotations[i].startswith(','):
			variant_type = annotations[i+1].split('&')
			for var_type in variant_type:
				if var_type != '':
					variant_types.setdefault(var_type, 0)
					variant_types[var_type] += 1

allele_counts = [round(allele_freqs[i]*allele_nums[i]) for i in range(len(allele_freqs))]

stuff = sorted([(x[1],x[0]) for x in variant_types.items()])
stuff = stuff[::-1]

variant_type_labels = [' '.join(item[1].split('_')) for item in stuff]
variant_type_counts = [item[0] for item in stuff]


fig, axs = plt.subplots(nrows=2,ncols=2, figsize=(10,10))
plt.title("Genotype Variation in Saccharomyces cerevisiae")

axs[0,0].hist(qualities, bins=100, range=(0,2000), edgecolor='#E6E6E6')
axs[0,0].set_xlabel('Genotype Quality')
axs[0,0].set_ylabel('Frequency of SNPs')
axs[0,0].set_facecolor('#E6E6E6')
axs[0,0].set_title('SNP Genotype Quality Distribution')
axs[1,0].hist(depths, bins=50, range=(0,150), color='green', edgecolor='#E6E6E6')
axs[1,0].set_xlabel('Read Depth')
axs[1,0].set_ylabel('Frequency of SNPs')
axs[1,0].set_facecolor('#E6E6E6')
axs[1,0].set_title('SNP Read Depth Distribution')
axs[0,1].hist(allele_counts, color='red', edgecolor='#E6E6E6')
axs[0,1].set_xlabel('Population Allele Frequency')
axs[0,1].set_ylabel('Frequency of SNPs')
axs[0,1].set_facecolor('#E6E6E6')
axs[0,1].set_title('SNP Allele Frequency Spectrum')
axs[1,1].bar(range(len(variant_type_labels)), variant_type_counts, color='yellow')
axs[1,1].set_xticks(range(len(variant_type_labels)))
axs[1,1].set_xticklabels(variant_type_labels, rotation=90)
axs[1,1].set_yscale('log')
axs[1,1].set_ylabel('Frequency of SNPs')
axs[1,1].set_facecolor('#E6E6E6')
axs[1,1].set_title('Predicted SNP Effects')


plt.tight_layout()
fig.savefig('variant_analysis.png')