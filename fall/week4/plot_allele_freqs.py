"""
Usage: python plot_allele_freqs.py <FILE_NAME.vcf> <out_file.png>
"""

import sys
import matplotlib.pyplot as plt

vcf_file = open(sys.argv[1])

allele_freqs = []

for line in vcf_file:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	allele_freq = fields[7][3:].split(',')
	for freq in allele_freq:
		allele_freqs.append(float(freq))

fig, ax = plt.subplots()
plt.hist(allele_freqs, bins=200)
ax.set_xlabel('Allele Frequency')
ax.set_ylabel('# Variants')
ax.set_title('Variant Allele Frequency Spectrum')
plt.tight_layout()
fig.savefig(sys.argv[2])