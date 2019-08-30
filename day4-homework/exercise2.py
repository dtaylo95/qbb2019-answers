#!/usr/bin/env Python3

"""Usage: ./exercise2.py <gene_name> <FPKMs.csv>

Boxplot all transcripts for a given gene
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


gene_name = sys.argv[1]
fpkm_file = sys.argv[2]

df = pd.read_csv(fpkm_file, sep=',', index_col='t_name')

goi = df.loc[:,"gene_name"] == gene_name

fpkms = df.drop(columns='gene_name')
filtered_fpkms = fpkms.loc[goi,:]

male_ids = [y for y in list(filtered_fpkms) if y.startswith('male')]
male_transcripts = np.log2(filtered_fpkms[male_ids] + 0.001)

female_ids = [x for x in list(filtered_fpkms) if x.startswith('female')]
female_transcripts = np.log2(filtered_fpkms[female_ids] + 0.001)

male_labels = [y[5:] for y in male_ids]
female_labels = [x[7:] for x in female_ids]

fig, (ax1, ax2) = plt.subplots(nrows=2)
ax1.boxplot(male_transcripts.T)
ax1.set_xticklabels(male_labels, fontsize=7)
ax1.set_title('Male Sxl Time Series')
ax1.set_xlabel('Developmental Stage')
ax1.set_ylabel('FPKM value of transcripts')
ax2.boxplot(female_transcripts.T)
ax2.set_xticklabels(female_labels, fontsize=7)
ax2.set_title('Female Sxl Time Series')
ax2.set_xlabel('Developmental Stage')
ax2.set_ylabel('FPKM value of transcripts')
plt.tight_layout()
fig.savefig('boxplot.png')
plt.close()
