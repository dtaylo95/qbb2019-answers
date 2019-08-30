#!/usr/bin/env Python3

"""
Usage: ./exercise4.py <all_FPKMs.csv> <SRR_sample1> <SRR_sample2>

Creates an MA plot for two selected data sets
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fpkms = pd.read_csv(sys.argv[1], index_col = 't_name')

transcripts_sample1 = list(np.log2(fpkms.loc[:,sys.argv[2]] + 0.01))
transcripts_sample2 = list(np.log2(fpkms.loc[:,sys.argv[3]] + 0.01))

my_data = []
for i in range(0, len(transcripts_sample1)):
	my_data.append((transcripts_sample1[i], transcripts_sample2[i]))

M_vals = []
A_vals = []

for point in my_data:
	R = point[0]
	G = point[1]
	M_vals.append(R - G)
	A_vals.append(0.5*(R +G))

fig, ax = plt.subplots()
ax.scatter(x=A_vals,y=M_vals, s=8, alpha=0.1, color='black')
ax.set_xlabel('A values')
ax.set_ylabel('M values')
fig.savefig('MA_plot.png')
plt.close(fig)