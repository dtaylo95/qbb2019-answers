#!/usr/bin/env Python3

"""
Usage: ./exercise3b.py <t_name> <samples.csv> <replicates.csv> <all(replicates)FPKMs_by_SRR.csv>

Create a timecourse of a given transcript for males and females
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Specify transcript of interest
t_name = sys.argv[1]

# Load metadata

SRR_to_dev_map = {}

for i, line in enumerate(open(sys.argv[2])):
	if i == 0:
		continue
	fields = line.strip().split(',')
	SRR_to_dev_map[fields[0]] = '_'.join(fields[1:])

for i, line in enumerate(open(sys.argv[3])):
	if i == 0:
		continue
	fields = line.strip().split(',')
	SRR_to_dev_map[fields[0]] = '_'.join(fields[1:])


def sex_sorter(sex):
	# Load FPKMs
	fpkms = pd.read_csv(sys.argv[4], index_col='t_name')
	#Extract 
	my_data = []
	identifiers = []
	for srr_id in SRR_to_dev_map:
		if SRR_to_dev_map[srr_id].startswith(sex):
			identifiers.append(SRR_to_dev_map[srr_id].split('_')[1])
			my_data.append(fpkms.loc[t_name,srr_id])
	return my_data, identifiers


male_data, all_lables = sex_sorter('male')
female_data, all_lables = sex_sorter('female')
male_sample = male_data[:-4]
male_replicates = male_data[-4:]
female_sample = female_data[:-4]
female_replicates = female_data[-4:]
sample_labels = all_lables[:-4]
replicates_lables = all_lables[-4:]


fig, ax = plt.subplots()
ax.plot(male_sample, color='blue', label='Male')
ax.plot(female_sample, color='red', label='Female')
ax.plot([5,6,7,8], male_replicates, color='cyan', label='Male Replicates')
ax.plot([5,6,7,8], female_replicates, color='pink', label='Female Replicates')
plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', ncol=1)
ax.set_xticks(np.arange(len(sample_labels)))
ax.set_xticklabels(sample_labels, fontsize=7, rotation=89.5)

ax.set_title('Transcript %s Expression over fly Development' %(t_name))
ax.set_xlabel('Developmental Stage')
ax.set_ylabel('FPKM value')

plt.tight_layout()

fig.savefig('timecourse.png')
plt.close(fig)