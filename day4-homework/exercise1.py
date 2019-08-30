#!/usr/bin/env Python3 

"""
Usage: ./exercise1.py <metadata.csv> <ctab_dir>

	<ctab_dir> e.g. ~/qbb2019-answers/results/stringtie

Create all.csv with FPKMs

	t_name, gene_name, sample1, ..., samplen
"""

import sys
import pandas as pd
import os

metadata = sys.argv[1]
ctab_dir = sys.argv[2]

fpkms = {}

for i, line in enumerate(open(metadata)):
	if i == 0:
		continue
	fields = line.strip().split(',')
	srr_id = fields[0]
	info = '_'.join(fields[1:])
	ctab_path = os.path.join(ctab_dir, srr_id, 't_data.ctab')

	df = pd.read_csv(ctab_path, sep='\t', index_col='t_name')
	fpkms['gene_name'] = df.loc[:,'gene_name']
	fpkms[info] = df.loc[:,'FPKM']

df_fpkms = pd.DataFrame(fpkms)
pd.DataFrame.to_csv(df_fpkms, 'all.csv')