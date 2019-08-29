#!/usr/bin/env Python3

"""
Usage: ./exercise3.py <threshold> <A|O> <ctab_file1> <ctab_file2> ... <ctab_file3>

Based on a threshold FPKM value, returns the set of transcripts whose FPKM value exceeds
that value in one or all input samples.
"""

import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

samples = {}

for i in range(3, len(sys.argv)):
	sample_name = sys.argv[i].split(os.sep)[-2]
	ctab = pd.read_csv(sys.argv[i], sep='\t', index_col='t_name')
	samples[sample_name] = ctab.loc[:,'FPKM']

new_samples = pd.DataFrame(samples)

if sys.argv[2] == 'A':
	sum_test = new_samples.sum(axis=1) > float(sys.argv[1])
	threshold_covered_transcripts = list(new_samples.loc[sum_test,:].index)


elif sys.argv[2] == 'O':
	any_test = new_samples.max(axis=1) > float(sys.argv[1])
	threshold_covered_transcripts = list(new_samples.loc[any_test, :].index)

print(len(threshold_covered_transcripts))