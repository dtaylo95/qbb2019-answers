#!/usr/bin/env Python3

"""
Usage: ./lin_reg_resids.py <ctab> [mark_FILE.tab ...] > [out_stats]
"""

import sys
import pandas as pd
import os
import statsmodels.api as sm
import matplotlib.pyplot as plt
import numpy as np

fpkms = pd.DataFrame(pd.read_csv(sys.argv[1], sep = '\t', index_col='t_name' )['FPKM'])

marks = []

for i in range(2, len(sys.argv)):
	mark_id = os.path.split(sys.argv[i])[-1][:-4]
	marks.append(mark_id)
	mark_df = pd.DataFrame(pd.read_csv(sys.argv[i], sep='\t', header=None, names=['t_name','size','covered','sum','mean0','mean']).loc[:,['t_name','mean']])
	mark_df.set_index('t_name', inplace=True)
	fpkms[mark_id] = mark_df['mean']

model = sm.formula.ols(formula='FPKM ~ %s' %(' + '.join(marks)), data=fpkms)
ols_results = model.fit()

residuals = ols_results.resid

fig, ax = plt.subplots()

ax.hist(residuals, bins=1000, range=(-100,100))
ax.set_xlim((-100,100))
ax.set_xlabel('Residual value')
ax.set_ylabel('Counts')
ax.set_title('Distribution of Regression Residuals')
fig.savefig('resid_hist.png')
plt.close(fig)
