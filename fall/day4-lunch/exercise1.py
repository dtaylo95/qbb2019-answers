#!/usr/bin/env Python3

"""
Usage: ../exercise1.py <ctab> <a> <loc> <scale>

Plot FPKM values as a histogram
"""

import sys
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

ctab = pd.read_csv(sys.argv[1], sep='\t')
goi = ctab.loc[:,'FPKM'] > 0

my_data = np.log2(ctab.loc[goi,'FPKM'])

title = sys.argv[1].split(os.sep)[-2]

mu = my_data.mean()
std = my_data.std()

a = float(sys.argv[2])
loc = float(sys.argv[3])
scale = float(sys.argv[4])
x = np.linspace(my_data.min(), my_data.max(),100)
y = stats.skewnorm.pdf(x, a, loc, scale)

y2 = stats.norm.pdf(x, mu, std)

fig, ax = plt.subplots()
ax.hist(my_data, bins=100, density=True)
ax.plot(x, y, color='red', label='Skewed Normal Distribution\na = %s\nloc = %s\nscale = %s' %(str(a), str(loc), str(scale)))
ax.plot(x, y2, color='green', label='Normal Distribution\nmu = %s\nstd_dev = %s' %(str(mu)[:5], str(std)[:5]))
ax.legend(loc='upper left', fontsize=8)
ax.set_xlabel('log2(FPKM)')
ax.set_ylabel('frequency')
ax.set_title('Distribution of %s non-zero FPKM values' %(title))
fig.savefig('exercise1.png')
plt.close(fig)