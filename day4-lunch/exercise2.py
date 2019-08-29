#!/usr/bin/env Python3

"""
Usage: ../exercise1.py <ctab> <a> <loc> <scale>

Plot FPKM values as a histogram
"""

import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

data_name1 = sys.argv[1].split(os.sep)[-2]
data_name2 = sys.argv[2].split(os.sep)[-2]
data_ctab1 = pd.read_csv(sys.argv[1], sep='\t', index_col='t_name')
data_ctab2 = pd.read_csv(sys.argv[2], sep='\t', index_col='t_name')

FPKMs_data1 = np.log2(data_ctab1.loc[:,'FPKM'] + 0.01) 
FPKMs_data2 = np.log2(data_ctab2.loc[:,'FPKM'] + 0.01)

m, b = np.polyfit(FPKMs_data1, FPKMs_data2, 1)

x = [FPKMs_data1.min(), FPKMs_data1.max()]

fig, ax = plt.subplots()
ax.scatter(x=FPKMs_data1, y=FPKMs_data2, s=3, alpha=0.1)
ax.plot(x, [(m*x1)+b for x1 in x], color='red', label='y = %sx + %s' %(str(m)[:5], str(b)[:5]))
ax.legend(loc='upper left', fontsize=8)
ax.set_xlabel('%s log2(FPKM)' %(data_name1))
ax.set_ylabel('%s log2(FPKM)' %(data_name2))
ax.set_title('Transcript FPKM distribution of %s vs %s' %(data_name1, data_name2))
fig.savefig('%s_FPKM-v-%s_FPKM.png' %(data_name1, data_name2))
plt.close(fig)