"""
Usage: python diff_exp.py <gene_expression_by_celltypes.txt>
"""

import pandas as pd
import numpy as np
import sys
from scipy import stats


df_data = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
diff_exp_high = ((df_data['CFU'] + df_data['unk'])/2)/((df_data['poly'] + df_data['int'])/2) >= 2
diff_exp_low = ((df_data['CFU'] + df_data['unk'])/2)/((df_data['poly'] + df_data['int'])/2) <= 0.5

diff_exp_genes = df_data[diff_exp_high | diff_exp_low]

out_file = open('diff_expressed_genes.txt','w')
out_file.write('gene\tCFU\tply\tunk\tint\tmys\tmid\tpval\n')

for gene_name, row in diff_exp_genes.iterrows():
	sample1 = [row['CFU'], row['unk']]
	sample2 = [row['poly'], row['int']]
	pvalue = stats.ttest_rel(sample1, sample2).pvalue
	if pvalue < 0.05:
		out_file.write(gene_name + '\t' + '\t'.join([str(x) for x in list(df_data.loc[gene_name])]) + '\t' + str(pvalue) + '\n')

