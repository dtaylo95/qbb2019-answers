"""
Usage: python cluster_genes.py <gene_expression_by_celltypes.txt>

"""

import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hac
import sys
import matplotlib.pyplot as plt
import seaborn as sns

df_data = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
all_data = np.array(df_data)

transposed_data = all_data.transpose()

cell_linkage_matrix = hac.linkage(transposed_data, method='average')
sorted_indices = hac.leaves_list(cell_linkage_matrix)
print(df_data.columns)
print(sorted_indices)

cell_types = df_data.columns

sorted_cell_types = []
sorted_data = []
for index in sorted_indices:
	sorted_data.append(transposed_data[index])
	sorted_cell_types.append(cell_types[index])
sorted_data = np.array(sorted_data)
sorted_data = sorted_data.transpose()

gene_linkage_matrix = hac.linkage(sorted_data, method='average')
double_sorted_indices = hac.leaves_list(gene_linkage_matrix)

double_sorted_data = []
for index in double_sorted_indices:
	double_sorted_data.append(sorted_data[index])

double_sorted_data = np.array(double_sorted_data)

fig, ax = plt.subplots( figsize=(4,8))
sns.heatmap(double_sorted_data, cbar_kws={'label': 'FPKM'})
ax.set_xticks([0.5,1.5,2.5,3.5,4.5,5.5])
ax.set_xticklabels(sorted_cell_types)
ax.set_yticks([])
ax.set_xlabel('Cell Type')
ax.set_ylabel('Genes')
fig.savefig('hema_heatmap.png')
plt.close()

fig, ax = plt.subplots()
hac.dendrogram(cell_linkage_matrix, labels=sorted_cell_types)
# ax.set_xticklabels(sorted_cell_types)
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig('hema_dendrogram.png')

