"""
Usage: python find_similar.py <gene_expression_by_celltypes.txt> <goi>
"""

import pandas as pd
import sys
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

df_data = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
rel_data = df_data[['CFU','poly']]

kmeans = KMeans(n_clusters=7).fit(rel_data)
labels = list(kmeans.labels_)
genes = list(rel_data.index.values)

goi_index = genes.index(sys.argv[2])
goi_cluster = labels[goi_index]

related_genes = []
for i, gene in enumerate(genes):
	if labels[i] == goi_cluster:
		related_genes.append(gene)

out_file = open('clustered_diff_exp_genes.txt','w')
for gene_name in related_genes:
	out_file.write(gene_name + '\n')