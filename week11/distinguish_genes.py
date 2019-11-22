"""
Usage: python distinguish_genes.py <RNAseq_data.h5>
"""

import scanpy as sc
import sys
from copy import deepcopy
import matplotlib.pyplot as plt

adata = sc.read_10x_h5(sys.argv[1])
adata.var_names_make_unique()
sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.louvain(adata, resolution=0.3)

sc.tl.rank_genes_groups(adata, groupby='louvain', method='t-test')
sc.pl.rank_genes_groups(adata, groupby='louvain', method='t-test', show=False, save='_ttest.png')
# print(adata.uns['rank_genes_groups']['names'][0])


sc.tl.rank_genes_groups(adata, groupby='louvain', method='logreg')
sc.pl.rank_genes_groups(adata, groupby='louvain', method='logreg', show=False, save='_logreg.png')
# print(adata.uns['rank_genes_groups']['names'][0])
