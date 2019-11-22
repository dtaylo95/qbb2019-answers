"""
Usage: python UMAP_by_gene.py <RNAseq_data.h5> <output_all_UMAPs.png>
"""

import scanpy as sc
import sys
import matplotlib.pyplot as plt

adata = sc.read_10x_h5(sys.argv[1])
adata.var_names_make_unique()
sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.louvain(adata, resolution=0.3)

sc.tl.umap(adata)
axes = sc.pl.umap(adata, color=['louvain','Malat1', 'Ccna2', 'mt-Atp6', 'Lhx6', 'mt-Cytb', 'Zbtb20', 'Npas1', 'Ftl1', 'Trem2', 'Reln', 'Col3a1'], show=False, legend_loc=None)

plt.tight_layout()
plt.show()