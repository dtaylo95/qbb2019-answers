"""
Usage: python UMAP.py <RNAseq_data.h5> <output_UMAP.png>
"""

import scanpy as sc
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import seaborn as sns

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16,8))

adata = sc.read_10x_h5(sys.argv[1])
adata.var_names_make_unique()
sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.louvain(adata, resolution=0.3)
colors = adata.obs['louvain']
print(colors)

sc.tl.tsne(adata)
sc.pl.tsne(adata, color='louvain', ax=ax1, show=False, legend_loc='on data')
plt.tight_layout()

sc.tl.umap(adata)
sc.pl.umap(adata, color='louvain', ax=ax2, show=False, legend_loc='on data')

ax1.set_title('tSNE Analysis')
ax2.set_title('UMAP Analysis')
fig.savefig('tSNE_UMAP_filtered_data.png')