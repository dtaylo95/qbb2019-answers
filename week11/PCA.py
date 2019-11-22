"""
Usage: python clustering.py <RNAseq_data.h5> <output_PCA.png>
"""

import scanpy as sc
import sys
from copy import deepcopy
import matplotlib.pyplot as plt


fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16,8))

adata = sc.read_10x_h5(sys.argv[1])
adata.var_names_make_unique()

unfiltered_adata = deepcopy(adata)

sc.tl.pca(unfiltered_adata, n_comps=2, return_info=True)
sc.pl.pca(unfiltered_adata, ax=ax1, show=False)
unfiltered_variance = unfiltered_adata.uns['pca']['variance_ratio']

sc.pp.recipe_zheng17(adata, n_top_genes=1000, log=True, plot=False, copy=False)

sc.tl.pca(adata, n_comps=2, return_info=True)
sc.pl.pca(adata, ax=ax2, show=False)
filtered_variance = adata.uns['pca']['variance_ratio']

ax1.set_title('Unfiltered')
ax1.set_xlabel('PC1 ({}%)'.format(str(unfiltered_variance[0]*100)[:5]))
ax1.set_ylabel('PC2 ({}%)'.format(str(unfiltered_variance[1]*100)[:5]))
ax2.set_xlabel('PC1 ({}%)'.format(str(filtered_variance[0]*100)[:5]))
ax2.set_ylabel('PC2 ({}%)'.format(str(filtered_variance[1]*100)[:5]))
ax2.set_title('Filtered')
plt.tight_layout()
fig.savefig(sys.argv[2])