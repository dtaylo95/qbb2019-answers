"""
Usage: python k_means.py <gene_expression_by_celltypes.txt>
"""

import pandas as pd
import sys
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def find_closest(point, centroids):
	closest_centroid = None
	closest_distance = 0
	for centroid in centroids:
		distance = (((centroid[0] - point[0])**2) + ((centroid[1] - point[1])**2))**0.5
		if closest_distance == 0 or distance < closest_distance:
			closest_centroid = centroid
			closest_distance = distance
	return closest_centroid

df_data = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
rel_data = df_data[['CFU','poly']]
# rel_data = np.array(rel_data)

# whitened = whiten(rel_data)
# codebook, distortion = kmeans(whitened, 5)
# codebook = [tuple(x) for x in codebook]

# colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'pink']

# cmap = {}
# for i in range(len(codebook)):
# 	cmap[codebook[i]] = colors[i]

# fig, ax = plt.subplots()
# for point in rel_data:
# 	closest_centroid = find_closest(point, codebook)
# 	ax.scatter(point[0], point[1], c=cmap[closest_centroid])

# plt.show()

kmeans = KMeans(n_clusters=7).fit(rel_data)
centroids = kmeans.cluster_centers_

fig, ax = plt.subplots()
ax.scatter(rel_data['CFU'], rel_data['poly'], c=kmeans.labels_.astype(float), alpha=0.5)
ax.scatter(centroids[:, 0], centroids[:, 1], c='red', marker='x')
ax.set_xlabel('CFU Expression (FPKM)')
ax.set_ylabel('Poly Expression (FPKM)')
plt.tight_layout()
fig.savefig('kmeans_clustering.png')
plt.close()