"""
Usage: python plot_PCA.py <FILE_NAME.eigenvec> <out_file.png>
"""

import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

pca_data = open(sys.argv[1])

all_data = []

families = set()

for line in pca_data:
	fields = line.strip().split()
	family = fields[0]
	indiv_id = fields[1]
	pca1 = float(fields[2])
	pca2 = float(fields[3])
	all_data.append((family, pca1, pca2))
	families.add(family)

colors = ['#800000', '#9A6324', '#e6194B', '#808000', '#ffe119', '#469990', '#000075', '#000000', '#f032e6', '#aaffc3', '#a9a9a9']
families = list(families)

fam_color_map = {}
for i in range(len(families)):
	fam_color_map[families[i]] = colors[i]

patches = [mpatches.Patch(color=color, label=fam) for fam,color in fam_color_map.items()]


fig, ax = plt.subplots()
for data_point in all_data:
	ax.scatter(data_point[1], data_point[2], color=fam_color_map[data_point[0]])

ax.set_xlabel('PCA 1')
ax.set_ylabel('PCA 2')
ax.legend(handles=patches, loc='center left', bbox_to_anchor=(1.0,0.5))
ax.set_title('2D Principle Component Analysis')
plt.tight_layout()
fig.savefig(sys.argv[2])
