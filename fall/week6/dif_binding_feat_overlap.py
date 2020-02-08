"""
Usage: python dif_binding_feat_overlap.py <dif_binding_ER4.bed> <dif_binding_G1E.bed> <features_ER4.bed> <features_G1E.bed>
"""

import sys
import matplotlib.pyplot as plt
import numpy as np


ER4_dif_binding = open(sys.argv[1])
G1E_dif_binding = open(sys.argv[2])

pre_dif_CTCF_counts = 0
post_dif_CTCF_counts = 0

for line in ER4_dif_binding:
	fields = line.strip().split()
	if len(fields) == 6:
		post_dif_CTCF_counts += 1

for line in G1E_dif_binding:
	fields = line.strip().split()
	if len(fields) == 6:
		pre_dif_CTCF_counts +=1


ER4_features = open(sys.argv[3])
G1E_features = open(sys.argv[4])

pre_dif_features = {}
post_dif_features = {}
for line in ER4_features:
	fields = line.strip().split()
	feature = fields[3]
	post_dif_features.setdefault(feature,0)
	post_dif_features[feature] += 1
for line in G1E_features:
	fields = line.strip().split()
	feature = fields[3]
	pre_dif_features.setdefault(feature,0)
	pre_dif_features[feature] += 1


plt.style.use('dark_background')

fig, (ax1, ax2) = plt.subplots(ncols=2,figsize=(8,6))

ax2.bar(x=['Lost during\nDifferentiation','Gained during\nDifferentiation'], height=[pre_dif_CTCF_counts,post_dif_CTCF_counts], color=['#1f78b3','#fe7f02'])
ax2.set_ylabel('# of CTCF sites')
ax2.set_title('Differential CTCF Binding\nin Erythroid Differentiation')

labels = list(pre_dif_features.keys())
lost_sites = [pre_dif_features[label] for label in labels]
gained_sites = [post_dif_features[label] for label in labels]

x = np.arange(len(labels))
width = 0.3

ax1.bar(x=x - width/2 - 0.05, height=lost_sites, width=width, color='#1f78b3', label='Pre-differentiation')
ax1.bar(x=x + width/2 + 0.05, height=gained_sites, width=width, color='#fe7f02', label = 'Post-differentiation')
ax1.set_xticks(x)
ax1.set_xticklabels(labels)
ax1.set_ylabel('# of CTCF sites')
ax1.legend()
ax1.set_title('Differential Feature Binding of\nCTCF in Erythroid Differentiation')

plt.tight_layout()
fig.savefig('Differential_CTCF_binding.png')