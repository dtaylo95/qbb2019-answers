"""
Usage: python motif_location_density.py <fimo.gff> <fimo_ref_seqs.bed> <out.png>
"""

import sys
import matplotlib.pyplot as plt
import seaborn as sns


motif_locs = {}
for line in open(sys.argv[1]):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample = int(fields[0])
	start = int(fields[3])
	end = int(fields[4])
	strand = fields[6]
	middle = (start + end)/2
	motif_locs.setdefault(sample, [])
	motif_locs[sample].append((middle,strand))


rel_locs = []

for i, line in enumerate(open(sys.argv[2])):
	if i not in motif_locs:
		continue
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	seq_len = int(fields[2]) - int(fields[1])
	for hit in motif_locs[i]:
		if hit[1] == '+':
			rel_locs.append(hit[0]/seq_len)
		else:
			rel_locs.append(1-(hit[0]/seq_len))

fig, ax = plt.subplots()
sns.distplot(rel_locs, bins=20, ax=ax)
ax.set_xlabel('Relative motif location within ER4 binding site')
ax.set_ylabel('Probability Density')
ax.set_title('Distribution of relative ER4 motif location')
fig.savefig(sys.argv[3])
