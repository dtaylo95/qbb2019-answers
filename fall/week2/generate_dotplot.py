"""
Usage: python generate_dotplot.py <lastz_alignment_file> <output_file> <plot_title>
"""

import sys
import matplotlib.pyplot as plt

alignment_file = open(sys.argv[1])
output_file = sys.argv[2]

alignment_dict = {}

for line in alignment_file:
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	name = fields[6].split('_')[1]
	ref_start = int(fields[4])
	ref_end = int(fields[5])
	if fields[7] == '+':
		cont_start = int(fields[9])
		cont_end = int(fields[10])
	elif fields[7] == '-':
		cont_start = int(fields[10])
		cont_end = int(fields[9])
	alignment_dict.setdefault(name,[])
	alignment_dict[name].append((ref_start, ref_end, cont_start, cont_end))

sorted_contigs = sorted(alignment_dict.keys(), key=lambda x: min([y[0] for y in alignment_dict[x]]))
ref_start = min([min([y[0] for y in alignment_dict[x]]) for x in alignment_dict])
ref_end = max([max([y[1] for y in alignment_dict[x]]) for x in alignment_dict])

fig, ax = plt.subplots()

x_start = 0 
for contig in sorted_contigs:
	contig_length = max([abs(alignment[3] - alignment[2]) for alignment in alignment_dict[contig]])
	for alignment in alignment_dict[contig]:
		ax.plot([x_start+min(alignment[2:]), x_start+max(alignment[2:])], [alignment[0], alignment[1]], '-', color='black')
	x_start += contig_length


ax.set_ylim(0, 110000)
ax.set_xlim(0)
ax.set_xlabel('Assembled Contigs')
ax.set_ylabel('Reference Genome (to 110kb)')
ax.set_yticks([0,100000])
ax.set_title(sys.argv[3])
plt.tight_layout()
fig.savefig(output_file)
