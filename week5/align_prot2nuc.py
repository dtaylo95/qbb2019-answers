"""
Usage: python align_prot2nuc.py <mafft_alignment_file.out> <nuc_alignment_file.fa>
"""

import sys
import statistics as stats
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

class FASTAReader(object):

	def __init__(self, fh):
		self.fh = fh
		self.last_ident = None
		self.eof = False

	def next(self):
		if self.eof:
			return None, None
		elif self.last_ident is None:
			line = self.fh.readline()
			assert line.startswith('>'), "Not a FASTA file"
			ident = line[1:].rstrip('\n')
		else:
			ident = self.last_ident
		sequence = ''
		while True:
			line = self.fh.readline()
			if line == '':
				self.eof = True
				break
			elif line.startswith('>'):
				self.last_ident = line[1:].rstrip('\n')
				break
			else:
				sequence += line.strip()
		return ident, sequence

table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    	}


protein_alignment = open(sys.argv[1],'r')
nucleotide_hits = open(sys.argv[2],'r')

protein_reader = FASTAReader(protein_alignment)
nucleotide_reader = FASTAReader(nucleotide_hits)

prot_seqs = {}
nucl_seqs = {}

while True:
	prot_id, prot_seq = protein_reader.next()
	if prot_id is None:
		break
	prot_id = prot_id.split('_')[0]
	prot_seqs[prot_id] = prot_seq
	nucl_id, nucl_seq = nucleotide_reader.next()
	nucl_seqs[nucl_id] = nucl_seq

aligned_nucl_seqs = {}

for seq_id in prot_seqs:
	nucl_start = 0
	output = ''
	prot_seq = prot_seqs[seq_id]
	nucl_seq = nucl_seqs[seq_id]
	for char in prot_seq:
		if char == '-':
			output += '---'
		else:
			output += nucl_seq[nucl_start:nucl_start+3]
			nucl_start += 3
	aligned_nucl_seqs[seq_id] = output

Ds = []
log_odds_list = []
for i in range(0, len(aligned_nucl_seqs['query']), 3):
	q_codon = aligned_nucl_seqs['query'][i:i+3]
	if q_codon not in table:
		continue
	dS = 0
	dN = 0
	for seq_id in aligned_nucl_seqs:
		if seq_id == 'query':
			continue
		seq_codon = aligned_nucl_seqs[seq_id][i:i+3]
		if seq_codon in table:
			if seq_codon != q_codon:
				if table[seq_codon] == table[q_codon]:
					dS += 1

				else:
					dN += 1
	log_odds = np.log2((dN+0.0001)/(dS+0.0001))
	log_odds_list.append(log_odds)
	D = dN - dS
	Ds.append(D)

D_SD = stats.stdev(Ds)

zscores = []

for D in Ds:
	zscore = D/D_SD
	zscores.append(zscore)

colors = ['#fe7f02' if zscore > 1.645 or zscore < -1.645 else '#1f78b3' for zscore in zscores]



plt.style.use('dark_background')
fig,ax = plt.subplots()
ax.scatter([x for x in range(len(zscores))], log_odds_list, color=colors, s=3)
legend_elements = [Line2D([0], [0], marker='o', color='black', label='Significant at a=0.1', markerfacecolor='#fe7f02' ,markersize=5),
                   Line2D([0], [0], marker='o', color='black', label='Not significant', markerfacecolor='#1f78b3' ,markersize=5)]
ax.legend(handles=legend_elements, bbox_to_anchor=[1,0.5],loc='center left')
ax.set_xticklabels([])
ax.set_xlabel('Codons')
ax.set_ylabel('log2(dN/dS)')
plt.tight_layout()
fig.savefig('codon_selective_pressure.png')
plt.close()








