"""
Usage: python get_assembly_stats.py <assembly_contigs.fa> > <out_file>
"""

import sys

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


def avg(lst):
	return float(sum(lst))/len(lst)


def N50(lst):
	lst.sort
	halfway = float(sum(lst))/2
	counter = 0
	for contig_length in lst:
		counter += contig_length
		if counter >= halfway:
			return contig_length


for i in range(1,len(sys.argv)):
	contig_file = open(sys.argv[i])
	contig_reader = FASTAReader(contig_file)

	contig_counter = 0
	contig_lengths = []

	while True:
		ident, sequence = contig_reader.next()
		if not ident is None: 
			contig_counter += 1
			contig_lengths.append(len(sequence))
		else:
			break
	
	num_contigs = contig_counter
	min_contig_length = min(contig_lengths)
	max_contig_length = max(contig_lengths)
	avg_contig_length = avg(contig_lengths)
	N50_contig_length = N50(contig_lengths)

	print('For contigs in %s' %(sys.argv[i]))
	print('\tNumber of Contigs: %d' %(num_contigs))
	print('\tMinimum Contig Length: %d' %(min_contig_length))
	print('\tMaximum Contig Length: %d' %(max_contig_length))
	print('\tAverage Contig Length: %d' %(avg_contig_length))
	print('\tN50 Contig Length: %d\n' %(N50_contig_length))

