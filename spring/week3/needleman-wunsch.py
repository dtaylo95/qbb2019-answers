import numpy as np

seq1 = 'CATAAACCCTGGCGCGCTCGCGGCCCGGCACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCACCATGGTGCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGCGGAGGCCCTGGAGAGGATGTTCCTGTCCTTCCCCACCACCAAGACCTACTTCCCGCACTTCGACCTGAGCCACGGCTCTGCCCAGGTTAAGGGCCACGGCAAGAAGGTGGCCGACGCGCTGACCAACGCCGTGGCGCACGTGGACGACATGCCCAACGCGCTGTCCGCCCTGAGCGACCTGCACGCGCACAAGCTTCGGGTGGACCCGGTCAACTTCAAGCTCCTAAGCCACTGCCTGCTGGTGACCCTGGCCGCCCACCTCCCCGCCGAGTTCACCCCTGCGGTGCACGCCTCCCTGGACAAGTTCCTGGCTTCTGTGAGCACCGTGCTGACCTCCAAATACCGTTAAGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGCAAAAAAAAAAAAAAAAAAAAAA'
seq2 = 'GGGGCTGCCAACACAGAGGTGCAACCATGGTGCTGTCCGCTGCTGACAAGAACAACGTCAAGGGCATCTTCACCAAAATCGCCGGCCATGCTGAGGAGTATGGCGCCGAGACCCTGGAAAGGATGTTCACCACCTACCCCCCAACCAAGACCTACTTCCCCCACTTCGATCTGTCACACGGCTCCGCTCAGATCAAGGGGCACGGCAAGAAGGTAGTGGCTGCCTTGATCGAGGCTGCCAACCACATTGATGACATCGCCGGCACCCTCTCCAAGCTCAGCGACCTCCATGCCCACAAGCTCCGCGTGGACCCTGTCAACTTCAAACTCCTGGGCCAATGCTTCCTGGTGGTGGTGGCCATCCACCACCCTGCTGCCCTGACCCCGGAGGTCCATGCTTCCCTGGACAAGTTCTTGTGCGCCGTGGGCACTGTGCTGACCGCCAAGTACCGTTAAGACGGCACGGTGGCTAGAGCTGGGGCCAACCCATCGCCAGCCCTCCGACAGCGAGCAGCCAAATGAGATGAAATAAAATCTGTTGCATTTGTGCTCCAG'

base_index = {
	'A' : 0,
	'C' : 1,
	'G' : 2,
	'T' : 3
}

scoring_matrix = np.array([
		[   91, -114,  -31, -123 ],
        [ -114,  100, -125,  -31 ],
        [  -31, -125,  100, -114 ],
        [ -123,  -31, -114,   91 ] ])

gap_penalty = 300


def get_align(align1, align2):
	out_align = ''
	for i in range(len(align1)):
		if align1[i] == align2[i]:
			out_align += '|'
		else:
			out_align += ' '
	return out_align


alignment_matrix = np.zeros((len(seq1)+1, len(seq2)+1))
traceback_matrix = np.empty((len(seq1),len(seq2)), dtype='str')


for i in range(len(seq1)+1):
	alignment_matrix[i,0] = -1*i*gap_penalty

for j in range(len(seq2)+1):
	alignment_matrix[0,j] = -1*j*gap_penalty

for i in range(len(seq1)):
	for j in range(len(seq2)):
		possible_scores = {
			'v' : alignment_matrix[i,j+1] - gap_penalty, # vertical
			'h' : alignment_matrix[i+1,j] - gap_penalty, # horizontal
			'd' : alignment_matrix[i,j] + scoring_matrix[base_index[seq1[i]],base_index[seq2[j]]] # diagonal
		}
		
		maxes = None
		for traceback in possible_scores:
			if maxes == None or possible_scores[traceback] > possible_scores[maxes]:
				maxes = traceback
			elif possible_scores[traceback] == possible_scores[maxes[0]]:
				maxes = min(traceback, maxes)
		traceback_matrix[i,j] = maxes
		alignment_matrix[i+1,j+1] = possible_scores[maxes]

edit_distance = alignment_matrix[len(seq1),len(seq2)]



i = len(seq1)-1
j = len(seq2)-1


seq1_align = ''
seq2_align = ''

while i != -1 and j!= -1:
	direction = traceback_matrix[i,j]
	if direction == 'd':
		seq1_align = seq1[i] + seq1_align
		seq2_align = seq2[j] + seq2_align
		i -= 1
		j -= 1
	if direction == 'v':
		seq1_align = seq1[i] + seq1_align
		seq2_align = '-' + seq2_align
		i -= 1
	if direction == 'h':
		seq1_align = '-' + seq1_align
		seq2_align = seq2[j] + seq2_align
		j -= 1

if i == -1 and j != -1:
	seq1_align = '-'*(j+1) + seq1_align
	seq2_align = seq2[:j+1] + seq2_align
elif i != -1 and j == -1:
	seq1_align = seq1[:i+1] + seq1_align
	seq2_align = '-'*(i+1) + seq2_align


# seq1_check = ''
# for base in seq1_align:
# 	if base != '-':
# 		seq1_check += base

# if seq1_check == seq1:
# 	print('seq1 alignment is good')

# seq2_check = ''
# for base in seq2_align:
# 	if base != '-':
# 		seq2_check += base

# if seq2_check == seq2:
# 	print('seq2 alignment is good')

print('Alignment Score: {}'.format(edit_distance))

outfile = 'alignment.txt'
output = open(outfile,'w')
output.write('Alignment Score: {}\n\n'.format(edit_distance))

counter = 0

while counter < len(seq1_align):
	output.write('seq1\t' + seq1_align[counter:counter+60] + '\n    \t' + get_align(seq1_align[counter:counter+60], seq2_align[counter:counter+60]) + '\nseq2\t' + seq2_align[counter:counter+60] + '\n\n')
	counter += 60		