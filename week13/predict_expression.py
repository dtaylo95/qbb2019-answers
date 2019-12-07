"""
Usage: python predict_expression.py <PROJECT_FNAME.h5> <start_loc> <end_loc> <bin_size> <ACTVITY_FNAME.bed> <RNA_FNAME.bed>
"""

import sys
import hifive
import numpy
import copy

start_coord = int(sys.argv[2])
end_coord = int(sys.argv[3])
binsize = int(sys.argv[4])

hic = hifive.HiC(sys.argv[1], 'r')
data = hic.cis_heatmap(chrom='chr10', start=start_coord, stop=end_coord, binsize=binsize, datatype='fend', arraytype='full')
where = numpy.where(data[:, :, 1] > 0)
data[where[0], where[1], 0] /= data[where[0], where[1], 1]
data = numpy.log(data[:, :, 0] + 0.1)
data -= numpy.amin(data)

data_copy = copy.deepcopy(data)

activity_data = {}
expression_data = {}

activity_data_file = open(sys.argv[5])
for i, line in enumerate(activity_data_file):
	if i == 0:
		continue
	fields = line.strip().split()
	chromosome = fields[0]
	start = int(fields[1])
	activity = float(fields[4])
	if start >= start_coord and start < end_coord:
		index = (start - start_coord) / binsize
		activity_data[index] = activity

rna_data_file = open(sys.argv[6])
for i, line in enumerate(rna_data_file):
	if i == 0:
		continue
	fields = line.strip().split()
	chromosome = fields[0]
	start = int(fields[1])
	rna_level = float(fields[4])
	if start >= start_coord and start < end_coord:
		index = (start - start_coord) / binsize
		expression_data[index] = rna_level

for index in activity_data:
	data[:,index] *= activity_data[index]


selected_data = data[sorted(expression_data.keys()), :][:, sorted(activity_data.keys())]

#Predicted expression of genes (for which RNAseq data exists)
predicted_expression = numpy.sum(selected_data, axis=1)

gene_indices = sorted(expression_data.keys())
actual_expression = [expression_data[x] for x in gene_indices]

R = numpy.corrcoef(actual_expression, predicted_expression)[0,1]
print('Activity-by-Contact R = {}'.format(R))


outfile = open('expression_predictions.tsv', 'w')
outfile.write('bin_index\tactual_expression\tactivity_by_contact_pred_expression\tdistance_decay_pred_expression')

N = data_copy.shape[0]
triu = numpy.triu_indices(N, 0)
distance_signal = numpy.bincount(numpy.abs(triu[1] - triu[0]), weights=data[triu])
distance_signal /= numpy.bincount(numpy.abs(triu[1] - triu[0]))
distance_matrix = numpy.zeros((N, N), dtype=numpy.float32)
distance_matrix[triu] = distance_signal[numpy.abs(triu[1] - triu[0])]
distance_matrix[triu[1], triu[0]] = distance_matrix[triu]


distance_decay_data = distance_matrix[sorted(expression_data.keys()), :][:, sorted(activity_data.keys())]
predicted_decay_expression = numpy.sum(distance_decay_data, axis=1)

R_decay = numpy.corrcoef(actual_expression, predicted_decay_expression)[0, 1]
print('Distance Decay R = {}'.format(R_decay))



for i in range(len(gene_indices)):
	outfile.write('{}\t{}\t{}\t{}\n'.format(gene_indices[i], actual_expression[i], predicted_expression[i], predicted_decay_expression[i]))
outfile.close()

r_out = open('r_coeffs.txt','w')
r_out.write('Activity-by-Contact Expression Prediction\nR = {}\n\n'.format(R))
r_out.write('Distance Decay Expression Prediction\nR = {}'.format(R_decay))
r_out.close()
