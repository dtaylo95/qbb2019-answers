import matplotlib.pyplot as plt
from fractions import Fraction
from scipy.stats import binom
import seaborn as sns
import numpy as np

def model(start_freq, population_size):
	gen_count = 0
	while True:
		end_freq = Fraction(binom.rvs(2*population_size, start_freq), 2*population_size)
		gen_count += 1
		if end_freq == 0 or end_freq == 1:
			return gen_count
		else:
			start_freq = float(end_freq)




fixation_times = []
for i in range(1000):
	fixation_times.append(model(0.5, 100))

fig, ax = plt.subplots()
ax = sns.distplot(fixation_times, bins=100)
ax.set_xlabel('Generation Time to Fixation')
ax.set_ylabel('Density')
ax.text(0.65, 0.95, 'allele freq. = 0.5\npop. size = 100', transform=ax.transAxes, verticalalignment = 'top', fontsize=12)
ax.set_title('Wright-Fisher time to fixation')
fig.savefig('time_to_fix_hist.png')




pop_sizes = [10, 20, 100, 200, 1000, 2000, 10000, 20000]
pop_size_fix_times = {}

for pop_size in pop_sizes:
	pop_size_fix_times[pop_size] = []
	for i in range(10):
		pop_size_fix_times[pop_size].append(model(0.5, pop_size))

means = []
stds = []

for pop_size in pop_sizes:
	means.append(np.mean(pop_size_fix_times[pop_size]))
	stds.append(np.std(pop_size_fix_times[pop_size]))

fig, ax = plt.subplots()
ax.errorbar(pop_sizes, means, yerr=stds, fmt='o', c='black')
ax.set_xscale('log')
ax.set_xlabel('Population Size')
ax.set_ylabel('Generations to Fixation')
ax.set_title('Fixation vs. Population Size')
ax.text(0.05, 0.95, '5 trials at each pop. size\nStarting allele freq. = 0.5', transform=ax.transAxes, verticalalignment='top', fontsize=12)
fig.savefig('time_to_fix_pop_size.png')




allele_freqs = [x/20 for x in range(1,20)]
allele_freq_fix_times = {}
for freq in allele_freqs:
	allele_freq_fix_times[freq] = []
	for i in range(100):
		allele_freq_fix_times[freq].append(model(freq, 100))

means = []
stds = []

for freq in allele_freqs:
	means.append(np.mean(allele_freq_fix_times[freq]))
	stds.append(np.std(allele_freq_fix_times[freq]))

fig, ax = plt.subplots()
ax.scatter(allele_freqs, means)
ax.errorbar(allele_freqs, means, yerr=stds, linestyle="None", c='black')
ax.set_xlabel('Starting Allele Frequency')
ax.set_ylabel('Generations to Fixation')
ax.set_title('Fixation vs. starting allele frequency')
ax.text(0.05, 0.95, '100 trial at each freq.\npop. size = 100', transform=ax.transAxes, verticalalignment='top', fontsize=8)
fig.savefig('time_to_fix_allele_freq.png')




def selection_model(start_freq, population_size, selection_coeff):
	gen_count = 0
	while True:
		allele_count = start_freq*population_size*2
		prob = (allele_count * (1 + selection_coeff))/((2*population_size) - allele_count + allele_count*(1+selection_coeff))
		end_freq = Fraction(binom.rvs(2*population_size, prob), 2*population_size)
		gen_count += 1
		if end_freq == 0 or end_freq == 1:
			return gen_count
		else:
			start_freq = float(end_freq)

sel_coeffs = [-1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
sel_coeff_fix_times = {}
for s in sel_coeffs:
	sel_coeff_fix_times[s] = []
	for i in range(1000):
		sel_coeff_fix_times[s].append(selection_model(0.5, 100, s))

means = []
stds = []

for s in sel_coeffs:
	means.append(np.mean(sel_coeff_fix_times[s]))
	stds.append(np.std(sel_coeff_fix_times[s]))

fig, ax = plt.subplots()
ax.scatter(sel_coeffs, means)
ax.errorbar(sel_coeffs, means, yerr=stds, linestyle="None", c='black')
ax.set_xlabel('Selection Coefficient')
ax.set_ylabel('Generations to Fixation')
ax.set_title('Fixation vs. selection coefficient')
ax.text(0.05, 0.95, '100 trials at each sel. coeff.\npop. size = 100', transform=ax.transAxes, verticalalignment='top', fontsize=10)
fig.savefig('time_to_fix_sel_coeff.png')


