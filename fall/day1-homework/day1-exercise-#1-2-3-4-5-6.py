#!/usr/bin/env Python3

import sys

if len(sys.argv) > 1:
	f = open(sys.argv[1])
else:
	f = sys.stdin


align_count = 0
perfect_match_count =0
single_align_count = 0
MAPQsum = 0
num6_count = 0
print('Exercise #4')
for line in f:
	fields = line.rstrip('\n').split('\t')
	if fields[2] != '*':
		align_count += 1
		if align_count <= 10:
			print(fields[2])
		for field in fields:
			if field.startswith('NM'):
				score = field.split(':')[2]
				if score == '0':
					perfect_match_count += 1
			elif field.startswith('NH'):
				cur_align_count = field.split(':')[2]
				if cur_align_count == '1':
					single_align_count += 1
		MAPQsum += int(fields[4])
	if fields[2] == '2L' and int(fields[3]) in range(10000,20001):
		num6_count += 1

MAPQ_avg = float(MAPQsum)/align_count

print('\nExercise #1\n', align_count)
print('Exercise #2\n', perfect_match_count)
print('Exercise #3\n', single_align_count)
print('Exercise #5\n', MAPQ_avg)
print('Exercise #6\n', num6_count)