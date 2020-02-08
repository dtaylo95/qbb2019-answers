#!/usr/bin/env Python3

import sys

if len(sys.argv) > 1:
	f = open(sys.argv[1])
else:
	f = sys.stdin

for i, line in enumerate(f):
	line = line.rstrip('\n')
	if not line.startswith('#'):
		fields = line.split('\t')
		if fields[2] == 'gene':
			print(int(fields[4]) - int(fields[3]))

#for i, line in enumerate(f):
#	print(line.rstrip('\n'))
#	if i >= 10:
#		break

# print(f.readline().rstrip('\n'))