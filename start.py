#!/usr/bin/env python

import sys

with open(sys.argv[1], 'r') as forward_reads:
	print(forward_reads.read())

with open(sys.argv[2], 'r') as reverse_reads:
	print(reverse_reads.read())

