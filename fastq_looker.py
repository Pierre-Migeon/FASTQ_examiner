#!/usr/bin/env python

"""
I have just now begun to learn python. It seems fun and straightforward. 

Friday: implemented parameters and usage statement etc. Wrote function to check to see if the fastq files are truncated and for general malformation. Good start.


"""

#############################################
#
#  A program to assess FASTQ format
#  correctness and assess fastq quality. 
#
#  Also written in order to learn Python.
#  Created by Pierre Migeon
#
#  Usage: fastq_looker.py 
#
#############################################

import 	sys
import	re
import	argparse


def check_correct_nucleotides(line):
	acceptable = ['A', 'C', 'T', 'G', 'N']
	for i in line:
		if i not in acceptable:
			print i
			return (False)
	return (True)

def check_truncated(f_file):
	i = 0;
	line_1 = re.compile('^@.*')
	line_3 = re.compile('^\+')
	for line in f_file:
		if i == 0 and not line_1.match(line):
			return (True)
		if i == 1:
			if not check_correct_nucleotides(line):
				return (True)
			length_seq = len(line)
		if i == 2 and not line_3.match(line):
			return (True)
		if i == 3:
			if len(line) != length_seq:
				return (True)
			i = -1
		i += 1
	return (False)


def main():
#################################
#  Parse command line arguements 
#################################
	parser = argparse.ArgumentParser(
					description='This is a script to check \
						the validity of FASTQ files, to \
						possibly correct problems with \
						formatting found, and subsequently \
						to produce summary statistics about \
						the FASTQ files (for instance, user \
						may be interested in examining files \
						before and after cleaning...)',
					epilog='Perhaps you should try to run it using these pointers?\n\n'
					)
	parser.add_argument('-c', help='correct invalid files')
	parser.add_argument('-v', help='check validity of files and then exit')
	parser.add_argument('-f1', '--forward', dest='fastq_1', help='The path to the first fastq', required=True)
	parser.add_argument('-f2', '--reverse', dest='fastq_2', help='The path to the second fastq', required=False)
	args = parser.parse_args()

	forward_file = open(args.fastq_1, 'r')
	if args.fastq_2:
		reverse_file = open(args.fastq_2, 'r')

	run_checks(forward_file, reverse_file)
	#run_graphs(forward_file, reverse_file)



def run_checks(forward_file, reverse_file):
	if check_truncated(forward_file):
		print ("The first file was truncated in some way!!!")
	if (reverse_file):
		if check_truncated(reverse_file):
			print ("The first second file was truncated in some way!!!")


#def run_graphs():



if __name__ == '__main__':
	main()
