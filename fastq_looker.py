#!/usr/bin/env python3

"""
I have just now begun to learn python. It seems fun and straightforward. 

Friday: implemented parameters and usage statement etc. Wrote function to check to see if the fastq files are truncated and for general malformation. Good start.

Saturday: corrected the error in the check_correct_nucleotides script, started making file input more flexible. Started writing some of the section to summarize quality in graph form, starting with what will probably the easiest one.

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
import  numpy as np
#############################################




#############################################
##	  Summary statistics go here!
#############################################
#This is to summarize the bases where Ns were found.
#I am assuming that you aren't getting reads that are longer than 500 bp.
def summarize_ns(file_name):
	N_count = np.zeros((500,), dtype=int)
	Max_length = 0;
	i = 0
	file = open(file_name, 'r')
	for line in file:
		if i == 1:
			j = 0
			for char in line:
				if char == "N":
					N_count[j] += 1
				j += 1
			if j > Max_length:
				Max_length = j
		if i == 3:
			i = -1
		i += 1	
	for p in N_count:
		print(p)



#############################################
## 	      quality check here
#############################################
def check_correct_nucleotides(line):
	acceptable = ['A', 'C', 'T', 'G', 'N']
	line = line.rstrip()
	for i in line:
		if i not in acceptable:
			return (False)
	return (True)

#this function checks to see if there is any sort of truncation:
#makes sure that you have lines with standard format, that qual line == seq line
#and that the sequence is only nucleotides or N.
def check_truncated(file_name):
	file = open(file_name, 'r')
	i = 0;
	line_1 = re.compile('^@.*')
	line_3 = re.compile('^\+')
	for line in file:
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

	files = []
	forward_file = args.fastq_1 #open(args.fastq_1, 'r')
	files.append(forward_file)
	if args.fastq_2:
		reverse_file = args.fastq_2 #open(args.fastq_2, 'r')
		files.append(reverse_file)
	run_checks(files)
	run_graphs(files)


def run_checks(files):
	for file in files:
		if check_truncated(file):
			print ("%s was truncated in some way!!!" % file)


def run_graphs(files):
	for file in files:
		summarize_ns(file)



######################################
#  Run Main!
######################################
if __name__ == '__main__':
	main()
