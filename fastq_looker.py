#!/usr/bin/env python3

"""
I have just now begun to learn python. It seems fun and straightforward. 

Friday: implemented parameters and usage statement etc. Wrote function to check to see if the fastq files are truncated and for general malformation. Good start.

Saturday: corrected the error in the check_correct_nucleotides script, started making file input more flexible. Started writing some of the section to summarize quality in graph form, starting with what will probably the easiest one.

Sunday: Added check for wrapped lines.

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
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
#############################################
##     Summary stats and plots go here!
#############################################
#This is to summarize the bases where Ns were found.
#I am assuming that you aren't getting reads that are longer than 750 bp.
def summarize_ns(file_name, plot_number):
	N_count = np.zeros((1500,), dtype=float)
	max_length = 0;
	i = 0
	file = open(file_name, 'r')
	for line in file:
		if i == 1:
			j = 0
			for char in line:
				N_count[j + 750] += 1
				if char == "N":
					N_count[j] += 1
				j += 1
			if j > max_length:
				max_length = j
		if i == 3:
			i = -1
		i += 1
	if plot_number:
		plt.plot(100 * (N_count[0:max_length] / N_count[750:max_length + 750]), color = 'green', linestyle='dashed')
		plt.xlabel("Position in Read")
		plt.ylabel("Percent N")
		plt.title("%%N by length of read for %s" % os.path.basename(file_name))
		plt.show()
	return(N_count)

def plot_total_ns(N_count):
	max_length = N_count.argmin()
	plt.plot(100 * (N_count[0:max_length] / N_count[750:max_length + 750]), color = 'green', linestyle='dashed')
	plt.xlabel("Position in Read")
	plt.ylabel("Percent N")
	plt.title("%N by length of read cumulative for all files")
	plt.show()


######################################
#  	Run Quality Graphs
######################################
def run_graphs(files, print_num):
        total_ns = np.zeros((1500,), dtype=float)
        for file in files:
                total_ns += summarize_ns(file, print_num)
        plot_total_ns(total_ns)


#############################################
## 	      quality check here
#############################################
def check_correct_nucleotides(line):
	acceptable = ['A', 'C', 'T', 'G', 'N', 'U', 'a', 't', 'c', 'g', 'n', 'u']
	line = line.rstrip()
	for i in line:
		if i not in acceptable:
			return (False)
	return (True)

'''
 The following is only relevant in some rare 
 cases (Sanger style fastq, for instance)
 But will still be useful to detect in such 
 a case. Only checks sequence lines assuming
 that these are the only things that are 
 wrapped. Wrapped seq ID lines will throw a
 trunctation error.
'''
#############################################
#  Unwrap the files 
#############################################
'''
def unwrap(file_name):
	file = open(file_name, 'r')
	growing_line = ""
	i = 0
	line_0 = re.compile('^@')
	line_2 = re.compile('^\+')
	for line in file:
		if i == 0:
		if i == 3:
			#cat together the line.
		if i == 2 and line_3.match(line):
			i += 1
'''		
				

#############################################
#  Check to see if the file is wrapped
#############################################
def check_wrapped(file_name):
	file = open(file_name, 'r')
	i = 0;
	for line in file:
		if check_correct_nucleotides(line):
			i += 1
		else :
			i = 0	
		if i > 1:
			return (True)
	return (False)

#############################################
#  Check truncation
#############################################
#this function checks to see if there is any sort of truncation:
#makes sure that you have lines with standard format, that qual line == seq line
#and that the sequence is only nucleotides or N.
def check_truncated(file_name):
	file = open(file_name, 'r')
	i = 0;
	line_0 = re.compile('^@.*')
	line_2 = re.compile('^\+')
	for line in file:
		if i == 0 and not line_0.match(line):
			return (True)
		if i == 1:
			if not check_correct_nucleotides(line):
				return (True)
			length_seq = len(line)
		if i == 2 and not line_2.match(line):
			return (True)
		if i == 3:
			if len(line) != length_seq:
				return (True)
			i = -1
		i += 1
	return (False)

#################################
# Run QC checks:
#################################
def run_checks(files):
        for file in files:
		if check_wrapped(file):
			print ("%s included wrapped text!!!" % file)
			#unwrap(file)
                elif check_truncated(file):
                        print ("%s was truncated in some way!!!" % file)

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
					)
	parser.add_argument('-c', help='correct invalid files')
	parser.add_argument('-v', help='check validity of files and then exit')
	parser.add_argument('-f1', '--forward', dest='fastq_1', help='The path to the first fastq', required=True)
	parser.add_argument('-f2', '--reverse', dest='fastq_2', help='The path to the second fastq', required=False)
	parser.add_argument('-ip', '--plot_individual', dest='plot_num', help='produce individual summary plots', required=False, action='store_true')
	args = parser.parse_args()

######################################
#  Open files
######################################
	files = []
	forward_file = args.fastq_1 
	files.append(forward_file)
	if args.fastq_2:
		reverse_file = args.fastq_2 
		files.append(reverse_file)
	#Run QC checks
	run_checks(files)
	#Run graphs and summary statistics
	run_graphs(files, args.plot_num)


######################################
#  Run Main!
######################################
if __name__ == '__main__':
	main()
