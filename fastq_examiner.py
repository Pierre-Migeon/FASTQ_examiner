#!/usr/bin/env python3

"""
I have just now begun to learn python. It seems fun and straightforward. 

"""

#############################################
#
#  A program to assess FASTQ format
#  correctness and assess fastq quality. 
#
#  Also written in order to learn Python.
#  Created by Pierre Migeon
#
#  Usage: fastq_looker.py -f1 file_1.fq [-f2 file_2.fq] 
#
#############################################
import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
#############################################
##     Summary stats and plots go here
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
		plt.plot(100 * (N_count[0:max_length] / N_count[750:max_length + 750]), color = 'blue') #linestyle='dashed')
		plt.xlabel("Position in Read")
		plt.ylabel("Percent N")
		plt.title("%%N by length of read for %s" % os.path.basename(file_name))
		plt.grid()
		plt.show()
	return(N_count)

def plot_total_ns(N_count):
	max_length = N_count.argmin()
	plt.plot(100 * (N_count[0:max_length] / N_count[750:max_length + 750]), color = 'blue') #, linestyle='dashed')
	plt.xlabel("Position in Read")
	plt.ylabel("Percent N")
	plt.title("%N by length of read cumulative for all files")
	plt.grid()
	plt.show()

def plot_number_of_x_length(lens, max_len, file_name):
	plt.xlabel("Length of Read")
        plt.ylabel("Number of Reads")
	plt.plot(lens[0:max_len + 2], color = 'blue') #, linestyle='dash_dot')
	plt.title("Read Length Distribution for %s" % os.path.basename(file_name))
	plt.grid()
	plt.show()

#This function plots the length distribution for the reads.
def number_of_x_length(seqs, file_plots):
	t_lens = np.zeros((750,), dtype=int)
	max_len = 0;
	for i in range(len(seqs)):
		lens = np.zeros((750,), dtype=int)
		for j in range(len(seqs[i])):
			length = len(seqs[i][j]["seq"])
			lens[length] += 1
			if length > max_len:
				max_len = length
		t_lens += lens
		if (file_plots):
			plot_number_of_x_length(lens, max_len, seqs[i][0]["filename"])
	plot_number_of_x_length(t_lens, max_len, "All Files")

def polish_arrays(nucleotides, max_len):
	for letter in nucleotides:
		while len(nucleotides[letter]) < max_len:
			nucleotides[letter] = np.append(nucleotides[letter], 0)

def plot_percent_gc(nucleotides, max_len, file_name):
	sum = np.zeros((max_len,), dtype=float)
	polish_arrays(nucleotides, max_len)
	sum = nucleotides["A"] + nucleotides["T"] + nucleotides["C"] + nucleotides["G"] + nucleotides["N"]
	for letter in nucleotides:
		nucleotides[letter] /= sum
	plt.xlabel("Position in Read")
	plt.ylabel("Proportion of Base")
	plt.title("Proportion of Each Base by Position in Read in %s" % os.path.basename(file_name))
	plt.plot(nucleotides["A"][0:max_len], color = 'red')
	plt.plot(nucleotides["T"][0:max_len], color = 'blue')
	plt.plot(nucleotides["C"][0:max_len], color = 'yellow')
	plt.plot(nucleotides["G"][0:max_len], color = 'green')
	plt.legend(["A", "T", "C", "G"])
	plt.grid()
	plt.show()

# Might consider converting all characters to 
# uppercase during intake of files, in case you
# run into fastq files with lowercase, would
# break this graph. Unlikely edge case but still
# You need to do one for the files cumulatively here as well.
def percent_gc(seqs, file_plots):
	nucleotides = []
	max_len = 0
	for i in (range(len(seqs))):
		nucleotides.append({"A" : np.array([], 'f'), "T" : np.array([], 'f'), "C" : np.array([], 'f'), "G" : np.array([], 'f'), "N" : np.array([], 'f')})
		for j in (range(len(seqs[i]))):
			x = 0
			for letter in seqs[i][j]["seq"]:
				while x not in range(len(nucleotides[i][letter])):
					nucleotides[i][letter] = np.append(nucleotides[i][letter], 0)
				nucleotides[i][letter][x] += 1
				x += 1
				if x > max_len:
					max_len = x
		if (file_plots):
			plot_percent_gc(nucleotides[i], max_len, seqs[i][0]["filename"])


def plot_quality_by_base(sum):
	plt.plot(sum)


def average_qual(qual, encoding):
	len = len(qual)
	qual = 0
	if encoding == 33:
		base = ord('!')
	else:
		base = ord('@')
	for char in qual:
		qual += ord(char) - ord('!')
	return (qual / len)

def get_encoding(seqs):
                min = '~'
                max = '!'
                for entry in range(len(seqs)):
                        for char in seqs[entry]["qual"]:
                                if char < min:
                                        min = char
                                if char > max:
                                        max = char
                if min < ':' or max < k:
			return (33)
                elif max > k:
                        return (64)

def quality_by_base(seqs, print_num):
	sum = np.zeros((45,), dtype=int)
	for file in range(len(seqs)):
		encoding = get_encoding(seqs["file"])
		for entry in range(len(seqs[file])):
			sum[average_qual(seqs[file][entry]["qual"]), encoding] += 1
		plot_quality_by_base(sum)

######################################
#  	Run Quality Graphs
######################################
def run_graphs(files, print_num, seqs):
	total_ns = np.zeros((1500,), dtype=float)
	for file in files:
		total_ns += summarize_ns(file, print_num)
        plot_total_ns(total_ns)
	#rename this function- percent by base
	percent_gc(seqs, print_num)
	number_of_x_length(seqs, print_num)
	quality_by_base(seqs, print_num)

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
#need to update this so that it crawls through a set of 4 guarenteed non-truncated entries.
#Otherwise, trunctation early in the file will throw everything off.
def get_header(file_name):
	header = re.compile('^@')
	plus = re.compile('^\+')
	file = open(file_name, 'r')
	first_line = file.readline()
	if len(first_line) < 2:
		print "Minimal information stored in the header line... Possible bugs ahead. You may want to reformat headerlines."
	seq_line = file.readline()
	x = 0
	if len(seq_line) == len(first_line):
		x = 1
	i = 2
	last_line = seq_line
	for line in file:
		if header.match(line) and not plus.match(last_line):
			if len(line) != len(seq_line) or x:
				if len(os.path.commonprefix([first_line, line])) > 1: 
					first_line = os.path.commonprefix([first_line, line])
		last_line = line
		i += 1
	file.close()
	return (first_line)

'''
def unwrap(file_name):
	header = get_header(file_name)
	line_2 = re.compile('^\+|^\+header')
	header = re.compile(header)
	file = open(file_name, 'r')
	out_path = "./out/" + os.path.splitext(os.path.basename(file_name))[0] + ".unwrapped.fastq"
	out = open(out_path, 'w')
	growing_line = ""
	i = 0
	for line in file:
		if header.match(line):
			if i == 0:
				i += 1
				out.write(line)
			else:
				 out.write("\n" + line)
		elif line_2.match(line):
			out.write("\n" + line)
		else :
			out.write(line.rstrip())
	out.write("\n")
	file.close()
	out.close()
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
'''
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
'''

def as_read(dict_entry):
	str = ""
	if "header" in dict_entry:
		str += dict_entry["header"]
	if "seq" in dict_entry:
		str += dict_entry["seq"]
		str += "\n"
	if "plus" in dict_entry:
		str += dict_entry["plus"]
        if "qual" in dict_entry:
                str += dict_entry["qual"]
		str += "\n"
	return (str)

def check_truncated(seqs):
	line_0 = re.compile('^@.*')
	line_2 = re.compile('^\+')
	outfile = open("./out/truncated_reads.fastq", 'w')
	removed = 0
	for i in range(len(seqs)):
		for j in range(len(seqs[i])):
			if j > len(seqs[i]) - 1:
				break
			if "header" not in seqs[i][j]:
				removed += 1
				outfile.write(as_read(seqs[i][j]))
				seqs[i].remove(seqs[i][j])
			elif "seq" not in seqs[i][j] or not check_correct_nucleotides(seqs[i][j]["seq"]):
				removed += 1
				outfile.write(as_read(seqs[i][j]))
				seqs[i].remove(seqs[i][j])
			elif "plus" not in seqs[i][j] or not line_2.match(seqs[i][j]["plus"]):
				removed += 1
				outfile.write(as_read(seqs[i][j]))
				seqs[i].remove(seqs[i][j])
			elif "qual" not in seqs[i][j]:
				removed += 1
				outfile.write(as_read(seqs[i][j]))
				seqs[i].remove(seqs[i][j])
			elif not len(seqs[i][j]["seq"]) == len(seqs[i][j]["qual"]):
				removed += 1
				outfile.write(as_read(seqs[i][j]))
				seqs[i].remove(seqs[i][j])
	return (removed);

#################################
# Run QC checks:
#################################
def run_checks(files, seqs):
	for file in files:
		if check_wrapped(file):
			print ("%s included wrapped text. The file will automatically be unwrapped" % file)
		i = check_truncated(seqs)
		if i:
                	print ("%s was truncated in some way. %i Truncated entries were automatically removed" % (file, i))


#############################################
# Place each file into array of dictionaries
#############################################
def put_in_struct(file_name):
	header = get_header(file_name)
	header = re.compile(header)
	non_seq = re.compile('^\+')
	lines = []
	j = -1
	file = open(file_name, 'r')
	for line in file:
		if header.match(line):
			lines.append({"header" : line})
			if j == -1:
				lines[0]["filename"] = file_name
			i = 1
			j += 1
		if i == 1 and check_correct_nucleotides(line):
			if "seq" in lines[j]:
				lines[j]["seq"] += line.rstrip()
			else:
				lines[j]["seq"] = line.rstrip()
		if i == 2:
			if "qual" in lines[j]:
				lines[j]["qual"] += line.rstrip()
			else:
				lines[j]["qual"] = line.rstrip()
		if i == 1 and non_seq.match(line):
			lines[j]["plus"] = line;
			i += 1
	return (lines)

####################################
# Check basic format correctness
####################################
def third_line(file_name):
	line_3 = re.compile('^\+')
	w = check_wrapped(file_name)
	file = open(file_name, 'r')
	for i in range(0, 3):
        	line = file.readline()
	if w:
		while check_correct_nucleotides(line):
			line = file.readline()
	file.close()
	return (line_3.match(line))

def legal_seq(file_name):
	file = open(file_name, 'r')
	for i in range(0, 2):
		line = file.readline()
	file.close()
	return (check_correct_nucleotides(line))

def starts_at(file_name):
	file = open(file_name)
	header = re.compile('^@')
	line = file.readline()
	file.close()
	return (header.match(line))

#run the above functions, check basic format correct.
def is_it_fastq(files):
	for file in files:
		if not starts_at(file):
			print("First entry in %s doesn't start with a @. %s might not actually be FASTQ format. Exiting." % (os.path.basename(file), os.path.basename(file)))
			sys.exit(0)
		if not legal_seq(file):
			print("The first entry in %s includes non-standard bases (something besides ATCGNU). %s might not actyally be FASTQ format. Exiting." % (os.path.basename(file), os.path.basename(file)))
			sys.exit(0)
		if not third_line(file):
			print("The first entry in %s is missing the '+' line. %s might not actyally be FASTQ format. Exiting." % (os.path.basename(file), os.path.basename(file)))
			sys.exit(0)

def main():
#################################
#  Parse command line arguements 
#################################
	parser = argparse.ArgumentParser(
					description='This is a script to check \
						the validity of FASTQ files, to \
						possibly correct formating errors, \
						and subsequently to produce summary \
						statistics about the FASTQ files\
						(for instance, the user may be\
						interested in examining files \
						before and after cleaning...)',
					)
	parser.add_argument('-c', help='correct invalid files')
	parser.add_argument('-v', help='check validity of files and then exit')
	parser.add_argument('-f1', '--forward', dest='fastq_1', help='The path to the first fastq', required=True)
	parser.add_argument('-f2', '--reverse', dest='fastq_2', help='The path to the second fastq', required=False)
	parser.add_argument('-i', '--plot_individual', dest='plot_num', help='produce individual summary plots', required=False, action='store_true')
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
	#Verify correct format:
	is_it_fastq(files)
	#if correct, read into array of dictionaries:
	seqs = []
	for file in files:
		seqs.append(put_in_struct(file))
	#Run QC checks. Truncated entries removed and files unwrapped.
	run_checks(files, seqs)
	#print summary table... Number of reads, quality encoding
	#Summary_table(seqs)
	#Run graphs and summary statistics
	run_graphs(files, args.plot_num, seqs)

######################################
#  Run Main!
######################################
if __name__ == '__main__':
	main()
