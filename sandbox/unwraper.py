

import os
import sys
import re

def get_header(file_name):
	file = open(file_name, 'r')
	header_line = re.compile('^@.*:.*:.*:')
	for line in file:
		if header_line.match(line):
			header = line[0:line.index(':')]
			file.close()
			return (header)

def unwrap(file_name):
	header = get_header(file_name)
	line_2 = re.compile('^\+|^\+header')
	header = re.compile(header)
	file = open(file_name, 'r')
	out_path = "./out/" + os.path.splitext(os.path.basename(file_name))[0] + "unwrapped.fastq"
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


unwrap(sys.argv[1])
