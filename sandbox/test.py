


'''
This is a big 'ol
multi-line comment yeah!
'''



def check_nucleotides(line):
        acceptable = ['A', 'C', 'T', 'G', 'N']
	line = line.rstrip()
        for i in line:
		if i not in acceptable:
			print ("ERROR!")


check_nucleotides("ATC GN\n")
