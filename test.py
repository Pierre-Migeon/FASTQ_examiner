


'''
This is a big 'ol
multi-line comment yeah!
'''



def check_nucleotides(line):
        acceptable = ['A', 'C', 'T', 'G', 'N']
        for i in line:
		if i not in acceptable:
			print ("ERROR!")



