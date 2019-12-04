
import matplotlib.pyplot as plt
import numpy as np



#def print_stuff(seqs):
#        for i in range(len(seqs)):
#		for j in range(len(seqs[i])):
#			print seqs[i][j]["header"]

#seqs = []
#seqs1 = []
#seqs2 = []
#seqs1.append({"header" : "headline"})
#seqs2.append({"header" : "Foobar!"})
#seqs.append(seqs1)
#seqs.append(seqs2)
#print_stuff(seqs)

list = [1,2,3,4,5,6,7]

#x = 0
#for i in list:
#	print i, x
#

#dict = {"a" : []}


#dict["a"].append(1)
#dict["a"].append(2)
#dict["a"].append(3)

#plt.plot(dict["a"])
#plt.show()

#array = np.array([], 'i')
#array = np.append(array, 0)
#array = np.append(array, 1)
#print array


qrray = np.array([2,4,6,8,10,12,14], 'i')
array2 = np.array([1,2,3,4,5,6,7], 'i')
array3 = np.array([1,2,3,4,5,6,7], 'i')

array /= array2 + array3
print array







