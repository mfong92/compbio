#! /usr/bin/python
import sys
import random
import re
import time

#timeStart = time.time();  #time check statement

'''
Implementation of the Nussinov Algorithm given an input of RNA as either upper or lower case.
Returns the number of maximum basepairs, assuming that there are no pseudoknots.

Input format: python nussinov.py augcaugc 
(augcaugc is an arbitrary RNA sequence)
'''

#check to see if an RNA sequence is input, assign seq and lenSeq to the sequence
if len(sys.argv) < 2:
	print "No sequence provided!"
	sys.exit()
seq = str(sys.argv[1]).upper()
lenSeq = len(seq)



#check of input RNA sequence is valid
def isValid(seq):
	disallowedBases = re.search(r"[^acguACGU]", seq)
	if disallowedBases:
		print "\nNot valid sequence!"
		sys.exit()
	else:
		return 1

#function that is called to give any WC match a score of 1

def wc(pair):
	if ((pair == 'AU') or (pair == 'UA') or (pair == 'GC') or (pair == 'CG')):
		return 1
	else:
		return 0

def nussinov(seq):
	#dictionary that corresponds to the score of each base pair - either 1 for the WC basepairs or 0 for all else.
	arr = []
	
	#initialize the 2d matrix as a array of arrays
	for i in xrange(lenSeq):
		arr.append([])
		for j in xrange(lenSeq):
			arr[i].append([])
			
	#initialize entire array as zeroes rather than just diagonals to improve performance 
	for a in range(lenSeq):
		for i in range(lenSeq):
			arr[a][i] = 0

	for x in range (1,lenSeq): #one dimension of the dynamic programming table
		for j in range(x,lenSeq): 
		
			i=j-x #need to fill out table moving along diagonal
			
			temp1=arr[i][j-1] #internal bp with i
			temp2=arr[i+1][j] #internal bp with j
			temp3=arr[i+1][j-1]+wc(seq[i]+seq[j]) #ends base pair, check if WC pair
			
			bifurcation = [] #temp array to hold possible bifurcation values
			temp4=0 #in case the following loop is not entered
			if i+1 < j: #don't enter loop if the for loop below will never execute
				for k in range(i+1,j): #keep track of all possible bifurcations
					bifurcation.append(arr[i][k]+arr[k+1][j])
				temp4=max(bifurcation)
				
				
			arr[i][j]=max(temp1, temp2, temp3, temp4) #return max to arr[i][j]
	return arr[0][lenSeq-1]


#actually run nussinov on input	
if isValid(seq):
	print "\nMaximum number of basepairs is " +  str(nussinov(seq))


#time check statement
#timeEnd = time.time()-timeStart
#print "\n"+ str(timeEnd) + " seconds"

		