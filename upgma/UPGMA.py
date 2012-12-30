	#! /usr/bin/python
import sys
import string
import textwrap
import math

'''
Input: UPGMA.py [distmat file]
Output: Newick format tree generated with UPGMA algorithm.
'''


def importDist(fileName):
	
	infile = open(fileName)
	N = 0
	listTable=[]
	distMatrix=[]
	for line in infile:
		N = N+1
		line = line.strip()
		listTable.append(line)
	
	arr = []
	Name = []
	for i in xrange(N):
		arr.append([])
		for j in xrange(N):
			arr[i].append([])
			
	for a in range(N):
		for i in range(N):
			arr[a][i] = 0

	a = 0
	infile = open(fileName)

	for line in infile:
		line = line.strip().split()
		for i in range(N):
			arr[a][i] = float(line[i+1])

		Name.append(line[0])

		a+=1
	return (arr, N, Name)
			
			
def errorCheck():
	#check if more than one argument was provided
	if len(sys.argv) < 2:
		print "No file name provided!"
		sys.exit()
	
	if len(sys.argv) > 3:
		print "Too many arguments!"
		sys.exit()
		
def upgma(scoreMatrix, N, Name):

	#intialize a height matrix
	H = []
	for i in range(N):
		H.append(0)
		
	
	while (len(scoreMatrix)>1):
		#initialize current to be a very large number
		current = 100000000
		
		#find the lowest score in the entire matrix
		for i in range(len(scoreMatrix)):	
			sortedMatrix = sorted(scoreMatrix[i])
			lowScore = sortedMatrix[1]	
			if (lowScore < current):
				current = lowScore
				currentIndex = i
		
		#find the column number that the lowest score occurs at
		for i in range(len(scoreMatrix[currentIndex])):
			if (scoreMatrix[currentIndex][i] == current):
				currentIndex2 = i
		
		#delete the two columns correlating to those two leaves
		for i in range (len(scoreMatrix)):
			del scoreMatrix[i][currentIndex]
			del scoreMatrix[i][currentIndex2-1] #has -1 factor because the first column has already been deleted
			
		Dkn = []
		
		#store the list of distances for the rows that correspond to nodes to find new distances
		Din = scoreMatrix[currentIndex]
		Dij = scoreMatrix[currentIndex2]
		
		#delete those two rows
		del scoreMatrix[currentIndex]
		del scoreMatrix[currentIndex2-1] #to account for the first one being already deleted


		#append a new column to correspond to new node
		for i in range (len(scoreMatrix)):
			scoreMatrix[i].append((Din[i]+Dij[i])*0.5)
			
		#bottom right corner of matrix will be a 0
		Din.append(0)
		Dij.append(0)
		
		#construct a new row to correspond to new node
		for i in range(len(Din)):
			Dkn.append((Din[i] + Dij[i]) * 0.5)

		#append this new row (end of list)
		scoreMatrix.append(Dkn)
		
		#update the list of heights based on formula
		Hnew = (H[currentIndex]+H[currentIndex2]+current)*0.5
		H.append(Hnew)
		
		#name the new node
		newName = "(" + Name[currentIndex] + ":" + str(Hnew-H[currentIndex]) + "," + Name[currentIndex2] + ":" + str(Hnew-H[currentIndex2]) + ")"
		Name.append(newName)
		
		
		#these deletions help keep all three matrices consistent in terms of size: scoreMatrix, Name, H
		
		#delete names of nodes that were combined
		del Name[currentIndex]
		del Name[currentIndex2-1]
		
		#delete heights of nodes that were combined
		del H[currentIndex]
		del H[currentIndex2-1]
		
		
		
	return Name[0]
		
def main():

	errorCheck()
	fileName = sys.argv[1]	

	try:
		infile = open(fileName)
	except IOError:
		print "This file does not exist!"
		sys.exit()

	#import matrix
	(scoreMatrix, N, Name) = importDist(fileName)

	#find final newick tree
	newick = upgma(scoreMatrix,N,Name)
	
	print newick
		
if __name__ == "__main__":
	main()