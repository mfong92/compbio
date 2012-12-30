#! /usr/bin/python
import sys
import string
import textwrap
import random
import re

# input from command line looks like this:
# $ python hw6.py matrixfile sequencefile

'''
Implementation of Needleman-Wunsch sequence alignment algorithm
Input: alignment.py [matrixFile] [sequenceFile]

matrixFile format: first line - list of N characters, separated by a space
next N lines: N x N matrix of values, all separated by a space, indicating gap penalty

sequenceFile format: fasta format, contains at least 2 sequences.

Returns: highest scoring alignment of the first 2 sequences, along with the alignment score.
'''


def readFasta(infile):
	name, seq = None, []
	for line in infile:
		if line.startswith(">"):
			if name:
				yield(name, ''.join(seq))
			name, seq = line.strip(), []
		else:
			seq.append(line.strip())
	if name:
		yield (name, ''.join(seq))

def prepMatrix(firststr,secondstr,scorelist,gap):

	len1 = len(firststr)
	len2 = len(secondstr)

	array=[]
	for i in range(len2+2):
		array.append([])
		for j in range(len1+2):
			array[i].append([])
			
	#initialize top row and left column of matrix
	array[0][0] = "*"
	array[0][1] = "*"
	array[1][0] = "*"

	for x in range (len1):
		array[0][x+2] = firststr[x]

	for x in range (len2):
		array[x+2][0] = secondstr[x]

	#initialize gap scores
	array[1][1] = 0
	for x in range(len1):
		array[1][x+2] = gap * (x+1)

	for x in range(len2):
		array[x+2][1] = gap * (x+1)

	for i in range(2, len2+2):
		for j in range(2, len1+2):
			temp1 = array[i-1][j-1] + scorelist[firststr[j-2] + secondstr[i-2]]
			temp2 = array[i][j-1] + gap
			temp3 = array[i-1][j] + gap
			array[i][j] = max(temp1, temp2, temp3)
			
	return array

def printAlignment(name1,firststr,name2,secondstr,array,gap):

	len1 = len(firststr)
	len2 = len(secondstr)

	arrString = ""	
	for x in range(len2+2):
		for y in range(len1+2):
			arrString += ("|" + str(array[x][y]))
		arrString += ("|" + '\n')

	a = len2+1
	b = len1+1
	align1 = ""
	align2 = ""

	while ((a > 1) or ( b > 1)):

		currentScore = array[a][b]

		left = array[a][b-1]
		up = array[a-1][b]
		diag = array[a-1][b-1]
		#print currentScore
		if (currentScore == left + gap):
		#	print "1"
			b = b-1
			align2 += "-"
			align1 += (str(firststr[b-1]))
		elif (currentScore == up + gap):
		#	print "2"
			a = a-1
			align1 += "-"
			align2 += (str(secondstr[a-1]))

		else:
		#	print "3"
			a = a-1
			b = b-1
			if ((b > 0) & (a > 0)):
				align2 += (str(secondstr[a-1]))
				align1 += (str(firststr[b-1]))

	align1 = align1[::-1]
	align1 = textwrap.wrap(align1, 80)
	align1 = '\n'.join(align1)

	align2 = align2[::-1]
	align2 = textwrap.wrap(align2, 80)
	align2 = '\n'.join(align2)
	
	print '\n'
	print arrString

	print name1 + " aligned to " + name2
	print align1
	print align2
	print '\n'


def main():

	# Read from command line and check for valid input files
	if len(sys.argv) < 3:
		print "Not enough arguments! Need input file and output file."
		sys.exit()
		
	input1File = sys.argv[1]
	input2File = sys.argv[2]
		
	try:
		infile = open(input1File)
	except IOError:
		print "The input file does not exist!"
		sys.exit()

	try:
		infile = open(input2File)
	except IOError:
		print "The output file does not exist!"
		sys.exit()

	# Parse the matrix file
	matrixFile = open(input1File)
	characters = matrixFile.readline().upper()
	
	
	#ERROR CHECKING BEGINS HERE ----------------
	nonchars = re.search(r"^[a-zA-Z]", characters)
	
	if not nonchars:
		raise AttributeError("Invalid characters in first line of matrix file!")
	
	characters = characters.split()
	num = len(characters)

	scores = {}
	for i in range (num):
	
		currentLine = matrixFile.readline()
		currentLine = currentLine.split()
		
		if len(currentLine) != num:
			raise TypeError("Incorrect Matrix format!")
		
		for j in range (num): 
			stringPair = characters[i] + characters[j]
			stringPair = stringPair.upper()
			
			try:
				scores[stringPair] = int(currentLine[j])
			except ValueError:
				print "Incorrect Matrix line format!"
				sys.exit()
				
	gapScore = matrixFile.readline()
	
	lastlinecheck = re.search(r"^-[0-9]", gapScore)
	
	if not lastlinecheck:
		raise AttributeError("Invalid Gap Score line format!")
	
	gapScore = gapScore.split()
	
	if len(gapScore) > 1:
		raise TypeError("Too many numbers on gap score line!")
	
	gapScore = int(gapScore[0])
	
	if matrixFile.readline():
		raise TypeError("Too many lines in matrix file!")
	
	matrixFile.close()

	# Parse the sequence file
	stringFile = open(input2File)
	
	seqList = []
	for name,seq in readFasta(stringFile):
		seqList.append([name,seq])
		
	stringFile.close()
	
	if len(seqList) < 2:
		raise IOError("Not enough input sequences!")
	
	SEQ1LOC		= 0
	SEQ2LOC		= 1
	NAMELOC 	= 0
	SEQLOC		= 1
	
	string1name = seqList[SEQ1LOC][NAMELOC][1:].strip()
	string1 	= seqList[SEQ1LOC][SEQLOC]
	string2name	= seqList[SEQ2LOC][NAMELOC][1:].strip()
	string2		= seqList[SEQ2LOC][SEQLOC]
	
	len1 = len(string1)	
	len2 = len(string2)
	
	for char in string1:
		if not char in characters:
			raise AttributeError("Invalid character in sequence file!")
	
	for char in string2:
		if not char in characters:
			raise AttributeError("Invalid character in sequence file!")
			
	#ERROR CHECKING ENDS HERE ---------------

	# Computes the alignment score matrix
	arr = prepMatrix(string1,string2,scores,gapScore)
	
	# Prints out the alignment
	printAlignment(string1name,string1,string2name,string2,arr,gapScore)

if __name__ == "__main__":
	main()
