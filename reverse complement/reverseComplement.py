#! /usr/bin/python
import sys
import string
import textwrap

'''
Returns the reverse complement (as DNA or RNA) of a specified FASTA file
Input: reverseComplement.py fileName <rna>
'''

#check if more than one argument was provided
if len(sys.argv) < 2:
	print "No file name provided!"
	sys.exit()
	
if len(sys.argv) > 3:
	print "Too many arguments!"
	sys.exit()

fileName = sys.argv[1]

#check to see if file could be opened
try:
	infile = open(fileName)
except IOError:
	print "This file does not exist!"
	sys.exit()

individualLines = [] #array to store individual lines of sequence
counter = 0
lineList = [] #list of all sequences
headerList = [] #list of all headers

DNADict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'U':'A'}

for line in infile:

	line = line.rstrip()
	
	#search for headers
	if (line.startswith(">")):
		headerList.append(line)
		counter = counter + 1
	
	#if hit new header, store all of the lines in the temporary storage (individualLines) into permanent storage(lineList)
	if (counter > 1):
		sequence = ""
		sequence = "".join(individualLines)
		counter = counter - 1
		individualLines = []
		lineList.append(sequence)

	#if not header, add new line of sequence
	if not (line.startswith(">") or line.startswith(";")):
		individualLines.append(line)

#finish storing last lines of sequence into the permanent storage		
if len(headerList) > len(lineList):
	sequence = ""
	sequence = "".join(individualLines)
	lineList.append(sequence)
	
n = len(headerList)
m = len(lineList)

	
#this for loop prints each sequence header and reverse sequence in sequential order
for x in range(0,n):
	sequence = lineList[x]
	header = headerList[x]
	reverseSequence = sequence[::-1] #reverses the sequence of the list
	
	
	#base pair complementation
	reverseSequence = reverseSequence.upper() #convert all bases to upper case
	transcript = ""
	for y in range(len(reverseSequence)):
		currentBase = reverseSequence[y]
		transcript += DNADict[currentBase]
	
	#check to see if user wants rna
	if len(sys.argv) == 3:
		if ((sys.argv[2] == 'rna') | (sys.argv[2] == 'RNA')):
			transcript = transcript.replace('T','U')
		
	
	#add number of bp to header
	print headerList[x] + ", " + str(len(transcript)) + " bp"
	text = textwrap.wrap(transcript, 80) #keep max of 80 bases/line
	finalText = '\n'.join(text)
	print finalText
	print '\n'
	
 