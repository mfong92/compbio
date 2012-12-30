#! /usr/bin/python
import sys
import string
import textwrap
import math

'''
Input: bayesian.py [fastaFile]
Output: log-ratio of posterior probabilities of whether the particular sequence 
		was randomly generated while excluding stop codons, or not excluding stop codons.
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

def errorCheck():
	#check if more than one argument was provided
	if len(sys.argv) < 2:
		print "No file name provided!"
		sys.exit()
	
	if len(sys.argv) > 3:
		print "Too many arguments!"
		sys.exit()


	
def bayesian(runningSequence):	
	numberSequences = 0
	ACount = 0
	TCount = 0
	CCount = 0
	GCount = 0
	totalCount = 0


	#keeps running counter of the four bases		
	runningSequence = runningSequence.upper()			
	for base in runningSequence:
		if base == "A":
			ACount=ACount+1
		if base == "T":
			TCount=TCount+1
		if base == "C":
			CCount=CCount+1
		if base == "G":
			GCount=GCount+1
	
	
	totalCount = ACount+TCount+CCount+GCount #total number of bases

	#Perc of each base
	APerc = float(ACount) / float(totalCount)
	CPerc = float(CCount) / float(totalCount)
	GPerc = float(GCount) / float(totalCount)
	TPerc = float(TCount) / float(totalCount)
	


	prob = {}
	prob['A'] = APerc
	prob['C'] = CPerc
	prob['G'] = GPerc
	prob['T'] = TPerc
	
	#initialize probabilities
	prob_S_G = 1
	prob_S_NotG =1

	#probability of seeing a stop codon if just IID of nucleotides
	stopPerc = float(TPerc*APerc*GPerc + TPerc*APerc*APerc + TPerc*GPerc*APerc)
	
	for i in range(totalCount/3):
		codon = (runningSequence[i*3:i*3+3])

		if (codon == 'TAA') or (codon == 'TAG') or (codon == 'TGA'):
			#cannot be (G=1) if a stop codon is ever encountered
			prob_S_G = 0
		else:
			#probability, assuming that there are no stop codons

			prob_S_G = float(prob_S_G*prob[codon[0]]*prob[codon[1]]*prob[codon[2]]/(1-stopPerc))

			#probability, assuming that there are no limitations on nucleotides
			prob_S_NotG = float(prob_S_NotG*prob[codon[0]]*prob[codon[1]]*prob[codon[2]])

	#overall probability of S occuring, given that the a priori probabilities are 0.5
	prob_S = (prob_S_G+prob_S_NotG)/2



	#bayesian inference formula
	probG_S = prob_S_G*.5/prob_S
	probNotG_S = prob_S_NotG*.5/prob_S
	
	

	print "P(G=1|S) = " + str(probG_S)
	if (probG_S):
		print "log-ratio of posterior probabilities = " + str(float(math.log(probG_S/probNotG_S,2)))
	else:
		#log(0) does not exist
		print "The log-ratio of posterior probabilities is undefined."

def main():

	errorCheck()
	fileName = sys.argv[1]	

	try:
		infile = open(fileName)
	except IOError:
		print "This file does not exist!"
		sys.exit()

	infile = open(fileName)
	seqList = []

	#run bayesian for each sequence in the fasta file
	for name,seq in readFasta(infile):

		print name
		bayesian(seq)
		print "\n"



if __name__ == "__main__":
	main()
