#! /usr/bin/python
import sys
import string
import textwrap
import random
import optparse
import re
from optparse import OptionParser

'''
Input: sequenceSimulator.py [option]

Generates a random FASTA file according to the nucleotide composition in the input file, or random (25%/nucleotide if no file is given).
--calc takes an input FASTA file and prints a randomly generated FASTA file with the same length/nucleotide characteristics
--load loads a parameter file that must be in the format: '%A: float %C: float %G: float %T: float avgLength: float numberSequences: int',
	where float and int are the data types that follow each label (must keep exact spacing).
--save saves the parameters of input FASTA file, or input param file with the last input FASTA file appended on the bottom.
--stats takes a polynucleotide sequence to calculate the statistics for in a given FASTA file: ie. --stats ACGTC
--statsFile saves the polynucleotide statistics to the given file (saves the statistics of the searched sequence in both the 
input file and the randomly generated file)

Returns: randomly generated FASTA file, meeting all FASTA specifications. Stop codons are removed from the randomly generated sequence.
'''




myparser = OptionParser()
myparser.add_option("-c", "--calc", action ="store", dest = "calc", help = "calculates nucleotide composition and length statistics of input FASTA file; prints a randomly generated FASTA file with the same characteristics   format: --calc fasta.file" ) 
myparser.add_option("-l", "--load", action ="store", dest = "load", help = "loads text file that contains the following parameters, all separated by spaces: '%A: float %C: float %G: float %T: float avgLength: float numberSequences: int'   format: --load myparams" ) 
myparser.add_option("-s", "--save", action ="store", dest = "save", help = "saves parameters of input param file or FASTA file into target param file. also outputs randomly generated FASTA file \n\n format: --save myparams" )
myparser.add_option("--stats", action ="store", dest = "stats", help = "input polynucleotide sequence to calculate statistics for (in the FASTA file specified by --calc)")
myparser.add_option("--statsFile", action = "store", dest = "statsFile", help= "saves statistics of --stats polynucleotide sequence to given file" ) #can compare with purely random seq probability


(options, args) = myparser.parse_args()





#generates a "random" base determined by the input percentages of each base
def generateBase(APercentage, CPercentage, GPercentage, TPercentage):
	x = random.random()
	if x < APercentage:
		base = 'A'
	elif ((x >= APercentage) & (x < APercentage+TPercentage)):
		base = 'T'
	elif ((x >= APercentage+TPercentage) & (x < APercentage+TPercentage+CPercentage)):
		base = 'C'
	else:
		base = 'G'
	return base
	

#print FASTA sequence file of uniform nucleotide composition (no input file)
if ((options.calc == None) & (options.load == None)):
	sequence = ""
	APercentage = 0.25
	CPercentage = 0.25
	GPercentage = 0.25
	TPercentage = 0.25
	sequence = 'ATG'
	for i in range (997):
		base = generateBase(APercentage, CPercentage, GPercentage, TPercentage)
		sequence = sequence+base
	averageLength = len(sequence)
	numberSequences = 1

		
#calculating statistics of input FASTA file
if options.calc != None:
	numberSequences = 0
	ACounter = 0
	TCounter = 0
	CCounter = 0
	GCounter = 0
	runningSequence = ""
	try:
		infile = open(options.calc)
	except IOError:
		print "This file does not exist!"
		sys.exit()
		
	#count the number of sequences
	for line in infile:
		line = line.rstrip()
		if (line.startswith(">")):
			numberSequences = numberSequences + 1
		if not (line.startswith(">")):
			runningSequence = runningSequence + line
	for base in runningSequence:
		if base == "A":
			ACounter=ACounter+1
		if base == "T":
			TCounter=TCounter+1
		if base == "C":
			CCounter=CCounter+1
		if base == "G":
			GCounter=GCounter+1

	
	totalCounter = ACounter+TCounter+CCounter+GCounter #total number of bases

	#percentage of each base
	APercentage = float(ACounter) / float(totalCounter)
	CPercentage = float(CCounter) / float(totalCounter)
	GPercentage = float(GCounter) / float(totalCounter)
	TPercentage = float(TCounter) / float(totalCounter)

	averageLength = float(totalCounter) / float(numberSequences)
	
	
	sequence = "ATG" #gives the sequence a start codon

	#generates random sequences weighted by the percentages in given file
	for i in range (int(averageLength-3)):
		base = generateBase(APercentage, CPercentage, GPercentage, TPercentage)
		sequence = sequence+base

#loading params from parameter file
if options.load != None:
	try:
		infile = open(options.load)
	except IOError:
		print "This file does not exist!"
		sys.exit()

	#parsing file for pertinent data
	for line in infile:
		line = line.rstrip()
		entries = line.split(' ')
		APercentage = float(entries[1]) / 100
		CPercentage = float(entries[3]) / 100
		GPercentage = float(entries[5]) / 100
		TPercentage = float(entries[7]) / 100
		averageLength = float(entries[9])
		numberSequences = int(entries[11])
		
	sequence = "ATG"

		#generates random sequences weighted by the percentages in given file
	for i in range (int(averageLength-3)):
		base = generateBase(APercentage, CPercentage, GPercentage, TPercentage)
		sequence = sequence+base

#if there is a save file and NO load file, save the calculated percentages 
if ((options.save != None) & (options.load == None)):
	outfile = open(options.save, "a")
	outfile.write("%A: " + str(APercentage*100) + " %C: " + str(CPercentage*100) + " %G: " +str(GPercentage*100) +" %T: "+str(TPercentage*100) +" avgLength: "+str(averageLength)+" numberSequences: " + str(numberSequences) + '\n')
	outfile.close()

#if there is both a load file and save file, copy the load file params into save file
if ((options.save != None) & (options.load != None)):
	try:
		infile = open(options.load)
	except IOError:
		print "This file does not exist!"
		sys.exit()

	outfile = open(options.save, "w")
	for line in infile:
		outfile.write(line)
	outfile.close()
	infile.close()

#error correcting for stop codons: if a stop codon is found, replace one of the nucleotides with another, at random
for i in range (int(len(sequence)/3)):
	if ((options.calc == None) & (options.load == None)):
		APercentage = 0.25
		GPercentage = 0.25
		CPercentage = 0.25
		TPercentage = 0.25
	sequenceList = list(sequence)
	if ((sequence[3*i]=='T') & (sequence[3*i+1] == 'G') & (sequence[3*i+2] == 'A')):
		sequenceList[3*i+1] = generateBase(APercentage + GPercentage/3, CPercentage + GPercentage/3, 0, TPercentage+GPercentage/3)
	if ((sequence[3*i]=='T') & (sequence[3*i+1] == 'A') & (sequence[3*i+2] == 'G')):
		sequenceList[3*i+1] = generateBase(0, CPercentage + APercentage/3, GPercentage+APercentage/3, TPercentage+APercentage/3 )
	if ((sequence[3*i]=='T') & (sequence[3*i+1] == 'A') & (sequence[3*i+2] == 'A')):
		sequenceList[3*i] = generateBase(APercentage + TPercentage/3, CPercentage + TPercentage/3, GPercentage+TPercentage/3, 0)
	sequence = "".join(sequenceList)
	

#looks for requested polynucleotide string within original sequence and generated sequence
if (options.stats != None):
	seqString = options.stats
	pattern = re.compile(seqString)

	if (options.calc != None):
		matches = runningSequence.count(seqString)	
	generatedMatches = sequence.count(seqString)

	if(options.statsFile != None):
		outfile = open(options.statsFile, "a")
		if (options.calc != None):
			outfile.write('\n'+ "Seq: " + seqString + ".  Original file. Number of Matches: " + str(matches) + " Percentage: " + str(100*float(matches)/float(totalCounter- len(seqString) + 1))+ '\n' )
			outfile.write("Seq: " + seqString + ". Generated file. Number of Matches: " + str(generatedMatches) + " Percentage: " + str(100*float(generatedMatches)/float((len(sequence)-len(seqString) +1))) + '\n' )
		else:
			outfile.write('\n'+ "Seq: " + seqString + ". Generated file. Number of Matches: " + str(generatedMatches) + " Percentage: " + str(100*float(generatedMatches)/float((len(sequence)-len(seqString) +1))) + '\n' )

		outfile.close()
	
			
			
			
textSequence = textwrap.wrap(sequence, 80)
finalText = '\n'.join(textSequence)

print finalText	



			
			
			
			
			
			
			
