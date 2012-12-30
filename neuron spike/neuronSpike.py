#! /usr/bin/python
import sys
import string
import textwrap
import math

# This programs takes the input of three files:
# The first represents the spike times of a resting neuron
# The second represents the spike times of an active neuron
# The third represents the spike times of an unknown neuron

# ASSUMPTION: This program assumes that the neural spikes follow a Poisson distribution,
# where lambda is represented by the expected number of spikes over a given time period
# and NOT the expected frequency of spikes.

'''
Input: neuronSpike.py [resting] [active] [unknown] 
(files are described above)
Output: posterior probability that the unknown file is that of an active neuron.
'''

def errorCheck():
	#check if more than one argument was provided
	if len(sys.argv) < 4:
		print "Not enough file names provided!"
		sys.exit()
	
	if len(sys.argv) > 4:
		print "Too many arguments!"
		sys.exit()


def retrieve(fileName, type):
	try:
		infile = open(fileName)
	except IOError:
		print "This file does not exist!"
		sys.exit()
	
	counter = 0
	arrayTime=[]

	#counting the number of spikes in each file
	for line in infile:
		counter = counter + 1
		time = line.rstrip()
		arrayTime.append(time)
	counter = counter - 1

	#finding the total time that each experiment was run for
	timeTotal = float(arrayTime[counter])

	#finds the rate of spikes
	rate = float(counter / timeTotal)
	if (type == "rate"):
		return rate
	elif (type == "counter"):
		return counter
	elif (type == "time"):
		return timeTotal


def main():

	#defining a factorial function because python 2.4 doesn't support it.
	fac = lambda n:[1,0][n>0] or fac(n-1)*n

	#check for proper input format
	errorCheck()	
	
	#input files E1, E2, E3 (resting, active, unknown)
	fileName1 = sys.argv[1]	
	fileName2 = sys.argv[2]	
	fileName3 = sys.argv[3]	

	#call retrieve function to parse data
	restingRate = retrieve(fileName1,"rate")
	activeRate = retrieve(fileName2,"rate")
	unknownCounter = retrieve(fileName3,"counter")
	unknownTime = retrieve(fileName3,"time")
	

	print "Rate of Spikes_resting (/s): " + str(restingRate)
	print "Rate of Spikes_active (/s): " + str(activeRate)
	

	#lambda for poisson distribution (resting, active)
	restingLambda = float(restingRate * unknownTime)
	activeLambda = float(activeRate * unknownTime)

	print "Resting lambda - expected number of spikes over time period of unknown, if resting: " + str(restingLambda)
	print "Active lambda - expected number of spikes over time period of unknown, if active: " + str(activeLambda)


	#k for poisson distribution: actual number of spikes for neuron
	k = unknownCounter
	print "k - actual number of spikes over time period of unknown: " + str(k)
	
	#poisson: f = lambda^k * e^(-lambda) / k!
	#probabilities represent P(E3|resting), P(E3|active)
	probKResting = float((math.pow(restingLambda,k)) * (math.exp(-1*restingLambda)) / fac(k))
	probKActive = float((math.pow(activeLambda,k)) * (math.exp(-1*activeLambda)) / fac(k))


	pActive = 0.01
	pResting = 0.99

	#Bayesian inference formula
	#P(E3) = P(E3|active) + P(E3|resting)
	pE3 = float(probKActive*pActive + probKResting * pResting)

	#P(active|E3) = P(E3|active) * P(active) / P(E3)
	pActiveK = float(probKActive*pActive/(pE3))
	print "\n"
	print "Posterior probability that neuron was activated during E3: " + str(pActiveK)


if __name__ == "__main__":
	main()
