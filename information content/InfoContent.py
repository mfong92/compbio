#! /usr/bin/python

import sys
import re
import math

'''
Input: InfoContent.py [stockholmFile]
Output: Can be configured to print the mutual information between any two pairs of sequences, as well as the entropy of any sequence
        Currently prints the 10 lowest entropy sequences, along with the top 50 mutual information pairs.
'''        



# Extracts all sequences from a Stockholm file
# Input:
#   infile - sequence alignment in Stockholm format
# Output:
#   namelist - ordered list of sequence names
#   seqlist - ordered list of sequences
#   * namelist and seqlist have a 1-to-1 correspondence based on indices

def processStockholm(infile):
    namelist = []
    seqlist = []
    NAME = 0
    SEQUENCE = 1
    for line in infile:
        if line.startswith("//"):
            break
        elif line.startswith("#"):
            continue
        else:
            seqline = line.strip().split()
            namelist.append(seqline[NAME])
            seqline[SEQUENCE] = seqline[SEQUENCE].upper()
            charcheck = re.search(r"[^ACGU.-]", seqline[SEQUENCE])
            if charcheck:
                raise TypeError("Invalid character in alignment! Please use an RNA alignment.")
            seqlist.append(seqline[SEQUENCE])
    if len(namelist) != len(seqlist):
        raise IOError("Incorrect file read/format!")
    return namelist, seqlist

# Returns a list of probability distributions for each column in an alignment
# Input:
#   seqlist - a list of sequence data in the alignment, which can include
#       gap characters, represented by a '.' or '-'
# Output:
#   colDistrib - a list of dictionaries holding probability distributions for
#       characters at a specified column

def genPD(seqlist):
    colDistrib = []
    totalchars = len(seqlist[0])
    seqnum = len(seqlist)
    for pos in range(totalchars):
        chardict = { 'A':0,
                     'U':0,
                     'G':0,
                     'C':0,
                     '-':0
                     }
        for seq in seqlist:
            char = seq[pos]
            if (char == '-') or (char == '.'):
                chardict['-'] = chardict['-'] + 1
            else:
                chardict[char] = chardict[char] + 1
        for i in chardict:
            chardict[i] = chardict[i] * 1.0 / seqnum
        colDistrib.append(chardict)
    return colDistrib

def testGenPD():
    
    testlist1 = ['AUCG-', 'AUCG-']
    distrib1 = genPD(testlist1)
    assert distrib1[0].get('A') == 1.0, "test1 'A' failed"
    assert distrib1[1].get('U') == 1.0, "test1 'U' failed"
    assert distrib1[2].get('C') == 1.0, "test1 'C' failed"
    assert distrib1[3].get('G') == 1.0, "test1 'G' failed"
    assert distrib1[4].get('-') == 1.0, "test1 '-' failed"

    testlist2 = ['AAUU', 'ACUG']
    distrib2 = genPD(testlist2)
    assert distrib2[0].get('A') == 1.0, "test2 col0 failed"
    assert distrib2[1].get('A') == 0.5, "test2 col1-1 failed"
    assert distrib2[1].get('C') == 0.5, "test2 col1-2 failed"
    assert distrib2[2].get('U') == 1.0, "test2 col2 failed"
    assert distrib2[3].get('U') == 0.5, "test2 col3-1 failed"
    assert distrib2[3].get('G') == 0.5, "test2 col3-2 failed"

    testlist3 = ['AA', 'UU', 'G-', 'C-']
    distrib3 = genPD(testlist3)
    assert distrib3[0].get('A') == 0.25, "test3 col0 failed"
    assert distrib3[0].get('U') == 0.25, "test3 col0 failed"
    assert distrib3[0].get('G') == 0.25, "test3 col0 failed"
    assert distrib3[0].get('C') == 0.25, "test3 col0 failed"
    assert distrib3[0].get('-') == 0.0, "test3 col0 failed"
    assert distrib3[1].get('A') == 0.25, "test3 col1 failed"
    assert distrib3[1].get('U') == 0.25, "test3 col1 failed"
    assert distrib3[1].get('G') == 0.0, "test3 col1 failed"
    assert distrib3[1].get('C') == 0.0, "test3 col1 failed"
    assert distrib3[1].get('-') == 0.5, "test3 col1 failed"
    
    print "genPD() tests passed!"

# Calculates the entropy of a probability distribution
# Input:
#   pdist - a dictionary holding the probability distribution of a column
# Output:
#   entropy - the entropy of a column

def calcEntropy(pdist):
    entropy = 0
    for char in pdist:
        prob = pdist[char]
        if prob == 0.0:
            continue
        else:
            entropy = entropy + prob*math.log(prob, 2)
    return entropy

# Calculates the joint probability distribution for each pair of columns
# Input:
#   seqlist - a list of sequences from the alignment
# Output:
#   jpmatrix - a matrix holding the joint probability distributions of column pairs
#       with format p(i,j) = jpmatrix[i][j-i-1]

def calcJP(seqlist):
    seqnum = len(seqlist)
    columns = len(seqlist[0])
    jpmatrix = []
    for i in range(seqnum):
        jpmatrix.append([])
    for i in range(columns):
        for j in range(i+1, columns):
            jpdict = {'AA':0,
                      'AC':0,
                      'AG':0,
                      'AU':0,
                      'A-':0,
                      'CA':0,
                      'CC':0,
                      'CG':0,
                      'CU':0,
                      'C-':0,
                      'GA':0,
                      'GC':0,
                      'GG':0,
                      'GU':0,
                      'G-':0,
                      'UA':0,
                      'UC':0,
                      'UG':0,
                      'UU':0,
                      'U-':0,
                      '-A':0,
                      '-C':0,
                      '-G':0,
                      '-U':0,
                      '--':0,
                      }
            for row in seqlist:
                letter1 = row[i]
                if letter1 == '.':
                    letter1 = '-'
                letter2 = row[j]
                if letter2 == '.':
                    letter2 = '-'
                pairname = letter1 + letter2
                jpdict[pairname] = jpdict[pairname] + 1
            for pair in jpdict:
                jpdict[pair] = jpdict[pair] * 1.0 / seqnum
            jpmatrix[i].append(jpdict)
    return jpmatrix

def testCalcJP():
    list1 = ['AAG','GGA']
    jpmat1 = calcJP(list1)
    assert jpmat1[0][0]['AA'] == 0.5, "AA test failed"
    assert jpmat1[0][0]['GG'] == 0.5, "GG test failed"
    assert jpmat1[0][1]['AG'] == 0.5, "AG test failed"
    assert jpmat1[0][1]['GA'] == 0.5, "GA test failed"
    assert jpmat1[1][0]['AG'] == 0.5, "GG test failed"
    assert jpmat1[1][0]['GA'] == 0.5, "GA test failed"
    print "calcJP() tests passed!"

# Calculates the mutual information content
# Input:
#   pdlist - a list of the probability distributions
#   jpmatrix - matrix containing the joint probability distributions
# Output:
#   mimatrix - matrix containing the mutual information values

def calcMutInfo(pdlist, jpmatrix):
    columns = len(pdlist)
    mimatrix = []
    for i in range(len(jpmatrix)):
        mimatrix.append([])
    for i in range(columns):
        for j in range(i+1, columns):
            mival = 0.0
            midict = {'AA':0,
                      'AC':0,
                      'AG':0,
                      'AU':0,
                      'A-':0,
                      'CA':0,
                      'CC':0,
                      'CG':0,
                      'CU':0,
                      'C-':0,
                      'GA':0,
                      'GC':0,
                      'GG':0,
                      'GU':0,
                      'G-':0,
                      'UA':0,
                      'UC':0,
                      'UG':0,
                      'UU':0,
                      'U-':0,
                      '-A':0,
                      '-C':0,
                      '-G':0,
                      '-U':0,
                      '--':0,
                      }
            for charpair in midict:
                first = charpair[0]
                second = charpair[1]
                joint = jpmatrix[i][j-i-1][charpair]
                pi = pdlist[i][first]
                pj = pdlist[j][second]
                if joint == 0.0 or pi == 0.0 or pj == 0.0:
                    continue
                else:
                    logpart = joint*math.log(joint * 1.0 / (pi * pj))
                    midict[charpair] = midict[charpair] + logpart
            for charpair in midict:
                mival = mival + midict[charpair]
            mimatrix[i].append(mival)
    return mimatrix

# Handles sorting of entropy values from low to high and prints the lowest
# 10 values while preserving the format of the input list
# Input:
#   entlist - list of entropy values, unsorted
# Output:
#   Prints the ten lowest entropy values (by magnitude) along with their
#   column indices in the alignment

def entropyLow10(entlist):
    print "Lowest 10 entropy values (lowest magnitude = closest to 0):"
    printlist = []
    for entval in entlist:
        printlist.append(entval)
    printlist.sort()
    for i in range(10):
	currentEntry = printlist[i]
	currentEntry = currentEntry.split()
	col = str(currentEntry[1])
        val = str(currentEntry[0])
        print "Col(i): (%s) %s" % (col, val)

# Handles pre-sorting of mutual information values from high to low
# while also preserving the format of the input matrix
# Input:
#   mimatrix - matrix of mutual information values
# Output:
#   milist - a linear sorted list of mutual information values

def condenseMI(mimatrix):
    milist = []
    for group in mimatrix:
        for val in group:
            milist.append(val)
    milist.sort()
    milist = milist[::-1]
    return milist

# Returns the column indices for a MI value
# Input:
#   mival - a floating point mutual information value
#   mimatrix - matrix of mutual information values
# Output:
#   coli - the position of column i used to calculate mival
#   colj - the position of column j used to calculate mival

def findMIPos(mival, mimatrix):
    for i in range(len(mimatrix)):
        for j in range(len(mimatrix[i])):
            if mival == mimatrix[i][j]:
                coli = i
                colj = j+i+1
                return (coli, colj)

# Prints the top 20 Mutual Information Values
# Input:
#   mimatrix - matrix of mutual information values
# Output:
#   Prints the top 20 mutual information values

def MItop20(mimatrix):
    print "Top 20 mutual information values:"
    printlist = condenseMI(mimatrix)
    for i in range(20):
        colstr = findMIPos(printlist[i], mimatrix)
        print "Col(i,j): %s %.12f" % (colstr, printlist[i])
        
# Prints the top 50 Mutual Information values
# Input:
#   mimatrix - matrix of mutual information values
# Output:
#   Prints the top 50 mutual information values

def MItop50(mimatrix):
    print "Top 50 mutual information values:"
    printlist = condenseMI(mimatrix)
    for i in range(50):
        colstr = findMIPos(printlist[i], mimatrix)
        print "Col(i,j): %s %.12f" % (colstr, printlist[i])

# Executes all test functions
# This is only used for testing!

def testAll():
    testGenPD()
    testCalcJP()

# Executes the program when run from the command line
# $ python InfoContent.py STOCKHOLMFILE    

def main():

    # read command-line arguments of the form $ InfoContent.py STOCKHOLMFILE
    if (len(sys.argv) == 1):
        raise IOError("Please enter an alignment file in Stockholm format")
    elif (len(sys.argv) == 2):
        try:
            alignfile = open(sys.argv[1])
        except IOError:
            raise IOError("Invalid file name!")
        alignfile = open(sys.argv[1])
    else:
        raise IOError("Too many arguments!")

    print "Now reading alignment data..."

    # Read alignment data from Stockholm file
    names, sequences = processStockholm(alignfile)
    alignfile.close()

    print "Now calculating probability distributions..."
    
    # Calculate the probability distribution that a randomly selected row
    # contains character x at column i
    pdlist = genPD(sequences)

    print "Now calculating column entropy..."

    # Calculate the entropy of each column in the alignment, storing
    # entropy values with their column indices
    entlist = []
    i = 0
    for pd in pdlist:
        entr = calcEntropy(pd)
        strEntr = str(entr) + " " + str(i)
	i += 1
	entlist.append(strEntr)

    print "Now calculating joint probability distributions..."
    print "(this one's a little slow, so hang on tight!)"

    # Calculate the joint probability distribution that a randomly selected row
    # contains character x in column i and character y in column j
    jpmatrix = calcJP(sequences)

    print "Now calculating mutual information values..."

    # Calculate the mutual information
    mimatrix = calcMutInfo(pdlist, jpmatrix)

    print "Now printing output values..."

    # Print out the values
    print ""
    entropyLow10(entlist)
    print ""
    MItop20(mimatrix)
    print ""
    MItop50(mimatrix)

if __name__ == "__main__":
    main()