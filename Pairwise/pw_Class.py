##########################################
### Pairwise Test
### class
### Marcus Blaisdell
##########################################

import random
from functions import *

### Create class pw_Class for pairwise:
class pw_Class ():
    w = {}
    scoreList = []
    eta = 1.0 # eta is learning rate
    T = 3
    k = 1
    b = 0.0 # default biasVariable
    trainMistakes = 0
    trainTotal = 0
    trainGood = 0
    trainnpr = 0 # numerator for trainPrecision and trainRecall
    traindp = 0 # denominator for trainPrecision
    trainPrecision = 0.0
    trainRecall = 0.0
    trainF1 = 0.0
    testMistakes = 0
    testTotal = 0
    testGood = 0
    testnpr = 0 # numerator for trainPrecision and trainRecall
    testdp = 0 # denominator for trainPrecision
    testPrecision = 0.0
    testRecall = 0.0
    testF1 = 0.0
    #fileLength = 38000000 # temporary placeholder, replace by actual file size (computed)
    fileLength = 10000

    def __init__(self):
        pass

    ### pairwise function:

    def pairwise(self, sampleSizeK, thePath, gene, kmer):
        fileName = thePath + gene + '-' + str(kmer)

        # 1. sample K, k-mers:
        # 2. score k-mers using current weight
        # 3. sort them according to scores
        # 4. identify mistaken pairs (good, bad)
        # update weight: w_t+1 = w_t + eta * (x_good - x_bad)

        #for x in range(sampleSizeK):
        ### create a list of size sampleSizeK of random line numbers
        ### from the file:
        sampleList = self.getRandomList(sampleSizeK)

        ### read the selected samples from the file into a list:
        readFile = open (fileName, 'r')

        r = 0
        self.scoreList = []
        for s in sampleList:
            self.trainTotal += 1
            t = s - r
            ### iterate to the line before the one in s:
            for x in range(t):
                tempString = readFile.readline()
            sampleX = processLine(readFile, 'bogus')
            theScore = self.calcScore (sampleX)
            scoreTuple = (theScore, sampleX)
            self.scoreList.append(scoreTuple) # **replace with insertSort
            r = s

        readFile.close()

        ### sort the scores: (This could be replaced by insertSort)
        self.scoreList.sort()

        ### check each pair, if we have an unmatched pair, (good,bad), (bad,good),
        ### update the weight with their difference:

        for i in range(len(self.scoreList) - 1):
            if self.scoreList[i][1][1] == 1:
                self.trainnpr +=1
                self.trainGood += 1
            else:
                self.trainMistakes += 1

            for j in range (i + 1, len(self.scoreList)):
                if (self.scoreList[i][1][1] != self.scoreList[j][1][1]):
                    #self.trainMistakes += 1
                    self.updateWeight (self.scoreList[i], self.scoreList[j])
                else:
                    self.trainnpr += 1

    ### end pairwise function:

    ### saveWeight function, save the weight vector to file:

    def saveWeight (self):
        writeFile = open ('weight.txt', 'w')
        weightList = self.w.keys()
        weightList.sort()
        for index in weightList:
            writeFile.write ( str(index) + ', ' + str(self.w.get(index)) + '\n')
        writeFile.close()

    ### end saveWeight function

    ### get a list of random numbers of size sampleSizeK

    def getRandomList (self, sampleSizeK):
        randomList = []
        for s in range (sampleSizeK):
            randomList.append( random.randrange(0, self.fileLength) )
        randomList.sort()
        return randomList

    ### end getRandomList function

    ### calcScore function, calculate the score of a sample k-mer:

    def calcScore (self, sampleX):
        xit = sampleX[0]
        yit = sampleX[1]
        theScore = 0.0

        ### The score is the dot-product of the weight and the
        ### feature vector:

        for index in xit:
            if self.w.get(index[0], '--') == '--':
                self.w[index[0]] = 0.0
            theScore += self.w.get(index[0]) * index[1]

        return theScore

    ### end calcScore function

    ### update the weight:

    def updateWeight (self, i, j):
        ii = 0
        jj = 0

        while ii < len(i) and jj < len(j):
            if i[1][0][ii][0] < j[1][0][jj][0]:
                self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * i[1][0][ii][1] * i[1][1]
                ii += 1
            if i[1][0][ii][0] == j[1][0][jj][0]:
                self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * (i[1][0][ii][0] * i[1][1] + j[1][0][jj][0] * j[1][1])
                ii += 1
                jj += 1
            if i[1][0][ii][0] > j[1][0][jj][0]:
                self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                jj += 1

    ### end updateWeight function

### End of file
