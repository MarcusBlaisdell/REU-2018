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
        self.trainTotal = sampleSizeK
        self.trainMistakes = 0

        ### read the selected samples from the file into a list:
        readFile = open (fileName, 'r')

        r = 0
        self.scoreList = []

        for s in sampleList:
            t = s - r - 1
            ### iterate to the line before the one in s:
            for x in range(t):
                tempString = readFile.readline()
            ### read the sample:
            sampleX = processLine(readFile, 'bogus')
            if sampleX == 0:
                continue
            theScore = self.calcScore (sampleX)
            scoreTuple = (theScore, sampleX)
            self.scoreList.append(scoreTuple) # **replace with insertSort
            r = s

        readFile.close()

        ### sort the scores: (This could be replaced by insertSort)
        ### Put the highest scores at the top
        ### to rank them from highest to lowest
        self.scoreList.sort(reverse=True)

        ### For the top 20% of scores,
        ### check each pair, if we have an unmatched pair, (good,bad), (bad,good),
        ### update the weight with their difference:

        ### get the total number of good from the list of samples:

        for k in range(len(self.scoreList)):
            if self.scoreList[k][1][1] == 1:
                self.trainGood += 1

        ### a mistake is an irrelevant record scored higher than a relevant record:

        k = 0
        flag = 0

        while flag == 0:
            count = 1
            if self.scoreList[k][1][1] == -1:
                l = k + 1
                while (l < len(self.scoreList)):
                    if (self.scoreList[l][1][1] == 1):
                        self.trainMistakes += count
                        k = l
                        if k == len(self.scoreList) - 1:
                            flag = 1
                        l = len(self.scoreList)
                    else:
                        l += 1
                        if l == len(self.scoreList):
                            flag = 1
                        count += 1
            else:
                k += 1
                if k == (len(self.scoreList) - 1):
                    flag = 1

        '''
        print 'ScoreList:'
        for x in self.scoreList:
            print x[0], ' : ', x[1][1]
        print 'trainGood: ', self.trainGood
        print 'trainMistakes: ', self.trainMistakes
        '''

        ### evaluate pairs, if mismatch, update weight:

        self.traindp = (len(self.scoreList) / 5)
        self.trainnpr = 0

        for i in range((len(self.scoreList) / 5) - 1):
            if self.scoreList[i][1][1] == 1:
                #self.trainnpr += 1
                ### Number of good in the top 20%
                self.trainnpr += 1
            else:
                self.trainMistakes += 1

            for j in range (i + 1, (len(self.scoreList) / 5)):
                if (self.scoreList[i][1][1] != self.scoreList[j][1][1]):
                    self.updateWeight (self.scoreList[i], self.scoreList[j])
                #else:
                    #self.trainnpr += 1

    ### end pairwise function:

    #def testWeight (self, trainDataSet, testDataSet, outFile, t):
    def testWeight (self, testDataSet, t):
        # reset mistakes count so each iteration starts at 0

        xit = [] # x-sub-i-sub-t, the training vector for the current iterations
        yit = 0 # y-sub-i-sub-t, the label for the current training vector (yStar)
        yHat = 0.0

        ### evaluate all test samples:
        self.testTotal = 0
        self.testMistakes = 0

        for i in range(len(testDataSet)):
            self.testTotal += 1
            xit = testDataSet[i][0]
            yit = testDataSet[i][1]

            ### make the prediction:
            #yHat = yit * self.dotProd(xit) # Method A
            #yHat = yit * (self.dotProd(xit) + self.b) # Method B
            #yHat = sign(yit * (self.dotProd(xit) + self.b)) # Method C
            yHat = sign(yit * self.dotProd(xit)) # Method D

            ### if the value is actually good,
            ### and we predicted good, increment npr
            ### which is the # lines predicted as good that are actually good
            if yit == 1:
                self.testGood += 1
                if yHat > 0:
                    self.testnpr += 1
                    self.testdp += 1
                else:
                    self.testMistakes += 1
            if yit == -1:
                if yHat > 0:
                    self.testMistakes += 1
                    self.testdp += 1
        print 'testMistakes: ', self.testMistakes
        print 'testTotal: ', self.testTotal

    ### end testWeight function

    ### function dotProd(), dot product by index
    def dotProd(self, xArray):
        result = 0.0

        for element in xArray:
            if self.w.get(element[0], '--') == '--':
                #self.w[element[0]] = 0
                continue
            else:
                result += self.w[element[0]] * element[1]

        return result

    ### end dotProd() function

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

        ### Check each index in each feature set,
        ### If they match, use the difference to update,
        ### otherwise, use the value of the unmatched index

        iFlag = 0
        jFlag = 0
        while iFlag * jFlag != 1:

            if iFlag == 0:
                if jFlag == 0:
                    ### if the current index of i is lower than the current index of j:
                    if (i[1][0][ii][0] < j[1][0][jj][0]):
                        self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * i[1][0][ii][1] * i[1][1]
                        ii += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                    ### if the current index of i matches the current index of j:
                    elif (i[1][0][ii][0] == j[1][0][jj][0]):
                        self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * ((i[1][0][ii][0] * i[1][1]) + (j[1][0][jj][0] * j[1][1]))
                        ii += 1
                        jj += 1
                        if ii == len(i[1][0]):
                            iFlag = 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                    ### if the current index of i is higher than the current index of j:
                    elif (i[1][0][ii][0] > j[1][0][jj][0]):
                        self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                        jj += 1
                        if jj == len(j[1][0]):
                            jFlag = 1
                ### jFlag == 1, there are no more indexes in j:
                else:
                    self.w[i[1][0][ii][0]] = self.w[i[1][0][ii][0]] + self.eta * i[1][0][ii][1] * i[1][1]
                    ii += 1
                    if ii == len(i[1][0]):
                        iFlag = 1

            ### iFlag == 1, no more indexes in i, iterate through remainder of j:
            else:
                if (jFlag == 0):
                    self.w[j[1][0][jj][0]] = self.w[j[1][0][jj][0]] + self.eta * j[1][0][jj][1] * j[1][1]
                    jj += 1
                    if jj == len(j[1][0]):
                        jFlag = 1

    ### end updateWeight function

### End of file
