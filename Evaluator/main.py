##########################################
### PEGASOS versus Perceptron Test
### for genomic data
### Marcus Blaisdell
### Parallel Processing
###
##########################################

import multiprocessing as mp
from pg_Class import *
from pr_Class import *
from functions import *
import logging
import sys
import time

### create array of the values to label the training, testing, results files:

#geneList = ['0157', '2016C', 'CH611', 'Co6114', 'ED1a', 'EDL933-1', 'FAP1', '_isolate102', 'RS76', 'UMN026']
#geneList = ['2016C', 'CH611', 'Co6114']
geneList = ['Co6114', 'CH611', '2016C']
#geneList = ['2016C']
geneTest = '0157'
#kmerList = [11, 13, 15, 17]
kmerList = [11]
#subList = ['a', 'b', 'c', 'd']
subList = ['a']
#biasVariableList = [1, 10, 50, 100]
biasVariableList = [0]

maxReadSize = 2100

### run the comparison for each training/testing set:

#thePath = "/home/marcus/Data/Genomic-Data/Testing/"
#thePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/"
thePath = "/Users/MarcusBlaisdell/Documents/LinuxShare/tenK/"
#thePath = ""
#resultsStringPre = "/home/marcus/Documents/LinuxShare/Results/resultszU-A-b-1-100_"
resultsStringPre = "results_" + geneTest + "_"
resultsStringPost = ".csv"

### create empty log file:
logging.basicConfig(filename='run.log',level=logging.DEBUG)
logging.info("Log File created")

def individualProcess(kmer, letter):
    myPGClass = pg_Class()
    myPRClass = pr_Class()

    ### The datasets are the same for each function, create them once to save memmory:
    trainDataSet = []
    testDataSet = []

    myPGClass.kmer = kmer
    myPRClass.kmer = kmer

    resultsFile = resultsStringPre + str(kmer) + resultsStringPost
    print 'resultsFile: ', resultsFile
    logging.info ('results File: ' + resultsFile)

    print 'set: ', str(kmer)
    logging.info ('set: ' + str(str(kmer)))

    print 'loading data:'
    logging.info('loading data')

    trainDataSet = []
    testDataSet = []

    print 'opening results file:'
    logging.info ('opening results file:' + '\n')
    outFile = open(resultsFile, "w")

    outFile.write ('Perceptron'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')
    logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running perceptron . . . '
    logging.info ('running perceptron . . . ' + '\n')

    startTime = time.time()
    for biasVariable in biasVariableList:
        myPRClass.b = biasVariable
        #####
        for t in range(myPRClass.T):
            #print 'iteration: ', t
            #####
            for gene in geneList:
                readGene (thePath, gene, kmer, maxReadSize, "perceptron", myPRClass, t)

                ### report results of training

                myPRClass.trainAccuracy = 100 - (100 * (myPRClass.trainMistakes / float(myPRClass.trainTotal)) )
                #print 'trainMistakes = ', myPRClass.trainMistakes, '   : train success = ', \
                      #myPRClass.trainAccuracy, '%'
                if myPRClass.traindp > 0:
                    myPRClass.trainPrecision = myPRClass.trainnpr / float(myPRClass.traindp)
                if myPRClass.trainGood > 0:
                    myPRClass.trainRecall = myPRClass.trainnpr / float(myPRClass.trainGood)
                #print 'trainPrecision: ', myPRClass.trainPrecision
                #print 'trainRecall: ', myPRClass.trainRecall
                if (myPRClass.trainPrecision + myPRClass.trainRecall) > 0:
                    myPRClass.trainF1 = 2 * myPRClass.trainPrecision * myPRClass.trainRecall / (myPRClass.trainPrecision + myPRClass.trainRecall)
                #print 'F1: ', myPRClass.trainF1
                outFile.write('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                              + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                              + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                              + str(myPRClass.b) + ','\
                              + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')
                logging.info ('train' + ',' + str(myPRClass.trainPrecision) + ',' + str(myPRClass.trainRecall) \
                              + ',' + str(myPRClass.trainF1) + ',' + str(myPRClass.trainAccuracy)\
                              + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                              + str(myPRClass.b) + ','\
                              + str(myPRClass.trainMistakes) + ',' + str(t + 1) + '\n')
                print len(myPRClass.w)

            logging.info('training Complete')

            readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPRClass, t)
            ### report results of testing

            myPRClass.testAccuracy = 100 - (100 * (myPRClass.testMistakes / float(myPRClass.testTotal)) )
            #print 'testMistakes = ', myPRClass.testMistakes, '   : test success = ', \
                  #myPRClass.testAccuracy, '%'

            if myPRClass.testdp > 0:
                myPRClass.testPrecision = myPRClass.testnpr / float(myPRClass.testdp)
            myPRClass.testRecall = myPRClass.testnpr / float(myPRClass.testGood)
            #print 'testPrecision: ', myPRClass.testPrecision
            #print 'testRecall: ', myPRClass.testRecall
            if (myPRClass.testPrecision + myPRClass.testRecall) > 0:
                myPRClass.testF1 = 2 * myPRClass.testPrecision * myPRClass.testRecall / (myPRClass.testPrecision + myPRClass.testRecall)

            outFile.write('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                          + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                          + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                          + str(myPRClass.b) + ','\
                          + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')
            logging.info('test' + ',' + str(myPRClass.testPrecision) + ',' + str(myPRClass.testRecall) \
                          + ',' + str(myPRClass.testF1) + ',' + str(myPRClass.testAccuracy)\
                          + ',' + str(myPRClass.T) + ',' + str(myPRClass.k) + ',' \
                          + str(myPRClass.b) + ','\
                          + str(myPRClass.testMistakes) + ',' + str(t + 1) + '\n')

    endTime = time.time()
    print 'Perceptron runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    logging.info ('Perceptron runTime: ' + str(endTime - startTime) + '\n')


    ### end test perceptron
    '''

    ### Add a header for Pegasos:
    outFile.write ('\n' + '\n' + 'Pegasos' + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running Pegasos . . .'
    logging.info('running Pegasos . . .' + '\n')

    ### record start time to determine runtime of algorithm
    startTime = time.time()

    ### run for each value of T:
    for i in myPGClass.Tlist:
        myPGClass.T = i
        ### run for each value of k:
        for j in myPGClass.klist:
            myPGClass.k = j
            for l in myPGClass.lamList:
                myPGClass.lam = l

                ### repeat for T iterations:
                for t in range(myPGClass.T):
                    #myPGClass.pegasosBatch(trainDataSet, testDataSet, t)

                    outFile.write('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                                  + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')

                    #myPGClass.testWeight(testDataSet)
                    outFile.write('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                                  + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.testMistakes) + '\n')


    endTime = time.time()
    print 'Pegasos runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))
    logging.info('Pegasos runTime: ' + str(endTime - startTime) + '\n')
    '''
    print 'closing output file . . .'
    logging.info ('closing output file . . .' + '\n')

    outFile.close()

### end individualProcess function

### Run each N in parallel:

processes = [mp.Process(target=individualProcess, args=(k, l)) for k in kmerList for l in subList]

for p in processes:
    p.start()

for p in processes:
    p.join()

print '*** program complete ***'
logging.info ('*** program complete ***' + '\n')
