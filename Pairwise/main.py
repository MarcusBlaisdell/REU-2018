##########################################
### Pairwise Test
### for genomic data
### Marcus Blaisdell
### Parallel Processing
###
##########################################

import multiprocessing as mp
from pg_Class import *
from pr_Class import *
from pw_Class import *
from functions import *
import logging
import sys
import time

### create array of the values to label the training, testing, results files:

#geneList = ['0157', '2016C', 'CH611', 'Co6114', 'ED1a', 'EDL933-1', 'FAP1', '_isolate102', 'RS76', 'UMN026']
#geneList = ['2016C', 'CH611', 'Co6114']
geneList = ['2016C']
geneTest = '0157'
#geneTest = '2016C'
#kmerList = [11, 13, 15, 17]
kmerList = [11]
#subList = ['a', 'b', 'c', 'd']
subList = ['a']
#biasVariableList = [1, 10, 50, 100]
biasVariableList = [0]

maxReadSize = 500000
sampleSizeK = 1000

### run the comparison for each training/testing set:

#thePath = "/home/marcus/Data/Genomic-Data/Testing/"
#thePath = "/data/doppa/users/mblaisdell/genome-assembly-datasets/"
thePath = "/home/marcus/Testing/LOO2/tenK/"
#thePath = ""
#resultsStringPre = "/home/marcus/Documents/LinuxShare/Results/resultszU-A-b-1-100_"
resultsStringPre = "results_" + geneTest + "_"
resultsStringPost = ".csv"

### create empty log file:
logging.basicConfig(filename='run.log',level=logging.DEBUG)
logging.info("Log File created")

def individualProcess(kmer, letter):
    myPGClass = pg_Class ()
    myPRClass = pr_Class ()
    myPWClass = pw_Class ()

    myPGClass.kmer = kmer
    myPRClass.kmer = kmer
    myPWClass.k = sampleSizeK

    resultsFile = resultsStringPre + str(kmer) + resultsStringPost
    print 'resultsFile: ', resultsFile
    logging.info ('results File: ' + resultsFile)

    print 'set: ', str(kmer)
    logging.info ('set: ' + str(str(kmer)))

    print 'loading data:'
    logging.info('loading data')

    print 'opening results file:'
    logging.info ('opening results file:' + '\n')
    outFile = open(resultsFile, "w")

    outFile.write ('Pairwise'  + '\n' + '\n')
    outFile.write ('Gene left out: ' + geneTest + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')
    logging.info('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Bias' + ','\
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running pairwise . . . '
    logging.info ('running pairwise . . . ' + '\n')

    startTime = time.time()
    for biasVariable in biasVariableList:
        myPRClass.b = biasVariable
        #####
        for t in range(myPWClass.T):
            #print 'iteration: ', t
            #####
            for gene in geneList:
                myPWClass.pairwise(sampleSizeK, thePath, gene, kmer)

                ### report results of training

                myPWClass.trainAccuracy = 100 - (100 * (myPWClass.trainMistakes / float(myPWClass.trainTotal)) )
                if myPWClass.traindp > 0:
                    myPWClass.trainPrecision = myPWClass.trainnpr / float(myPWClass.traindp)
                if myPWClass.trainGood > 0:
                    myPWClass.trainRecall = myPWClass.trainnpr / float(myPWClass.trainGood)
                #print 'trainPrecision: ', myPWClass.trainPrecision
                #print 'trainRecall: ', myPWClass.trainRecall
                if (myPWClass.trainPrecision + myPWClass.trainRecall) > 0:
                    myPWClass.trainF1 = 2 * myPWClass.trainPrecision * myPWClass.trainRecall / (myPWClass.trainPrecision + myPWClass.trainRecall)
                #print 'F1: ', myPWClass.trainF1
                outFile.write('train' + ',' + str(myPWClass.trainPrecision) + ',' + str(myPWClass.trainRecall) \
                              + ',' + str(myPWClass.trainF1) + ',' + str(myPWClass.trainAccuracy)\
                              + ',' + str(myPWClass.T) + ',' + str(myPWClass.k) + ',' \
                              + str(myPWClass.b) + ','\
                              + str(myPWClass.trainMistakes) + ',' + str(t + 1) + '\n')
                logging.info ('train' + ',' + str(myPWClass.trainPrecision) + ',' + str(myPWClass.trainRecall) \
                              + ',' + str(myPWClass.trainF1) + ',' + str(myPWClass.trainAccuracy)\
                              + ',' + str(myPWClass.T) + ',' + str(myPWClass.k) + ',' \
                              + str(myPWClass.b) + ','\
                              + str(myPWClass.trainMistakes) + ',' + str(t + 1) + '\n')

                '''
                readGene (thePath, gene, kmer, maxReadSize, "pairwise", myPRClass, t)

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
            '''

            logging.info('training Complete')

            myPWClass.saveWeight()

            readGene (thePath, geneTest, kmer, maxReadSize, "testWeight", myPWClass, t, logging)
            ### report results of testing

            myPWClass.testAccuracy = 100 - (100 * (myPWClass.testMistakes / float(myPWClass.testTotal)) )
            #print 'testMistakes = ', myPWClass.testMistakes, '   : test success = ', \
                  #myPWClass.testAccuracy, '%'

            if myPWClass.testdp > 0:
                myPWClass.testPrecision = myPWClass.testnpr / float(myPWClass.testdp)
            myPWClass.testRecall = myPWClass.testnpr / float(myPWClass.testGood)
            #print 'testPrecision: ', myPWClass.testPrecision
            #print 'testRecall: ', myPWClass.testRecall
            if (myPWClass.testPrecision + myPWClass.testRecall) > 0:
                myPWClass.testF1 = 2 * myPWClass.testPrecision * myPWClass.testRecall / (myPWClass.testPrecision + myPWClass.testRecall)

            outFile.write('test' + ',' + str(myPWClass.testPrecision) + ',' + str(myPWClass.testRecall) \
                          + ',' + str(myPWClass.testF1) + ',' + str(myPWClass.testAccuracy)\
                          + ',' + str(myPWClass.T) + ',' + str(myPWClass.k) + ',' \
                          + str(myPWClass.b) + ','\
                          + str(myPWClass.testMistakes) + ',' + str(t + 1) + '\n')
            logging.info('test' + ',' + str(myPWClass.testPrecision) + ',' + str(myPWClass.testRecall) \
                          + ',' + str(myPWClass.testF1) + ',' + str(myPWClass.testAccuracy)\
                          + ',' + str(myPWClass.T) + ',' + str(myPWClass.k) + ',' \
                          + str(myPWClass.b) + ','\
                          + str(myPWClass.testMistakes) + ',' + str(t + 1) + '\n')

            '''
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
            '''

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
