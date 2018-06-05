##########################################
### PEGASOS versus Perceptron Test
### for genomic data
### Marcus Blaisdell
##########################################

from pg_Class import *
from pr_Class import *
from functions import *

myPGClass = pg_Class()
myPRClass = pr_Class()

### create arrays with all the training, testing files:

trainArray = ["11a-train", "11b-train", "11c-train", "11d-train", \
                "13a-train", "13b-train", "13c-train", "13d-train", \
                "15a-train", "15b-train", "15c-train", "15d-train", \
                "17a-train", "17b-train", "17c-train", "17d-train" ]

testArray = ["11a-test", "11b-test", "11c-test", "11d-test", \
                "13a-test", "13b-test", "13c-test", "13d-test", \
                "15a-test", "15b-test", "15c-test", "15d-test", \
                "17a-test", "17b-test", "17c-test", "17d-test" ]

resultsArray = ["11a", "11b", "11c", "11d", \
                "13a", "13b", "13c", "13d", \
                "15a", "15b", "15c", "15d", \
                "17a", "17b", "17c", "17d" ]

### run the comparison for each training/testing set:

thePath = "/home/marcus/Data/Genomic-Data/Testing/"
resultsStringPre = "results_"
resultsStringPost = ".csv"

for index in range(len(trainArray)):
    trainData = thePath + trainArray[index]
    testData = thePath + testArray[index]
    resultsFile = resultsStringPre + resultsArray[index] + resultsStringPost

    print 'set: ', resultsArray[index]

    print 'loading data for perceptron:'
    loadData(trainData, myPRClass.trainData)
    loadData(testData, myPRClass.testData)
    print 'load complete'
    myPRClass.countGood()

    print 'opening results file:'
    outFile = open(resultsFile, "w")

    outFile.write ('Perceptron'  + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running perceptron . . . '

    startTime = time.time()
    myPRClass.perceptron(outFile)
    endTime = time.time()
    print 'Perceptron runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))

    ### end test perceptron

    print 'loading data for Pegasos . . . '

    loadData(trainData, myPGClass.trainData)
    loadData(testData, myPGClass.testData)

    print 'load complete'

    myPGClass.countGood()


    ### Add a header for Pegasos:
    outFile.write ('\n' + '\n' + 'Pegasos' + '\n' + '\n')
    outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
                  'Accuracy' + ',' + 'T' + ',' + 'K' + ',' + 'Lamda' + ',' \
                  + 'Mistakes' + ',' + 'Iteration' + '\n')

    print 'running Pegasos . . .'

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

                # initialize weight vector, w to zero:

                kmerSize = ((4**myPGClass.kmer) + 1)
                myPGClass.w = np.zeros(kmerSize, dtype=float)
                #myPGClass.w = np.zeros(10000000, dtype = float)

                ### repeat for T iterations:
                for t in range(myPGClass.T):
                    myPGClass.pegasosBatch(t)

                    outFile.write('train' + ',' + str(myPGClass.trainPrecision) + ',' + str(myPGClass.trainRecall) \
                                  + ',' + str(myPGClass.trainF1) + ',' + str(myPGClass.trainAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.trainMistakes) + ',' + str(t + 1) + '\n')

                    myPGClass.testWeight()
                    outFile.write('test' + ',' + str(myPGClass.testPrecision) + ',' + str(myPGClass.testRecall) \
                                  + ',' + str(myPGClass.testF1) + ',' + str(myPGClass.testAccuracy)\
                                  + ',' + str(myPGClass.T) + ',' + str(myPGClass.k) + ',' \
                                  + str(myPGClass.lam) + ',' + str(myPGClass.testMistakes) + '\n')


    endTime = time.time()
    print 'Pegasos runTime: ', endTime - startTime
    outFile.write (str(endTime - startTime))

    print 'closing output file . . .'

    outFile.close()

print '*** program complete ***'
