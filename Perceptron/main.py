##########################################
### Perceptron Test
### Marcus Blaisdell
##########################################

from p_Class import *
#import time

myClass = p_Class()

myClass.loadData("/home/marcus/Data/Genomic-Data/training_Set-3", myClass.trainData)
#myClass.loadData("/home/marcus/Data/Genomic-Data/shortEDL-test.txt", myClass.validationData)
myClass.loadData("/home/marcus/Data/Genomic-Data/testing_Set-3", myClass.testData)
'''
myClass.loadData("/home/marcus/Data/Genomic-Data/training", myClass.trainData)
myClass.loadData("/home/marcus/Data/Genomic-Data/testing", myClass.testData)
'''
myClass.countGood()

outFile = open("/home/marcus/Documents/LinuxShare/Results/Perceptron/results_Set-3.csv", "w")

outFile.write('type' + ',' + 'Precision' + ',' + 'Recall' + ',' + 'F1' + ',' + \
              'Accuracy' + ',' + 'T' + ',' + 'K' + ',' \
              + 'Mistakes' + ',' + 'Iteration' + '\n')

startTime = time.time()
myClass.perceptron(outFile)
endTime = time.time()
print 'runTime: ', endTime - startTime

outFile.close()
