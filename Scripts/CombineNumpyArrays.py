# Numpy array combiner

import os, sys
import numpy as np 

combineList = sys.argv[2:]
print "Processing array %s..." % sys.argv[1]
combinedArray = np.load(sys.argv[1])
#combinedArray = []
for item in combineList:
	print "Appending array %s..." % item
	newData = np.load(item)
	combinedArray = np.concatenate((combinedArray, newData))
	#combinedArray.append(np.load(item))

print "Writing combined array..."
np.save("CombinedArray.npy", combinedArray)