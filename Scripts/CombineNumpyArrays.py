# Numpy array combiner

import os, sys
import numpy as np 

# combineList = sys.argv[2:]
combineDir = sys.argv[1]
path = os.path.expanduser(combineDir)
combineList = [path+f for f in os.listdir(path) if str(f).endswith(".npy")]
combinedArray = np.load(combineList[0])
for i in range(1, len(combineList)):
	print "Appending array %s..." % combineList[i]
	newData = np.load(combineList[i])
	if newData.size > 0:
		combinedArray = np.concatenate((combinedArray, newData))
	else:
		print "Array is empty, skipping."

print "Writing combined array..."
np.save("CombinedArray.npy", combinedArray) 