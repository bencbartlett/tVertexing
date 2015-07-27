# Quick and dirty multiprocessing queue to process batch jobs by energy on lxplus
# Ben Bartlett
# 23 Jul 2015

import os, sys
from multiprocessing import Pool

# Change these as needed
processString = "/store/cmst3/group/hgcal/CMSSW/GluGluHtoGG_CMSSW_6_2_0_SLHC26_patch2/RECO-PU0-HLLHC_Fix-FTnosmear_v2/Events_211_25_"
#outdirectory = "~/work/public/Events_22_AllEnergies/Events_22_"
outdirectory = "/tmp/Hgg_All/"
outfile = "Events_211_"
cuts = [1,2,3,4,5,6,7,8,9] # Fill with energy levels to be processed

def process(i):
	print "Processing group: %i" % i
	command =  "(cmsRun ../EDMProcessor.py " + processString + str(i) + "* " + outdirectory + outfile + str(i) + \
					".root) 2>&1|tee /tmp/Hgg_All/output"+str(i)+".txt"
	print command
	os.system(command)

pool = Pool(processes=len(cuts))
pool.map(process, cuts)


