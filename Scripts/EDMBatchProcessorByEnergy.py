# Quick and dirty multiprocessing queue to process batch jobs by energy on lxplus
# Ben Bartlett
# 23 Jul 2015

import os, sys
from multiprocessing import Pool

# Change these as needed
processString = "/store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC26_patch2/RECO-PU0-FTnosmear_v2/Events_22_"
#outdirectory = "~/work/public/Events_22_AllEnergies/Events_22_"
outdirectory = "/tmp/Events_22_All/"
outfile = "Events_22_"
energies = [100,10,125,175,20,250,2,3,400,40,500,50,5,75,8] # Fill with energy levels to be processed

# Saves in the form:
# /store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC26_patch2/RECO-PU0-FTnosmear_v2/Events_22_500_1.root
# /store/cmst3/group/hgcal/CMSSW/Single22_CMSSW_6_2_0_SLHC26_patch2/RECO-PU0-FTnosmear_v2/Events_22_500_2.root
# etc.
# goes into:
# ~/work/public/Events_22_AllEnergies/Events_22_500.root

def process(i):
	print "Processing energy level: %i" % i
	os.system("cmsRun ../EDMProcessor.py " + processString + str(i) + "_* " + outdirectory + outfile + str(i) + ".root")

pool = Pool(processes=len(energies))
pool.map(process, energies)


