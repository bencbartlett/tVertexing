#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #
#                                                                                                   #
#     Ben Bartlett                                     File description:                            #
#     California Institute of Technology                 ROOT file processor: takes .root files     #
#     bartlett@caltech.edu                               and converts them to a saved numpy array.  #
#     ben.bartlett@cern.ch                               (since numpy is easier to work with)       #
#     benjamincbartlett@gmail.com                      Notes:                                       #
#                                                        After some fidgeting, must be run on       #
# Created:       6 June 2015                             a CERN lxplus machine, as it requires      #
# Last modified: 8 June 2015 			     			 libFWCoreFWLite.so and AutoLibraryLoader   #
#####################################################################################################

'''
Usage (in cmsenv): RootProcessor.py <[file or directory]> <[options]>

See parser for more details.
'''
c = 29.9792458                                                                                      #|Lightspeed in cm/ns
# Argument parser
import argparse
parser = argparse.ArgumentParser(description = 
	'''Converts HGCROIAnalysis-format .root files into RecHit numpy arrays.
	If converting all files in a directory using -f, make sure all files are HGCROI-formatted.''',\
	epilog = '''NOTE: This program must be run while cmsenv is active.

	''')
parser.add_argument("-f", "--folder", \
	help="Process all files in this directory.", action="store_true")
parser.add_argument("input", \
	help="File to process. If -f is used, full directory to process.")
parser.add_argument("-o", "--outdir", \
	help="Directory to save output to. Defaults to Data/.", action="store")
args = parser.parse_args()

# Importations
print "Importing libraries..."
import ROOT
import os, sys, platform
from Misc import progressbar
try: 
	import numpy as np
except ImportError:
	print "Numpy not found. Reverting to known path..."
	sys.path.append("/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/")
	import numpy as np

# Load necessary libraries
if platform.uname()[0] == "Darwin":
	ROOT.gSystem.Load("libFWCoreFWLite.dylib")                                                      #|Provides dictionaries for opening the various slimmed vertex types and such.
elif platform.uname()[0] == "Linux":
	ROOT.gSystem.Load("libFWCoreFWLite.so")
else:
	print "This program requires FWLite running on OSX or Linux." 
	sys.exit()                                                            
ROOT.AutoLibraryLoader.enable()                                                                     

# Main processing function
def process(f, outdir):
	print "Processing file: " + str(f) + "..."                                                                                                                
	outArray   = []                                           
	fIn        = ROOT.TFile.Open(f)
	tree       = fIn.Get('analysis/HGC') 
	numentries = tree.GetEntries() 
	pbar       = progressbar("Processing &count& events:", numentries) 
	pbar.start()                                                
	for i in xrange(0, numentries): 
		# print "-> Processing event "+str(i)+"/"+str(numentries)+"..."
		tree.GetEntry(i) 
		x         = []
		y         = []
		z         = []
		t         = []
		correctT  = []
		en        = []
		clusterID = []                                                          
		for hit in tree.RecHits:                                                                    #|Modify properties you want to extract at will.
			x.append(hit.x_)                                                                        #|Loops over files, extracting the
			y.append(hit.y_)                                                                        #|xyzt data from rechits in each file
			z.append(hit.z_)                                                                        #|and saving each rechit array
			t.append(hit.t_)                                                                        #|to a .npy file for later use.
			correctT.append(hit.t_ + np.sqrt(hit.x_**2 + hit.y_**2 + hit.z_**2)/c)
			en.append(hit.en_)
			clusterID.append(hit.clustId_)
		outArray.append(np.core.records.fromarrays([x,y,z,t,en,correctT,clusterID],
						names = 'x,y,z,t,en,correctT,clusterID'))                                   #|Converts to a 2D numpy structured array indexed by event number
		pbar.update(i)
	filename = str(f[:-5]) + ".npy"
	pbar.finish()
	filename = filename.split("/")[-1] #|Removes directory prefixes
	print "Writing file " + filename + "..."
	filepath = outdir+filename
	np.save(filepath, outArray)


if __name__ == "__main__":
	# Process arguments and do the stuff
	if args.outdir == None:
		outdir = "Data/"
	else:
		outdir = args.outdir
	if args.folder:
		directory = args.input
		files = [f for f in os.listdir(directory) if f.endswith(".root")]                           #|Get root files in directory
		os.chdir(directory)
		for f in files:
			process(f, outdir)
	else:
		process(args.input, outdir)

  








