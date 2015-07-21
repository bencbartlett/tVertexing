#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #
#                                                                                                   #
#     Ben Bartlett                                     File description:                            #
#     California Institute of Technology                 ROOT file processor: takes .root files     #
#     bartlett@caltech.edu                               and converts them to a saved numpy array.  #
#     ben.bartlett@cern.ch                               (since numpy is easier to work with)       #
#     benjamincbartlett@gmail.com                      Notes:                                       #
#                                                        After some fidgeting, must be run on       #
# Created:       6 July 2015                             a CERN lxplus machine, as it requires      #
# Last modified: 8 July 2015 			     			 libFWCoreFWLite.so and AutoLibraryLoader   #
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
from Libraries.FastProgressBar import progressbar 
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
	'''
	Usage: process(fileName, writeDirectory)

	Takes a HGCROI-formatted .root file and converts it to a list of structured numpy arrays.

	The format is as follows:
		List of events
			For each event, array for various data structures (RecHits, Clusters, etc.).
			Unfortunately, there is no way to do recursive recarrays, so we use the following
			organisation method:
			eventArray[0]=RecHits
			eventArray[1]=Clusters
			eventArray[2]=ROIs
			eventArray[2]=Vertices
			eventArray[3]=GenVertex
				For each data structure array, recarray of properties

	For example:
		(Index 0 = Event 0)
			(Index 0 = RecHits)
				'x' -> [1.23, 4.25, ...]
				'y' -> [5.24, 6.42, ...]
				...
				'clusterID' -> [1, 2, 1, 1, ...]
			(Index 1 = Clusters)
				'centerX' -> [3.21, 2.56, ...]
				...
			...
		(Index 1 = Event 1)
			(Index 0 = RecHits)
				...
			(Index 1 = Clusters)
				...
			...
		...

	So, for example, to access the array of x positions in the RecHits of the 7th event,
	you would use:
		xpos = Data[7][0]['x']
	'''
	print "Processing file: " + str(f) + "...\n"                                                                                                                
	outArray   = []                                           
	fIn        = ROOT.TFile.Open(f)
	tree       = fIn.Get('analysis/HGC') 
	numentries = tree.GetEntries() 
	pbar       = progressbar("Processing &count& events:", numentries) 
	pbar.start()                                                
	for i in xrange(0, numentries): 
		tree.GetEntry(i) 
		eventArray = []
		names      = ""
		# RecHits
		RecHits    = True                                                                           #|Store rechits
		x          = []                                                                             #|Position of the rechit
		y          = []
		z          = []
		t          = []
		tofT   	   = []                                                                             #|Time-of-flight corrected time data
		en         = []                                                                             #|Energy 
		clusterID  = []                                                                             #|What cluster the hit belongs to
		detID      = []
		layerID    = [] 
		isIn3x3    = []                                                                             #|Is in 3x3 grid from center of energy
		isIn5x5    = []                                                                             #|Is in 5x5 grid from center of energy
		isIn7x7    = []                                                                             #|Is in 7x7 grid from center of energy

		# Clusters
		Clusters   = True                                                                           #|Store clusters
		center_x_  = []                                                                             #|Center of the cluster (energy weighted)
		center_y_  = []
		center_z_  = []
		axis_x_    = []                                                                             #|Direction the cluster is pointing (looking back at the beamline from the cluster)
		axis_y_    = []
		axis_z_    = []
		ev_1_      = []                                                                             #|Eigenvalues from principal component analysis
		ev_2_      = []
		ev_3_      = []
		clusteren  = []                                                                             #|Energy
		clustereta = []                                                                             #|Eta
		clusterphi = []                                                                             #|Phi
		slfClustID = []                                                                             #|Self-referencing cluster ID
		clusterroi = []                                                                             #|ROI ID for the cluster
		# ROIs
		ROIs       = True 
		roiID      = []                                                                             #|Self-referencing ROI ID
		roieta     = []                                                                             #|Energy-weighted eta
		roiphi     = []                                                                             #|Energy-weighted phi
		roipt      = []                                                                             #|Energy-weighted pt
		roimass    = []
		roiarea    = []
		roigenpt   = []
		roigeneta  = []
		roigenphi  = []
		roigenmass = []
		roigenarea = []
		roistablex = []
		roistabley = []
		roistablez = []
		roistablID = []
		# Vertices
		Vertices   = False                                                                          #|Store vertices, won't work for photon gun
		vertex_x_  = []                                                                             #|Reconstructed vertex location using tracker 
		vertex_y_  = []
		vertex_z_  = []
		# Generated vertices ("true" vertices)
		GenVert    = True                                                                           #|Store generated vertices
		gen_x_     = []                                                                             #|Actual vertex location from simulated event
		gen_y_     = []
		gen_z_     = []

		if RecHits:
			for hit in tree.RecHits:                                                                #|Modify properties you want to extract at will.
				x        .append(hit.x_)                                                            #|Loops over files, extracting the xyzt data from rechits in each file and saving each rechit array to a .npy file for later use.
				y        .append(hit.y_)                                                                        
				z        .append(hit.z_)                                                                        
				t        .append(hit.t_)                                                                        
				tofT     .append(hit.t_ + np.sqrt(hit.x_**2 + hit.y_**2 + hit.z_**2)/c)
				en       .append(hit.en_)
				clusterID.append(hit.clustId_)
				detID    .append(hit.detId_)
				layerID  .append(hit.layerId_)
				isIn3x3  .append(hit.isIn3x3_)
				isIn5x5  .append(hit.isIn5x5_)
				isIn7x7  .append(hit.isIn7x7_)
			recHitsArray = np.core.records.fromarrays([x, y, z, t, en, tofT, clusterID,
													   detID, layerID, isIn3x3, isIn5x5, isIn7x7],
											  names = 'x,y,z,t,en,tofT,clusterID,\
											  		   detID,layerID,isIn3x3,isIn5x5,isIn7x7')      #|Form rechit array
			eventArray.append(recHitsArray)                                                         #|Append to event array
			names += 'RecHits'                                                                      #|Add to names list
		else:
			eventArray.append([])                                                                   #|This is to keep the index of the arrays the same

		clusterindex = 0
		if Clusters:
			for cluster in tree.Clusters:
				center_x_ .append(cluster.center_x_)
				center_y_ .append(cluster.center_y_)
				center_z_ .append(cluster.center_z_)
				axis_x_   .append(cluster.axis_x_)
				axis_y_   .append(cluster.axis_y_)
				axis_z_   .append(cluster.axis_z_)
				ev_1_     .append(cluster.ev_1_)
				ev_2_     .append(cluster.ev_2_)
				ev_3_     .append(cluster.ev_3_)
				clusteren .append(cluster.en_)
				clustereta.append(cluster.eta_)
				clusterphi.append(cluster.phi_)
				slfClustID.append(clusterindex)
				clusterroi.append(cluster.roiidx_)
				clusterindex += 1
			clusterArray = np.core.records.fromarrays([center_x_, center_y_, center_z_,
													   axis_x_, axis_y_, axis_z_,
													   ev_1_, ev_2_, ev_3_, clusteren, clustereta,
													   clusterphi, slfClustID, clusterroi],
													   names = 'centerX,centerY,centerZ,\
													   			axisX,axisY,axisZ,ev1,ev2,ev3,en,\
													   			eta,phi,clusterID,ROI')             #|Form array for clusters
			eventArray.append(clusterArray)                                                         #|Append to event array
			names += ',Clusters'
		else:
			eventArray.append([])                                                                   #|This is to keep the index of the arrays the same

		ROIindex = 0
		if ROIs:
			for ROI in tree.ROIs:
				roiID     .append(ROIindex)
				roipt     .append(ROI.pt_)
				roieta    .append(ROI.eta_)
				roiphi    .append(ROI.phi_)
				roimass   .append(ROI.mass_)
				roiarea   .append(ROI.area_)
				roigenpt  .append(ROI.genpt_)
				roigeneta .append(ROI.geneta_)
				roigenphi .append(ROI.genphi_)
				roigenmass.append(ROI.genmass_)
				roigenarea.append(ROI.genarea_)
				roistablex.append(ROI.stablex_)
				roistabley.append(ROI.stabley_)
				roistablez.append(ROI.stablez_)
				roistablID.append(ROI.stableid_)
				ROIindex +=1
			ROIArray = np.core.records.fromarrays([roiID, roipt, roieta, roiphi, roimass, roiarea,
												   roigenpt, roigeneta, roigenphi, roigenmass,
												   roigenarea, roistablex, roistabley, roistablez,
												   roistablID],
												   names = 'roiID,pt,eta,phi,mass,area,getpt,\
												   			geneta,getphi,genarea,stablex,stabley,\
												   			stablez,stableID')
			eventArray.append(ROIArray)
			names += ',ROIs'
		else:
			eventArray.append([])

		if Vertices:
			for vertex in tree.Vertices:
				vertex_x_.append(vertex.x_)
				vertex_y_.append(vertex.y_)
				vertex_z_.append(vertex.z_)
			vertexArray = np.core.records.fromarrays([vertex_x_, vertex_y_, vertex_z_],
													  names='x,y,z')
			eventArray.append(vertexArray)                                                          #|Vertices array
			names += ',Vertices'
		else:
			eventArray.append([])                                                                   #|This is to keep the index of the arrays the same

		if GenVert:
			gen_x_.append(tree.GenVertex.x())                                                       #|GenVertex is not iterable like the other classes, since there is only one per event.
			gen_y_.append(tree.GenVertex.y())
			gen_z_.append(tree.GenVertex.z())
			genVertexArray = np.core.records.fromarrays([gen_x_, gen_y_, gen_z_],
													  	 names='x,y,z')
			eventArray.append(genVertexArray)                                                       #|Generated vertices array
			names += ',GenVertex'
		else:
			eventArray.append([])                                                                   #|This is to keep the index of the arrays the same

		# Combine arrays for single event and append to outArray
		outArray.append(eventArray)                                                                 #|Converts to a 2D numpy structured array indexed by event number
		pbar.update(i)

	# Finish up and save array to file
	pbar.finish()
	filename = str(f[:-5]) + ".npy"                                                                 #|Replace .root with .npy
	filename = filename.split("/")[-1]                                                              #|Removes directory prefixes
	filepath = outdir+filename
	print "\nWriting file " + os.path.abspath(filepath) + "..."
	np.save(filepath, outArray)
	print "Processing complete.\n"


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

  








