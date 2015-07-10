#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #|Comments Section:
#                                                                                                   #|
#     Ben Bartlett                                     File description:                            #|For the sake of any poor person having to read my code, comments starting with "#|" will be aligned.
#     California Institute of Technology                 Master file. Implements the                #|
#     bartlett@caltech.edu                               actual algorithm.                          #|
#     ben.bartlett@cern.ch                                                                          #|
#     benjamincbartlett@gmail.com                      Notes:                                       #|
#                                                        See dedicated notes section.               #|
# Created:       3 June 2015                                                                        #|
# Last modified: 3 June 2015                                                                        #|
#####################################################################################################|
# Notes:
# ->If numpy does not import (and you have it installed) it is likely that you are still in the
#	reassignment of python required by CMSSW locally. Quit the terminal and rerun.
#


# Importation
import Plotter
import numpy as np


# Variable declarations:
c          = 29.9792458                                                                             #|Lightspeed in cm/ns
ECalRadius = 129.0                                                                                 #|Radius of ECal in cm
ECalZ      = 317.0                                                                                 #|Half-length of ECal in cm (from - to + this value)

def timingCorrection():
	'''similar to cedric's timing correction stuff for time of arrival. Might not be necessary.'''
	pass

def minimizeVertex():
	'''Runs a minimization procedure similar to the one in Sepehr's algorithm on a set of two
	clusters to try to estimate the interaction vertex. This takes the two arrival times of the
	clusters as arguments as computed in the getArrival() function.'''
	pass

def twoVertex(hits):
	'''Simple vertex finder for a two-photon system. Inherits first interaction xyzt location.'''
	x1,y1,z1,t1 = hits[0]
	x2,y2,z2,t2 = hits[1]
	r1          = np.sqrt(x1**2 + y1**2)                                                            #|Super-explicit form rendered from Mathematica's FullSimplify@Solve[] function
	r2          = np.sqrt(x2**2 + y2**2)
	soln1 = -1 * ((np.sqrt(c**2 * (t1-t2)**2 * (r1**4-2 * r1**2 * (r2**2 + c**2 * (t1-t2)**2 - \
		(z1-z2)**2) + (r2**2 - c**2 * (t1-t2)**2 + (z1-z2)**2)**2)) - c**2 * (t1-t2)**2 * (z1+z2) \
		+ (z1-z2) * (r1**2 - r2**2 + z1**2 - z2**2)) / (2 * (c * (t1-t2) + z1 - z2) * \
		(c * (t1-t2) - z1 + z2)))
	soln2 = ((np.sqrt(c**2 * (t1-t2)**2 * (r1**4 - 2*r1**2 * (r2**2 + c**2 * (t1-t2)**2 - \
		(z1-z2)**2) + (r2**2 - c**2 * (t1-t2)**2 + (z1-z2)**2)**2)) + c**2 * (t1-t2)**2 * (z1+z2) \
		- (z1-z2) * (r1**2 - r2**2 + z1**2 - z2**2)) / (2 * (c * (t1-t2) + z1 - z2) * \
		(c * (t1-t2) - z1 + z2)))
	solnZ = min(soln1, soln2, key=abs)                                                              #|Return the solution closer to the middle of the beam
	solnT = np.mean([t1 - np.sqrt(x1**2 + y1**2 + (z1-solnZ)**2)/c,
					 t2 - np.sqrt(x2**2 + y2**2 + (z2-solnZ)**2)/c])
	return [solnZ, solnT]

def getArrival():
	'''Computes the xyzt coordinates representing the locaiton of first interaction in a cluster.'''
	pass




if __name__ == '__main__':
    # hits = [[ECalRadius+1, 0.0, 2500.0, 0.0],
    # 		[0.0, ECalRadius+1, 1500.0, -2.5]]
    # Plotter.vertexPlot(twoVertex(hits), hits)
    data = np.load("Data/GammaGunProcessed22_100_11.npy")
    Plotter.showerAnimator(data[2], "GammaGun", 2, clusterID = 1)











