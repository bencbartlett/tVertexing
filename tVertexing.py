#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #|Comments Section:
#                                                                                                   #|
#     Ben Bartlett                                     File description:                            #|For the sake of any poor person having to read my code, comments starting with "#|" will be aligned.
#     California Institute of Technology                 Master file. Implements the                #|
#     bartlett@caltech.edu                               actual algorithm.                          #|
#     ben.bartlett@cern.ch                                                                          #|
#     benjamincbartlett@gmail.com                      Notes:                                       #|
#                                                        See dedicated notes section.               #|
# Created:       3 July 2015                                                                        #|
# Last modified: 3 July 2015                                                                        #|
#####################################################################################################|
# Notes:
# ->If numpy does not import (and you have it installed) it is likely that you are still in the
#   reassignment of python required by CMSSW locally. Quit the terminal and rerun.
#


# Importation
import Plotter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import warnings
from scipy.stats import norm

warnings.simplefilter("error")


# Variable declarations:
c          = 29.9792458                                                                             #|Lightspeed in cm/ns
ECalRadius = 129.0                                                                                  #|Radius of ECal in cm
ECalZ      = 317.0                                                                                  #|Half-length of ECal in cm (from - to + this value)

def timingCorrection():
    '''similar to cedric's timing correction stuff for time of arrival. Might not be necessary.'''
    pass

def minimizeVertex():
    '''Runs a minimization procedure similar to the one in Sepehr's algorithm on a set of two
    clusters to try to estimate the interaction vertex. This takes the two arrival times of the
    clusters as arguments as computed in the getArrival() function.'''
    pass

def twoVertex(hit1, hit2):
    '''Simple vertex finder for a two-photon system. Inherits first interaction xyzt location.'''
    x1, x2 = hit1[0], hit2[0]
    y1, y2 = hit1[1], hit2[1]
    z1, z2 = hit1[2], hit2[2]
    t1, t2 = hit1[3], hit2[3]
    r1     = np.sqrt(x1**2 + y1**2)                                                                 #|Super-explicit form rendered from Mathematica's FullSimplify@Solve[] function
    r2     = np.sqrt(x2**2 + y2**2)
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

def getClusterArrivalTimes(RecHits, clusterID):
    '''
    Computes the energy-weighted arrival times of a given cluster (not given in tree).
    Inefficient, I know (you should make it compute all at once instead of calling it several
    times), but it keeps the code simple, which is more important.
    '''
    # This should be energy weighted once Lindsey fixes the energy on the datasets
    arrivalTimes = np.extract(np.logical_and(RecHits['t']>0, 
                              RecHits['clusterID']==clusterID), RecHits)['t']
    return np.mean(arrivalTimes)

def clusterFilters(clusters, quiet=False):
    pass
    

if __name__ == '__main__':
    # hits = np.core.records.fromarrays(np.transpose([[ECalRadius+1, 0.0, 250.0, 0.0],
    #                                 [0.0, ECalRadius+1, -150.0, -2.5]]),
    #                                 names = 'x,y,z,t')
    # Plotter.vertexPlot(twoVertex(hits), hits)

    # data = np.load("Data/GammaGunProcessed22_100_11.npy")
    # Plotter.showerAnimator(data[2][0], "GammaGun", 2, projections=False, delete=True)

    # eventArray[0]=RecHits
    # eventArray[1]=Clusters
    # eventArray[2]=ROIs
    # eventArray[2]=Vertices
    # eventArray[3]=GenVertex


    data           = np.load("Data/GammaGunProcessed22_100_11.npy")
    i              = 0
    tVertexZList   = []
    genVertexZList = []
    for event in data:
        try: 
            clust        = np.sort(event[1], order='en')[::-1]                                      #|Use two most energetic clusters
            c0ID         = clust['clusterID'][0]
            c1ID         = clust['clusterID'][1] 
            if clust['ROI'][0] == clust['ROI'][1]:
                print ">>Most energetic clusters in same ROI! Skipping this set..."
                continue   
            #if np.absolute(np.absolute(clust['phi'][0]-clust['phi'][1])-np.pi)
            c0t          = getClusterArrivalTimes(event[0], c0ID)
            c1t          = getClusterArrivalTimes(event[0], c1ID)
            cluster0Hits = [clust['centerX'][0], clust['centerY'][0], clust['centerZ'][0], c0t]
            cluster1Hits = [clust['centerX'][1], clust['centerY'][1], clust['centerZ'][1], c1t]
            tVertexZ     = twoVertex(cluster0Hits, cluster1Hits)[0]
            genVertexZ   = event[0,4]['z'][0]
            tVertexZList.append(tVertexZ)
            genVertexZList.append(genVertexZ)
            errorZ = np.absolute(tVertexZ-genVertexZ)
            string = "Event: %3i   NumClusters: %2i  NumROIs: " % (i,len(event[1]))
            string += "tVertexZ: {0:12.5f}   genVertexZ:{1:12.5f}   error{2:>12.5f}"\
                                                    .format(tVertexZ, genVertexZ, errorZ)
            if np.absolute(tVertexZ-genVertexZ) < .1:
                string += "<"+"-"*15
            elif np.absolute(tVertexZ-genVertexZ) < .5:
                string += "<"+"-"*7
            elif np.absolute(tVertexZ-genVertexZ) < 1:
                string += "<"+"-"*2

        except RuntimeWarning:
            print ">>Runtime warning, skipping this event..."

        except IndexError:
            string = "Event: %3i   "%i + "-"*15 + "Two clusters not found" + "-"*15

        print string
        i += 1

    diffs = np.subtract(tVertexZList, genVertexZList)
    diffs = np.extract(diffs < 20, diffs)
    diffs = np.extract(diffs > -20, diffs)
    n, bins, patches = plt.hist(diffs, normed = True, bins = 150)
    (mu, sigma) = norm.fit(diffs)
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, "r--", linewidth = 2)
    plt.xlabel("Error (cm)") 
    plt.ylabel("Counts (normalized)")
    plt.show()
    diffs = np.absolute(diffs)
    print "Mu = %f,  sigma = %f" % (mu, sigma)
    print "Average error magnitude: %.4f \n" % np.mean(diffs)








