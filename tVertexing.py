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
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import warnings
from Libraries.FastProgressBar import progressbar
from scipy.stats import norm

warnings.simplefilter("error") 

# Variable declarations:
c          = 29.9792458                                                                             #|Lightspeed in cm/ns
ECalRadius = 129.0                                                                                  #|Radius of ECal in cm
ECalZ      = 317.0                                                                                  #|Half-length of ECal in cm (from - to + this value)
etaHGC     = 1.479                                                                                  #|Eta where the HGC starts

def timingCorrection():
    '''Similar to cedric's timing correction stuff for time of arrival. Might not be necessary.'''
    pass

def minimizeVertex():
    '''Runs a minimization procedure similar to the one in Sepehr's algorithm on a set of two
    clusters to try to estimate the interaction vertex. This takes the two arrival times of the
    clusters as arguments as computed in the getArrival() function.'''
    pass

def timeSmearing(data, sigma):
    '''Simulates time smearing in the arrival times of the data. Takes data and smearing sigma
    in picoseconds (not nanoseconds); returns smeared data set.'''
    # Check for sigma=0 and do nothing
    if sigma == 0: return data
    # Loop over events
    for event in data:
        event[0]['tofT'] = np.random.normal(event[0]['tofT'], float(sigma)/1000)                    #|Smear time-of-flight corrected data along a Gaussian
    return data

def twoVertex(hit1, hit2):
    '''Simple vertex finder for a two-photon system. Inherits first interaction xyzt location.'''
    x1, x2 = hit1[0], hit2[0]
    y1, y2 = hit1[1], hit2[1]
    z1, z2 = hit1[2], hit2[2]
    t1, t2 = hit1[3], hit2[3]
    r1     = np.sqrt(x1**2 + y1**2)                                                                 #|Super-explicit form rendered from Mathematica's FullSimplify@Solve[] function
    r2     = np.sqrt(x2**2 + y2**2)
    soln1  = -1 * ((np.sqrt(c**2 * (t1-t2)**2 * (r1**4 - 2 * r1**2 * (r2**2 + c**2 * (t1-t2)**2 - \
        (z1-z2)**2) + (r2**2 - c**2 * (t1-t2)**2 + (z1-z2)**2)**2)) - c**2 * (t1-t2)**2 * (z1+z2) \
        + (z1-z2) * (r1**2 - r2**2 + z1**2 - z2**2)) / (2 * (c * (t1-t2) + z1 - z2) * \
        (c * (t1-t2) - z1 + z2)))

    soln2  = ((np.sqrt(c**2 * (t1-t2)**2 * (r1**4 - 2*r1**2 * (r2**2 + c**2 * (t1-t2)**2 - \
        (z1-z2)**2) + (r2**2 - c**2 * (t1-t2)**2 + (z1-z2)**2)**2)) + c**2 * (t1-t2)**2 * (z1+z2) \
        - (z1-z2) * (r1**2 - r2**2 + z1**2 - z2**2)) / (2 * (c * (t1-t2) + z1 - z2) * \
        (c * (t1-t2) - z1 + z2)))

    error1 = np.absolute(c*(t1-t2) - np.sqrt(r1**2 + (z1-soln1)**2) \
                                   + np.sqrt(r2**2 + (z2-soln1)**2))                                #|This correctness function compares the goodness-of-fit of each solution with the actual measurements and picks the better one.   
    error2 = np.absolute(c*(t1-t2) - np.sqrt(r1**2 + (z1-soln2)**2) \
                                   + np.sqrt(r2**2 + (z2-soln2)**2))
    if error1 <= error2:
        solnZ = soln1 
        error = error1
    else:
        solnZ = soln2
        error = error2
    solnT = np.mean([t1 - np.sqrt(x1**2 + y1**2 + (z1-solnZ)**2)/c,
                     t2 - np.sqrt(x2**2 + y2**2 + (z2-solnZ)**2)/c])
    return [solnZ, solnT, soln1, soln2, error]

def getClusterArrivalTimes(RecHits, clusterID):
    '''
    Computes the energy-weighted arrival times of a given cluster (not given in tree).
    Inefficient, I know (you should make it compute all at once instead of calling it several
    times), but it keeps the code simple, which is more important.
    '''
    # This should be energy weighted once Lindsey fixes the energy on the datasets
    extractMask     = np.logical_and.reduce((RecHits['t']>0, 
                                             RecHits['clusterID']==clusterID,
                                             RecHits['isIn3x3']))                                   #|Extraction mask for cluster t and xyz processing; decides initially which hits to be counted
    arrivalTimes    = np.extract(extractMask, RecHits['tofT'])
    arrivalEnergies = np.extract(extractMask, RecHits['en'])
    return np.mean(arrivalTimes)
    #return np.dot(arrivalTimes,arrivalEnergies)/np.sum(arrivalEnergies)

def getClusterXYZ(RecHits, clusterID):
    '''
    Computes the log-energy-weighted xyzt coordinates of a given cluster. The energy is initially
    reduced by a threshold amount to discard noisy hits, and then log'd to account for energy
    collection fluctuations in the material.
    '''
    extractMask = np.logical_and.reduce((RecHits['t'] > 0, 
                                         RecHits['clusterID'] == clusterID,
                                         RecHits['z'] < 350,
                                         RecHits['isIn3x3']))                                       #|Extraction mask for cluster t and xyz processing; decides initially which hits to be counted
    newRecHits = np.compress(extractMask, RecHits)
    # minLayer = np.min(newRecHits['layerID'])
    # newRecHits = np.compress(newRecHits['layerID'] == minLayer, newRecHits)
    # x = np.extract(extractMask, RecHits['x'])
    # y = np.extract(extractMask, RecHits['y'])
    # z = np.extract(extractMask, RecHits['z'])
    # E = np.extract(extractMask, RecHits['en'])
    x = newRecHits['x']
    y = newRecHits['y']
    z = newRecHits['z']
    E = newRecHits['en']

    return np.mean(x), np.mean(y), np.mean(z)
    logE = E/250.0                                                                                  #|Set energy threshold, in this case, 30MIPs
    x = np.extract(E > 250, x)
    y = np.extract(E > 250, y)
    z = np.extract(E > 250, z)
    E = np.extract(E > 250, E)
    logE = np.log(np.extract(logE > 1, logE))
    # a = zip(logE,x,y,z)
    # a.sort(key=lambda val: abs(val[3]))
    # for thing in a:
    #     print thing
    return np.dot(x,logE)/np.sum(logE), np.dot(y,logE)/np.sum(logE), np.mean(z)
    #return np.dot(x,logE)/np.sum(logE), np.dot(y,logE)/np.sum(logE), np.dot(z,logE)/np.sum(logE)


def clusterFilters(clusters, quiet=False):
    pass

# def gammaGunTests(clust):
#     '''
#     Tests for gamma gun data set. Returns False if failed, else True.
#     '''
#     if clust['ROI'][0] == clust['ROI'][1]:                                                        #|Skips if the two most energetic clusters are in the same ROI, can cause problems
#         print ">>Most energetic clusters in same ROI! Skipping this set..."
#         return False
#     else:
#         return True

class fourVector(object):                                                                           #|Energy-momentum four-vector for Higgs filter
    def __init__(self,E,px,py,pz):
        self.__E  = E
        self.__px = px
        self.__py = py
        self.__pz = pz
    def dot(self, other):
        return self.__E*other.__E - (self.__px*other.__px+self.__py*other.__py+self.__pz*other.__pz)
    def __add__(self, other):
        self.__E  += other.__E
        self.__px += other.__px
        self.__py += other.__py
        self.__pz += other.__pz
        return self

def HiggsInvariantMassFilter(event):
    # Extract data
    eventpts    = np.sort(event[2], order='pt')
    pt1, pt2    = eventpts['pt'][0:2]
    eta1, eta2  = eventpts['eta'][0:2]
    phi1, phi2  = eventpts['phi'][0:2]
    # Calculate needed parameters
    E1, E2      = pt1*np.cosh(eta1), pt2*np.cosh(eta2)
    px1, px2    = pt1*np.cos(phi1), pt2*np.cos(phi2)
    py1, py2    = pt1*np.sin(phi1), pt2*np.sin(phi2)
    pz1, pz2    = E1*np.tanh(eta1), E2*np.tanh(eta2)
    # Calculate invariant mass
    eventVector = fourVector(E1,px1,py1,pz1) + fourVector(E2,px2,py2,pz2)
    invmass     = np.sqrt(eventVector.dot(eventVector))
    # Return if within 10GeV of Higgs mass
    return invmass > 115 and invmass < 135                                                          #|Filters for Higgs mass range


def HggFilter(data, quiet=True):
    '''
    Filters Hgg data set to clean it up some, removing some of the crap.
    '''
    if not quiet: print "Filtering"
    eventpTs = [event[2]['pt'] for event in data]                                                   #|Get ROI pT in descending order    
    twoROIs  = np.compress([len(pT) >= 2 for pT in eventpTs], data, axis=0)                         #|Get at least two clusters    
    pTs      = [np.sort(event[2]['pt'])[::-1] for event in twoROIs]                                 #|Get new sorted ROI pT from refined list
    data     = np.compress([pT[1] > 20 for pT in pTs], twoROIs, axis=0)                             #|Second biggest cluster energy must also have 20GeV 
    data     = np.compress([HiggsInvariantMassFilter(event) for event in data], data, axis=0)       #|Pick events with invariant mass near 125GeV

    return data


def tVertexData(data, pbar=False, runNumber=False, plotData=False, quiet=False, returnDiffs=True):
    '''
    Script for running through a large number of GammaGun events.
    '''
    i                 = 0
    incorrectCount    = 0                                                                           #|Number of incorrect vertex predictions
    tVertexZList      = []                                                                          #|List of calculated vertices
    genVertexZList    = []                                                                          #|List of actual vertices
    correctVertexList = []
    ldata             = len(data)
    nAttempted        = len(data)
    for event in data:
        try: 
            # Get data, run tests
            ROIs         = np.sort(event[2], order='pt')[::-1]                                      #|Sort ROIs by decreasing pt
            cluster0     = np.sort(np.compress(event[1]['ROI'] == ROIs['roiID'][0], event[1]),
                                   order='en')[::-1]                                                #|Gets highest energy cluster of highest-pT ROI
            cluster1     = np.sort(np.compress(event[1]['ROI'] == ROIs['roiID'][1], event[1]),
                                   order='en')[::-1]

            # Parse data 
            c0ID         = cluster0['clusterID'][0] 
            c1ID         = cluster1['clusterID'][0] 
            c0t          = getClusterArrivalTimes(event[0], c0ID)                                   #|Get flat-averaged arrival times
            c1t          = getClusterArrivalTimes(event[0], c1ID)
            # c0x,c0y,c0z  = cluster0['centerX'][0], cluster0['centerY'][0], cluster0['centerZ'][0]
            # c1x,c1y,c1z  = cluster1['centerX'][0], cluster1['centerY'][0], cluster1['centerZ'][0]
            c0x,c0y,c0z  = getClusterXYZ(event[0], c0ID)                                            #|Get log-energy weighted arrival times
            c1x,c1y,c1z  = getClusterXYZ(event[0], c1ID)
            cluster0Hits = [c0x, c0y, c0z, c0t]                                                     #|Package together for twoVertex() function
            cluster1Hits = [c1x, c1y, c1z, c1t]
            genVertexZ   = event[4]['z'][0]

            # Do calculations
            vertices     = twoVertex(cluster0Hits, cluster1Hits)
            tVertexZ     = vertices[0]                                                              #|Optimal solution 
            soln1, soln2 = vertices[2:4]                                                            #|Pass both solutions to compare
            error        = vertices[4] 
            if np.absolute(tVertexZ) > 15.0:                                                        #|Automatically discard any fishy looking solutions
                nAttempted -= 1
                if not quiet: print ">> Skipped!"
                continue                                                                            #|Get rid of the ~3sigma solutions that are probably errors.

            # Append to list and compare
            errorZ = np.absolute(tVertexZ-genVertexZ)
            tVertexZList.append(tVertexZ)
            genVertexZList.append(genVertexZ)

            # Build output string
            if not quiet:
                string = "Event: %3i  |  NumClusters: %2i  |  NumROIs: %2i  |  " \
                        % (i,len(event[1]), len(event[2])) 
                string += "tVertexZ: {0:10.5f}cm  |  genVertexZ:{1:10.5f}cm  |  Error:{2:>9.5f}cm "\
                        .format(tVertexZ, genVertexZ, errorZ)
                if np.absolute(tVertexZ-genVertexZ) < .001: string += "<"+"-"*15                    #|Points out very good fits
                elif np.absolute(tVertexZ-genVertexZ) < .01: string += "<"+"-"*7

            # Test for incorrectly identified vertices
            if np.absolute(soln1-genVertexZ) < errorZ:
                incorrectCount += 1
                if not quiet:
                    string += "<-- Incorrect choice of point. Correct: %.5fcm with error %.5fcm." \
                            % (soln1, np.absolute(soln1-genVertexZ))
            elif np.absolute(soln2-genVertexZ) < errorZ:
                incorrectCount +=1
                if not quiet:
                    string += "<-- Incorrect choice of point. Correct: %.5fcm with error %.5fcm." \
                            % (soln2, np.absolute(soln2-genVertexZ))
            else:
                correctVertexList.append(tVertexZ-genVertexZ)                                       #|Pick out vertices the algorithm correctly identified

        except RuntimeWarning:            
            string     = ">> Runtime warning, skipping this event..."
            nAttempted -= 1

        except IndexError:
            string     = "Event: %3i   "%i + "-"*15 + "Two clusters not found" + "-"*15
            nAttempted -= 1

        if not quiet: print string
        i += 1

        if runNumber:
            pbar.update(ldata*runNumber + i)

    # Do some manipulation to the lists and possibly clean it up to get a fit
    diffs        = np.subtract(tVertexZList, genVertexZList)
    oldlen       = len(diffs)
    worstPercent = np.percentile(np.absolute(diffs),100-3)                                          #|Extract magnitude of worst estimate
    diffs        = np.extract(np.absolute(diffs) < worstPercent, diffs)                             #|Pull out all the stuff past 2mm error
    
    if not quiet: 
        print "Correctly identified %i/%i events." % (nAttempted-incorrectCount, nAttempted)
        print "%s events trimmed" % str(oldlen-len(diffs))

    if plotData:
        Plotter.tVertexErrorHist(10*diffs, len(diffs), ranges=[-5,5], quiet=False)

    if returnDiffs:
        return diffs 



if __name__ == '__main__':
    print "Loading data..."
    # hits = np.core.records.fromarrays(np.transpose([[ECalRadius+1, 0.0, 250.0, 0.0],
    #                                 [0.0, ECalRadius+1, -150.0, -2.5]]),
    #                                 names = 'x,y,z,t')
    # Plotter.vertexPlot(twoVertex(hits), hits)

    # data = np.load("Data/GammaGunProcessed22_100_11.npy")
    # Plotter.showerAnimator(data[2][0], 2, "GammaGun 2", clusterID=1, projections=True, frames=50,
    #                                                     transparency=False, endLoop = 30, 
    #                                                     delete=False)
    
    data  = np.load("Data/500GeVPhoton.npy")
    res = 50
    diffsList = np.array([])
    smearingValues = np.array([])
    num = 0
    pbar = progressbar("Computing &count&:", res*len(data))
    pbar.start()
    for smearingVal in np.linspace(0,50,res):
        data  = np.load("Data/500GeVPhoton.npy")                                                    #|Quick and dirty reload, numpy was doing funky shit
        diffs = tVertexData(timeSmearing(data, smearingVal), runNumber=num, pbar=pbar, quiet=True)
        diffsList = np.append(diffsList, diffs)
        smearingValues = np.append(smearingValues, np.repeat(smearingVal, len(diffs)))    
        num += 1
    pbar.finish()
    np.save("diffs.npy", np.multiply(10,diffsList))
    np.save("vals.npy", smearingValues)
    Plotter.tVertexErrorHist2D(np.multiply(10,diffsList), smearingValues)









