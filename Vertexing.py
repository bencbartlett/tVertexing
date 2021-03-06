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
# ->If numpy does not import locally (and you have it installed) it is likely that you are still in 
#   the reassignment of python required by CMSSW locally. Quit the terminal and rerun.
# ->This program does some funky shit while trying to multiprocess in cmsenv


# Importation
import Plotter
import numpy as np
import warnings
import os, sys
from Libraries.FastProgressBar import progressbar
from multiprocessing import Pool
from ExternalLibraries.blessings import Terminal
from scipy.optimize import minimize
term = Terminal()

warnings.simplefilter("error") 

# Variable declarations:
c          = 29.9792458                                                                             #|Lightspeed in cm/ns
ECalRadius = 129.0                                                                                  #|Radius of ECal in cm
ECalZ      = 317.0                                                                                  #|Half-length of ECal in cm (from - to + this value)
etaHGC     = 1.478                                                                                  #|Eta where the HGC starts
HiggsMass  = 125.09                                                                                 #|Higgs mass in GeV/c^2

class Writer(object):
    """Helper class for making multiprocessing progress bars"""
    def __init__(self, location):
        """Input: location - tuple of ints (x, y), the position of the bar in the terminal"""
        self.location = location
    def write(self, string):
        with term.location(*self.location):
            print(string)


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


def XYZtoEtaPhi(xyz):
    '''Convert xyz coordinates to an eta-phi map'''
    x, y, z = xyz 
    phi = np.arctan2(y, x)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    eta = -np.log(np.tan(theta/2)) 
    return eta, phi


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
    # arrivalEnergies = np.extract(extractMask, RecHits['en'])
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
                                         RecHits['isIn3x3']))                                       #|Extraction mask for cluster t and xyz processing; decides initially which hits to be counted
    #minLayer = np.min(newRecHits['layerID'])
    #newRecHits = np.compress(newRecHits['layerID'] == minLayer, newRecHits)
    x = np.extract(extractMask, RecHits['x'])
    y = np.extract(extractMask, RecHits['y'])
    z = np.extract(extractMask, RecHits['z'])
    E = np.extract(extractMask, RecHits['en'])

    return (np.mean(x), np.mean(y), np.mean(z))                                                     #|Flat average

    # Energy weight
    w0 = 7                                                                                          #|Arbitrary weighting term that works well, as found by Geoffrey Monet
    Eweight = np.maximum(np.log(E/np.sum(E)) + w0, 0)


    xw, yw, zw = (np.dot(x,Eweight)/np.sum(Eweight), 
                  np.dot(y,Eweight)/np.sum(Eweight), 
                  np.dot(z,Eweight)/np.sum(Eweight))                                                #|Weighted average

    return xw, yw, zw


def getClusterPointer(RecHits, center, clusterID, quiet=True):
    '''
    Computes the energy-weighted pointing direction of the cluster, using a different method than
    PCA analysis, which in some sets seems to be inaccurate.
    '''
    extractMask   = np.logical_and.reduce((RecHits['t'] > 0, 
                                           RecHits['clusterID'] == clusterID,
                                           RecHits['isIn3x3']))                                     #|Extraction mask for cluster t and xyz processing; decides initially which hits to be counted
    points        = np.vstack((np.extract(extractMask, RecHits['x']),
                               np.extract(extractMask, RecHits['y']),
                               np.extract(extractMask, RecHits['z']))).T
    centerDist    = np.linalg.norm(center)
    furtherPoints = np.compress([np.linalg.norm(point) >= centerDist for point in points], 
                                 points - center, axis=0)
    closerPoints  = np.compress([np.linalg.norm(point) < centerDist for point in points], 
                                 center - points, axis=0) 
    pointerVector = np.mean(np.concatenate((furtherPoints, closerPoints)), axis=0)
    pointerVector /= np.linalg.norm(pointerVector)
    # "Vertex" the point using the pointing data
    if not quiet:
        pass
    return pointerVector/np.linalg.norm(pointerVector)
    

def pVertexPCA(event, cID):
    '''
    Gets cluster pointing vector from internal axis variables and estimates the vertex.
    '''
    pVertexX = event[1]['centerX'][cID] / event[1]['axisX'][cID] * event[1]['axisZ'][cID] - \
               event[1]['centerZ'][cID] 
    pVertexY = event[1]['centerY'][cID] / event[1]['axisY'][cID] * event[1]['axisZ'][cID] - \
               event[1]['centerZ'][cID]
    pVertexZ = np.average([pVertexX, pVertexY], 
               weights=np.absolute([event[1]['centerX'][cID], event[1]['centerY'][cID]]))
    return pVertexZ, pVertexY, pVertexX


def pVertex(x, y, z, t): 
    '''
    Regresses x, y, z, and t using a nonlinear least squares method to solve for the interaction 
    vertex of a single cluster.
    ''' 
    xfit = lambda p, t: c * np.sin(p[0]) * (t-1.0) * np.sin(p[1])                                   #|Offset of 1 nanosecond due to digitizer reasons I don't fully understand, but whatevs
    yfit = lambda p, t: c * np.sin(p[0]) * (t-1.0) * np.cos(p[1])
    zfit = lambda p, t: c * np.cos(p[0]) * (t-1.0) + p[2]
    def err(p, x, y, z, t):                                                                         #|Error function for least squares
        return np.sum(np.square(xfit(p,t) - x) + np.square(yfit(p,t) - y) + \
               np.square(zfit(p,t) - z))
    def errWrap(p):                                                                                 #|Wrapper function for error, used to keep things consistent
        return err(p,x,y,z,t)
    theta = np.mean(np.arctan2(np.sqrt(np.square(x) + np.square(y)), z))                            #|Initial guess for theta is where theta for the mean of the clusters is
    phi = np.mean(np.arctan2(y, x))                                                                 #|Guess phi from the mean
    z0i = 0.0
    # t0i = 0.0
    params = np.array((theta, phi, z0i))#, t0i)) #| Use (t-1.0-p[3]) if including t0
    result = minimize(errWrap, params.copy(), method='BFGS')
    return result


class fourVector(object):                                                                           #|Energy-momentum four-vector for Higgs filter
    '''Energy-momentum four vector for a physics object.'''
    def __init__(self,E,px,py,pz):
        self.E  = E
        self.px = px
        self.py = py
        self.pz = pz
    def dot(self, other):
        return self.E*other.E - (self.px*other.px + self.py*other.py + self.pz*other.pz)
    def __add__(self, other):
        self.E  += other.E
        self.px += other.px
        self.py += other.py
        self.pz += other.pz
        return self


def invariantMass(event):
    '''Get the invariant mass of a system. Used in a comparison of invariant mass vs error.'''
    # Extract data
    eventpts    = np.sort(event[2], order='pt')[::-1]
    pt1, pt2    = eventpts['pt'][0:2]
    eta1, eta2  = eventpts['eta'][0:2]
    phi1, phi2  = eventpts['phi'][0:2]
    # Calculate needed parameters
    E1,  E2     = pt1*np.cosh(eta1), pt2*np.cosh(eta2)
    px1, px2    = pt1*np.cos(phi1),  pt2*np.cos(phi2)
    py1, py2    = pt1*np.sin(phi1),  pt2*np.sin(phi2)
    pz1, pz2    = pt1*np.sinh(eta1), pt2*np.sinh(eta2)
    # Calculate invariant mass
    eventVector = fourVector(E1,px1,py1,pz1) + fourVector(E2,px2,py2,pz2)
    invmass     = np.sqrt(eventVector.dot(eventVector))
    return invmass


def HiggsInvariantMassFilter(event):
    # Return if within variance of Higgs mass
    variance = 10.0
    invmass = invariantMass(event)
    return invmass > HiggsMass - variance and invmass < HiggsMass + variance                        #|Filters for Higgs mass range


def ROIEnergySpreadFilter(event):
    '''
    Looks at how spread the energy is over clusters. Only keeps events with energies very 
    concentrated in the most energetic clusters.
    '''
    pass


def HggFilter(data, quiet=True):
    '''Filters Hgg data set to clean it up some, removing some of the crap.'''
    if not quiet: print "Filtering"   
    data     = np.compress([len(event[2]) >= 2 for event in data], data, axis=0)                    #|Get at least two ROI's    
    data     = np.compress([len(event[1]) <= 7 for event in data], data, axis=0)                    #|Use at events with at most 7 clusters to cut down on crap
    #pTs      = [np.sort(event[2]['pt'])[::-1] for event in data]                                   #|Get new sorted ROI pT from refined list
    #data     = np.compress([pT[1] > 20 for pT in pTs], data, axis=0)                               #|Second biggest cluster energy must also have 20GeV 
    data     = np.compress([HiggsInvariantMassFilter(event) for event in data], data, axis=0)       #|Pick events with invariant mass near 125GeV
    data     = np.compress([np.max(event[2]['eta'])*np.min(event[2]['eta']) < 0 for event in data], 
                           data, axis=0)                                                            #|Require an eta separation of opposite endcaps to prevent high errors
    return data


def gammaGunFilter(data, quiet=True):
    '''Filters gamma gun data set to clean it up some, removing some of the crap.'''
    if not quiet: print "Filtering"   
    data     = np.compress([len(event[2]) >= 2 for event in data], data, axis=0)                    #|Get at least two ROI's    
    data     = np.compress([np.max(event[2]['eta'])*np.min(event[2]['eta']) < 0 for event in data], 
                           data, axis=0)                                                            #|Require an eta separation of opposite endcaps to prevent high errors
    return data


def vertexData(data, pbar=None, currentFile="", runNumber=None, dataLength=None,
                      errorPlot=False, plotTitle=None, quiet=False, returnDiffs=True, 
                      invMassVsError=False):
    '''Apply the tVertexing/pVertexing algorithm to a large number of events.'''
    if invMassVsError: returnDiffs = False
    i                 = 0
    skippedCount      = 0                                                                           #|Number of skipped events
    incorrectCount    = 0                                                                           #|Number of incorrect vertex predictions
    tVertexZList      = []                                                                          #|List of tVertexed vertices
    pVertexZList      = []                                                                          #|List of pVertexed vertices
    pVertex0List      = []
    pVertex1List      = []
    genVertexZList    = []                                                                          #|List of actual vertices
    correctVertexList = []
    invariantMassList = []
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
            genVertexZ   = event[4]['z'][0] #|Generator level "true" vertices

            # tVertexing
            vertices     = twoVertex(cluster0Hits, cluster1Hits)
            tVertexZ     = vertices[0]                                                              #|Optimal solution 
            soln1, soln2 = vertices[2:4]                                                            #|Pass both solutions to compare
            # error        = vertices[4] 
            if np.absolute(tVertexZ) > 15.0:                                                        #|Automatically discard any fishy looking solutions
                nAttempted   -= 1
                i            += 1
                skippedCount += 1
                if not quiet: print ">> Skipped!"
                continue                                                                            #|Get rid of the ~3sigma solutions that are probably errors.

            # pVertexing
            c0hits       = np.compress(np.logical_and(
                                event[0]['clusterID']==c0ID, event[0]['t']>0), event[0], axis=0)
            c1hits       = np.compress(np.logical_and(
                                event[0]['clusterID']==c1ID, event[0]['t']>0), event[0], axis=0)
            x0,y0,z0,t0  = c0hits['x'], c0hits['y'], c0hits['z'], c0hits['tofT']
            x1,y1,z1,t1  = c1hits['x'], c1hits['y'], c1hits['z'], c1hits['tofT']
            # pVertex0, pVertex0t = pVertex(x0,y0,z0,t0)['x'][2:]
            # pVertex1, pVertex1t = pVertex(x1,y1,z1,t1)['x'][2:]
            pVertex0 = pVertex(x0,y0,z0,t0)['x'][2]
            pVertex1 = pVertex(x1,y1,z1,t1)['x'][2]
            pVertexZ = np.mean((pVertex0, pVertex1))

            # Append to lists
            tVertexError = np.absolute(tVertexZ-genVertexZ)
            pVertexError = np.absolute(pVertexZ - genVertexZ)
            tVertexZList.append(tVertexZ)
            pVertex0List.append(pVertex0)
            pVertex1List.append(pVertex1)
            pVertexZList.append(pVertexZ)
            genVertexZList.append(genVertexZ)
            invariantMassList.append(invariantMass(event))

            # Build output string
            if not quiet:
                string = "Event: %3i  |  NumClusters: %2i  |  NumROIs: %2i  |  " \
                        % (i,len(event[1]), len(event[2])) 
                string += "tVertexZ: {0:9.5f}cm  |  genVertexZ:{1:9.5f}cm  |  Error:{2:>9.5f}cm"\
                        .format(tVertexZ, genVertexZ, tVertexError)
                string += "  |  pV0:{0:>9.5f}cm  |  pV1:{1:>9.5f}cm  |  pVertex Error:{2:>9.5f}cm"\
                        .format(pVertex0, pVertex1, pVertexError)
                string += "  |  tV-pV:{0:>9.5f}cm"\
                        .format(np.absolute(tVertexZ-pVertexZ))
                # string += "  |  pV0t:{0:>9.2f}ns  |  pV1t:{1:>9.2f}ns"\
                #         .format(pVertex0t, pVertex1t)        
                if np.absolute(tVertexZ-genVertexZ) < .001:  string += "<"+"-"*7                     #|Points out very good fits
                elif np.absolute(tVertexZ-genVertexZ) < .01: string += "<"+"-"*3

            # Test for incorrectly identified vertices
            if np.absolute(soln1-genVertexZ) < tVertexError:
                incorrectCount += 1
                if not quiet:
                    string += "<-- Incorrect choice of point. Correct: %.5fcm with error %.5fcm." \
                            % (soln1, np.absolute(soln1-genVertexZ))
            elif np.absolute(soln2-genVertexZ) < tVertexError:
                incorrectCount +=1
                if not quiet:
                    string += "<-- Incorrect choice of point. Correct: %.5fcm with error %.5fcm." \
                            % (soln2, np.absolute(soln2-genVertexZ))
            else:
                correctVertexList.append(tVertexZ-genVertexZ)                                       #|Pick out vertices the algorithm correctly identified

        except RuntimeWarning:            
            string        = ">> Runtime warning, skipping this event..."
            nAttempted   -= 1
            skippedCount += 1
        # except IndexError:
        #     string     = "Event: %3i   "%i + "-"*15 + "Two clusters/ROIs not found" + "-"*15
        #     nAttempted -= 1
        #     skippedCount += 1

        if not quiet: print string
        i += 1

        # Update progressbar
        if runNumber != None and pbar != None:
            pbar.update(ldata*runNumber + i) 
        elif pbar!=None:
            pbar.update(i)

    # Do some manipulation to the lists and possibly clean it up to get a fit
    diffs    = np.subtract(tVertexZList, genVertexZList)
    pDiffs   = np.subtract(pVertexZList, genVertexZList)
    pDiffs0  = np.subtract(pVertex0List, genVertexZList)
    pDiffs1  = np.subtract(pVertex1List, genVertexZList)
    oldlen   = len(diffs)
    if oldlen == 0: return np.array([0])                                                            #|-1 is returned for empty array
    # Extract worst x% of events to improve fitting
    diffs    = np.extract(np.absolute(diffs) < np.percentile(np.absolute(diffs), 99.0), diffs)  #|np.percentile() doesn't exist in lxplus numpy. Remove this if using there.
    pDiffs   = np.extract(np.absolute(pDiffs) < np.percentile(np.absolute(pDiffs), 99.0), pDiffs)
    pDiffs0  = np.extract(np.absolute(pDiffs0) < np.percentile(np.absolute(pDiffs0), 99.0), pDiffs0)
    pDiffs1  = np.extract(np.absolute(pDiffs1) < np.percentile(np.absolute(pDiffs1), 99.0), pDiffs1)
    #invariantMassList = np.extract(np.absolute(diffs) < 1.5, invariantMassList)

    if not quiet: 
        print "Correctly identified %i/%i events in %s." % \
                (nAttempted-incorrectCount, nAttempted, currentFile)
        print "%s events trimmed" % str(oldlen-len(diffs))
        print "%i events skipped." % skippedCount
    if dataLength:
        print "Processed %i/%i events (%.2f%%)." % \
                (nAttempted, dataLength, 100.0*float(nAttempted)/dataLength)

    if errorPlot:
        Plotter.tVertexErrorHist(10*diffs, len(diffs), title=plotTitle, 
                                 ranges=[-7,7], quiet=False)
        Plotter.tVertexErrorHist(10*pDiffs, len(pDiffs), 
                                 title="pVertexed $z$ - genVertex $z$ for 500GeV $\gamma$-gun", 
                                 ranges=[-300,300], quiet=False)
        Plotter.tVertexErrorHist(10*pDiffs0, len(pDiffs0), 
                                 title="pVertexed $z$ - genVertex $z$ for 500GeV $\gamma$-gun (C0)", 
                                 ranges=[-300,300], quiet=False)
        Plotter.tVertexErrorHist(10*pDiffs1, len(pDiffs1), 
                                 title="pVertexed $z$ - genVertex $z$ for 500GeV $\gamma$-gun (C1)", 
                                 ranges=[-300,300], quiet=False)

    if invMassVsError:
        Plotter.invariantMassErrorPlot(np.absolute(diffs), invariantMassList)

    if returnDiffs:
        return diffs 


def timeSmearingTest():
    '''Iterator function to run a resolution analysis on all of the data available.'''
    data           = np.load("Data/500GeVPhoton.npy")
    res            = 50
    diffsList      = np.array([])
    smearingValues = np.array([])
    num            = 0
    pbar           = progressbar("Computing &count&:", res*len(data))
    pbar.start()
    for smearingVal in np.linspace(0,50,res):
        data           = np.load("Data/500GeVPhoton.npy")                                           #|Quick and dirty reload, numpy was doing funky shit
        diffs          = vertexData(timeSmearing(data, smearingVal), runNumber=num, pbar=pbar, quiet=True)
        diffsList      = np.append(diffsList, diffs)
        smearingValues = np.append(smearingValues, np.repeat(smearingVal, len(diffs)))    
        num            += 1
    pbar.finish()
    np.save("diffs.npy", np.multiply(10,diffsList))
    np.save("vals.npy", smearingValues)
    Plotter.tVertexErrorHist2D(np.multiply(10,diffsList), smearingValues)


def energyResolution((directory, f, threadID)):
    '''Multiprocessed function for analysing resolution as a function of energy.'''
    smearingValue = 50
    writer        = Writer((0,2*threadID+2))
    pbar          = progressbar("Loading %s" % f, 1, fd=writer)
    pbar.start()
    pbar.update(0)
    data          = np.load(os.path.expanduser(directory+f))
    data          = timeSmearing(data, smearingValue)
    pbar.finish()
    name          = "{0:22s} &count&: ".format(str(f))
    # en            = f.split("_")[-1].split(".")[0]
    pbar          = progressbar(name, len(data) + 1, fd=writer)
    pbar.start()
    diffs         = vertexData(data, pbar=pbar, quiet=True)
    pbar.finish()
    outfile       = os.path.expanduser("Data/Processed/50ps/50ps"+str(f[:-4]) + "_tVertexed.npy")
    np.save(outfile, np.multiply(10,diffs))
    # print "Done."


def singleThreadEnergyResolution(directory, files):
    '''Singly-processed function for analysing resolution as a function of energy.'''
    dataSets = []
    for f in files:
        print "Loading %s..." % f 
        dataSets.append(np.load(os.path.expanduser(directory+f)))
    totalLength = 0
    for entry in dataSets:
        totalLength += len(entry)

    for data in dataSets: 
        print "Processing %i of %i files..." % (dataSets.index(data), len(dataSets))
        pbar  = progressbar("Computing &count&:", totalLength + 1)
        pbar.start()                         
        diffs = vertexData(data, pbar=pbar, quiet=True)  
        np.save(directory+str(f[:-5]) + "_tVertexed.npy", np.multiply(10,diffs))
        pbar.finish() 
        print "Finished; closing file %s.\n" % f


def roiEnergyAnalysis(data):
    '''Troubleshooting function, compares the observed sum energy in an ROI to the genEnergy''' 
    genEnergies = [] 
    sumEnergies = []
    pbar = progressbar("Processing event &count&:", len(data)+1)
    pbar.start()
    count = 0
    for event in data:
        genEnergy = event[2]['getpt'] * np.cosh(event[2]['geneta'])         
        for i in range(len(genEnergy)):
            clustersIndices = np.compress(event[1]['ROI'] == i, event[1]['clusterID'], axis=0)      #|Only take clusters corresponding to right ROI
            clusterEnergies = []
            for clusterID in clustersIndices:                                                       #|Only take hits corresponding to correct cluster
                hits = np.compress(event[0]['clusterID'] == clusterID, event[0], axis=0) 
                energies = hits['en'] 
                for energy in energies: 
                    clusterEnergies.append(energy)                                                  #|Add the energy to the cluster energies
            ROIEnergy = np.sum(clusterEnergies)
            # Append to original lists
            genEnergies.append(genEnergy[i])
            sumEnergies.append(ROIEnergy)
        pbar.update(count)
        count += 1
    pbar.finish()
    # np.save("sums.npy", sumEnergies)
    # np.save("gens.npy", genEnergies)
    # Plot it
    Plotter.sumEnergyVsGenEnergy(sumEnergies, genEnergies) 


if __name__ == '__main__':  
    print "Loading a lot of shit, may take a minute..."
    #timeSmearingTest()

    #data0 = np.load(os.path.expanduser("Data/Hgg_All.npy"))
    #data = HggFilter(data0)

    #Plotter.invariantMassDistribution(data)


    data = gammaGunFilter(np.load("Data/500GeVPhoton.npy"))

    # data = gammaGunFilter(np.load(os.path.expanduser(
    #            "~/work/public/Events_22_All_Processed/CombinedArray.npy")))
    # roiEnergyAnalysis(data)

    #vertexData(data, plotData=True)
    #data = np.load(os.path.expanduser("Data/Events_211_2.npy"))

    vertexData(data, errorPlot = True, 
        plotTitle="tVertexed $z$ - genVertex $z$ for 500GeV $\gamma$-gun", 
        dataLength=len(data), invMassVsError=False)
    
    #Plotter.showerAnimator(data[6][0], 1, "", frames=5, clusterID=0)

    #for i in range(4):
        #Plotter.showerAnimator(data[72][0], 72+i/10.0, "Hgg Clusters", clusterID=i)
    #    print data[72,1]['en']

    # directory = "~/work/public/Events_22_All_Processed/"
    # #os.chdir(directory)
    # files = [os.path.expanduser(f) for f in os.listdir(os.path.expanduser(directory)) if f.endswith(".npy")]
    # #files = [os.path.absolute(directory+f) for f in files]

    # #singleThreadEnergyResolution(directory, files)
    # with term.fullscreen():
    #     os.system("clear")
    #     pool = Pool(processes = len(files))
    #     print "Processing %i files in %s..." % (len(files), directory)
    #     pool.map(energyResolution, [(directory, files[i], i) for i in range(len(files))])
    #     # for i in range(len(files)):
    #     #     print "using file %s"%files[i]
    #     #     energyResolution((directory, files[i], i))

    #raw_input("Operation completed\nPress enter to return to the terminal.")

    # gens = np.load("Data/gens.npy")
    # sums = np.load("Data/sums.npy")
    # Plotter.sumEnergyVsGenEnergy(sums, gens)


################################################################################################
# Todo:
# Make invariant mass filter that looks over all clusters and picks the one with
#   invmass closest to Higgs mass.










