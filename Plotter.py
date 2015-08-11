#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #|Comments Section:
#                                                                                                   #|
#     Ben Bartlett                                     File description:                            #|For the sake of any poor person having to read my code,
#     California Institute of Technology                 Plotter file containing various            #|comments starting with "#|" will be aligned.
#     bartlett@caltech.edu                               functions for plotting results.            #|
#     ben.bartlett@cern.ch                                                                          #|
#     benjamincbartlett@gmail.com                      Notes:                                       #|
#                                                        None.                                      #|
# Created:       3 July 2015                                                                        #|
# Last modified: 3 July 2015                                                                        #|
#####################################################################################################|

# Importations
import os, shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pylab
import warnings
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
from Libraries.FastProgressBar import progressbar
from scipy.stats import norm
from ExternalLibraries.blessings import Terminal
term = Terminal()

# Constants
ECalRadius = 129.0                                                                                  #|Radius of ECal in cm
ECalZ      = 317.0                                                                                  #|Half-length of ECal in cm (from - to + this value)

def vertexPlot(vertex, hits):
    '''Plots a 3D representation of the vertex in comparison to the CMS detector.'''
    vertexZ, vertexT = vertex
    # Set up figure
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    # Plot CMS cylinder
    cylinderRes = 60                                                                                #|Controls how finely the cylinder is drawn. Bigger values increase the number of lines though.
    x           = np.linspace(-ECalRadius, ECalRadius, cylinderRes)
    z           = np.linspace(-ECalZ, ECalZ, cylinderRes)
    Xc, Zc      = np.meshgrid(x, z)
    Yc          = np.sqrt(ECalRadius**2-Xc**2)
    # Draw parameters
    rstride = 20
    cstride = 10
    ax.plot_surface(Zc, Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, edgecolors='#000066')  #|Plot (z,y,x) so that z is facing horizontally along beamline
    ax.plot_surface(Zc, -Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, edgecolors='#000066')
    # Draw hits and lines
    ax.plot([-500,500],[0,0],[0,0],c='m')                                                           #|Simulated beamline
    xh, yh, zh = hits['x'], hits['y'], hits['z']                                                    #|Separating hits
    ax.scatter(zh, yh, xh, c='r', s=120)                                                            #|Plot hits
    ax.scatter(vertexZ, 0, 0, c='r', s=150, marker='*')                                             #|Plot reconstructed vertex
    for i in range(len(xh)):
        ax.plot([vertexZ,zh[i]], [0,yh[i]], [0,xh[i]])                                              #|Plot lines to hits from vertex
    # Set labels and ranges
    ax.set_xlabel("Beamline (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("x (cm)")
    ax.set_xlim(-ECalZ, ECalZ)
    ax.set_ylim(-ECalZ, ECalZ) 
    ax.set_zlim(-ECalZ, ECalZ)

    plt.show()


def showerAnimator(hits, eventNo, title, clusterID=-1, delete=False, frames=30, endLoop = 0,
                                         projections=True, transparency=False):
    '''
    Usage: showerAnimator(recordedHits, eventToAnalyse, plotTitle, clusterToAnalyse=None, 
                          deleteFramesOnFinish=False, enableMatplotlibTransparency=False)

    Plots an animated 5D (x,y,z,t,energy) graph over time and renders it to an animated gif.
    If you want to make plots like this, I've also created a cleaned-up version of this code 
    and put it on GitHub as an easily-callable package. (github.com/bencbartlett/Animator5D)
    '''
    # Set up figure
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    import matplotlib.cm as cm

    # Parse data for proper cluster ID and existence of timing data
    if clusterID >= 0:
        hits = np.extract(np.logical_and(hits['clusterID']==clusterID, hits['t']>0), hits)          #|Select only the rechits corresponding to the given cluster and with time data
    else:
        hits = np.extract(hits['t']>0, hits)                                                        #|Only include timing data points

    # Further data parsing
    hits               = np.sort(hits, order = ['tofT'])
    xh,yh,zh,cth,Eh = hits['x'], hits['y'], hits['z'], hits['tofT'], hits['en']       #|Separating hits
    maxEh              = np.max(Eh)                                                     
    xmin, xmax         = min(xh), max(xh)                                                           #|Get and set limits
    ymin, ymax         = min(yh), max(yh)
    zmin, zmax         = min(zh), max(zh)
    xc                 = np.mean([xmin, xmax])                                                      #|x and y centroids, so we can keep scale proportional on all axes
    yc                 = np.mean([ymin, ymax])
    zd                 = zmax-zmin

    # Set limits and labels
    ax.set_zlim(xc-zd/2, xc+zd/2)                                                                   #|Not a mistake, the reversal of the xyz pairing is so that z appears horizontally on the plot.
    ax.set_ylim(yc-zd/2, yc+zd/2)
    ax.set_xlim(zmin, zmax)
    ax.set_zlabel('x (cm)')                                                                         #|Label stuff
    ax.set_ylabel('y (cm)')
    ax.set_xlabel('Beamline (cm)')

    # Clear out existing crap
    path = "Plots/ShowerAnimations/"+title+"/Event"+str(eventNo)
    if os.path.exists(path):                                                                        #|Remove previous crap
        shutil.rmtree(path)
        os.makedirs(path+"/frames")
    else:
        os.makedirs(path+"/frames")

    # Set up animation
    count    = 1
    t = t0   = np.min(cth)
    maxt     = np.max(cth)
    tstep    = (maxt - t)/float(frames)                                                             #|Time step in ns
    title    = ax.text2D(.3,1.0,'Shower simulation', transform=ax.transAxes, size='large')
    colorcal = ax.scatter([0,0],[0,0],[0,0], c=[0,maxEh + 1], cmap=cm.jet)
    cbar     = fig.colorbar(colorcal, shrink=.7)
    pbar     = progressbar("Rendering &count& frames:", int((maxt-t0)/tstep)+1) 
    cbar.set_label("Energy (MIPs)") 
    pbar.start()
    
    # Render frames
    while t <= np.max(maxt):
        mask = np.logical_and(t<cth, cth<=(t+tstep))                                                #|What to plot in this time step
        xplt = np.extract(mask, xh) 
        yplt = np.extract(mask, yh)
        zplt = np.extract(mask, zh)
        Eplt = np.extract(mask, Eh)
        txt  = ax.text2D(0.1, 0.9,'$t=%.3f$ns'%t, transform=ax.transAxes)
        cx   = np.ones_like(xplt) * ax.get_zlim3d()[0]                                              #|Again, not a typo with mixing x and z
        cy   = np.ones_like(yplt) * ax.get_ylim3d()[1]
        cz   = np.ones_like(zplt) * ax.get_xlim3d()[0]
        mark = ax.scatter(zplt, yplt, xplt, c=Eplt, cmap=cm.jet, vmin=0, vmax=maxEh, \
                                            s=100*Eplt/maxEh, marker=',', lw=1)
        if projections:
            ax.scatter(zplt, yplt, cx, c='#444444', marker=',', lw=0, s=100*Eplt/maxEh, alpha=0.3)  #|Plot the projections
            ax.scatter(zplt, cy, xplt, c='#444444', marker=',', lw=0, s=100*Eplt/maxEh, alpha=0.3)
            ax.scatter(cz, yplt, xplt, c='#444444', marker=',', lw=0, s=100*Eplt/maxEh, alpha=0.3)
        if not transparency: mark.set_edgecolors = mark.set_facecolors = lambda *args:None          #|Super-hacky way to disable transparency in the 3D plot, makes it cleaner to read.

        plt.draw()
        filename = path + "/frames/"+str(count).zfill(3)+".gif"
        plt.savefig(filename)                                                                       #|Save the frame
        txt.remove()
        pbar.update(count)
        count += 1
        t     += tstep

    pbar.finish()
    
    if endLoop: print "Copying tail frames..."
    for i in xrange(endLoop):
        shutil.copyfile(filename, filename[:-7]+str(count).zfill(3)+".gif")
        count += 1

    # Combine frames to animated gif
    print "Combining frames..."
    import subprocess
    args = (['convert', '-delay', '.1', '-loop', '0', 
             path+"/frames/*.gif", path+"/Shower.gif"])                                             #|This part requires ImageMagick to function
    print "Saved to path "+str(os.path.abspath(path))+"/Shower.gif"
    subprocess.check_call(args)

    # Delete frames if told to do so
    if delete:
        shutil.rmtree(path+"/frames")

def XYZtoEtaPhi(xyz):
    '''Convert xyz coordinates to an eta-phi map'''
    x, y, z = xyz 
    phi     = np.arctan2(y, x)
    theta   = np.arctan2(np.sqrt(x**2 + y**2), z)
    eta     = -np.log(np.tan(theta/2)) 
    return eta, phi

def layerVarianceFrame(enWeightedEtaPhi, etaPhi, frameNo, center, rng=None, maxE=1.0,
                       xbins=30, ybins=30, numEvents=1, saveAs=None):
    '''Plots individual frames of eta-phi 2D histogram showing energy variance'''
    from matplotlib.ticker import NullFormatter, MaxNLocator
    from matplotlib.colors import LogNorm
    def ellipse(ra,rb,ang,x0,y0,Nb=100):
        '''Plots 1,2,3 sigma rings'''
        xpos,ypos = x0,y0
        radm,radn = ra,rb
        an        = ang
        co,si     = np.cos(an),np.sin(an)
        the       = np.linspace(0,2*np.pi,Nb)
        X         = radm*np.cos(the)*co - si*radn*np.sin(the) + xpos
        Y         = radm*np.cos(the)*si + co*radn*np.sin(the) + ypos
        return X,Y     
    # Unpack data
    x, y = enWeightedEtaPhi
    eta, phi = etaPhi
    if rng != None:
        xlim, ylim = rng
    else:
        xlim = [np.min(x), np.max(x)]
        ylim = [np.min(y), np.max(y)]
    # Format axes
    nullfmt        = NullFormatter() 
    left, width    = 0.12, 0.54                                                                     #|These determine where the subplots go
    bottom, height = 0.12, 0.54
    bottom_h       = left_h = left+width+0.02
    rectscatter    = [left, bottom, width, height]
    recthistx      = [left, bottom_h, width, 0.2]
    recthisty      = [left_h, bottom, 0.2, height]
    rectcbar       = [left_h+0.22, bottom, 0.05, height]
    # Set up figure and axes
    plt.figure(1, figsize=(8,8))
    axHist2D = plt.axes(rectscatter)
    axHistx  = plt.axes(recthistx)
    axHisty  = plt.axes(recthisty)
    axCbar   = plt.axes(rectcbar)
    # Remove labels from supplementary plots
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    axCbar.xaxis.set_major_formatter(nullfmt)
    axCbar.yaxis.set_major_formatter(nullfmt)
    axCbar.set_visible(False)
    # Configure and plot center and sigma bounds
    contourcolor = 'black'
    xcenter      = np.mean(x)
    ycenter      = np.mean(y)
    ra           = np.std(x)
    rb           = np.std(y)
    ang          = 0
    axHist2D.plot(center[0], center[1], '+w', markersize=50, mew=3.0)
    X,Y = ellipse(ra,rb,ang,xcenter,ycenter)
    axHist2D.plot(X,Y,"k:",ms=1,linewidth=2.0)
    axHist2D.annotate('$1\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
                                    textcoords='offset points', horizontalalignment='right',
                                    verticalalignment='bottom',fontsize=25)    
    X,Y = ellipse(2*ra,2*rb,ang,xcenter,ycenter)
    axHist2D.plot(X,Y,"k:",color = contourcolor,ms=1,linewidth=2.0)
    axHist2D.annotate('$2\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
                                    textcoords='offset points',horizontalalignment='right',
                                    verticalalignment='bottom',fontsize=25, color=contourcolor)     
    X,Y = ellipse(3*ra,3*rb,ang,xcenter,ycenter)
    axHist2D.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
    axHist2D.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
                                    textcoords='offset points',horizontalalignment='right',
                                    verticalalignment='bottom',fontsize=25, color=contourcolor)   
    # Plot the main 2D histogram
    H, xedges, yedges, img = axHist2D.hist2d(x, y, bins=(xbins,ybins), range=rng,
                                                   norm=LogNorm(vmax=numEvents*maxE/50.0))
    cbar = plt.colorbar(img, ax=axCbar, fraction=1.0)                                               #|Fraction makes the colorbar fill the entire axis
    cbar.set_label("$\log$ Energy ($\log$ MIPs) (total for $n$ events)", rotation=90)
    # Set main limits
    axHist2D.set_xlim(xlim)
    axHist2D.set_ylim(ylim)
    axHist2D.set_xlabel("$\eta$",fontsize=25)
    axHist2D.set_ylabel("$\phi$",fontsize=25)
    axHist2D.set_axis_bgcolor(pylab.cm.jet(0.0))
    axHist2D.text(1.05, 1.20, "Layer %s/30"%str(frameNo),transform=axHist2D.transAxes,fontsize=25)
    # Plot two side histograms
    axHistx.hist(eta, range = rng[0], bins=xbins, color='blue')
    axHisty.hist(phi, range = rng[1], bins=ybins, orientation='horizontal', color='red')
    # Set those limits
    axHistx.set_xlim(axHist2D.get_xlim())
    axHistx.set_ylim([0,numEvents*5])
    axHisty.set_ylim(axHist2D.get_ylim())
    axHisty.set_xlim([0,numEvents*5])
    axHistx.set_ylabel("Energy-unweighted hits", rotation=90)
    axHisty.set_xlabel("Energy-unweighted hits")
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    # Finish everything up
    if saveAs:
        plt.savefig(saveAs)
    else:
        plt.show()
    plt.clf()

def layerVarianceAnalysis(data, eventNo, title, mode=None, clusterID=0, rng=None, endLoop=0,
                                                numEvents=1, delete=False, cumulative=False):
    '''
    Plots an evolving 2D histogram by layer of a hit.
    Warning: running this on cumulative mode for very large datasets can make you run out of RAM.
    '''
    # Parse data for proper cluster ID and existence of timing data
    if not cumulative:
        hits        = np.extract(data[eventNo][0]['clusterID']==clusterID, data[eventNo][0])        #|Select only the rechits corresponding to the given cluster
        # Further data parsing
        xh,yh,zh,Eh = hits['x'], hits['y'], hits['z'], hits['en']                                    #|Separating hits
        center      = [data[eventNo][1]['eta'][clusterID], data[eventNo][1]['phi'][clusterID]]
        eta, phi    = XYZtoEtaPhi([xh,yh,zh])
        layers      = hits['layerID']
        if mode == "log":
            repeats = np.ceil(np.log(Eh)).astype(int)
            etaw    = np.repeat(eta, repeats)                                                       #|Energy-weighted eta and phi values
            phiw    = np.repeat(phi, repeats)
            layersw = np.repeat(layers, repeats)
        elif mode == "linear":
            repeats = np.ceil(Eh).astype(int)
            etaw    = np.repeat(eta, repeats)                                                       #|Energy-weighted eta and phi values
            phiw    = np.repeat(phi, repeats)
            layersw = np.repeat(layers, repeats)
        elif mode == "flat":
            etaw, phiw, layersw = eta, phi, layers
    else: 
        center  = [0,0]
        eta     = np.array([])
        phi     = np.array([])
        etaw    = np.array([])
        phiw    = np.array([])
        layers  = np.array([])
        layersw = np.array([])
        Eh      = np.array([])
        k       = 0
        pbar    = progressbar("Processing &count& events:", len(data)+1)
        pbar.start()
        for event in data: 
            hits          = np.extract(event[0]['clusterID']==clusterID, event[0])                  #|Select only the rechits corresponding to the given cluster
            # Further data parsing
            xh,yh,zh,EhE  = hits['x'], hits['y'], hits['z'], hits['en']                             #|Separating hits
            etaE, phiE    = XYZtoEtaPhi([xh,yh,zh])
            layersE       = hits['layerID']
            for i in range(len(hits)):
                clusterID = hits[i]['clusterID']
                etaE[i]  -= event[1]['eta'][clusterID]
                phiE[i]  -= event[1]['phi'][clusterID]
                if phiE[i] > 3.0:
                    phiE[i] -= 2*np.pi                                                              #|This fixes some modular issues we were having, with -epsilon = 2pi-epsilon
                elif phiE[i] < -3.0:
                    phiE[i] += 2*np.pi
            if mode == "log":
                repeats  = np.ceil(np.log(EhE)).astype(int)
                etawE    = np.repeat(etaE, repeats)                                                 #|Energy-weighted eta and phi values
                phiwE    = np.repeat(phiE, repeats)
                layerswE = np.repeat(layersE, repeats)
            elif mode == "linear":
                repeats  = np.ceil(EhE).astype(int)
                etawE    = np.repeat(etaE, repeats)                                                 #|Energy-weighted eta and phi values
                phiwE    = np.repeat(phiE, repeats)
                layerswE = np.repeat(layersE, repeats)
            elif mode == "flat":
                etawE, phiwE, layerswE = etaE, phiE, layersE
            eta      = np.append(eta, etaE)
            phi      = np.append(phi, phiE)
            etaw     = np.append(etaw, etawE)
            phiw     = np.append(phiw, phiwE)
            layers   = np.append(layers, layersE)
            layersw  = np.append(layersw, layerswE)
            Eh       = np.append(Eh, EhE)
            k       += 1
            pbar.update(k)
        pbar.finish()
    #     print "Saving array..."
    #     np.save("Data/etaPhiProcessed.npy", [eta, phi, etaw, phiw, layers, layersw, Eh])
    # eta, phi, etaw, phiw, layers, layersw, Eh = np.load("Data/etaPhiProcessed.npy")
    # center = [0,0]
    # Set plot ranges
    if rng == None:
        pltrange = np.array([(np.min(eta), np.max(eta)), (np.min(phi), np.max(phi))])
    else:
        pltrange = rng
    # Clear out existing crap
    path = "Plots/LayerHitAnimations/"+title+"/Event"+str(eventNo)
    if os.path.exists(path):                                                                        #|Remove previous crap
        shutil.rmtree(path)
        os.makedirs(path+"/frames")
    else:
        os.makedirs(path+"/frames")
    # Set up animation
    minlayer = 1                                                                                    #|Minimum layer
    maxlayer = 30                                                                                   #|Maximum layer
    count    = minlayer
    pbar     = progressbar("Rendering &count& frames:", maxlayer-minlayer+1)
    pbar.start()
    # Render frames
    while count <= maxlayer:
        etapltw = np.extract(layersw == count, etaw) 
        phipltw = np.extract(layersw == count, phiw)
        etaplt = np.extract(layers == count, eta)
        phiplt = np.extract(layers == count, phi)
        # Eplt = np.extract(layers == count, Eh)
        filename = path + "/frames/"+str(count).zfill(3)+".gif"
        if len(etaplt) > 0:
            layerVarianceFrame([etapltw, phipltw], [etaplt, phiplt], count, center, rng=pltrange,
                                numEvents=numEvents, maxE=np.max(Eh), xbins=50, ybins=50, 
                                saveAs=filename)
        pbar.update(count)
        count += 1
    pbar.finish()
    # Render tail if needed:
    if endLoop: print "Copying tail frames..."
    for i in xrange(endLoop):
        shutil.copyfile(filename, filename[:-7]+str(count).zfill(3)+".gif")
        count += 1
    # Combine frames to animated gif
    print "Combining frames..."
    import subprocess
    args = (['convert', '-delay', '25', '-loop', '0', path+"/frames/*.gif", path+"/Shower.gif"])    #|This part requires ImageMagick to function
    print "Saved to path "+str(os.path.abspath(path))+"/Shower.gif"
    subprocess.check_call(args)
    # Delete frames if told to do so
    if delete:
        shutil.rmtree(path+"/frames")



def tVertexErrorHist(diffs, nEvents, title=None, ranges=None, quiet=True):
    '''
    Plots an error histogram for tVertex-genVertex z values.
    Usage: tVertexErrorHist(differences, number of counts, quiet=True)
    '''
    absDiffs         = np.absolute(diffs)
    fig, ax          = plt.subplots()                                                               #|Set up graph
    n, bins, patches = ax.hist(diffs, normed = False, range=ranges, bins = 100)
    (mu, sigma)      = norm.fit(diffs)                                                              #|Fit curve
    muerr            = sigma/np.sqrt(nEvents)
    dx               = bins[1]-bins[0]                                                              #|Get bin width
    scale            = dx * len(absDiffs)                                                           #|Scale histogram
    fitline          = mlab.normpdf(bins, mu, sigma) * scale
    ax.plot(bins, fitline, "r--", linewidth = 2) #|Plot fit line
    ax.set_xlabel("Error (mm)") 
    ax.set_ylabel("Counts ($\Sigma=%i$)" % nEvents)
    if title == None:
        ax.set_title("tVertexed $z$ - genVertex $z$ for 500GeV $\gamma$-gun")
    else: 
        ax.set_title(title)
    # Build output and figure text string
    string =  "$\mu$ = %.3f$\pm$%.3fmm, $\sigma_{fitted}$ = %.3fmm \n" % (mu, muerr, sigma)           #|Print out info box 
    string += "68% error magnitude:      {0:>6.3f}mm \n".format(np.percentile(absDiffs, 68.27)) 
    string += "Median error magnitude:   {0:>6.3f}mm \n".format(np.median(absDiffs))
    string += "Mean error magnitude:     {0:>6.3f}mm \n".format(np.mean(absDiffs)) 
    # Changing font to look nicer and line up colons
    font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 9}
    matplotlib.rc('font', **font)                                                                   #|Changes fond to monospace 
    # Annotate plot
    ax.text(.02, .99, string, transform=ax.transAxes, verticalalignment='top')
    if not quiet: print string 
    # Show it
    plt.show()


def tVertexErrorHist2D(diffsList, values):
    '''
    Plots a 2D error histogram for the tVertex-genVertex z values along some other axis,
    such as time smearing.
    Usage: tVertexErrorHist2D(differencesList)
    Differences are a 2D list of differences for each time smearing value.
    '''
    pylab.hist2d(diffsList, values, bins=(100, len(np.unique(values))))
    pylab.xlabel("tVertexed $z$ - genVertex $z$ (mm)")
    pylab.ylabel("Smearing $\sigma$ (ps)")
    pylab.title("tVertexing Precision Over Time Smearing")
    cbar = pylab.colorbar()
    cbar.set_label("Counts (normalized per $\sigma$ value)")
    plt.show()


def smearingBarPlot(diffsList, values):
    '''Plots a bar chart for median vertexing errors over time smearing values'''
    medians = [np.mean(np.extract(values == i, np.absolute(diffsList))) for i in np.unique(values)]
    occurrences = [len(np.extract(values == i, values)) for i in np.unique(values)]
    pylab.bar(np.unique(values), medians, color='#000088', yerr=medians/np.sqrt(occurrences))
    pylab.xlim(np.min(values), np.max(values))
    pylab.xlabel("Smearing $\sigma$ (ps)")
    pylab.ylabel("Median error (mm)")
    pylab.title("tVertexing Precision Over Time Smearing")
    plt.show()


def energyResolutionBarPlot(en, medians, occurrences):
    pylab.bar(en, medians, width=np.multiply(0.2,en), color='#000088', 
                           yerr=medians/np.sqrt(occurrences))
    pylab.xlim(9.0, np.max(en)+300)
    pylab.xscale("log")
    pylab.xlabel("Photon energy (GeV - log scale)")
    pylab.ylabel("Median error (mm)")
    pylab.title("tVertexing Precision Over $\gamma$ Energy - 50ps Time Smearing")
    comment = "Clustering algorithm\nfails at very low\n" + \
              "energies, likely\ndue to energy\nthresholding."
    pylab.text(2.2,0.25,comment)
    tickLabels = [str(eng)+"GeV" for eng in en]
    pylab.xticks(en, tickLabels, rotation="vertical")
    plt.show()


def energyResolutionLinePlot(en, medians, occurrences):
    pylab.errorbar(en, medians, color='#000088', lw=2.0, 
                                yerr=medians/np.sqrt(occurrences), ecolor="r")
    pylab.fill_between(en, medians + medians/np.sqrt(occurrences), 
                           medians - medians/np.sqrt(occurrences), color='#BBBBBB')
    pylab.xlim(9.0, np.max(en)+100)
    pylab.ylim(0.0,4.0)
    pylab.xscale("log")
    pylab.xlabel("Photon energy (GeV - log scale)")
    pylab.ylabel("Median error (mm)")
    pylab.title("tVertexing Precision Over $\gamma$ Energy - No Time Smearing")
    #pylab.text(2.2,0.25,"Clustering algorithm\nfails at very low\nenergies, likely\ndue to energy\nthresholding.")
    tickLabels = [str(eng)+"GeV" for eng in en]
    pylab.xticks(en, tickLabels, rotation="vertical")
    plt.savefig("LinePlot0ps.png")
    plt.show()


def energySpectrum(data):
    '''Quick, sloppy funciton. Plots pT spectrum for a distribution.'''
    print "Getting energies, etas, pTs..."  
    eventpTs = [np.sort(event[2]['pt'])[::-1] for event in data]                                    #|Get ROI pT in descending order    
    pTs      = np.extract([len(event) >= 2 for event in eventpTs], eventpTs)                        #|Get at least two clusters    
    pTs      = np.extract([pT[1] > 20 for pT in pTs], pTs)                                          #|Second biggest cluster energy must also have 20GeV 
    flatpTs  = np.hstack(pTs)
    fig, ax  = plt.subplots()
    ax.hist(flatpTs, bins = 100, range = [0,300])
    ax.set_xlabel("ROI $p_T \cdot c$ (GeV)")
    ax.set_ylabel("Count")
    ax.set_title("$H\\rightarrow\gamma\gamma$ cluster $p_T$ spectrum")
    plt.show()

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

def invariantMassDistribution(data):
    invMasses = []
    for event in data:
        invMasses.append(invariantMass(event))
    fig, ax  = plt.subplots()
    ax.hist(invMasses, bins = 70, range = [100,150])
    ax.set_xlabel("Invariant mass (GeV/$c^2$)")
    ax.set_ylabel("Count")
    ax.set_title("Invariant mass of filtered $H\\rightarrow\gamma\gamma$ events")
    plt.show()

def invariantMassErrorPlot(absDiffs, masses):
    '''Plots error magnitude against invariant mass'''
    plt.scatter(masses, absDiffs)
    plt.xlabel("Invariant mass (GeV/$c^2$)")
    plt.ylabel("Error magnitude")
    plt.title("Error vs Invariant Mass for $H\\rightarrow\gamma\gamma$")
    # plt.xlim([0,1000])
    # plt.ylim([0,1.6])
    plt.show()

def sumEnergyVsGenEnergy(sums, gens):
    '''Plots the sum of energies in a cluster vs the generator level energies of the cluster.'''
    plt.scatter(gens, sums)
    plt.xlabel("Generator level energy")
    plt.ylabel("Total ROI energy")
    plt.title("Observed vs Generated Energy")
    #plt.xlim([499,501])
    # plt.ylim([0,1.6])
    plt.show()

def sumEnergyHist(sums):
    fig, ax  = plt.subplots()
    ax.hist(sums, bins = 70)
    ax.set_xlabel("Observed energy (MIPs)")
    ax.set_ylabel("Count")
    ax.set_title("Observed energy for given emitted energy")
    plt.show()


if __name__ == '__main__':                                                                       
    # Test stuff
    # smearingBarPlot(np.load("diffs.npy"), np.load("vals.npy"))


    # directory    = "Data/Processed/"
    # files        = [f for f in os.listdir(directory) if f.endswith("_tVertexed.npy")]
    # en           = [] 
    # medians      = []
    # occurrences  = []
    # for f in files:        
    #     absdiffs = np.absolute(np.load(directory+f))
    #     if len(absdiffs) <= 1:
    #         continue
    #     en.append(int(f.split("_")[-2]))
    #     medians.append(np.median(absdiffs))
    #     occurrences.append(len(absdiffs))
    # # Sort by energy
    # order       = np.argsort(en)
    # en          = np.array(en)[order]
    # medians     = np.array(medians)[order]
    # occurrences = np.array(occurrences)[order]
    # energyResolutionLinePlot(np.array(en), medians, occurrences)

    # data = [np.random.randn(1000), np.random.randn(1000)]
    # layerVarianceFrame(data,1)

    data = np.load("Data/500GeVPhoton.npy")
    #data = data[:50]
    #hit = data
    layerVarianceAnalysis(data, 0, "500GeVPhoton", mode="linear",numEvents=len(data),
                          rng=[[-.2,.2],[-.2,.2]], endLoop=10, cumulative=True)


