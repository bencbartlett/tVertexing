#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #|Comments Section:
#                                                                                                   #|
#     Ben Bartlett                                     File description:                            #|For the sake of any poor person having to read my code,
#     California Institute of Technology                 Plotter file containing various            #|comments starting with "#|" will be aligned.
#     bartlett@caltech.edu                               functions for plotting results.            #|
#     ben.bartlett@cern.ch                                                                          #|
#     benjamincbartlett@gmail.com                      Notes:                                       #|
#                                                        None.                                      #|
# Created:       3 June 2015                                                                        #|
# Last modified: 3 June 2015                                                                        #|
#####################################################################################################|

# Importations
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d
from Misc import progressbar
import os

# Constants
ECalRadius = 129.0                                                                                 #|Radius of ECal in cm
ECalZ      = 317.0                                                                                 #|Half-length of ECal in cm (from - to + this value)

def vertexPlot(vertex, hits):
    vertexZ, vertexT = vertex

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')

    # Plot CMS cylinder
    cylinderRes = 100                                                                               #|Controls how finely the cylinder is drawn. Bigger values increase the number of lines though.
    x           = np.linspace(-ECalRadius, ECalRadius, cylinderRes)
    z           = np.linspace(-ECalZ, ECalZ, cylinderRes)
    Xc, Zc      = np.meshgrid(x, z)
    Yc          = np.sqrt(ECalRadius**2-Xc**2)

    # Draw parameters
    rstride = 20
    cstride = 10
    ax.plot_surface(Zc, Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, edgecolors='#0000AA')  #|Plot (z,y,x) so that z is facing horizontally along beamline
    ax.plot_surface(Zc, -Yc, Xc, alpha=0.1, rstride=rstride, cstride=cstride, edgecolors='#0000AA')

    # Draw hits and lines
    ax.plot([-500,500],[0,0],[0,0],c='m')                                                         #|Simulated beamline
    xh, yh, zh = hits['x'], hits['y'], hits['z']                                                    #|Separating hits
    ax.scatter(zh, yh, xh, c='r', s=100)                                                            #|Plot hits
    ax.scatter(vertexZ, 0, 0, c='g', s=120, marker='*')                                             #|Plot reconstructed vertex
    for i in range(len(xh)):
        ax.plot([vertexZ,zh[i]], [0,yh[i]], [0,xh[i]])                                              #|Plot lines to hits from vertex

    # Set labels and ranges
    ax.set_xlabel("Beamline (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("x (mm)")
    ax.set_xlim(-ECalZ, ECalZ)
    ax.set_ylim(-ECalZ, ECalZ) 
    ax.set_zlim(-ECalZ, ECalZ)

    plt.show()


def showerAnimator(hits, title, eventNo, clusterID=0):
    fig = plt.figure()
    #ax = Axes3D(fig)
    ax  = fig.add_subplot(111, projection='3d')
    #cm  = plt.get_cmap("gist_rainbow")                                                             #|Colormap for 3D plot
    import matplotlib.cm as cm

    if clusterID > 0:
        hits = np.extract(np.logical_and(hits['clusterID']==clusterID, hits['t']>0), hits)          #|Select only the rechits corresponding to the given cluster and with time data
    else:
        hits = np.extract(hits['t']>0, hits)                                                        #|Only include timing data points

    hits               = np.sort(hits, order = ['correctT'])
    xh,yh,zh,th,cth,Eh = hits['x'], hits['y'], hits['z'], hits['t'], hits['correctT'], hits['en']   #|Separating hits
    maxEh              = np.max(Eh)
    normEh             = Eh / maxEh                                                        
    xmin, xmax         = min(xh), max(xh)                                                           #|Get and set limits
    ymin, ymax         = min(yh), max(yh)
    zmin, zmax         = min(zh), max(zh)
    xc                 = np.mean([xmin, xmax])  #|x and y centroids, so we can keep scale proportional on all axes
    yc                 = np.mean([ymin, ymax])
    zd                 = zmax-zmin
    ax.set_zlim(xc-zd/2, xc+zd/2)                                                                   #|Not a mistake, the reversal of the xyz pairing is so that z appears horizontally on the plot.
    ax.set_ylim(yc-zd/2, yc+zd/2)
    ax.set_xlim(zmin, zmax)
    ax.set_zlabel('x (cm)')                                                                         #|Label stuff
    ax.set_ylabel('y (cm)')
    ax.set_xlabel('Beamline (cm)')

    path = "Plots/ShowerAnimations/"+title+"/Event"+str(eventNo)
    if os.path.exists(path):                                                                        #|Remove previous crap
        import shutil
        shutil.rmtree(path)
        os.makedirs(path)
    else:
        os.makedirs(path)


    # Render frames
    count    = 1
    tstep    = .005                                                                                 #|Time step in ns
    t = t0   = np.min(cth)
    maxt     = np.max(cth)
    title    = ax.text2D(.3,1.0,'Shower simulation', transform=ax.transAxes, size='large')
    colorcal = ax.scatter([0,0],[0,0],[0,0], c=[0,maxEh + 1], cmap=cm.rainbow)
    cbar     = fig.colorbar(colorcal)
    pbar     = progressbar("Rendering &count& frames:", int((maxt-t0)/tstep)+1) 
    cbar.set_label("Energy (keV)") 
    pbar.start()
    while t <= np.max(maxt):
        mask = np.logical_and(t<cth, cth<=(t+tstep))                                                #|What to plot in this time step
        xplt = np.extract(mask, xh) 
        yplt = np.extract(mask, yh)
        zplt = np.extract(mask, zh)
        Eplt = np.extract(mask, Eh)
        txt  = ax.text2D(0.1, 0.9,'$t=%.3f$ns'%t, transform=ax.transAxes)
        mark = ax.scatter(zplt, yplt, xplt, c=Eplt, cmap=cm.rainbow, vmin=0, vmax=maxEh, \
                                            s=100*Eplt/maxEh, marker='o', lw=1)
        mark.set_edgecolors = mark.set_facecolors = lambda *args:None                               #|Super-hacky way to disable transparency in the 3D plot, makes it cleaner to read.
        
        # print "Render frame "+str(count)+" with t="+str(t)+"/"+str(maxt)+"..."
        plt.draw()
        plt.savefig(path + "/"+str(count).zfill(3)+".gif")                                          #|Save the frame
        txt.remove()
        pbar.update(count)
        count += 1
        t     += tstep
    pbar.finish()
    # Combine frames to animated gif
    print "Combining frames..."
    import subprocess
    args = (['convert', '-delay', '.1', '-loop', '0', path+"/*.gif", path+"/Shower.gif"])           #|This part requires ImageMagick to function
    subprocess.check_call(args)


if __name__ == '__main__':                                                                          #| Test stuff
    vertexZ = 500
    hits    = [[ECalRadius,0.0,2000.0,0],[0.0,ECalRadius,2000.0,0]]
    vertexPlot(vertexZ, hits)






