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
import sys
import time

from ExternalLibraries.progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
                                          FileTransferSpeed, FormatLabel, Percentage, \
                                          ProgressBar, ReverseBar, RotatingMarker, \
                                          SimpleProgress, Timer

def progressbar(title, max):
    '''Usage: pbar = progressbar("firstPart &count& secondPart", count)
                or
              pbar = progressbar("Title", count)
              pbar.update(count)
       Looks like:
       Rendering 33 of 33 frames: 99% [#################################] ETA: 0:00:03
       '''
    if "&count&" in title: #|&count& can be used to keep an updating queue list
        title = title.split("&count&")
        widgets = [title[0], SimpleProgress(), title[1], Percentage(), ' ',
                   Bar(marker='#',left='[',right=']'), ' ', ETA()]
    else:
        widgets = [title, ' ', Percentage(), ' ',
                   Bar(marker='#',left='[',right=']'), ' ', ETA()]
    pbar = ProgressBar(widgets=widgets, maxval= max)
    return pbar

# Doesn't work yet
def timer(title):
    widgets = [title, Timer()]
    pbar = ProgressBar(widgets=widgets)

if __name__ == '__main__':
    pbar = progressbar("test")
    pbar.start()
    for i in range(100,500+1,50):
        time.sleep(0.2)
        pbar.update(i)
    pbar.finish()





