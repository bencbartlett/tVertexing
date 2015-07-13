import ROOT
import sys
ROOT.gSystem.Load("libFWCoreFWLite.dylib")
ROOT.AutoLibraryLoader.enable() 
sys.path.append("/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/")
import numpy as np
fIn = ROOT.TFile.Open("ProcessedEDM/GammaGunProcessed22_100_11.root")
tree = fIn.Get('analysis/HGC')
tree.GetEntry(0)