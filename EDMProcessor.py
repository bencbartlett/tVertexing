#####################################################################################################|
# Time-Based Vertex Reconstruction Algorithm for the Compact Muon Solenoid                          #
#                                                                                                   #
#     Ben Bartlett                                     File description:                            #
#     California Institute of Technology                 EDM processor: takes full EDM              #
#     bartlett@caltech.edu                               files and converts them to a more          #
#     ben.bartlett@cern.ch                               usable file format.                        #
#     benjamincbartlett@gmail.com                      Notes:                                       #
#                                                        Modified from a previously existing file:  #
# Created:       4 June 2015                             runHGCROIAnalyzer_cfg.py. This does not    #
# Last modified: 4 June 2015                             run well locally; only use lxplus.         #
#####################################################################################################

'''
Usage in cmsenv: cmsRun EDMProcessor.py <directory or file> [<output file name>] [<files to ignore>]
If output file name not provided, uses "ProcessedEDM/ProcessedEDM.root"
For EOS inputs, use /store/...., without the /eos/cms/
Filenames to ignore can be partial - the full file path is not needed.
'''

import FWCore.ParameterSet.Config as cms

# Process Loading
process = cms.Process("ROIAnalysis")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')    
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")                                               
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False),
                                       SkipEvent = cms.untracked.vstring('ProductNotFound')) 

# Getting arguments from terminal
import os,sys                                                                                       #|Configure from command line
from fnmatch import fnmatch                                                                         #|Used for wildcard matching
if(len(sys.argv)<2):                                                                                #|cmsRun test/runHGCSimHitsAnalyzer_cfg.py tag
    print '\n Usage: cmsRun EDMProcessor.py input [outputfile]\n'                                   #|where tag can be any sub-directory under /store/cmst3/group/hgcal/CMSSW
    sys.exit()                                                                                      #|or any upgrade relval sample (may need tweaking for new releases...)
input2process=sys.argv[2]
outputName='../ProcessedEDM/ProcessedEDM.root'
if(len(sys.argv)==4):
    outputName=sys.argv[3]
elif(len(sys.argv)>4):
    outputName=sys.argv[3]
    ignore=sys.argv[4:]                                                                             #|Files to ignore

if "*" in input2process:
    origInput = "*"+input2process                                                                   #|Adding the initial * allows for relative file locations.
    useWildCards = True
    input2process = input2process.rsplit("/",1)[0]                                                  #|Gets folder above directory
else: 
    useWildCards = False

print '[EDMProcessor] processing from %s and output name is %s'%(input2process,outputName)

from UserCode.HGCanalysis.storeTools_cff import fillFromStore                                       #|Configure the source (list all files in directory within range [ffile,ffile+step]
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
if input2process.find('file') >= 0:
    process.source.fileNames = cms.untracked.vstring(input2process)
else :                       
    process.source.fileNames = fillFromStore(input2process)

#print process.source.fileNames

# Postprocessing - use wildcards
fileList = []
if useWildCards:
    for filename in process.source.fileNames:
        #print "Filename: %s \norigInput: %s" % (filename, origInput)
        if fnmatch(filename, origInput):                                                            #|Only use files matching wildcard
            fileList.append(filename)
    process.source.fileNames = cms.untracked.vstring(fileList)


# Postprocessing - remove ignore files
if(len(sys.argv)>4):
    fileList = [filename for filename in process.source.fileNames]
    ignoreList = []
    for filename in fileList:
        for ignorefile in ignore:
            if ignorefile in filename:
                fileList.remove(filename)
                ignoreList.append(filename)
    process.source.fileNames = cms.untracked.vstring(fileList)


print "Processing files:"
print process.source.fileNames
print "Ignoring files: "
if(len(sys.argv)>4):
    for i in ignoreList:
        print i 

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#from RecoJets.Configuration.RecoPFJets_cff import *                                                #|Prepare alternative jet collections
#process.ak3PFJetsPandora = ak4PFJets.clone( src=cms.InputTag('pandorapfanew'), rParam = 0.3 )
#process.load('RecoJets.Configuration.GenJetParticles_cff')
#from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
#process.ak3GenJets = ak5GenJets.clone(rParam = 0.3)

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputName))               #|Load the analyzer
process.load('UserCode.HGCanalysis.hgcROIAnalyzer_cfi')

process.p = cms.Path(process.analysis)                                                              #|Run it




