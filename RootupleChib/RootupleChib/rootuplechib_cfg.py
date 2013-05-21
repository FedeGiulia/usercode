import FWCore.ParameterSet.Config as cms
import os
import fnmatch

process = cms.Process("Rootuple")

#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("debug_rootupleChib.txt")
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("INFO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# input dir
#dir = '/data1/gdujany/chiB/CMSSW_5_3_3/src/SimChiB/crab/ChiB_1P_12_UpsilonPt_10M/res/'
dir = '/data1/gdujany/chiB/CMSSW_5_3_3/src/2012TotalData/chib-Run2012A1/res/'

inputfiles = []
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, 'testdimuon*.root'):
        filename = 'file:'+dir+str(file)
        inputfiles.append(filename)
inputfiles = inputfiles[:1]
        
process.source = cms.Source("PoolSource",
    	fileNames = cms.untracked.vstring( inputfiles )
)
process.TFileService = cms.Service("TFileService", 
	fileName = cms.string("ChiB_1P_12_UpsilonPt.root"),
	closeFileFast = cms.untracked.bool(True)
)

process.load('RootupleChib.RootupleChib.rootuplechib_cfi')
#rootuple.isMC = True

process.p = cms.Path(process.rootuple)
