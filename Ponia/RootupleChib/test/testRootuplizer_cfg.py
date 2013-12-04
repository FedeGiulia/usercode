import FWCore.ParameterSet.Config as cms
import os
import fnmatch

process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.destinations = cms.untracked.vstring("debug_rootupleChib.txt")
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("INFO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# input dir
#dir = '/data1/gdujany/chiB/CMSSW_5_3_3/src/SimChiB/crab/Multicrab_ChiB_1P_1_2_UpsilonPt_10M_each/ChiB_1P_1_UpsilonPt_10M/res/'

inputfiles = ['file:sel_chib_10_1_gQw.root']


#inputfiles = []
#for file in os.listdir(dir):
#    if fnmatch.fnmatch(file, 'testdimuon*.root'):
#        filename = 'file:'+dir+str(file)
#        inputfiles.append(filename)

process.source = cms.Source("PoolSource",
    	fileNames = cms.untracked.vstring( inputfiles )
)
process.TFileService = cms.Service("TFileService", 
	fileName = cms.string("rootuplaTest.root"),
	closeFileFast = cms.untracked.bool(True)
)


process.load('RootupleChib.RootupleChib.rootuplechib_cfi')

process.p = cms.Path(process.rootuple)
