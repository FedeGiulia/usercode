import FWCore.ParameterSet.Config as cms
import os
import fnmatch


# input dir: full path of the dir in which the files are
dir = 'chic_flatPt_20M/res/'

# regex of files: indicate the regular expression to distinguish the files you want (eg. gen_*.root, *.root...)
regName = 'gen_*.root'

# output name: chose the file name of the output produced
outputname = 'chic_rootuple_MC_10M_gen.root'

# genparticle: does the files contains genparticles from MC? 
genpart = 1
#genpart = 0

process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.destinations = cms.untracked.vstring("debug_rootupleChic.txt")
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("INFO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

inputfiles = []
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, regName):
        filename = 'file:'+dir+str(file)
        inputfiles.append(filename)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( inputfiles )
)
process.TFileService = cms.Service("TFileService", 
	fileName = cms.string(outputname),
	closeFileFast = cms.untracked.bool(True)
)

if(genpart):
	process.load('Ponia.RootupleChib.rootuplechicGen_cfi')
else:
	process.load('Ponia.RootupleChib.rootuplechic_cfi')

process.p = cms.Path(process.rootuple)
