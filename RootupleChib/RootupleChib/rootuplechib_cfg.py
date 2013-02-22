import FWCore.ParameterSet.Config as cms

process = cms.Process("Rootuple")

#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.destinations = cms.untracked.vstring("debug_rootupleChib.txt")
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("INFO")

from filelist import inputfiles

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

process.source = cms.Source("PoolSource",
    	fileNames = cms.untracked.vstring( inputfiles )
)
process.TFileService = cms.Service("TFileService", 
#	fileName = cms.string("/data1/degano/storage/ChiB/RooTuples/2012_ALL/2012_ALL_pi0Fixed.root"),
	fileName = cms.string("test_Y_sigmas.root"),
	closeFileFast = cms.untracked.bool(True)
)

process.rootuple = cms.EDAnalyzer('RootupleChib',
	chi_cand = cms.InputTag("chiCandProducer","chicand"),
	pi0_comb = cms.InputTag("chiCandProducer", "piZeroRejectCand"),
	ups_cand = cms.InputTag("dimuonProducer", "UpsilonCandLorentzVector")
)

process.p = cms.Path(process.rootuple)
