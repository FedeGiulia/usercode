debug_filename = 'debug.txt'

import FWCore.ParameterSet.Config as cms
process = cms.Process('test')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_E_V31::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.destinations = cms.untracked.vstring(debug_filename)
process.MessageLogger.cout = process.MessageLogger.cerr
process.MessageLogger.cout.threshold = cms.untracked.string("DEBUG")

input_filename = "file:9A42402F-DF68-E211-9E73-E41F131812CC.root"

input_filename = 'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/9A42402F-DF68-E211-9E73-E41F131812CC.root'

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(input_filename)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_Dimuon8_Jpsi_v*'),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "gtDigis" ),
                                        l1tIgnoreMask = cms.bool( False ),
                                        l1techIgnorePrescales = cms.bool( False ),
                                        daqPartitions = cms.uint32( 1 ),
                                        throw = cms.bool( True )
                                        )

process.load('Ponia.Modules.CHARM_chiCandProducer_cff')

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('sel_chic.root'),
    outputCommands =  cms.untracked.vstring('drop *',
					    'keep *_offlinePrimaryVertices_*_*',
#					    Temporary, testing if it solves excessive use of ram in chic analysis
#                                            'keep *_dimuonProducer_UpsilonCandLorentzVector_*',
                                            'keep *_chiCandProducer_chicand_*',
					    'keep *_chiCandProducer_piZeroRejectCand_*',
					    'keep *_refit_*_*'
                                            ),
    )

process.end = cms.EndPath(process.out)

#process.p = cms.Path(process.triggerSelection*process.chiSequence)
process.p = cms.Path(process.triggerSelection*process.chiSequence)
process.schedule = cms.Schedule(process.p,process.end)


# HACK ! Need to understand !
process.dimuonProducer.addMuonlessPrimaryVertex = False
