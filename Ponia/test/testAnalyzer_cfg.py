input_filename = 'file:/data1/vcogoni/CMSSW_5_3_6/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FE7159F0-33CF-E111-A303-002618FDA279.root'


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
#process.MessageLogger.destinations = cms.untracked.vstring(debug_filename)
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("DEBUG")

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(input_filename)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )

process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
                                         triggerConditions = cms.vstring(
        'HLT_Dimuon5_Upsilon_v*',
        'HLT_Dimuon7_Upsilon_v*',
        'HLT_Dimuon8_Upsilon_v*',
        'HLT_Dimuon11_Upsilon_v*',
        ),
                                         hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                         l1tResults = cms.InputTag( "gtDigis" ),
                                         l1tIgnoreMask = cms.bool( False ),
                                         l1techIgnorePrescales = cms.bool( False ),
                                         daqPartitions = cms.uint32( 1 ),
                                         throw = cms.bool( True )
                                         )

process.load('Ponia.Modules.chiCandProducer_cff')

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('sel_chib.root'),
    outputCommands =  cms.untracked.vstring('drop *',
					    'keep *_offlinePrimaryVertices_*_*',
                                            'keep *_dimuonProducer_UpsilonCandLorentzVector_*',
                                            'keep *_chiCandProducer_chicand_*',
					    'keep *_chiCandProducer_piZeroRejectCand_*',
					    'keep *_refit1S_*_*',
					    'keep *_refit2S_*_*',
					    'keep *_refit3S_*_*',
                                            ),
    )

process.end = cms.EndPath(process.out)

process.p = cms.Path(process.triggerSelection*process.chibSequence)
process.schedule = cms.Schedule(process.p,process.end)


# HACK ! Need to understand !
process.dimuonProducer.addMuonlessPrimaryVertex = False
