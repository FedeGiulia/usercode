import FWCore.ParameterSet.Config as cms


process = cms.Process('test')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Ponia.Configuration.MuonSelection')
process.load('Ponia.Modules.dimuonCandProducer_cfi')
process.load('MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_E_V31::All')

process.MessageLogger.destinations = cms.untracked.vstring("debug_ChiAnalyzer.txt")
process.MessageLogger.cout = process.MessageLogger.cerr
process.MessageLogger.cout.threshold = cms.untracked.string("DEBUG")

process.source = cms.Source(
    
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/data1/vcogoni/CMSSW_5_3_6/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FE7159F0-33CF-E111-A303-002618FDA279.root')
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.selMuons = process.muonSelection2012.clone()

process.diMuonCount = cms.EDFilter(
    'CandViewCountFilter',
    src = cms.InputTag('dimuon'),
    minNumber = cms.uint32(1),
    )

process.diMuonSelection = cms.EDFilter(
    'CandViewSelector',
    src = cms.InputTag('dimuon'),
    cut = cms.string('p4.M > 8.5 && p4.M < 11.0 && p4.Pt > 6.0')
)

process.PreDump = cms.EDAnalyzer(
    'DumpConversion',
    conversions = cms.InputTag('allConversions')
    )

process.conversionProducer = cms.EDProducer(
    'PhotonConversionCandProducer',
    conversions = cms.InputTag('allConversions'),
    dimuons     = cms.InputTag('dimuon'),
    pfphotons   = cms.InputTag('pfPhotons'),
    convAlgo    = cms.string('mixed'),
#    convQuality = cms.vstring('highPurity'),
    convQuality = cms.vstring(''),
    convSelection = cms.string('conversionVertex.position.rho>0.5 ')
    )

process.PostDump = cms.EDAnalyzer(
    'DumpConversion',
    conversions = cms.InputTag('conversionProducer','conversions')
    )

process.chiCandProducer = cms.EDProducer(
    'ChiCandProducer',
    conversions = cms.InputTag('conversionProducer','conversions'),
    dimuons     = cms.InputTag('dimuon'),
    pfphotons   = cms.InputTag('pfPhotons'),
    pfcandidates = cms.InputTag('particleFlow'),
    pi0SmallWindow   = cms.vdouble(0.130,0.140),
    pi0LargeWindow   = cms.vdouble(0.110,0.160),
    deltaMass   = cms.vdouble(0.0, 2.0),
    dzmax       = cms.double(0.5)
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('testdimuon.root'),
    outputCommands =  cms.untracked.vstring('drop *',
                                            'keep *_dimuon_*_*',
					    'keep *_*_conversions_*',
                                            'keep *_*_chicand_*',
					    'keep *_*_piZeroRejectCand_*',
                                            ),
    #SelectEvents = cms.untracked.PSet(
    #SelectEvents = cms.vstring('')
    #)
    )



process.load('PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi')

process.dimuonpath = cms.Path(process.patMuons*   #create PAT muons
                              process.selMuons*   #official muon selection
                              process.dimuon*     #dimuon producer
			      process.diMuonCount*
			      process.diMuonSelection*
#			      process.PreDump*
			      process.conversionProducer*
			      process.chiCandProducer)
#			      process.PostDump)

process.end        = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.dimuonpath,process.end)

# HACK ! Need to understand !
process.dimuon.addMuonlessPrimaryVertex = False
