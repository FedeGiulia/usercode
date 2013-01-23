####################################################################################
############################# CONFIGURATION PARAMETERS #############################
####################################################################################

debug_filename = 'debug_testAnalyzer.txt' # Name of the debug file
input_filename = 'file:/data1/vcogoni/CMSSW_5_3_6/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FE7159F0-33CF-E111-A303-002618FDA279.root' # Path and name of the file to be analyzed
tag_dimuon = 'dimuon' # Tag name of the dimuon collection as saved from process.muonSelection20xx
cut_dimuon_Mass_low = 8.5
cut_dimuon_Mass_high = 11.0
cut_dimuon_Pt_min = 6.0
tag_conversion = 'allConversions'
tag_pfPhotons = 'pfPhotons'
conv_algo = 'mixed'
conv_qual = 'highPurity'
tag_primary_vertex = 'offlinePrimaryVertices'
conv_vertex_rho = 1.5
conv_vtx_comp = True
conv_tk_vtx = 5
conv_inn_hits = True
conv_min_dof = 3
tag_chi_conv_prod = 'conversionProducer'
tag_chi_conv_lab = 'conversions'
tag_pfCandidates = 'particleFlow'
pi0_small_min = 0.130
pi0_small_max = 0.140
pi0_large_min = 0.110
pi0_large_max = 0.160
chi_deltaM_min = 0.0
chi_deltaM_max = 2.0
chi_dzMax = 0.5


############################# CONFIGURATION END ####################################
if(pi0_small_min < pi0_large_min or pi0_small_max > pi0_large_max or pi0_small_max < pi0_small_min or pi0_large_max < pi0_large_min):
	print "############################################################"
	print "Error! Check that pi0 rejection small windows is smaller than the large one and/or that the order of the window boundaries is correct!"
	print "############################################################"
	raise ValueError

import FWCore.ParameterSet.Config as cms
process = cms.Process('test')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Ponia.Configuration.MuonSelection')
process.load('MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_E_V31::All')

process.MessageLogger.destinations = cms.untracked.vstring(debug_filename)
process.MessageLogger.cout = process.MessageLogger.cerr
process.MessageLogger.cout.threshold = cms.untracked.string("DEBUG")

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(input_filename)
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

process.selMuons = process.muonSelection2012.clone()

process.dimuon = cms.EDProducer(
	'DiMuonCandProducer',
	muons=                    cms.InputTag('selMuons'),
	beamSpotTag =             cms.InputTag('offlineBeamSpot'),
	primaryVertexTag =        cms.InputTag('offlinePrimaryVertices'),
	addCommonVertex=          cms.bool(True),
	addMuonlessPrimaryVertex= cms.bool(True),
        dimuonSelection = cms.string( 'p4.M > {0} && p4.M < {1} && p4.Pt > {2}'.format(cut_dimuon_Mass_low, cut_dimuon_Mass_high, cut_dimuon_Pt_min) ),
	)

process.diMuonCount = cms.EDFilter(
    'CandViewCountFilter',
    src = cms.InputTag(tag_dimuon),
    minNumber = cms.uint32(1),
    filter = cms.bool(True)
    )

process.diMuonSelection = cms.EDFilter(
    'CandViewSelector',
    src = cms.InputTag(tag_dimuon),
    cut = cms.string('p4.M > {0} && p4.M < {1} && p4.Pt > {2}'.format(cut_dimuon_Mass_low, cut_dimuon_Mass_high, cut_dimuon_Pt_min) ),# Formatting for dimuon cut string
    filter = cms.bool(True)
    )

process.conversionProducer = cms.EDProducer(
    'PhotonConversionCandProducer',
    conversions = cms.InputTag(tag_conversion),
    dimuons     = cms.InputTag(tag_dimuon),
    pfphotons   = cms.InputTag(tag_pfPhotons),
    convAlgo    = cms.string(conv_algo),
    convQuality = cms.vstring(conv_qual),
    primaryVertexTag = cms.InputTag(tag_primary_vertex),
    convSelection = cms.string('conversionVertex.position.rho>{0}'.format(conv_vertex_rho) ),
    wantTkVtxCompatibility = cms.bool(conv_vtx_comp),
    sigmaTkVtxComp = cms.uint32(conv_tk_vtx),
    wantCompatibleInnerHits = cms.bool(conv_inn_hits),
    TkMinNumOfDOF = cms.uint32(conv_min_dof)
    )

process.chiCandProducer = cms.EDProducer(
    'ChiCandProducer',
    conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab),
    dimuons     = cms.InputTag(tag_dimuon),
    pfphotons   = cms.InputTag(tag_pfPhotons),
    pfcandidates = cms.InputTag(tag_pfCandidates),
    pi0SmallWindow   = cms.vdouble(pi0_small_min, pi0_small_max),
    pi0LargeWindow   = cms.vdouble(pi0_large_min, pi0_large_max),
    deltaMass   = cms.vdouble(chi_deltaM_min, chi_deltaM_max),
    dzmax       = cms.double(chi_dzMax)
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('testdimuon.root'),
    outputCommands =  cms.untracked.vstring('drop *',
                                            'keep *_dimuon_UpsilonCandLorentzVector_*',
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
			      process.conversionProducer*
			      process.chiCandProducer)

process.end        = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.dimuonpath,process.end)

# HACK ! Need to understand !
process.dimuon.addMuonlessPrimaryVertex = False
