####################################################################################
############################# CONFIGURATION PARAMETERS #############################
####################################################################################

#debug_filename = 'debug_testAnalyzer.txt' # Name of the debug file
input_filename = 'file:/data1/vcogoni/CMSSW_5_3_6/store/data/Run2012A/MuOnia/AOD/13Jul2012-v1/00000/FE7159F0-33CF-E111-A303-002618FDA279.root'
tag_dimuon = 'dimuonProducer' # Tag name of the dimuon collection as saved from process.dimuonProducer
cut_dimuon_Mass_low = 8.5
cut_dimuon_Mass_high = 11.0
cut_dimuon_Pt_min = 0.0
cut_dimuon_rapidity = 1.6
cut_dimuon_vprob = 0.01 # Minimum vertex probability for dimuon candidate
tag_conversion = 'allConversions'
tag_pfPhotons = 'pfPhotons'
conv_algo = 'undefined'
conv_qual = ['highPurity','generalTracksOnly']
tag_primary_vertex = 'offlinePrimaryVertices'
conv_vertex_rho = 1.5
conv_vtx_comp = False
conv_tk_vtx = 5
conv_inn_hits = True
conv_min_dof = 3
tag_chi_conv_prod = 'conversionProducer'
tag_chi_conv_lab = 'conversions'
tag_pfCandidates = 'particleFlow'
pi0_online_switch = True
pi0_small_min = 0.130
pi0_small_max = 0.140
pi0_large_min = 0.110
pi0_large_max = 0.160
chi_deltaM_min = 0.0    # This two values define the minimum and the maximum values 
chi_deltaM_max = 2.0    # required for the QValue of the Chi candidate
chi_dzMax = 0.5
upsilon_masses = [9.4603, 10.02326, 10.3552]

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
    input = cms.untracked.int32(1000)
    )

process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
      'HLT_Dimuon5_Upsilon_v3',
      'HLT_Dimuon7_Upsilon_v3',
      'HLT_Dimuon8_Upsilon_v3',
      'HLT_Dimuon11_Upsilon_v3',
),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    l1tResults = cms.InputTag( "gtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

process.selMuons = process.muonSelection2012.clone()

process.dimuonProducer = cms.EDProducer(
	'DiMuonCandProducer',
	muons=                    cms.InputTag('selMuons'),
	beamSpotTag =             cms.InputTag('offlineBeamSpot'),
	primaryVertexTag =        cms.InputTag('offlinePrimaryVertices'),
	addCommonVertex=          cms.bool(True),
	addMuonlessPrimaryVertex= cms.bool(True),
        dimuonSelection = cms.string( 'p4.M > {0} && p4.M < {1} && p4.Pt > {2} && abs(y) < {3} && userFloat("vProb") > {4}'.format(cut_dimuon_Mass_low, cut_dimuon_Mass_high, cut_dimuon_Pt_min, cut_dimuon_rapidity, cut_dimuon_vprob) ),
	)

process.diMuonCount = cms.EDFilter(
    'CandViewCountFilter',
    src = cms.InputTag('dimuonProducer'),
    minNumber = cms.uint32(1),
    filter = cms.bool(True)
    )

process.conversionProducer = cms.EDProducer(
    'PhotonConversionCandProducer',
    conversions = cms.InputTag(tag_conversion),
    pfphotons   = cms.InputTag(tag_pfPhotons),
    convAlgo    = cms.string(conv_algo),
    convQuality = cms.vstring(conv_qual),
    primaryVertexTag = cms.InputTag(tag_primary_vertex),
    convSelection = cms.string('conversionVertex.position.rho>{0}'.format(conv_vertex_rho) ),
    wantTkVtxCompatibility = cms.bool(conv_vtx_comp),
    sigmaTkVtxComp = cms.uint32(conv_tk_vtx),
    wantCompatibleInnerHits = cms.bool(conv_inn_hits),
    TkMinNumOfDOF = cms.uint32(conv_min_dof),
    wantHighpurity = cms.bool(False),
    vertexChi2ProbCut = cms.double(0.0005),
    trackchi2Cut = cms.double(10),
    minDistanceOfApproachMinCut = cms.double(-0.25),
    minDistanceOfApproachMaxCut = cms.double(1.00)
    )

process.chiCandProducer = cms.EDProducer(
    'ChiCandProducer',
    conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab),
    dimuons     = cms.InputTag('dimuonProducer'),
    pfphotons   = cms.InputTag(tag_pfPhotons),
    pfcandidates = cms.InputTag(tag_pfCandidates),
    pi0OnlineSwitch = cms.bool(pi0_online_switch),
    pi0SmallWindow   = cms.vdouble(pi0_small_min, pi0_small_max),
    pi0LargeWindow   = cms.vdouble(pi0_large_min, pi0_large_max),
    deltaMass   = cms.vdouble(chi_deltaM_min, chi_deltaM_max),
    dzmax       = cms.double(chi_dzMax)
    )



process.refit1S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[0]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand1S")
				 )
process.refit2S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[1]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand2S")
				 )
process.refit3S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[2]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand3S")
				 )




process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('testdimuon.root'),
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



process.load('PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi')

process.dimuonpath = cms.Path(process.triggerSelection*
			      process.patMuons*   #create PAT muons
                              process.selMuons*   #official muon selection
                              process.dimuonProducer*     #dimuon producer
			      process.diMuonCount*
			      process.conversionProducer*
			      process.chiCandProducer*
			      process.refit1S*
			      process.refit2S*
			      process.refit3S
			      )

process.end        = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.dimuonpath,process.end)

# HACK ! Need to understand !
process.dimuonProducer.addMuonlessPrimaryVertex = False
