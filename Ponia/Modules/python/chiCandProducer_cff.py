import FWCore.ParameterSet.Config as cms
####################################################################################
############################# CONFIGURATION PARAMETERS #############################
####################################################################################

tag_dimuon = 'dimuonProducer' # Tag name of the dimuon collection as saved from process.dimuonProducer
cut_dimuon_Mass_low = 8.5
cut_dimuon_Mass_high = 11.0
cut_dimuon_Pt_min = 0.0
cut_dimuon_rapidity = 1.6
cut_dimuon_vprob = 0.01 # Minimum vertex probability for dimuon candidate
tag_conversion = 'allConversions'
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
pi0_online_switch = False
pi0_small_min = 0.130
pi0_small_max = 0.140
pi0_large_min = 0.110
pi0_large_max = 0.160
chi_deltaM_min = 0.0    # This two values define the minimum and the maximum values 
chi_deltaM_max = 2.0    # required for the QValue of the Chi candidate
chi_dzMax = 0.5
upsilon_masses = [9.4603, 10.02326, 10.3552]

############################# CONFIGURATION END ####################################

triggerSelection = cms.EDFilter( "TriggerResultsFilter",
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

from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons

from Ponia.Configuration.MuonSelection import muonSelection2012
selMuons = muonSelection2012.clone()

dimuonProducer = cms.EDProducer(
	'DiMuonCandProducer',
	muons=                    cms.InputTag('selMuons'),
	beamSpotTag =             cms.InputTag('offlineBeamSpot'),
	primaryVertexTag =        cms.InputTag('offlinePrimaryVertices'),
	addCommonVertex=          cms.bool(True),
	addMuonlessPrimaryVertex= cms.bool(True),
        dimuonSelection = cms.string( 'p4.M > {0} && p4.M < {1} && p4.Pt > {2} && abs(y) < {3} && userFloat("vProb") > {4}'.format(cut_dimuon_Mass_low, cut_dimuon_Mass_high, cut_dimuon_Pt_min, cut_dimuon_rapidity, cut_dimuon_vprob) ),
	)

diMuonCount = cms.EDFilter(
    'CandViewCountFilter',
    src = cms.InputTag('dimuonProducer'),
    minNumber = cms.uint32(1),
    filter = cms.bool(True)
    )

conversionProducer = cms.EDProducer(
    'PhotonConversionCandProducer',
    conversions = cms.InputTag(tag_conversion),
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

chiCandProducer = cms.EDProducer(
    'ChiCandProducer',
    conversions = cms.InputTag(tag_chi_conv_prod, tag_chi_conv_lab),
    dimuons     = cms.InputTag('dimuonProducer'),
    pfcandidates = cms.InputTag(tag_pfCandidates),
    pi0OnlineSwitch = cms.bool(pi0_online_switch),
    pi0SmallWindow   = cms.vdouble(pi0_small_min, pi0_small_max),
    pi0LargeWindow   = cms.vdouble(pi0_large_min, pi0_large_max),
    deltaMass   = cms.vdouble(chi_deltaM_min, chi_deltaM_max),
    dzmax       = cms.double(chi_dzMax)
    )



refit1S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[0]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand1S")
				 )
refit2S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[1]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand2S")
				 )
refit3S = cms.EDProducer('ChibKinematicRefit',
				 chi_cand = cms.InputTag("chiCandProducer","chicand"),
				 conversions = cms.InputTag("chiCandProducer","chiConversions"),
				 upsilon_mass = cms.double(upsilon_masses[2]), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
				 product_name = cms.string("chiCand3S")
				 )


chiPath = cms.Path(triggerSelection*
			      patMuons*   #create PAT muons
                              selMuons*   #official muon selection
                              dimuonProducer*     #dimuon producer
			      diMuonCount*
			      conversionProducer*
			      chiCandProducer*
			      refit1S*
			      refit2S*
			      refit3S
			      )
