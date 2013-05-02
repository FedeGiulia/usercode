import FWCore.ParameterSet.Config as cms

process = cms.Process('chic')

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("Configuration.Generator.PythiaUESettings_cfi")

# Event output
process.load("Configuration.EventContent.EventContent_cff")

#process.load('Configuration.StandardSequences.GeometryDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('RecoTracker.TrackProducer.TrackRefitter_cfi')
process.load('GammaConversion.MCMatching.FakePrimaryVertex_cfi')

#material budget studies
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#HLT
process.load('HLTrigger.Configuration.HLT_quarkonium_1E33_3E33_cff')

process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    )

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#    generator = cms.PSet(
#        initialSeed = cms.untracked.uint32(123456789),
#        engineName = cms.untracked.string('HepJamesRandom')
#    )
#)



process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("Pythia6CustomPtGunTFunction",
    maxEventsToPrint = cms.untracked.int32(5),
    pythiaPylistVerbosity = cms.untracked.int32(0),   #set to 1 for printout
    pythiaHepMCVerbosity = cms.untracked.bool(False), #set to True for printout   
    PGunParameters = cms.PSet(
        ParticleID = cms.vint32(445,20443), #Chi2
        AddAntiParticle = cms.bool(False),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        MinPt = cms.double(0.0),
        MaxPt = cms.double(40.0),
        MinEta = cms.double(-1.25),
        MaxEta = cms.double(1.25),
        OnOdd  = cms.bool(True),
        TFunction_string = cms.string("x*pow(1.+1./({0}-2.)*x*x/{1},-{0})".format(3.7,19.36)),
        TFunction_min = cms.double("5"),
   	  TFunction_max = cms.double("40")
        ),
     PythiaParameters = cms.PSet(
        process.pythiaUESettingsBlock,
        pythiaJpsiDecays = cms.vstring(
            'MSEL=61                       ! Quarkonia', 
            'MDME(858,1)=0                 ! J/psi -> ee turned OFF', 
            'MDME(859,1)=1                 ! J/psi -> mumu turned ON', 
            'MDME(860,1)=0                 ! J/psi -> random turned OFF',
            'BRAT(861)=1.0                 ! chi_2c->J/psi gamma', 
            'BRAT(862)=0.0                 ! chi_2c->rndmflav rndmflavbar', 
            'BRAT(1501)=0.0                ! chi_0c->J/psi gamma', 
            'BRAT(1502)=0.0                ! chi_0c->rndmflav rndmflavbar', 
            'BRAT(1555)=1.0                ! chi_1c->J/psi gamma', 
            'BRAT(1556)=0.0                ! chi_1c->rndmflav rndmflavbar'
            ),
           
     # This is a vector of ParameterSet names to be read, in this order
     parameterSets = cms.vstring('pythiaUESettings', 
                                 'pythiaJpsiDecays'
                                 )
     )
                                   
)

print process.generator.PGunParameters.TFunction_string

process.gammaConversionFilter = cms.EDFilter("GammaConversionFilter",
                                              trackingTruth = cms.untracked.InputTag('mergedtruth', 'MergedTrackTruth'),
                                              ToBeSelected = cms.bool(True)
                                              )



import Onia.Chi2OniaGamma.VarOptions_cff as var
# get and parse the command line arguments
var.myOptions.parseArguments()

from Onia.Chi2OniaGamma.chic_step1_cfi import chicstep1
from Onia.Chi2OniaGamma.chic_step4_cfi import chicstep4
#from Onia.Chi2OniaGamma.chic_step5_cfi import chicstep5


process=chicstep1(process,var)
process=chicstep4(process,var)
#process=chicstep5(process,var)


from Onia.Chi2OniaGamma.chic_step5_selection_cfi import chic_step5_selection
process=chic_step5_selection(process)

process.chicstep5 = process.chicstep4

process.chic.fillJpsiVsNtracks = True
process.chic.maxAbsZ=15


outpfx = var.myOptions.outprefix

from Onia.Chi2OniaGamma.ChicEventContent_cff import *

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outpfx+'_Vertices.root')
                                   )

chicOutputCommands = cms.untracked.vstring('keep recoVertexs_*_*_*',             ## take all vertices
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep patMuons_patMuons_*_*',         ## All PAT muons 
                                           'keep patCompositeCandidates_*_*_*',  ## PAT di-muons
                                           'keep recoConversions_*_*_*',
                                           'keep recoCompositeCandidates_*_*_*', ## jpsi,chi, and gamma from chic analyzer
                                           'keep *_chic_*_*',                    ## well, everything from chic analyzer
                                           )

RAWSIMEventContent.outputCommands.extend(ChicAnalysisEventContent.outputCommands)

process.chicout = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(outpfx+'_SelEvts.root'),
#    outputCommands = ChicAnalysisEventContent.outputCommands,    
    outputCommands = RAWSIMEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('reco')
       )
    )


process.genout = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('genchi.root'),
    outputCommands = cms.untracked.vstring('drop *',
					    'keep *_genParticles_*_*'),
    SelectEvents = cms.untracked.PSet(
	SelectEvents = cms.vstring('gen')
	)
    )


process.simout = cms.OutputModule(
	"PoolOutputModule",
	fileName = cms.untracked.string('convertedchi.root'),
	outputCommands = cms.untracked.vstring('drop *',
					       'keep *_genParticles_*_*',
					       #'keep *_mergedtruth_*_*'
					       ),
	SelectEvents = cms.untracked.PSet(
	SelectEvents = cms.vstring('sim')
	)
	)


process.GlobalTag.globaltag = 'START42_V14A::All'
#process.GlobalTag.globaltag = 'MC_42_V15B::All'


#HLT Global Tag
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuGMTParametersRcd' ),
        tag     = cms.string( 'L1MuGMTParameters_synctf_10_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
    ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuDTTFParametersRcd' ),
        tag     = cms.string( 'L1MuDTTFParameters_dttf11_TSC_09_17_col_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1MuCSCTFConfigurationRcd' ),
        tag     = cms.string( 'L1MuCSCTFConfiguration_90511_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCBxOrConfigRcd' ),
        tag     = cms.string( 'L1RPCBxOrConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConeDefinitionRcd' ),
        tag     = cms.string( 'L1RPCConeDefinition_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCConfigRcd' ),
        tag     = cms.string( 'L1RPCConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        ))
process.GlobalTag.toGet.append(
    cms.PSet(
        record  = cms.string( 'L1RPCHsbConfigRcd' ),
        tag     = cms.string( 'L1RPCHsbConfig_LHC7_1EX_mc' ),
        label   = cms.untracked.string( '' ),
        connect = cms.untracked.string('frontier://FrontierProd/CMS_COND_31X_L1T')
        )
    )


#process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
#    triggerConditions = cms.vstring(
#      'HLT_Dimuon5_Upsilon_Barrel_v3',
#      'HLT_Dimuon10_Jpsi_Barrel_v3',
#      'HLT_Dimuon7_Jpsi_X_Barrel_v3',
#      'HLT_Dimuon0_Jpsi_Muon_v1',
#      'HLT_Dimuon0_Upsilon_Muon_v1',
#      'HLT_Dimuon0_Jpsi_v3',
#      'HLT_Dimuon0_Upsilon_v3',
#      'HLT_Dimuon5_Upsilon_Barrel_v5',
#      'HLT_Dimuon10_Jpsi_Barrel_v5',
#      'HLT_Dimuon7_Jpsi_X_Barrel_v5',
#      'HLT_Dimuon0_Jpsi_Muon_v6',
#      'HLT_Dimuon0_Upsilon_Muon_v6',
#      'HLT_Dimuon0_Jpsi_v5',
#      'HLT_Dimuon0_Upsilon_v5',
#      'HLT_Dimuon0_Jpsi_NoVertexing_v2',
#      'HLT_Dimuon7_Upsilon_Barrel_v1',
#      'HLT_Dimuon9_Upsilon_Barrel_v1',
#      'HLT_Dimuon10_Jpsi_Barrel_v6',
#      'HLT_Dimuon13_Jpsi_Barrel_v1',
#      'HLT_Dimuon0_Jpsi_Muon_v7',
#      'HLT_Dimuon0_Upsilon_Muon_v7',
#      'HLT_Dimuon0_Jpsi_v6',
#      'HLT_Dimuon0_Upsilon_v6',
#      'HLT_Dimuon0_Jpsi_NoVertexing_v3',
#      'HLT_Mu5_L2Mu2_Jpsi_v8',
#      'HLT_Mu5_L2Mu2_Jpsi_v9',
#      'HLT_Mu5_Track2_Jpsi_v9',
#      'HLT_Mu7_Track7_Jpsi_v10',
#      'HLT_Dimuon7_PsiPrime_v3',
#      'HLT_Dimuon7_PsiPrime_v5',
#      'HLT_Dimuon9_PsiPrime_v1',
#      'HLT_Dimuon11_PsiPrime_v1'),
#    hltResults = cms.InputTag("TriggerResults"),
#    l1tResults = cms.InputTag(""),
#    l1tIgnoreMask = cms.bool( False ),
#    l1techIgnorePrescales = cms.bool( False ),
#    daqPartitions = cms.uint32( 1 ),
#    throw = cms.bool( True )
#)


#process.mrtrigger = cms.EDAnalyzer('MrTrigger'
#)

process.endreco = cms.EndPath(process.chicout)
process.endsim = cms.EndPath(process.simout)
process.endgen  = cms.EndPath(process.genout)

#custom for particle gun

process.chic.vtxMinDoF = -1. # we have no real vertex with particle gun
process.onia2MuMuPatTrkTrk.addMuonlessPrimaryVertex = False


process.gen = cms.Path(process.generator*
		       process.VtxSmeared*
		       process.genParticles)


        #                      process.genFilterSummary*
	
process.sim = cms.Path(process.psim*
		       process.pdigi*
		       process.gammaConversionFilter
		       )
		       
process.reco= cms.Path(
                       process.psim*
		       process.pdigi*
		       process.gammaConversionFilter*
	               process.SimL1Emulator*
		       process.DigiToRaw*   
	               process.RawToDigi*
		       process.pixelVertices*   # this should be in std reco ??
		       process.L1Reco*
		       process.reconstruction*
		       process.TrackRefitter*    # why do we need this ?
		       process.chicstep1*
                       #process.triggerSelection
                       #process.mrtrigger
		       process.chicstep5
               #process.chicstep5        
		       )
process.PreHLT= cms.Path(
                       process.psim*
		       process.pdigi*
		       #process.gammaConversionFilter*
	               process.SimL1Emulator*
		       process.DigiToRaw)




process.schedule = cms.Schedule(process.gen,
				process.sim,
				process.PreHLT)

process.schedule.extend(process.HLTSchedule)
process.schedule.append(process.reco)
process.schedule.append(process.endgen)
process.schedule.append(process.endsim)
process.schedule.append(process.endreco)

def fakeVertexing(process):
        process.offlinePrimaryVertices = process.offlinePrimaryVerticesFake
        process.offlinePrimaryVerticesWithBS = process.offlinePrimaryVerticesWithBSFake
        process.pixelVertices = process.pixelVerticesFake
        process.recopixelvertexing =    process.recopixelvertexingFake
        return process



process = fakeVertexing(process)

# customisation of the process.
process.es_prefer_magneticfield = cms.ESPrefer("VolumeBasedMagneticFieldESProducer","")
process.es_prefer_calotowergeometry = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
process.es_prefer_castorgeometry = cms.ESPrefer("CastorGeometryFromDBEP","")
process.es_prefer_ecalbarrelgeometry = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
process.es_prefer_ecalendcapgeometry = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
process.es_prefer_ecalpreshowergeometry = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
process.es_prefer_hcalgeometry = cms.ESPrefer("HcalGeometryFromDBEP","")
process.es_prefer_zdcgeometry = cms.ESPrefer("ZdcGeometryFromDBEP","")


# Automatic addition of the customisation function from Geometry.TrackerGeometryBuilder.customTrackerLiMax
#from Geometry.TrackerGeometryBuilder.customTrackerLiMax import customise 

#call to customisation function customise imported from Geometry.TrackerGeometryBuilder.customTrackerLiMax
#process = customise(process)

# End of customisation functions
