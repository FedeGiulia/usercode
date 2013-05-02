import FWCore.ParameterSet.Config as cms

process = cms.Process('chisimreco')

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

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
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

process.generator = cms.EDProducer("Pythia6PtGun",
    maxEventsToPrint = cms.untracked.int32(5),
    pythiaPylistVerbosity = cms.untracked.int32(0),   #set to 1 for printout
    pythiaHepMCVerbosity = cms.untracked.bool(False), #set to True for printout   
    PGunParameters = cms.PSet(
        ParticleID = cms.vint32(445), #Chi2
        AddAntiParticle = cms.bool(False),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        MinPt = cms.double(50.0),
        MaxPt = cms.double(50.0001),
        MinEta = cms.double(-1.25),
        MaxEta = cms.double(1.25),
        OnOdd  = cms.bool(True),
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

process.gammaConversionFilter = cms.EDFilter("GammaConversionFilter",
                                              trackingTruth = cms.untracked.InputTag('mergedtruth', 'MergedTrackTruth'),
                                              ToBeSelected = cms.bool(True)
                                              )



import Onia.Chi2OniaGamma.VarOptions_cff as var
# get and parse the command line arguments
var.myOptions.parseArguments()

from Onia.Chi2OniaGamma.chic_step1_cfi import chicstep1
from Onia.Chi2OniaGamma.chic_step4_cfi import chicstep4

process=chicstep1(process,var)
process=chicstep4(process,var)

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

process.AODSIMEventContent.outputCommands.extend(chicOutputCommands)

process.chicout = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(outpfx+'_SelEvts.root'),
    #    outputCommands = ONIA2MUMUEventContent.outputCommands,
    outputCommands = process.AODSIMEventContent.outputCommands,
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('chisimreco')
       )
    )    

process.GlobalTag.globaltag = 'START42_V14A::All'
    
process.e = cms.EndPath(process.chicout)

process.chisimreco = cms.Path(process.generator*
                              process.genParticles*
        #                      process.genFilterSummary*
                              process.psim*
                              process.pdigi*
                              process.gammaConversionFilter*
                              process.SimL1Emulator*
                              process.DigiToRaw*
                              process.RawToDigi*
                              process.L1Reco*
                              process.reconstruction*
                              process.TrackRefitter*
                              process.chicstep1*
                              process.chicstep4    
                              )


process.schedule = cms.Schedule(process.chisimreco,process.e)
