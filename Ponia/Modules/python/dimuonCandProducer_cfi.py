import FWCore.ParameterSet.Config as cms

dimuon = cms.EDProducer('DiMuonCandProducer',
                        muons=                    cms.InputTag('selMuons'),
                        beamSpotTag =             cms.InputTag('offlineBeamSpot'),
                        primaryVertexTag =        cms.InputTag('offlinePrimaryVertices'),
                        addCommonVertex=          cms.bool(True),
                        addMuonlessPrimaryVertex= cms.bool(True)
                        )
