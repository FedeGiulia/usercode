import FWCore.ParameterSet.Config as cms

dimuon = cms.EDProducer('DiMuonCandProducer',
                        muons=                    cms.InputTag('selMuons'),
                        beamSpotTag =             cms.InputTag('offlineBeamSpot'),
                        primaryVertexTag =        cms.InputTag('offlinePrimaryVertices'),
                        addCommonVertex=          cms.bool(True),
                        addMuonlessPrimaryVertex= cms.bool(True),
                        
                        theTriggerNames = cms.vstring('HLT_Dimuon7_Upsilon_v3',
                                                      'HLT_Dimuon7_Upsilon_v4',
                                                      'HLT_Dimuon7_Upsilon_v5',
                                                      'HLT_Dimuon7_Upsilon_v6',
                                                      'HLT_Dimuon7_Upsilon_v7'),
                        
                        HLTLastFilters = cms.vstring('hltVertexmumuFilterDimuon7Upsilon',
                                                     'hltVertexmumuFilterDimuon7Upsilon',
                                                     'hltVertexmumuFilterDimuon7Upsilon',
                                                     'hltVertexmumuFilterDimuon7Upsilon',
                                                     'hltVertexmumuFilterDimuon7Upsilon',) #HLT_Dimuon7_Upsilon_v{3,4,5,6,7}
                        )                        
