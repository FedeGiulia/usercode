import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('RootupleChib',
                          chi_cand = cms.InputTag("chiCandProducer","chicand"),
                          pi0_comb = cms.InputTag("chiCandProducer", "piZeroRejectCand"),
                          ups_cand = cms.InputTag("dimuonProducer", "UpsilonCandLorentzVector"),
                          refit1S = cms.InputTag("refit1S", "chiCand1S"),
                          refit2S = cms.InputTag("refit2S", "chiCand2S"),
                          refit3S = cms.InputTag("refit3S", "chiCand3S"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
								  TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False)
                          )
