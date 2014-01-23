import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('RootupleChic',
                          chi_cand = cms.InputTag("chiCandProducer","chicand"),
                          pi0_comb = cms.InputTag("chiCandProducer", "piZeroRejectCand"),
                          ups_cand = cms.InputTag("dimuonProducer", "UpsilonCandLorentzVector"),
                          refit1S = cms.InputTag("refit", "chiCandRefit"),
                          refit2S = cms.InputTag("refit2S", "chiCand2S"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC = cms.bool(False)
                          )
