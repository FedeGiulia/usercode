import FWCore.ParameterSet.Config as cms

chibKinematicRefit = cms.EDProducer('ChibKinematicRefit',
                                    chi_cand = cms.InputTag("chiCandProducer","chicand"),
                                    conversions = cms.InputTag("chiCandProducer","chiConversions"),
                                    upsilon_mass = cms.double(9.4603), # GeV   1S = 9.4603   2S = 10.02326    3S = 10.3552
                                    product_name = cms.string("chiCand_1S")
                                    )
