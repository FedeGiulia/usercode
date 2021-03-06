#
# Standard muon selection from muon POG
# $Id$
#

import FWCore.ParameterSet.Config as cms

#2012
muonSelection2012 = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('patMuons'),        
   cut = cms.string('muonID(\"TMOneStationTight\")'                    
                    '&& abs(innerTrack.dxy) < 3.0'                            
                    '&& abs(innerTrack.dz)  < 30'
                    '&& track.hitPattern.trackerLayersWithMeasurement > 5'
                    '&& innerTrack.hitPattern.pixelLayersWithMeasurement > 1'
                    '&& innerTrack.normalizedChi2 < 1.8' 
                    ),
   filter = cms.bool(True),
     
)

#2011
muonSelection2011 = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('patMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    '&& abs(innerTrack.dxy) < 3.0'
                    '&& abs(innerTrack.dz)  < 30'
                    '&& innerTrack.hitPattern.numberOfValidTrackerHits> 10'
                    '&& innerTrack.hitPattern.pixelLayersWithMeasurement > 1'
                    '&& innerTrack.normalizedChi2 < 1.8' 
                    ),
   filter = cms.bool(True),
)
