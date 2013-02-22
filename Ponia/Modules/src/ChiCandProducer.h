/**
   \file
   Declaration of ChiCandProducer

   \author Stefano Argiro
   \version $Id: ChiCandProducer.h,v 1.4 2013/02/14 09:39:39 degano Exp $
   \date 18 Dec 2012
*/

#ifndef __ChiCandProducer_h_
#define __ChiCandProducer_h_


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include <TLorentzVector.h>
#include <vector>

///For kinematic fit:
/*
#include <DataFormats/PatCandidates/interface/Muon.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include <boost/foreach.hpp>
*/

static const char CVSId__ChiCandProducer[] = 
"$Id: ChiCandProducer.h,v 1.4 2013/02/14 09:39:39 degano Exp $";

/**
   Create a Chi(b,c) candidate by mathing dimuon and conversion
 */
class ChiCandProducer : public edm::EDProducer {

 public:
  explicit ChiCandProducer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);
  
  virtual void endJob();

  edm::InputTag convCollection_;
  edm::InputTag diMuonCollection_;
  edm::InputTag pfCandidateCollection_;


  // create a collection of PF photons
  const reco::PFCandidateCollection  selectPFPhotons(const reco::PFCandidateCollection& 
						     pfcandidates);
  
  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);

  const pat::CompositeCandidate makePhotonCandidate(const reco::Conversion& conv);

  const pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate& dimuon, 
						 const pat::CompositeCandidate& photon);

  // check if the mass difference is in desired range
  bool cutDeltaMass(const pat::CompositeCandidate& chiCand,
		  const pat::CompositeCandidate& dimuonCand); 

  bool cutdz(float dz){return dz<dzMax_; }


  /** return a vector containing the invariant mass combinations 
    * of the conversion with PF photons.
    * Only combinations within the given window are retained */
  std::vector<float> invmCombinations(const reco::Conversion& conv,
				      const reco::PFCandidateCollection& photons,
				      bool* _pizero_rejected );

  // low and high window limits
  std::vector<double> pi0SmallWindow_;
  std::vector<double> pi0LargeWindow_;
  
  // delta mass range
  std::vector<double> deltaMass_;
  double dzMax_;

  int candidates;
  int delta_mass_fail;
  int dz_cut_fail;
  int pizero_fail;

  typedef std::vector<float> fl_vect;
  typedef std::vector<fl_vect> PiZeroRejCollection;
};



#endif // __ChiCandProducer_h_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "cd .. ; scram b"
// End:
