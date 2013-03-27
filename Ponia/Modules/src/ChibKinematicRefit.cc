// -*- C++ -*-
//
// Package:    ChibKinematicRefit
// Class:      ChibKinematicRefit
// 
/**\class ChibKinematicRefit ChibKinematicRefit.cc Ponia/ChibKinematicRefit/src/ChibKinematicRefit.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giulio Dujany
//         Created:  Tue Mar 26 14:51:59 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/EgammaCandidates/interface/Conversion.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h> // GD aggiunto questo include

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
#include <string>

//
// class declaration
//

class ChibKinematicRefit : public edm::EDProducer {
   public:
      explicit ChibKinematicRefit(const edm::ParameterSet&);
      ~ChibKinematicRefit();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag chi_cand_, conversions_;
  double upsilon_mass_;
  std::string product_name_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
ChibKinematicRefit::ChibKinematicRefit(const edm::ParameterSet& iConfig)
{
  chi_cand_ = iConfig.getParameter<edm::InputTag>("chi_cand");
  conversions_ = iConfig.getParameter<edm::InputTag>("conversions");
  upsilon_mass_ = iConfig.getParameter<double>("upsilon_mass");
  product_name_ = iConfig.getParameter<std::string>("product_name");


  //kinematic refit collections
  produces<pat::CompositeCandidateCollection>(product_name_);
 
  //now do what ever other initialization is needed
  
}


ChibKinematicRefit::~ChibKinematicRefit()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ChibKinematicRefit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> chiCandHandle;
  iEvent.getByLabel(chi_cand_ ,chiCandHandle);
  
  edm::Handle<reco::ConversionCollection> conversionsHandle;
  iEvent.getByLabel(conversions_ ,conversionsHandle);
  
  //Kinemati refit collection
  std::auto_ptr< pat::CompositeCandidateCollection > chicCompCandRefitColl(new pat::CompositeCandidateCollection);
  
 
  // Kinematic fit
  
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); //GD esetup -> iSetup
  
  int indexConversion=-1;
  
  for(pat::CompositeCandidateCollection::const_iterator chiCand=chiCandHandle->begin();chiCand!=chiCandHandle->end();++chiCand)   
    { 
      indexConversion++;
      reco::TrackRef JpsiTk[2]={
	( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon1") ) )->innerTrack(),
	( dynamic_cast<const pat::Muon*>(chiCand->daughter("dimuon")->daughter("muon2") ) )->innerTrack()
      };
      
      reco::TrackCollection convTracks;
      reco::ConversionRef conv(conversionsHandle,indexConversion);
      BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv->tracks()){
	convTracks.push_back(*tk);
      }
      
      std::vector<reco::TransientTrack> MuMuTT;
      MuMuTT.push_back((*theB).build(&JpsiTk[0]));
      MuMuTT.push_back((*theB).build(&JpsiTk[1]));
      
      std::vector<reco::TransientTrack> EETT;
      EETT.push_back((*theB).build(convTracks[0]));
      EETT.push_back((*theB).build(convTracks[1]));
      
      const ParticleMass zero_mass(0);
      float zero_sigma = 1E-6;
      
      KinematicParticleFactoryFromTransientTrack pFactory;
      
      std::vector<RefCountedKinematicParticle> PhotonParticles;
      PhotonParticles.push_back(pFactory.particle(EETT[0],zero_mass,float(0),float(0),zero_sigma));
      PhotonParticles.push_back(pFactory.particle(EETT[1],zero_mass,float(0),float(0),zero_sigma));
      
      //KinematicParticleVertexFitter fitter(ps);
      KinematicParticleVertexFitter fitter;
      
      RefCountedKinematicTree photonVertexFitTree;
      photonVertexFitTree = fitter.fit(PhotonParticles);
      
      if (!photonVertexFitTree->isValid()) 
	{ 
	  edm::ParameterSet pSet;
	  
	  pSet.addParameter<double>("maxDistance", 3);
	  pSet.addParameter<int>("maxNbrOfIterations", 10000); //10
	  
	  KinematicParticleVertexFitter fitter2(pSet);
	  
	  //RefCountedKinematicTree photonVertexFitTree;
	  photonVertexFitTree = fitter2.fit(PhotonParticles);
	}
      
      if (photonVertexFitTree->isValid()) 
	{ 
	  
	  // now apply Photon mass constraint                                                                                                            
	  KinematicParticleFitter csFitterPhoton;
	  KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);
	  // add mass constraint to the photon fit to do a constrained fit:                                                                              
	  
	  photonVertexFitTree->movePointerToTheTop();
	  photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
	  
	  if (photonVertexFitTree->isValid()) 
	    {
	      const ParticleMass muonMass(0.1056583);
	      float muonSigma = muonMass*1E-6;
	      
	      photonVertexFitTree->movePointerToTheTop();
	      RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
	      
	      std::vector<RefCountedKinematicParticle> allChiBDaughters;
	      allChiBDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
	      allChiBDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
	      allChiBDaughters.push_back(fitPhoton);
	      
	      KinematicConstrainedVertexFitter constVertexFitter;
	      
	      MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
	      RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters,upsilon_mtc);
	      if(!ChiBTree->isEmpty()) 
		{
		  
		  ChiBTree->movePointerToTheTop();
		  RefCountedKinematicParticle fitChiB = ChiBTree->currentParticle();
		  RefCountedKinematicVertex ChiBDecayVertex = ChiBTree->currentDecayVertex();
		  
		  if (fitChiB->currentState().isValid()) 
		    {          
		      float ChiBM_fit = fitChiB->currentState().mass();
		      float ChiBPx_fit = fitChiB->currentState().kinematicParameters().momentum().x();
		      float ChiBPy_fit = fitChiB->currentState().kinematicParameters().momentum().y();
		      float ChiBPz_fit = fitChiB->currentState().kinematicParameters().momentum().z();
		      float ChiBVtxX_fit = ChiBDecayVertex->position().x();
		      float ChiBVtxY_fit = ChiBDecayVertex->position().y();
		      float ChiBVtxZ_fit = ChiBDecayVertex->position().z();
		      float ChiBVtxP_fit=ChiSquaredProbability((double)(ChiBDecayVertex->chiSquared()),(double)(ChiBDecayVertex->degreesOfFreedom()));
		      
		      reco::CompositeCandidate aChicCandRefit(0,math::XYZTLorentzVector(ChiBPx_fit,ChiBPy_fit,ChiBPz_fit,sqrt(ChiBM_fit*ChiBM_fit+ChiBPx_fit*ChiBPx_fit+ChiBPy_fit*ChiBPy_fit+ChiBPz_fit*ChiBPz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),50551);
		      
		      pat::CompositeCandidate patChicCandRefit(aChicCandRefit);
		      patChicCandRefit.addUserFloat("vProb",ChiBVtxP_fit);
		      chicCompCandRefitColl->push_back(patChicCandRefit);
		      
		    }
		}
	    }
	}
    }  
  // End kinematic fit
  
  iEvent.put(chicCompCandRefitColl,product_name_); //GD event -> iEvent
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
ChibKinematicRefit::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ChibKinematicRefit::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ChibKinematicRefit::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ChibKinematicRefit::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ChibKinematicRefit::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ChibKinematicRefit::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ChibKinematicRefit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ChibKinematicRefit);
