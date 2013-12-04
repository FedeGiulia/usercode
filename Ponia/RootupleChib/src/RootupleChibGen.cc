// -*- C++ -*-
//
// Package:    RootupleChibGen
// Class:      RootupleChibGen
// 
/**\class RootupleChibGen RootupleChibGen.cc RootupleChibGen/RootupleChibGen/src/RootupleChibGen.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giulio Dujany,32 1-C13,+41227678098
//         Created:  Fri Mar 08 16:32:37 CET 2013
// $Id: RootupleChibGen.cc,v 1.1 2013/07/05 08:38:18 gdujany Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>

//
// class declaration
//

class RootupleChibGen : public edm::EDAnalyzer 
{
public:
  explicit RootupleChibGen(const edm::ParameterSet&);
  ~RootupleChibGen();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // Tags  
  //std::string file_name;
  
  // Variables
  TLorentzVector chib_p4;
  Double_t chib_mass;
  Double_t chib_pt;
  Double_t chib_eta;
  Double_t chib_rapidity;
  Double_t chib_phi;
  Int_t chib_pdgId;
  TLorentzVector Upsilon_p4;
  Double_t Upsilon_mass;
  Double_t Upsilon_pt;
  Double_t Upsilon_eta;
  Double_t Upsilon_rapidity;
  Double_t Upsilon_phi;
  Int_t Upsilon_pdgId;
  TLorentzVector photon_p4;
  Double_t photon_pt;
  Double_t photon_eta;
  Double_t photon_phi;
  TLorentzVector muP_p4;
  TLorentzVector muM_p4;
  
  TTree* genParticles_tree;
  
};

//
// constants, enums and typedefs
//


//
// constructors and destructor
//
RootupleChibGen::RootupleChibGen(const edm::ParameterSet& iConfig)
{
	edm::Service<TFileService> fs;
	genParticles_tree = fs->make<TTree>("GenParticlesTree","Tree of genParticles");

	genParticles_tree->Branch("chib_p4", "TLorentzVector", &chib_p4);
	genParticles_tree->Branch("chib_mass", &chib_mass, "chib_mass/D");
	genParticles_tree->Branch("chib_pt", &chib_pt, "chib_pt/D");
	genParticles_tree->Branch("chib_eta", &chib_eta, "chib_eta/D");
	genParticles_tree->Branch("chib_rapidity", &chib_rapidity, "chib_rapidity/D");
	genParticles_tree->Branch("chib_phi", &chib_phi, "chib_phi/D");
	genParticles_tree->Branch("chib_pdgId", &chib_pdgId, "chib_pdgId/I");
	genParticles_tree->Branch("Upsilon_p4", "TLorentzVector", &Upsilon_p4);
	genParticles_tree->Branch("Upsilon_mass", &Upsilon_mass, "Upsilon_mass/D");
	genParticles_tree->Branch("Upsilon_pt", &Upsilon_pt, "Upsilon_pt/D");
	genParticles_tree->Branch("Upsilon_eta", &Upsilon_eta, "Upsilon_eta/D");
	genParticles_tree->Branch("Upsilon_rapidity", &Upsilon_rapidity, "Upsilon_rapidity/D");
	genParticles_tree->Branch("Upsilon_phi", &Upsilon_phi, "Upsilon_phi/D");
	genParticles_tree->Branch("Upsilon_pdgId", &Upsilon_pdgId, "Upsilon_pdgId/I");
	genParticles_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);
	genParticles_tree->Branch("photon_pt", &photon_pt, "photon_pt/D");
	genParticles_tree->Branch("photon_eta", &photon_eta, "photon_eta/D");
	genParticles_tree->Branch("photon_phi", &photon_phi, "photon_phi/D");
	genParticles_tree->Branch("muP_p4", "TLorentzVector", &muP_p4);
	genParticles_tree->Branch("muM_p4", "TLorentzVector", &muM_p4);
	
}


RootupleChibGen::~RootupleChibGen()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RootupleChibGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel("genParticles",GenParticles);   //Open GenParticles
  
  if(GenParticles.isValid() )
    {
      for(reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle)
	{
	  Int_t pdgId = itParticle->pdgId();
	  
	  if(pdgId == 10551 || pdgId == 110551 || pdgId == 210551 || //chib0 1,2,3 P
	     pdgId == 20553 || pdgId == 120553 || pdgId == 220553 || //chib1 1,2,3 P
	     pdgId ==   555 || pdgId == 100555 || pdgId == 200555 )  //chib2 1,2,3 P
	    {
	      chib_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
	      chib_pdgId = pdgId; 
	      chib_mass = chib_p4.M();
	      chib_pt = chib_p4.Pt();
	      chib_eta = chib_p4.Eta();
	      chib_rapidity = chib_p4.Rapidity();
	      chib_phi = chib_p4.Phi();
	      continue;
	    }
	  
	  if(pdgId ==    553 || // Upsilon 1S
	     pdgId == 100553 || // Upsilon 2S
	     pdgId == 200553 )  // Upsilon 3S
	    {
	      Upsilon_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
	      Upsilon_pdgId = pdgId; 
	      Upsilon_mass = Upsilon_p4.M();
	      Upsilon_pt = Upsilon_p4.Pt();
	      Upsilon_eta = Upsilon_p4.Eta();
	      Upsilon_rapidity = Upsilon_p4.Rapidity();
	      Upsilon_phi = Upsilon_p4.Phi();
	      continue;
	    }
	  
	  if( pdgId == 13 ) // mu-
	    {
	      muM_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
	      continue;
	    }

	  if( pdgId == -13 ) // mu+
	    {
	      muP_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
	      continue;
	    }

	  if( pdgId == 22 ) // photon
	    {
	      photon_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
	      photon_pt = photon_p4.Pt();
	      photon_eta = photon_p4.Eta();
	      photon_phi = photon_p4.Phi();
	    }

	}
    }
  
  genParticles_tree ->Fill();
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
RootupleChibGen::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RootupleChibGen::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
RootupleChibGen::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RootupleChibGen::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RootupleChibGen::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RootupleChibGen::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RootupleChibGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RootupleChibGen);
