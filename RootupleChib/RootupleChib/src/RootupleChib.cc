// -*- C++ -*-
//
// Package:    RootupleChib
// Class:      RootupleChib
// 
/**\class RootupleChib RootupleChib.cc RootupleChib/RootupleChib/src/RootupleChib.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alessandro Degano,32 1-C13,+41227678098
//         Created:  Tue Feb 19 14:32:37 CET 2013
// $Id$
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

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>

//
// class declaration
//

class RootupleChib : public edm::EDAnalyzer {
   public:
      explicit RootupleChib(const edm::ParameterSet&);
      ~RootupleChib();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
	std::string file_name;
	edm::InputTag chi_cand_Label;
	edm::InputTag pi0_comb_Label;
	edm::InputTag ups_Label;

	Double_t invm1S;
	Double_t invm2S;
	Double_t invm3S;
	Double_t dimuon_mass;
	Double_t dimuon_rapidity;
	Double_t dimuon_pt;
	Double_t photon_eta;
	Double_t photon_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t pi0_abs_mass;
	Double_t Y1S_nsigma;
	Double_t Y2S_nsigma;
	Double_t Y3S_nsigma;

	Double_t ups_mass;
	Double_t ups_rapidity;
	Double_t ups_pt;

	TTree* chib_tree;
	TTree* upsilon_tree;
	TLorentzVector lorVect;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass = 0.1349766;
static const Double_t Y1SMass =  9.46030;
static const Double_t Y2SMass = 10.02326;
static const Double_t Y3SMass = 10.3552;

static const double Y_sig_par_A = 0.058;
static const double Y_sig_par_B = 0.047;
static const double Y_sig_par_C = 0.22;


//
// constructors and destructor
//
RootupleChib::RootupleChib(const edm::ParameterSet& iConfig):
	chi_cand_Label(iConfig.getParameter<edm::InputTag>("chi_cand")),
	pi0_comb_Label(iConfig.getParameter<edm::InputTag>("pi0_comb")),
	ups_Label(iConfig.getParameter<edm::InputTag>("ups_cand"))
{
	edm::Service<TFileService> fs;
	chib_tree = fs->make<TTree>("chibTree","Tree of chib");

	chib_tree->Branch("invm1S", &invm1S, "invm1S/D");
	chib_tree->Branch("invm2S", &invm2S, "invm2S/D");
	chib_tree->Branch("invm3S", &invm3S, "invm3S/D");
	chib_tree->Branch("dimuon_mass", &dimuon_mass, "dimuon_mass/D");
	chib_tree->Branch("dimuon_rapidity", &dimuon_rapidity, "dimuon_rapidity/D");
	chib_tree->Branch("dimuon_pt", &dimuon_pt, "dimuon_pt/D");
	chib_tree->Branch("dimuon_pt", &dimuon_pt, "dimuon_pt/D");
	chib_tree->Branch("photon_eta", &photon_eta, "photon_eta/D");
	chib_tree->Branch("photon_pt", &photon_pt, "photon_pt/D");
	chib_tree->Branch("ctpv", &ctpv, "ctpv/D");
	chib_tree->Branch("ctpv_error", &ctpv_error, "ctpv_error/D");
	chib_tree->Branch("pi0_abs_mass", &pi0_abs_mass, "pi0_abs_mass/D");
	chib_tree->Branch("Y1S_nsigma", &Y1S_nsigma, "Y1S_nsigma/D");
	chib_tree->Branch("Y2S_nsigma", &Y2S_nsigma, "Y2S_nsigma/D");
	chib_tree->Branch("Y3S_nsigma", &Y3S_nsigma, "Y3S_nsigma/D");

	upsilon_tree = fs->make<TTree>("upsTree","Tree of Upsilon");
	upsilon_tree->Branch("ups_mass", &ups_mass, "ups_mass/D");
	upsilon_tree->Branch("ups_rapidity", &ups_rapidity, "ups_rapidity/D");
	upsilon_tree->Branch("ups_pt", &ups_pt, "ups_pt/D");
}


RootupleChib::~RootupleChib()
{

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RootupleChib::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   	using namespace edm;
	Handle<std::vector<std::vector<float > > > pi0_comb_handle;
	iEvent.getByLabel(pi0_comb_Label, pi0_comb_handle);

	Handle<std::vector<pat::CompositeCandidate> > chi_cand_handle;
	iEvent.getByLabel(chi_cand_Label, chi_cand_handle);

	Handle<std::vector<std::vector<float > > > ups_hand;
	iEvent.getByLabel(ups_Label, ups_hand);

	std::vector<double> pi0_abs_values;
	pat::CompositeCandidate chi_cand;
	std::vector<float > pi0_comb;

	double QValue = 0.0;

	if(chi_cand_handle.isValid() ){  
		for(unsigned int i=0; i< chi_cand_handle->size(); i++){
			chi_cand = chi_cand_handle->at(i);
			
			if( pi0_comb_handle.isValid() && pi0_comb_handle->at(i).size() > 0 ){
				pi0_comb = pi0_comb_handle->at(i);
				for(unsigned int k=0; k< pi0_comb.size(); k++){
					pi0_abs_values.push_back( fabs( pi0_comb[k] - pi0_mass) );
				}
				std::sort( pi0_abs_values.begin(), pi0_abs_values.end() );
			} else {
				pi0_abs_values.push_back( 1.0 );
			}
			
			QValue = chi_cand.mass() - chi_cand.daughter("dimuon")->mass();
		
			invm1S = QValue + Y1SMass;
			invm2S = QValue + Y2SMass;
			invm3S = QValue + Y3SMass;
			dimuon_mass = chi_cand.daughter("dimuon")->mass();
			dimuon_rapidity = chi_cand.daughter("dimuon")->y();
			dimuon_pt = chi_cand.daughter("dimuon")->pt();
			photon_eta = chi_cand.daughter("photon")->eta();
			photon_pt = chi_cand.daughter("photon")->pt();
			ctpv = ( dynamic_cast<pat::CompositeCandidate*>(chi_cand.daughter("dimuon")) )->userFloat("ppdlPV");
			ctpv_error = ( dynamic_cast<pat::CompositeCandidate*>(chi_cand.daughter("dimuon")) )->userFloat("ppdlErrPV");
			pi0_abs_mass = pi0_abs_values[0];

			double sigma = Y_sig_par_A+Y_sig_par_B*(fabs(dimuon_rapidity)-Y_sig_par_C);
			if (sigma < Y_sig_par_A) sigma = Y_sig_par_A;
			Y1S_nsigma = fabs(dimuon_mass - Y1SMass)/sigma;
			Y2S_nsigma = fabs(dimuon_mass - Y2SMass)/sigma;
			Y3S_nsigma = fabs(dimuon_mass - Y3SMass)/sigma;
			chib_tree->Fill();
		}
	}
	
	if(ups_hand.isValid() ){
		for(unsigned int i=0; i< ups_hand->size(); i++){
			lorVect.SetPxPyPzE(ups_hand->at(i)[0], ups_hand->at(i)[1], ups_hand->at(i)[2], ups_hand->at(i)[3]);
			ups_mass = lorVect.M();
			ups_rapidity = lorVect.Rapidity();
			ups_pt = lorVect.Pt();
			upsilon_tree->Fill();
		}
	}
}


// ------------ method called once each job just before starting event loop  ------------
void 
RootupleChib::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RootupleChib::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
RootupleChib::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RootupleChib::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RootupleChib::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RootupleChib::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RootupleChib::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RootupleChib);
