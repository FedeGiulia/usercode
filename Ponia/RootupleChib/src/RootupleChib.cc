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
// $Id: RootupleChib.cc,v 1.4 2013/06/21 12:49:49 gdujany Exp $
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


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
  edm::InputTag refit1S_Label;
  edm::InputTag refit2S_Label;
  edm::InputTag refit3S_Label;
  edm::InputTag primaryVertices_Label;
  edm::InputTag triggerResults_Label;
  bool isMC_;
  
  UInt_t   run;
  UInt_t   event;
  Double_t invm1S;
  Double_t invm2S;
  Double_t invm3S;
  Double_t chib_mass;
  Double_t chib_pt;
  Double_t chib_eta;
  Double_t chib_phi;
  Double_t dimuon_mass;
  Double_t dimuon_rapidity;
  Double_t dimuon_pt;
  Double_t photon_eta;
  Double_t photon_pt;
  Double_t ele_lowerPt_pt;
  Double_t ele_higherPt_pt;
  Double_t ctpv;
  Double_t ctpv_error;
  Double_t pi0_abs_mass;
  Double_t Y1S_nsigma;
  Double_t Y2S_nsigma;
  Double_t Y3S_nsigma;
  Double_t conv_vertex;
  Double_t dz;
  Double_t numPrimaryVertices;
  Int_t trigger;
  

  Double_t rf1S_chib_mass;
  Double_t probFit1S;
  Double_t rf1S_chib_pt;
  Double_t rf1S_chib_eta;
  Double_t rf1S_dimuon_mass;
  Double_t rf1S_dimuon_rapidity;
  Double_t rf1S_dimuon_pt;
  Double_t rf1S_photon_eta;
  Double_t rf1S_photon_pt;

  Double_t rf2S_chib_mass;
  Double_t probFit2S;
  Double_t rf3S_chib_mass;
  Double_t probFit3S;

  
  Double_t ups_mass;
  Double_t ups_rapidity;
  Double_t ups_pt;
  
  TTree* chib_tree;
  TTree* upsilon_tree;
  TLorentzVector lorVect;
  
  TLorentzVector chib_p4;
  Int_t chib_pdgId;
  TLorentzVector Upsilon_p4;
  Int_t Upsilon_pdgId;
  TLorentzVector photon_p4;
  TLorentzVector muP_p4;
  TLorentzVector muM_p4;

 
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

// 2011 par
//static const double Y_sig_par_A = 0.058;
//static const double Y_sig_par_B = 0.047;
//static const double Y_sig_par_C = 0.22;

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;


//
// constructors and destructor
//
RootupleChib::RootupleChib(const edm::ParameterSet& iConfig):
	chi_cand_Label(iConfig.getParameter<edm::InputTag>("chi_cand")),
	pi0_comb_Label(iConfig.getParameter<edm::InputTag>("pi0_comb")),
	ups_Label(iConfig.getParameter<edm::InputTag>("ups_cand")),
	refit1S_Label(iConfig.getParameter<edm::InputTag>("refit1S")),
	refit2S_Label(iConfig.getParameter<edm::InputTag>("refit2S")),
	refit3S_Label(iConfig.getParameter<edm::InputTag>("refit3S")),
	primaryVertices_Label(iConfig.getParameter<edm::InputTag>("primaryVertices")),
	triggerResults_Label(iConfig.getParameter<edm::InputTag>("TriggerResults")),
	isMC_(iConfig.getParameter<bool>("isMC"))
{
	edm::Service<TFileService> fs;
	chib_tree = fs->make<TTree>("chibTree","Tree of chib");

	
    chib_tree->Branch("run"  , &run,   "run/I");
    chib_tree->Branch("event", &event, "event/I");

    chib_tree->Branch("invm1S", &invm1S, "invm1S/D");
	chib_tree->Branch("invm2S", &invm2S, "invm2S/D");
	chib_tree->Branch("invm3S", &invm3S, "invm3S/D");
	chib_tree->Branch("chib_mass", &chib_mass, "chib_mass/D");
	chib_tree->Branch("chib_pt", &chib_pt, "chib_pt/D");
	chib_tree->Branch("chib_eta", &chib_eta, "chib_eta/D");
	chib_tree->Branch("chib_phi", &chib_phi, "chib_phi/D");
	chib_tree->Branch("dimuon_mass", &dimuon_mass, "dimuon_mass/D");
	chib_tree->Branch("dimuon_rapidity", &dimuon_rapidity, "dimuon_rapidity/D");
	chib_tree->Branch("dimuon_pt", &dimuon_pt, "dimuon_pt/D");
	chib_tree->Branch("dimuon_pt", &dimuon_pt, "dimuon_pt/D");
	chib_tree->Branch("photon_eta", &photon_eta, "photon_eta/D");
	chib_tree->Branch("photon_pt", &photon_pt, "photon_pt/D");
	chib_tree->Branch("ele_lowerPt_pt", &ele_lowerPt_pt, "ele_lowerPt_pt/D");
	chib_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");
	chib_tree->Branch("ctpv", &ctpv, "ctpv/D");
	chib_tree->Branch("ctpv_error", &ctpv_error, "ctpv_error/D");
	chib_tree->Branch("pi0_abs_mass", &pi0_abs_mass, "pi0_abs_mass/D");
	chib_tree->Branch("Y1S_nsigma", &Y1S_nsigma, "Y1S_nsigma/D");
	chib_tree->Branch("Y2S_nsigma", &Y2S_nsigma, "Y2S_nsigma/D");
	chib_tree->Branch("Y3S_nsigma", &Y3S_nsigma, "Y3S_nsigma/D");
	chib_tree->Branch("conv_vertex", &conv_vertex, "conv_vertex/D");
	chib_tree->Branch("dz", &dz, "dz/D");
	chib_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
	chib_tree->Branch("trigger", &trigger, "trigger/I");
	chib_tree->Branch("rf1S_chib_mass", &rf1S_chib_mass, "rf1S_chib_mass/D");
	chib_tree->Branch("probFit1S", &probFit1S, "probFit1S/D");

	chib_tree->Branch("rf1S_chib_pt", &rf1S_chib_pt, "rf1S_chib_pt/D");
	chib_tree->Branch("rf1S_chib_eta", &rf1S_chib_eta, "rf1S_chib_eta/D");
	chib_tree->Branch("rf1S_dimuon_mass", &rf1S_dimuon_mass, "rf1S_dimuon_mass/D");
	chib_tree->Branch("rf1S_dimuon_rapidity", &rf1S_dimuon_rapidity, "rf1S_dimuon_rapidity/D");
	chib_tree->Branch("rf1S_dimuon_pt", &rf1S_dimuon_pt, "rf1S_dimuon_pt/D");
	chib_tree->Branch("rf1S_dimuon_pt", &rf1S_dimuon_pt, "rf1S_dimuon_pt/D");
	chib_tree->Branch("rf1S_photon_eta", &rf1S_photon_eta, "rf1S_photon_eta/D");
	chib_tree->Branch("rf1S_photon_pt", &rf1S_photon_pt, "rf1S_photon_pt/D");
	chib_tree->Branch("rf2S_chib_mass", &rf2S_chib_mass, "rf2S_chib_mass/D");
	chib_tree->Branch("probFit2S", &probFit2S, "probFit2S/D");
    chib_tree->Branch("rf3S_chib_mass", &rf3S_chib_mass, "rf3S_chib_mass/D");
	chib_tree->Branch("probFit3S", &probFit3S, "probFit3S/D");

	if(isMC_)
	  {
	    chib_tree->Branch("chib_p4", "TLorentzVector", &chib_p4);
	    chib_tree->Branch("chib_pdgId", &chib_pdgId, "chib_pdgId/I");
	    chib_tree->Branch("Upsilon_p4", "TLorentzVector", &Upsilon_p4);
	    chib_tree->Branch("Upsilon_pdgId", &Upsilon_pdgId, "Upsilon_pdgId/I");
	    chib_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);
	    chib_tree->Branch("muP_p4", "TLorentzVector", &muP_p4);
	    chib_tree->Branch("muM_p4", "TLorentzVector", &muM_p4);
	  }

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
	using namespace std;
	Handle<std::vector<std::vector<float > > > pi0_comb_handle;
	iEvent.getByLabel(pi0_comb_Label, pi0_comb_handle);

	Handle<std::vector<pat::CompositeCandidate> > chi_cand_handle;
	iEvent.getByLabel(chi_cand_Label, chi_cand_handle);

	Handle<std::vector<std::vector<float > > > ups_hand;
	iEvent.getByLabel(ups_Label, ups_hand);

	Handle<std::vector<pat::CompositeCandidate> > refit1S_handle;
	iEvent.getByLabel(refit1S_Label, refit1S_handle);

	Handle<std::vector<pat::CompositeCandidate> > refit2S_handle;
	iEvent.getByLabel(refit2S_Label, refit2S_handle);

	Handle<std::vector<pat::CompositeCandidate> > refit3S_handle;
	iEvent.getByLabel(refit3S_Label, refit3S_handle);

	Handle<std::vector<reco::Vertex > > primaryVertices_handle;
	iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

	edm::Handle< edm::TriggerResults > triggerResults_handle;
	iEvent.getByLabel( triggerResults_Label , triggerResults_handle);
	//https://cmssdt.cern.ch/SDT/lxr/source/DQMOffline/JetMET/src/CaloTowerAnalyzer.cc#190

	numPrimaryVertices = primaryVertices_handle->size();

	std::vector<double> pi0_abs_values;
	pat::CompositeCandidate chi_cand;
	pat::CompositeCandidate refit1S;
	pat::CompositeCandidate refit2S;
	pat::CompositeCandidate refit3S;
	std::vector<float > pi0_comb;

	double QValue = 0.0;

	if(isMC_)
	  {
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
			continue;
		      }
		    
		    if(pdgId ==    553 || // Upsilon 1S
		       pdgId == 100553 || // Upsilon 2S
		       pdgId == 200553 )  // Upsilon 3S
		      {
			Upsilon_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass()); 
			Upsilon_pdgId = pdgId; 
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
		      }		    
		  }
	      }    
	  } // end if isMC

	//grab Trigger informations
	// save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
	// (pass 11)(pass 8)(pass 7)(pass 5)
	// es. 11 = pass 5, 7 and 11
	// es. 4 = pass only 8
	trigger = 0;
	if (triggerResults_handle.isValid())
	  {
	    const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	    
	    vector <unsigned int> bits_5;
	    for(int version = 3; version<7; version ++)
	      {
		//char str[30];
		//sprintf(str, "HLT_Dimuon5_Upsilon_v%i", version);
		//puts(str);
		stringstream ss;
		ss<<"HLT_Dimuon5_Upsilon_v"<<version;
		bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
	      }

	    vector <unsigned int> bits_7;
	    for(int version = 3; version<8; version ++)
	      {
		stringstream ss;
		ss<<"HLT_Dimuon7_Upsilon_v"<<version;
		bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
	      }

	    vector <unsigned int> bits_8;
	    for(int version = 3; version<7; version ++)
	      {
		stringstream ss;
		ss<<"HLT_Dimuon8_Upsilon_v"<<version;
		bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
	      }

	    vector <unsigned int> bits_11;
	    for(int version = 3; version<7; version ++)
	      {
		stringstream ss;
		ss<<"HLT_Dimuon11_Upsilon_v"<<version;
		bits_11.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
	      }
	    

	    for(unsigned int i=0; i<bits_5.size(); i++)
	      {
	    	unsigned int bit = bits_5[i];
	    	if( bit < triggerResults_handle->size() )
	    	  {
	    	    if( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) )
	    	      {
	    		//std::cout<<std::endl<<"Passed trigger 5"<<std::endl;
			trigger += 1;
	    	      }
	    	  }
	      }
	    for(unsigned int i=0; i<bits_7.size(); i++)
	      {
	    	unsigned int bit = bits_7[i];
	    	if( bit < triggerResults_handle->size() )
	    	  {
	    	    if( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) )
	    	      {
	    		//std::cout<<std::endl<<"Passed trigger 7"<<std::endl;
			trigger += 2;
	    	      }
	    	  }
	      }
	    for(unsigned int i=0; i<bits_8.size(); i++)
	      {
	    	unsigned int bit = bits_8[i];
	    	if( bit < triggerResults_handle->size() )
	    	  {
	    	    if( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) )
	    	      {
	    		//std::cout<<std::endl<<"Passed trigger 8"<<std::endl;
			trigger += 4;
	    	      }
	    	  }
	      }
	    for(unsigned int i=0; i<bits_11.size(); i++)
	      {
	    	unsigned int bit = bits_11[i];
	    	if( bit < triggerResults_handle->size() )
	    	  {
	    	    if( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) )
	    	      {
	    		//std::cout<<std::endl<<"Passed trigger 11"<<std::endl;
			trigger += 8;
	    	      }
	    	  }
	      }

	  }
	
	// grabbing chib inforamtion
	if(chi_cand_handle.isValid() &&  refit1S_handle.isValid() && refit2S_handle.isValid() && refit3S_handle.isValid())
	  {  
	    for(unsigned int i=0; i< chi_cand_handle->size(); i++)
	      {
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
		chib_mass = chi_cand.mass();
		chib_pt = chi_cand.pt();
		chib_eta = chi_cand.eta();
		chib_phi = chi_cand.phi();
		dimuon_mass = chi_cand.daughter("dimuon")->mass();
		dimuon_rapidity = chi_cand.daughter("dimuon")->y();
		dimuon_pt = chi_cand.daughter("dimuon")->pt();
		photon_eta = chi_cand.daughter("photon")->eta();
		photon_pt = chi_cand.daughter("photon")->pt();
		
		Double_t ele1_pt = chi_cand.daughter("photon")->daughter("ele1")->pt();
		Double_t ele2_pt = chi_cand.daughter("photon")->daughter("ele2")->pt();
		
		if(ele1_pt>ele2_pt)
		  {
		    ele_higherPt_pt = ele1_pt;
		    ele_lowerPt_pt = ele2_pt;
		  }
		else
		  {
		    ele_higherPt_pt = ele2_pt;
		    ele_lowerPt_pt = ele1_pt;
		  }
		
		ctpv = ( dynamic_cast<pat::CompositeCandidate*>(chi_cand.daughter("dimuon")) )->userFloat("ppdlPV");
		ctpv_error = ( dynamic_cast<pat::CompositeCandidate*>(chi_cand.daughter("dimuon")) )->userFloat("ppdlErrPV");
		pi0_abs_mass = pi0_abs_values[0];
		conv_vertex = chi_cand.daughter("photon")->vertex().rho();
		dz = chi_cand.userFloat("dz");
		
        // old (2011) parameterization
		//double sigma = Y_sig_par_A+Y_sig_par_B*(fabs(dimuon_rapidity)-Y_sig_par_C);
		//if (sigma < Y_sig_par_A) sigma = Y_sig_par_A;

        // 2012 parameterization
        double sigma = Y_sig_par_A + 
                       Y_sig_par_B * pow(fabs(dimuon_rapidity), 2) + 
                       Y_sig_par_C * pow(fabs(dimuon_rapidity), 3) ;

		Y1S_nsigma = fabs(dimuon_mass - Y1SMass)/sigma;
		Y2S_nsigma = fabs(dimuon_mass - Y2SMass)/sigma;
		Y3S_nsigma = fabs(dimuon_mass - Y3SMass)/sigma;
		
		if(i<refit1S_handle->size())
		  {
		    refit1S = refit1S_handle->at(i);
		    rf1S_chib_mass = refit1S.mass(); 
		    probFit1S = refit1S.userFloat("vProb");
		    rf1S_chib_pt = refit1S.pt(); 
		    rf1S_chib_eta = refit1S.eta(); 
		    rf1S_dimuon_mass = refit1S.daughter("dimuon")->mass();
		    rf1S_dimuon_rapidity = refit1S.daughter("dimuon")->y();
		    rf1S_dimuon_pt = refit1S.daughter("dimuon")->pt();
		    rf1S_photon_eta = refit1S.daughter("photon")->eta();
		    rf1S_photon_pt = refit1S.daughter("photon")->pt();
		  }
		else 
		  {
		    rf1S_chib_mass = invm1S; 
		    probFit1S = 0;
		  }
		
		if(i<refit2S_handle->size())
		  {
		    refit2S = refit2S_handle->at(i);
		    rf2S_chib_mass = refit2S.mass(); 
		    probFit2S = refit2S.userFloat("vProb");
		  }
		else 
		  {
		    rf2S_chib_mass = invm2S; 
		    probFit2S = 0;
		  }
		
		if(i<refit3S_handle->size())
		  {
		    refit3S = refit3S_handle->at(i);
		    rf3S_chib_mass = refit3S.mass(); 
		    probFit3S = refit3S.userFloat("vProb");
		  }
		else 
		  {
		    rf3S_chib_mass = invm3S; 
		    probFit3S = 0;
		  }
		
		run=   iEvent.id().run();
        event= iEvent.id().event();
        //std::cout << "Run : " << run <<  " Event : " << event<< std::endl;
		chib_tree->Fill();
	      } // for i 
	  } // if chi_cand_handle && refit1S &&
	
	
	
	if(ups_hand.isValid() )
	  {
	    for(unsigned int i=0; i< ups_hand->size(); i++)
	      {
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
