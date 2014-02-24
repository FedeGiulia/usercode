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
//         Created:  Tue Feb 19 14:32:37 CET 2013 //
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

class RootupleChic:public edm::EDAnalyzer {
      public:
	explicit RootupleChic(const edm::ParameterSet &);
	~RootupleChic();

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();

	virtual void beginRun(edm::Run const &, edm::EventSetup const &);
	virtual void endRun(edm::Run const &, edm::EventSetup const &);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

	// ----------member data ---------------------------
	std::string file_name;
	edm::InputTag chi_cand_Label;
	edm::InputTag pi0_comb_Label;
	edm::InputTag ups_Label;
	edm::InputTag refit1S_Label;
	edm::InputTag refit2S_Label;
	edm::InputTag primaryVertices_Label;
	edm::InputTag triggerResults_Label;
	bool isMC_;

	UInt_t run;
	UInt_t event;

	TLorentzVector chi_p4;
	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
	TLorentzVector photon_p4;

	TLorentzVector rf1S_chi_p4;
	TLorentzVector rf1S_dimuon_p4;
	TLorentzVector rf1S_muonP_p4;
	TLorentzVector rf1S_muonN_p4;
	TLorentzVector rf1S_photon_p4;

	TLorentzVector rf2S_chi_p4;
	TLorentzVector rf2S_dimuon_p4;
	TLorentzVector rf2S_muonP_p4;
	TLorentzVector rf2S_muonN_p4;
	TLorentzVector rf2S_photon_p4;

	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t pi0_abs_mass;
	Double_t psi1S_nsigma;
	Double_t psi2S_nsigma;
	Double_t conv_vertex;
	Double_t dz;
	Double_t numPrimaryVertices;
	Int_t trigger;

	Double_t probFit1S;
	Int_t rf1S_rank;
	Double_t probFit2S;

	TTree *chib_tree;

	TLorentzVector gen_chi_p4;
	Int_t gen_chi_pdgId;
	TLorentzVector gen_dimuon_p4;
	Int_t gen_dimuon_pdgId;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muP_p4;
	TLorentzVector gen_muM_p4;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass = 0.1349766;
static const Double_t psi1SMass = 3.09691;
static const Double_t psi2SMass = 3.68610;

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
RootupleChic::RootupleChic(const edm::ParameterSet & iConfig):chi_cand_Label(iConfig.getParameter < edm::InputTag > ("chi_cand")),
pi0_comb_Label(iConfig.getParameter < edm::InputTag > ("pi0_comb")),
ups_Label(iConfig.getParameter < edm::InputTag > ("ups_cand")),
refit1S_Label(iConfig.getParameter < edm::InputTag > ("refit1S")),
refit2S_Label(iConfig.getParameter < edm::InputTag > ("refit2S")),
primaryVertices_Label(iConfig.getParameter < edm::InputTag > ("primaryVertices")),
triggerResults_Label(iConfig.getParameter < edm::InputTag > ("TriggerResults")), isMC_(iConfig.getParameter < bool > ("isMC"))
{
	edm::Service < TFileService > fs;
	chib_tree = fs->make < TTree > ("chicTree", "Tree of chic");

	chib_tree->Branch("run", &run, "run/I");
	chib_tree->Branch("event", &event, "event/I");

	chib_tree->Branch("chi_p4", "TLorentzVector", &chi_p4);
	chib_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
	chib_tree->Branch("muonP_p4", "TLorentzVector", &muonP_p4);
	chib_tree->Branch("muonN_p4", "TLorentzVector", &muonN_p4);
	chib_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);

	chib_tree->Branch("rf1S_chi_p4", "TLorentzVector", &rf1S_chi_p4);
	chib_tree->Branch("rf1S_dimuon_p4", "TLorentzVector", &rf1S_dimuon_p4);
	chib_tree->Branch("rf1S_muonP_p4", "TLorentzVector", &rf1S_muonP_p4);
	chib_tree->Branch("rf1S_muonN_p4", "TLorentzVector", &rf1S_muonN_p4);
	chib_tree->Branch("rf1S_photon_p4", "TLorentzVector", &rf1S_photon_p4);

	chib_tree->Branch("rf2S_chi_p4", "TLorentzVector", &rf2S_chi_p4);
	chib_tree->Branch("rf2S_dimuon_p4", "TLorentzVector", &rf2S_dimuon_p4);
	chib_tree->Branch("rf2S_muonP_p4", "TLorentzVector", &rf2S_muonP_p4);
	chib_tree->Branch("rf2S_muonN_p4", "TLorentzVector", &rf2S_muonN_p4);
	chib_tree->Branch("rf2S_photon_p4", "TLorentzVector", &rf2S_photon_p4);

	chib_tree->Branch("ele_lowerPt_pt", &ele_lowerPt_pt, "ele_lowerPt_pt/D");
	chib_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");
	chib_tree->Branch("ctpv", &ctpv, "ctpv/D");
	chib_tree->Branch("ctpv_error", &ctpv_error, "ctpv_error/D");
	chib_tree->Branch("pi0_abs_mass", &pi0_abs_mass, "pi0_abs_mass/D");
	chib_tree->Branch("psi1S_nsigma", &psi1S_nsigma, "psi1S_nsigma/D");
	chib_tree->Branch("psi2S_nsigma", &psi2S_nsigma, "psi2S_nsigma/D");

	chib_tree->Branch("conv_vertex", &conv_vertex, "conv_vertex/D");
	chib_tree->Branch("dz", &dz, "dz/D");
	chib_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
	chib_tree->Branch("trigger", &trigger, "trigger/I");
	chib_tree->Branch("probFit1S", &probFit1S, "probFit1S/D");
	chib_tree->Branch("rf1S_rank", &rf1S_rank, "rf1S_rank/I");
	chib_tree->Branch("probFit2S", &probFit2S, "probFit2S/D");

	if (isMC_) {
		chib_tree->Branch("gen_chi_p4", "TLorentzVector", &gen_chi_p4);
		chib_tree->Branch("gen_chi_pdgId", &gen_chi_pdgId, "gen_chi_pdgId/I");
		chib_tree->Branch("gen_dimuon_p4", "TLorentzVector", &gen_dimuon_p4);
		chib_tree->Branch("gen_dimuon_pdgId", &gen_dimuon_pdgId, "gen_dimuon_pdgId/I");
		chib_tree->Branch("gen_photon_p4", "TLorentzVector", &gen_photon_p4);
		chib_tree->Branch("gen_muP_p4", "TLorentzVector", &gen_muP_p4);
		chib_tree->Branch("gen_muM_p4", "TLorentzVector", &gen_muM_p4);
	}

}

RootupleChic::~RootupleChic()
{

}

//
// member functions
//

// ------------ method called for each event  ------------
void
 RootupleChic::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
	using namespace edm;
	using namespace std;
	Handle < std::vector < std::vector < float >>>pi0_comb_handle;
	iEvent.getByLabel(pi0_comb_Label, pi0_comb_handle);

	Handle < std::vector < pat::CompositeCandidate > >chi_cand_handle;
	iEvent.getByLabel(chi_cand_Label, chi_cand_handle);

	Handle < std::vector < std::vector < float >>>ups_hand;
	iEvent.getByLabel(ups_Label, ups_hand);

	Handle < std::vector < pat::CompositeCandidate > >refit1S_handle;
	iEvent.getByLabel(refit1S_Label, refit1S_handle);

	Handle < std::vector < pat::CompositeCandidate > >refit2S_handle;
	iEvent.getByLabel(refit2S_Label, refit2S_handle);

	Handle < std::vector < reco::Vertex > >primaryVertices_handle;
	iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

	edm::Handle < edm::TriggerResults > triggerResults_handle;
	iEvent.getByLabel(triggerResults_Label, triggerResults_handle);
	//https://cmssdt.cern.ch/SDT/lxr/source/DQMOffline/JetMET/src/CaloTowerAnalyzer.cc#190

	numPrimaryVertices = primaryVertices_handle->size();

	std::vector < double >pi0_abs_values;
	pat::CompositeCandidate chi_cand;
	pat::CompositeCandidate refit1S;
	pat::CompositeCandidate refit2S;
	std::vector < float >pi0_comb;

	if (isMC_) {
		Handle < reco::GenParticleCollection > GenParticles;
		iEvent.getByLabel("genParticles", GenParticles);	//Open GenParticles

		if (GenParticles.isValid()) {
			for (reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle) {
				Int_t pdgId = itParticle->pdgId();

				if (pdgId == 10441 ||	//chic0
				    pdgId == 20443 ||	//chic1
				    pdgId == 445) {	//chic2 
					gen_chi_p4.SetPtEtaPhiM(itParticle->pt(), itParticle->eta(), itParticle->phi(), itParticle->mass());
					gen_chi_pdgId = pdgId;
					continue;
				}

				if (pdgId == 443 ||	// psi 1S
				    pdgId == 100553) {	// psi 2S  

					gen_dimuon_p4.SetPtEtaPhiM(itParticle->pt(), itParticle->eta(), itParticle->phi(), itParticle->mass());
					gen_dimuon_pdgId = pdgId;
					continue;
				}

				if (pdgId == 13) {	// mu-

					gen_muM_p4.SetPtEtaPhiM(itParticle->pt(), itParticle->eta(), itParticle->phi(), itParticle->mass());
					continue;
				}

				if (pdgId == -13) {	// mu+

					gen_muP_p4.SetPtEtaPhiM(itParticle->pt(), itParticle->eta(), itParticle->phi(), itParticle->mass());
					continue;
				}

				if (pdgId == 22) {	// photon

					gen_photon_p4.SetPtEtaPhiM(itParticle->pt(), itParticle->eta(), itParticle->phi(), itParticle->mass());
				}
			}
		}
	}			// end if isMC

	//grab Trigger informations
	// save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
	// (pass 11)(pass 8)(pass 7)(pass 5)
	// es. 11 = pass 5, 7 and 11
	// es. 4 = pass only 8
	trigger = 0;
	if (triggerResults_handle.isValid()) {
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);

		vector < unsigned int >bits_5;
		for (int version = 3; version < 8; version++) {
			//char str[30];
			//sprintf(str, "HLT_Dimuon5_Upsilon_v%i", version);
			//puts(str);
			stringstream ss;
			ss << "HLT_Dimuon8_Jpsi_v" << version;
			bits_5.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str()));
		}

		vector < unsigned int >bits_7;
		for (int version = 3; version < 8; version++) {
			stringstream ss;
			ss << "HLT_Dimuon10_Jpsi_v" << version;
			bits_7.push_back(TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str()));
		}

		for (unsigned int i = 0; i < bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if (bit < triggerResults_handle->size()) {
				if (triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					//std::cout<<std::endl<<"Passed trigger 5"<<std::endl;
					trigger += 1;
				}
			}
		}
		for (unsigned int i = 0; i < bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if (bit < triggerResults_handle->size()) {
				if (triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					//std::cout<<std::endl<<"Passed trigger 7"<<std::endl;
					trigger += 2;
				}
			}
		}

	}

	bool bestCandidateOnly_ = true;
	// grabbing chi inforamtion
	if (chi_cand_handle.isValid() && chi_cand_handle->size() > 0) {

		unsigned int csize = chi_cand_handle->size();
		if (bestCandidateOnly_)
			csize = 1;

		for (unsigned int i = 0; i < csize; i++) {
			chi_cand = chi_cand_handle->at(i);

			if (pi0_comb_handle.isValid() && pi0_comb_handle->at(i).size() > 0) {
				pi0_comb = pi0_comb_handle->at(i);
				for (unsigned int k = 0; k < pi0_comb.size(); k++) {
					pi0_abs_values.push_back(fabs(pi0_comb[k] - pi0_mass));
				}
				std::sort(pi0_abs_values.begin(), pi0_abs_values.end());
			} else {
				pi0_abs_values.push_back(1.0);
			}

			chi_p4.SetPtEtaPhiM(chi_cand.pt(), chi_cand.eta(), chi_cand.phi(), chi_cand.mass());

			dimuon_p4.SetPtEtaPhiM(chi_cand.daughter("dimuon")->pt(), chi_cand.daughter("dimuon")->eta(), chi_cand.daughter("dimuon")->phi(), chi_cand.daughter("dimuon")->mass());

			photon_p4.SetPtEtaPhiM(chi_cand.daughter("photon")->pt(), chi_cand.daughter("photon")->eta(), chi_cand.daughter("photon")->phi(), chi_cand.daughter("photon")->mass());

			reco::Candidate::LorentzVector vP = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			reco::Candidate::LorentzVector vM = chi_cand.daughter("dimuon")->daughter("muon2")->p4();

			if (chi_cand.daughter("dimuon")->daughter("muon1")->charge() < 0) {
				vP = chi_cand.daughter("dimuon")->daughter("muon2")->p4();
				vM = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			}

			muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
			muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

			Double_t ele1_pt = chi_cand.daughter("photon")->daughter("ele1")->pt();
			Double_t ele2_pt = chi_cand.daughter("photon")->daughter("ele2")->pt();

			if (ele1_pt > ele2_pt) {
				ele_higherPt_pt = ele1_pt;
				ele_lowerPt_pt = ele2_pt;
			} else {
				ele_higherPt_pt = ele2_pt;
				ele_lowerPt_pt = ele1_pt;
			}

			ctpv = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlPV");
			ctpv_error = (dynamic_cast < pat::CompositeCandidate * >(chi_cand.daughter("dimuon")))->userFloat("ppdlErrPV");
			pi0_abs_mass = pi0_abs_values[0];
			conv_vertex = chi_cand.daughter("photon")->vertex().rho();
			dz = chi_cand.userFloat("dz");

			// old (2011) parameterization
			//double sigma = Y_sig_par_A+Y_sig_par_B*(fabs(dimuon_rapidity)-Y_sig_par_C);
			//if (sigma < Y_sig_par_A) sigma = Y_sig_par_A;

			// 2012 parameterization
			double sigma = Y_sig_par_A + Y_sig_par_B * pow(fabs(dimuon_p4.Rapidity()), 2) + Y_sig_par_C * pow(fabs(dimuon_p4.Rapidity()), 3);

			psi1S_nsigma = fabs(dimuon_p4.M() - psi1SMass) / sigma;
			psi2S_nsigma = fabs(dimuon_p4.M() - psi2SMass) / sigma;

			if (refit1S_handle.isValid() && i < refit1S_handle->size()) {

				refit1S = refit1S_handle->at(i);

				rf1S_chi_p4.SetPtEtaPhiM(refit1S.pt(), refit1S.eta(), refit1S.phi(), refit1S.mass());

				rf1S_dimuon_p4.SetPtEtaPhiM(refit1S.daughter("dimuon")->pt(), refit1S.daughter("dimuon")->eta(), refit1S.daughter("dimuon")->phi(), refit1S.daughter("dimuon")->mass());

				rf1S_photon_p4.SetPtEtaPhiM(refit1S.daughter("photon")->pt(), refit1S.daughter("photon")->eta(), refit1S.daughter("photon")->phi(), refit1S.daughter("photon")->mass());

				reco::Candidate::LorentzVector vP = refit1S.daughter("dimuon")->daughter("muon1")->p4();
				reco::Candidate::LorentzVector vM = refit1S.daughter("dimuon")->daughter("muon2")->p4();

				if (refit1S.daughter("dimuon")->daughter("muon1")->charge() < 0) {
					vP = refit1S.daughter("dimuon")->daughter("muon2")->p4();
					vM = refit1S.daughter("dimuon")->daughter("muon1")->p4();
				}

				rf1S_muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
				rf1S_muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
				probFit1S = refit1S.userFloat("vProb");

				rf1S_rank = i;
			} else {
				rf1S_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf1S_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf1S_photon_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf1S_muonP_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf1S_muonN_p4.SetPtEtaPhiM(0, 0, 0, 0);
				probFit1S = 0;
				rf1S_rank = -1;
			}	// if rf1S is valid

			if (refit2S_handle.isValid() && i < refit1S_handle->size()) {
				refit2S = refit2S_handle->at(i);

				rf2S_chi_p4.SetPtEtaPhiM(refit2S.pt(), refit2S.eta(), refit2S.phi(), refit2S.mass());

				rf2S_dimuon_p4.SetPtEtaPhiM(refit2S.daughter("dimuon")->pt(), refit2S.daughter("dimuon")->eta(), refit2S.daughter("dimuon")->phi(), refit2S.daughter("dimuon")->mass());

				rf2S_photon_p4.SetPtEtaPhiM(refit2S.daughter("photon")->pt(), refit2S.daughter("photon")->eta(), refit2S.daughter("photon")->phi(), refit2S.daughter("photon")->mass());

				reco::Candidate::LorentzVector vP = refit2S.daughter("dimuon")->daughter("muon1")->p4();
				reco::Candidate::LorentzVector vM = refit2S.daughter("dimuon")->daughter("muon2")->p4();

				if (refit2S.daughter("dimuon")->daughter("muon1")->charge() < 0) {
					vP = refit2S.daughter("dimuon")->daughter("muon2")->p4();
					vM = refit2S.daughter("dimuon")->daughter("muon1")->p4();
				}

				rf2S_muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
				rf2S_muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

				probFit2S = refit2S.userFloat("vProb");

			} else {
				rf2S_chi_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf2S_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf2S_photon_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf2S_muonP_p4.SetPtEtaPhiM(0, 0, 0, 0);
				rf2S_muonN_p4.SetPtEtaPhiM(0, 0, 0, 0);
				probFit2S = 0;
			}	// if rf2S valid

			run = iEvent.id().run();
			event = iEvent.id().event();
			//std::cout << "Run : " << run <<  " Event : " << event<< " M() "<< rf1S_chi_p4.M()<<  " rank "<< rf1S_rank<< " p " << probFit1S << std::endl;
			chib_tree->Fill();
		}		// for i on chi_cand_handle

	}			// if chi handle valid
}

// ------------ method called once each job just before starting event loop  ------------
void RootupleChic::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void RootupleChic::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void RootupleChic::beginRun(edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a run  ------------
void RootupleChic::endRun(edm::Run const &, edm::EventSetup const &)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void RootupleChic::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void RootupleChic::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RootupleChic::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RootupleChic);
