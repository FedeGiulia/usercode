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

class RootupleChicGen : public edm::EDAnalyzer
{
public:
    explicit RootupleChicGen(const edm::ParameterSet&);
    ~RootupleChicGen();

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
    TLorentzVector chic_p4;
    Double_t chic_mass;
    Double_t chic_pt;
    Double_t chic_eta;
    Double_t chic_rapidity;
    Double_t chic_phi;
    Int_t chic_pdgId;
    TLorentzVector Jpsi_p4;
    Double_t Jpsi_mass;
    Double_t Jpsi_pt;
    Double_t Jpsi_eta;
    Double_t Jpsi_rapidity;
    Double_t Jpsi_phi;
    Int_t Jpsi_pdgId;
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
RootupleChicGen::RootupleChicGen(const edm::ParameterSet& iConfig)
{
    edm::Service<TFileService> fs;
    genParticles_tree = fs->make<TTree>("GenParticlesTree","Tree of genParticles");

    genParticles_tree->Branch("chic_p4", "TLorentzVector", &chic_p4);
    genParticles_tree->Branch("chic_mass", &chic_mass, "chic_mass/D");
    genParticles_tree->Branch("chic_pt", &chic_pt, "chic_pt/D");
    genParticles_tree->Branch("chic_eta", &chic_eta, "chic_eta/D");
    genParticles_tree->Branch("chic_rapidity", &chic_rapidity, "chic_rapidity/D");
    genParticles_tree->Branch("chic_phi", &chic_phi, "chic_phi/D");
    genParticles_tree->Branch("chic_pdgId", &chic_pdgId, "chic_pdgId/I");
    genParticles_tree->Branch("Jpsi_p4", "TLorentzVector", &Jpsi_p4);
    genParticles_tree->Branch("Jpsi_mass", &Jpsi_mass, "Jpsi_mass/D");
    genParticles_tree->Branch("Jpsi_pt", &Jpsi_pt, "Jpsi_pt/D");
    genParticles_tree->Branch("Jpsi_eta", &Jpsi_eta, "Jpsi_eta/D");
    genParticles_tree->Branch("Jpsi_rapidity", &Jpsi_rapidity, "Jpsi_rapidity/D");
    genParticles_tree->Branch("Jpsi_phi", &Jpsi_phi, "Jpsi_phi/D");
    genParticles_tree->Branch("Jpsi_pdgId", &Jpsi_pdgId, "Jpsi_pdgId/I");
    genParticles_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);
    genParticles_tree->Branch("photon_pt", &photon_pt, "photon_pt/D");
    genParticles_tree->Branch("photon_eta", &photon_eta, "photon_eta/D");
    genParticles_tree->Branch("photon_phi", &photon_phi, "photon_phi/D");
    genParticles_tree->Branch("muP_p4", "TLorentzVector", &muP_p4);
    genParticles_tree->Branch("muM_p4", "TLorentzVector", &muM_p4);

}


RootupleChicGen::~RootupleChicGen()
{

}

void
RootupleChicGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    Handle<reco::GenParticleCollection> GenParticles;
    iEvent.getByLabel("genParticles",GenParticles);   //Open GenParticles

    if(GenParticles.isValid() ) {
        for(reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle) {
            Int_t pdgId = itParticle->pdgId();

            if(pdgId == 20443 || pdgId == 445 ) { //chic1/chic2
                chic_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass());
                chic_pdgId = pdgId;
                chic_mass = chic_p4.M();
                chic_pt = chic_p4.Pt();
                chic_eta = chic_p4.Eta();
                chic_rapidity = chic_p4.Rapidity();
                chic_phi = chic_p4.Phi();
                continue;
            }

            if( pdgId == 443 ) { // Jpsi
                Jpsi_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass());
                Jpsi_pdgId = pdgId;
                Jpsi_mass = Jpsi_p4.M();
                Jpsi_pt = Jpsi_p4.Pt();
                Jpsi_eta = Jpsi_p4.Eta();
                Jpsi_rapidity = Jpsi_p4.Rapidity();
                Jpsi_phi = Jpsi_p4.Phi();
                continue;
            }

            if( pdgId == 13 ) { // mu-
                muM_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass());
                continue;
            }

            if( pdgId == -13 ) { // mu+
                muP_p4.SetPtEtaPhiM(itParticle->pt(),itParticle->eta(),itParticle->phi(),itParticle->mass());
                continue;
            }

            if( pdgId == 22 ) { // photon
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
RootupleChicGen::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
RootupleChicGen::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
RootupleChicGen::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
RootupleChicGen::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
RootupleChicGen::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
RootupleChicGen::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RootupleChicGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RootupleChicGen);
