
#include <Ponia/Modules/src/ChiCandProducer.h>
#include <DataFormats/EgammaCandidates/interface/ConversionFwd.h>
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

ChiCandProducer::ChiCandProducer(const edm::ParameterSet& ps):
  convCollection_(ps.getParameter<edm::InputTag>("conversions")),
  diMuonCollection_(ps.getParameter<edm::InputTag>("dimuons")),
  pfCandidateCollection_(ps.getParameter<edm::InputTag>("pfcandidates")),
  pi0SmallWindow_(ps.getParameter<std::vector<double> >("pi0SmallWindow")),
  pi0LargeWindow_(ps.getParameter<std::vector<double> >("pi0LargeWindow")),
  deltaMass_(ps.getParameter<std::vector<double> >("deltaMass")),
  dzMax_(ps.getParameter<double>("dzmax"))
{
  produces<pat::CompositeCandidateCollection>("chicand");
  produces<std::vector<std::vector<float> > >("piZeroRejectCand");
  candidates = 0;
  delta_mass_fail = 0;
  dz_cut_fail = 0;
}
 

void ChiCandProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::auto_ptr<pat::CompositeCandidateCollection> chiCandColl(new pat::CompositeCandidateCollection);
  std::auto_ptr<std::vector<std::vector<float> > > piZeroRejectCand( new std::vector<std::vector<float> > );

  bool pizero_rejected = false;
  bool* ptr_pizero_rejected;
  ptr_pizero_rejected = &pizero_rejected;

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  event.getByLabel(diMuonCollection_,dimuons);

  edm::Handle<reco::ConversionCollection> conversions;
  event.getByLabel(convCollection_,conversions);
  edm::Handle<reco::PFCandidateCollection> pfcandidates;
  event.getByLabel(pfCandidateCollection_,pfcandidates);
  
  // TODO: have a collection of refs
  const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);
 
  // the first candidate is the one 
  // with better vProb, take this one

  const pat::CompositeCandidate& dimuonCand = (*dimuons)[0]; 

  // loop on conversion candidates, make chi cand
  for (reco::ConversionCollection::const_iterator conv = conversions->begin(); conv!= conversions->end(); ++conv){

    *ptr_pizero_rejected = false;

    pat::CompositeCandidate photonCand = makePhotonCandidate(*conv);
    pat::CompositeCandidate chiCand = makeChiCandidate(dimuonCand,photonCand);
    
    if (!cutDeltaMass(chiCand,dimuonCand)){
	delta_mass_fail++;
	continue;
    }
    float dz = conv->dz(dimuonCand.vertex());
    chiCand.addUserFloat("dz",dz);

    if (!cutdz(dz)){
	dz_cut_fail++;	
	continue;
    }
    std::vector<float> invmc= invmCombinations(*conv, pfphotons, ptr_pizero_rejected);
    if (*ptr_pizero_rejected){
        pizero_fail++;
        continue;
    }
//    chiCand.addUserData("invmc",invmc);

    chiCandColl->push_back(chiCand);
    piZeroRejectCand->push_back(invmc);
    candidates++;    
  }
  event.put(chiCandColl,"chicand");
  event.put(piZeroRejectCand,"piZeroRejectCand");
}

void ChiCandProducer::endJob(){
  std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
  std::cout << "Dz fail: " << dz_cut_fail << std::endl;
  std::cout << "Found " << candidates << " candidates." << std::endl;
}
  

// create a collection of PF photons
const reco::PFCandidateCollection  
ChiCandProducer::selectPFPhotons(const reco::PFCandidateCollection& pfcandidates){
  
  reco::PFCandidateCollection pfphotons;

  for (reco::PFCandidateCollection::const_iterator cand =   pfcandidates.begin();
       cand != pfcandidates.end(); ++cand){
    if (cand->particleId() == reco::PFCandidate::gamma) pfphotons.push_back(*cand);
  }
  
  return  pfphotons;
}

const pat::CompositeCandidate 
ChiCandProducer::makePhotonCandidate(const reco::Conversion& conv){

  pat::CompositeCandidate photonCand;
  photonCand.setP4(convertVector(conv.refittedPair4Momentum()));
  photonCand.setVertex(conv.conversionVertex().position());

  // todo: add electrons as daughters

  return photonCand;
}

const pat::CompositeCandidate 
ChiCandProducer::makeChiCandidate(const pat::CompositeCandidate& dimuon, 
				  const pat::CompositeCandidate& photon){
  
  pat::CompositeCandidate chiCand;
  chiCand.addDaughter(dimuon,"dimuon");
  chiCand.addDaughter(photon,"photon");
  chiCand.setVertex(dimuon.vertex());
  
  reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
  chiCand.setP4(vChic);
  
  return chiCand;

}

// check if the mass difference is in desired range
bool ChiCandProducer::cutDeltaMass(const pat::CompositeCandidate& chiCand,
				   const pat::CompositeCandidate& dimuonCand){
  
  float deltam = chiCand.p4().M() - dimuonCand.p4().M();
  float m1     = deltaMass_[0];
  float m2     = deltaMass_[1];
  
  return (deltam > m1 && deltam<m2);
}



// return a vector containing the invariant mass combinations 
// of the conversion with PF photons
// Only combinations within the given window are retained

std::vector<float> 
ChiCandProducer::invmCombinations(const reco::Conversion& conv,
				  const reco::PFCandidateCollection& photons,
				  bool* _pizero_rejected ){
  
  std::vector<float> ret;
  float small1 = pi0SmallWindow_[0];
  float small2 = pi0SmallWindow_[1];
  float large1 = pi0LargeWindow_[0];
  float large2 = pi0LargeWindow_[1];

  for (reco::PFCandidateCollection::const_iterator photon = photons.begin();
       photon!=photons.end(); ++photon){
    
    float inv = (conv.refittedPair4Momentum() + photon->p4()).M(); 
    if (inv > small1 && inv < small2){
       *_pizero_rejected = true;
       return ret;
    }else if(inv > large1 && inv < large2){
       ret.push_back(inv);
    }else{
    }
  }

  return ret;

}

reco::Candidate::LorentzVector ChiCandProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(ChiCandProducer);
