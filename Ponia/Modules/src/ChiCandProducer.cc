
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
//kinematic refit collections
  produces<pat::CompositeCandidateCollection>("chicCompCandRefit1s");
  produces<pat::CompositeCandidateCollection>("chicCompCandRefit2s");
  produces<pat::CompositeCandidateCollection>("chicCompCandRefit3s");

  candidates = 0;
  delta_mass_fail = 0;
  dz_cut_fail = 0;
  pizero_fail = 0;
}
 

void ChiCandProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::auto_ptr<pat::CompositeCandidateCollection> chiCandColl(new pat::CompositeCandidateCollection);
  std::auto_ptr<PiZeroRejCollection> piZeroRejectCand( new PiZeroRejCollection );
//Kinemati refit collections
  std::auto_ptr< pat::CompositeCandidateCollection > chicCompCandRefitColl1s(new pat::CompositeCandidateCollection);
  std::auto_ptr< pat::CompositeCandidateCollection > chicCompCandRefitColl2s(new pat::CompositeCandidateCollection);
  std::auto_ptr< pat::CompositeCandidateCollection > chicCompCandRefitColl3s(new pat::CompositeCandidateCollection);

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

//  const pat::CompositeCandidate& dimuonCand = (*dimuons)[0]; 

  for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand){
     // loop on conversion candidates, make chi cand
     for (reco::ConversionCollection::const_iterator conv = conversions->begin(); conv!= conversions->end(); ++conv){

    	*ptr_pizero_rejected = false;

	pat::CompositeCandidate photonCand = makePhotonCandidate(*conv);
	pat::CompositeCandidate chiCand = makeChiCandidate(*dimuonCand,photonCand);
    
	if (!cutDeltaMass(chiCand,*dimuonCand)){
	   delta_mass_fail++;
	   continue;
	}
    	float dz = conv->dz(dimuonCand->vertex());
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


// Kinematic fit

    	edm::ESHandle<TransientTrackBuilder> theB; 
	esetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
    
	float ChiBM_fit = -1;
    float ChiBPx_fit = -1;
    float ChiBPy_fit = -1;
    float ChiBPz_fit = -1;
    float ChiBVtxX_fit = -1;
    float ChiBVtxY_fit = -1;
    float ChiBVtxZ_fit = -1;
    float ChiBVtxP_fit=-1;
    float ChiBM_fit2s = -1;
    float ChiBPx_fit2s = -1;
    float ChiBPy_fit2s = -1;
    float ChiBPz_fit2s = -1;
    float ChiBVtxX_fit2s = -1;
    float ChiBVtxY_fit2s = -1;
    float ChiBVtxZ_fit2s = -1;
    float ChiBVtxP_fit2s=-1;
    float ChiBM_fit3s = -1;
    float ChiBPx_fit3s = -1;
    float ChiBPy_fit3s = -1;
    float ChiBPz_fit3s = -1;
    float ChiBVtxX_fit3s = -1;
    float ChiBVtxY_fit3s = -1;
    float ChiBVtxZ_fit3s = -1;
    float ChiBVtxP_fit3s=-1;

    reco::TrackRef JpsiTk[2]={
      ( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2") ) )->innerTrack()
    };
    
    reco::TrackCollection convTracks;
        
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
    
    if (!photonVertexFitTree->isValid()) { 
      edm::ParameterSet pSet;
      
      pSet.addParameter<double>("maxDistance", 3);
      pSet.addParameter<int>("maxNbrOfIterations", 10000); //10

      KinematicParticleVertexFitter fitter2(pSet);
      
      //RefCountedKinematicTree photonVertexFitTree;
      photonVertexFitTree = fitter2.fit(PhotonParticles);
      }

    if (photonVertexFitTree->isValid()) { 
    
      // now apply Photon mass constraint                                                                                                            
      KinematicParticleFitter csFitterPhoton;
      KinematicConstraint * pho_c = new MassKinematicConstraint(zero_mass,zero_sigma);
      // add mass constraint to the photon fit to do a constrained fit:                                                                              

      photonVertexFitTree->movePointerToTheTop();
      photonVertexFitTree = csFitterPhoton.fit(pho_c,photonVertexFitTree);
     
    if (photonVertexFitTree->isValid()) {
      
	const ParticleMass muonMass(0.1056583);
	float muonSigma = muonMass*1E-6;

	photonVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle fitPhoton = photonVertexFitTree->currentParticle();
      
	std::vector<RefCountedKinematicParticle> allChiBDaughters;
	allChiBDaughters.push_back(pFactory.particle (MuMuTT[0], muonMass, float(0), float(0), muonSigma));
	allChiBDaughters.push_back(pFactory.particle (MuMuTT[1], muonMass, float(0), float(0), muonSigma));
	allChiBDaughters.push_back(fitPhoton);

	         
        KinematicConstrainedVertexFitter constVertexFitter;
	double nominalUpsilonMass = 9.46;
          
        MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(nominalUpsilonMass);
	RefCountedKinematicTree ChiBTree = constVertexFitter.fit(allChiBDaughters,upsilon_mtc);
	if(!ChiBTree->isEmpty()) {
	     
	  ChiBTree->movePointerToTheTop();
	  RefCountedKinematicParticle fitChiB = ChiBTree->currentParticle();
	  RefCountedKinematicVertex ChiBDecayVertex = ChiBTree->currentDecayVertex();
          
	  if (fitChiB->currentState().isValid()) {          
	    ChiBM_fit = fitChiB->currentState().mass();
	    ChiBPx_fit = fitChiB->currentState().kinematicParameters().momentum().x();
	    ChiBPy_fit = fitChiB->currentState().kinematicParameters().momentum().y();
	    ChiBPz_fit = fitChiB->currentState().kinematicParameters().momentum().z();
	    ChiBVtxX_fit = ChiBDecayVertex->position().x();
	    ChiBVtxY_fit = ChiBDecayVertex->position().y();
	    ChiBVtxZ_fit = ChiBDecayVertex->position().z();
	    ChiBVtxP_fit=ChiSquaredProbability((double)(ChiBDecayVertex->chiSquared()),(double)(ChiBDecayVertex->degreesOfFreedom()));            
	    reco::CompositeCandidate aChicCandRefit1s(0,math::XYZTLorentzVector(ChiBPx_fit,ChiBPy_fit,ChiBPz_fit,sqrt(ChiBM_fit*ChiBM_fit+ChiBPx_fit*ChiBPx_fit+ChiBPy_fit*ChiBPy_fit+ChiBPz_fit*ChiBPz_fit)),math::XYZPoint(ChiBVtxX_fit,ChiBVtxY_fit,ChiBVtxZ_fit),50551);
	    pat::CompositeCandidate patChicCandRefit1s(aChicCandRefit1s);
	    patChicCandRefit1s.addUserFloat("vProb",ChiBVtxP_fit);
	    chicCompCandRefitColl1s->push_back(patChicCandRefit1s);

	  }
        }

        KinematicConstrainedVertexFitter constVertexFitter2s;
	double nominalUpsilonMass2s = 10.023;
          
        MultiTrackKinematicConstraint *upsilon_mtc2s = new  TwoTrackMassKinematicConstraint(nominalUpsilonMass2s);
	RefCountedKinematicTree ChiBTree2s = constVertexFitter2s.fit(allChiBDaughters,upsilon_mtc2s);
	if(!ChiBTree2s->isEmpty()) {
	      
	  ChiBTree2s->movePointerToTheTop();
	  RefCountedKinematicParticle fitChiB2s = ChiBTree2s->currentParticle();
	  RefCountedKinematicVertex ChiBDecayVertex2s = ChiBTree2s->currentDecayVertex();
          
	  if (fitChiB2s->currentState().isValid()) {          
	    ChiBM_fit2s = fitChiB2s->currentState().mass();
	    ChiBPx_fit2s = fitChiB2s->currentState().kinematicParameters().momentum().x();
	    ChiBPy_fit2s = fitChiB2s->currentState().kinematicParameters().momentum().y();
	    ChiBPz_fit2s = fitChiB2s->currentState().kinematicParameters().momentum().z();
	    ChiBVtxX_fit2s = ChiBDecayVertex2s->position().x();
	    ChiBVtxY_fit2s = ChiBDecayVertex2s->position().y();
	    ChiBVtxZ_fit2s = ChiBDecayVertex2s->position().z();
	    ChiBVtxP_fit2s=ChiSquaredProbability((double)(ChiBDecayVertex2s->chiSquared()),(double)(ChiBDecayVertex2s->degreesOfFreedom()));       
	    reco::CompositeCandidate aChicCandRefit2s(0,math::XYZTLorentzVector(ChiBPx_fit2s,ChiBPy_fit2s,ChiBPz_fit2s,sqrt(ChiBM_fit2s*ChiBM_fit2s+ChiBPx_fit2s*ChiBPx_fit2s+ChiBPy_fit2s*ChiBPy_fit2s+ChiBPz_fit2s*ChiBPz_fit2s)),math::XYZPoint(ChiBVtxX_fit2s,ChiBVtxY_fit2s,ChiBVtxZ_fit2s),50551);
	    pat::CompositeCandidate patChicCandRefit2s(aChicCandRefit2s);
	    patChicCandRefit2s.addUserFloat("vProb",ChiBVtxP_fit2s);              
	    chicCompCandRefitColl2s->push_back(patChicCandRefit2s);
          }
	}
     
        KinematicConstrainedVertexFitter constVertexFitter3s;
	double nominalUpsilonMass3s = 10.355;
          
        MultiTrackKinematicConstraint *upsilon_mtc3s = new  TwoTrackMassKinematicConstraint(nominalUpsilonMass3s);
	RefCountedKinematicTree ChiBTree3s = constVertexFitter3s.fit(allChiBDaughters,upsilon_mtc3s);
	if(!ChiBTree3s->isEmpty()) {
	      
	  ChiBTree3s->movePointerToTheTop();
	  RefCountedKinematicParticle fitChiB3s = ChiBTree3s->currentParticle();
	  RefCountedKinematicVertex ChiBDecayVertex3s = ChiBTree3s->currentDecayVertex();
          
	  if (fitChiB3s->currentState().isValid()) {          
	    ChiBM_fit3s = fitChiB3s->currentState().mass();
	    ChiBPx_fit3s = fitChiB3s->currentState().kinematicParameters().momentum().x();
	    ChiBPy_fit3s = fitChiB3s->currentState().kinematicParameters().momentum().y();
	    ChiBPz_fit3s = fitChiB3s->currentState().kinematicParameters().momentum().z();
	    ChiBVtxX_fit3s = ChiBDecayVertex3s->position().x();
	    ChiBVtxY_fit3s = ChiBDecayVertex3s->position().y();
	    ChiBVtxZ_fit3s = ChiBDecayVertex3s->position().z();
	    ChiBVtxP_fit3s=ChiSquaredProbability((double)(ChiBDecayVertex3s->chiSquared()),(double)(ChiBDecayVertex3s->degreesOfFreedom()));       
	    reco::CompositeCandidate aChicCandRefit3s(0,math::XYZTLorentzVector(ChiBPx_fit3s,ChiBPy_fit3s,ChiBPz_fit3s,sqrt(ChiBM_fit3s*ChiBM_fit3s+ChiBPx_fit3s*ChiBPx_fit3s+ChiBPy_fit3s*ChiBPy_fit3s+ChiBPz_fit3s*ChiBPz_fit3s)),math::XYZPoint(ChiBVtxX_fit3s,ChiBVtxY_fit3s,ChiBVtxZ_fit3s),50551);
	    pat::CompositeCandidate patChicCandRefit3s(aChicCandRefit3s);
	    patChicCandRefit3s.addUserFloat("vProb",ChiBVtxP_fit3s);              
	    chicCompCandRefitColl3s->push_back(patChicCandRefit3s);
          }
	}
      }

    }

 // End kinematic fit


	chiCandColl->push_back(chiCand);
	piZeroRejectCand->push_back(invmc);
	candidates++;    
     }
  }
  event.put(chiCandColl,"chicand");
  event.put(piZeroRejectCand,"piZeroRejectCand");
  event.put(chicCompCandRefitColl1s,"chicCompCandRefit1s");
  event.put(chicCompCandRefitColl2s,"chicCompCandRefit2s");
  event.put(chicCompCandRefitColl3s,"chicCompCandRefit3s");
}

void ChiCandProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "Chi Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Delta mass fail: " << delta_mass_fail << std::endl;
  std::cout << "Dz fail: " << dz_cut_fail << std::endl;
  std::cout << "Pi0 fail: " << pizero_fail << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " Chi candidates." << std::endl;
  std::cout << "###########################" << std::endl;
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
  
  return (deltam > m1 && deltam < m2);
}



// return a vector containing the invariant mass combinations 
// of the conversion with PF photons
// Only combinations within the given window are retained

std::vector<float> 
ChiCandProducer::invmCombinations(const reco::Conversion& conv,
				  const reco::PFCandidateCollection& photons,
				  bool* _pizero_rejected ){

  // 2 windows are defined for Pi0 rejection, Conversions that, paired with others photons from the event, have an
  // invariant mass inside the "small" window will be istantly rejected. Those that falls in the large window will
  // be saved inside a vector for later refined offline selection. The rest are just supposed to have nothing to do 
  // with Pi0 decay and therefore nothing is done with them.

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
