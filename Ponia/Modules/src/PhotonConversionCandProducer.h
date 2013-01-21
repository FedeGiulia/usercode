/**
   \file
   Declaration of PhotonConversionCandProducer

   \author Stefano Argiro
   \version $Id: PhotonConversionCandProducer.h,v 1.1 2013/01/15 14:58:57 degano Exp $
   \date 18 Dec 2012
*/

#ifndef __PhotonConversionCandProducer_h_
#define __PhotonConversionCandProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

static const char CVSId__PhotonConversionCandProducer[] = 
"$Id: PhotonConversionCandProducer.h,v 1.1 2013/01/15 14:58:57 degano Exp $";

/**
   Select photon conversions and produce a conversion candidate collection

 */

class PhotonConversionCandProducer : public edm::EDProducer {

 public:
  explicit PhotonConversionCandProducer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);
  virtual void endJob() ;
  void removeDuplicates(reco::ConversionCollection& c);
  bool checkTkVtxCompatibility(const reco::Conversion&, const reco::VertexCollection&);
  bool foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB); 
 
  edm::InputTag convCollection_;
  edm::InputTag diMuonCollection_;
  edm::InputTag pfPhotonCollection_;
  edm::InputTag thePVs_;
  bool        wantTkVtxCompatibility_;
  uint32_t    sigmaTkVtxComp_;
  bool        wantCompatibleInnerHits_;
  uint32_t    TkMinNumOfDOF_;
  
  int convAlgo_;
  std::vector<int>   convQuality_;
  
  int selection_fail;
  int algo_fail;
  int flag_fail;
  int duplicates;
  int TkVtxC;
  int CInnerHits;

  std::string convSelectionCuts_;

};


#endif // __PhotonConversionCandProducer_h_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "cd .. ; scram b"
// End:
