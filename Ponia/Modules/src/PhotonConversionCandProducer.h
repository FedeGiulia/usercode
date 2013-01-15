/**
   \file
   Declaration of PhotonConversionCandProducer

   \author Stefano Argiro
   \version $Id$
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


static const char CVSId__PhotonConversionCandProducer[] = 
"$Id$";

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

  edm::InputTag convCollection_;
  edm::InputTag diMuonCollection_;
  edm::InputTag pfPhotonCollection_;

  
  int convAlgo_;
  std::vector<int>   convQuality_;
  
  int selection_fail;
  int algo_fail;
  int flag_fail;
  int duplicates;

  std::string convSelectionCuts_;

};


#endif // __PhotonConversionCandProducer_h_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "cd .. ; scram b"
// End:
