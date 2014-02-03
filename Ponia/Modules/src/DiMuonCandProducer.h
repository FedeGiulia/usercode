#ifndef ponia_modules_DiMuonCandProducer_h
#define ponia_modules_DiMuonCandProducer_h


// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include <vector>
#include "DataFormats/Math/interface/LorentzVector.h"

/** $Id: DiMuonCandProducer.h,v 1.1 2013/01/15 14:58:55 degano Exp $
 */

template<typename T>
struct GreaterByVProb {
  typedef T first_argument_type;
  typedef T second_argument_type;
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};


//
// class decleration
//

class DiMuonCandProducer : public edm::EDProducer {
 public:
  explicit DiMuonCandProducer(const edm::ParameterSet&);
 
 private:

  virtual void produce(edm::Event&, const edm::EventSetup&);
	int isTriggerMatched(pat::CompositeCandidate*);

  // ----------member data ---------------------------
 private:
  edm::InputTag muons_;
  edm::InputTag thePVs_;
  edm::InputTag thebeamspot_;
  bool addCommonVertex_, addMuonlessPrimaryVertex_;

  GreaterByPt<pat::CompositeCandidate> pTComparator_;
  GreaterByVProb<pat::CompositeCandidate> vPComparator_;

  InvariantMassFromVertex massCalculator;

  typedef std::vector<float> fl_vect;
  typedef std::vector<fl_vect> UpsilonCand;

  std::string dimuonSelectionCuts_;
  std::vector<std::string> HLTLastFilters;
  std::vector<std::string> theTriggerNames;


};



#endif
