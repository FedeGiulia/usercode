#include <Ponia/Modules/src/PhotonConversionCandProducer.h>
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include <TMath.h>
#include <boost/foreach.hpp>

// to order from high to low ProbChi2
bool ConversionLessByChi2(const reco::Conversion& c1, const reco::Conversion& c2){
  return TMath::Prob(c1.conversionVertex().chi2(),c1.conversionVertex().ndof()) > TMath::Prob(c2.conversionVertex().chi2(),c2.conversionVertex().ndof());
}


bool ConversionEqualByTrack(const reco::Conversion& c1, const reco::Conversion& c2){
  
  bool atLeastOneInCommon=false;
  BOOST_FOREACH(const edm::RefToBase<reco::Track> tk1, c1.tracks()){
    BOOST_FOREACH(const edm::RefToBase<reco::Track> tk2, c2.tracks()){
      if(tk1==tk2){
	atLeastOneInCommon=true;
	break;
      }
    } 
  } 
  return atLeastOneInCommon;
  

}

// define operator== for conversions, those with at least one track in common
namespace reco {

  bool operator==(const reco::Conversion& c1, const reco::Conversion& c2) {
    return c1.tracks()[0] == c2.tracks()[0] ||
           c1.tracks()[1] == c2.tracks()[1] ||
	   c1.tracks()[1] == c2.tracks()[0] ||
           c1.tracks()[0] == c2.tracks()[1] ;
  }
}

PhotonConversionCandProducer:: PhotonConversionCandProducer(const edm::ParameterSet& ps){

  convCollection_      = ps.getParameter<edm::InputTag>("conversions");
  diMuonCollection_    = ps.getParameter<edm::InputTag>("dimuons");
  pfPhotonCollection_  = ps.getParameter<edm::InputTag>("pfphotons");

  std::string algo = ps.getParameter<std::string>("convAlgo");
  convAlgo_ = StringToEnumValue<reco::Conversion::ConversionAlgorithm>(algo);

  std::vector<std::string> qual = 
    ps.getParameter<std::vector<std::string> >("convQuality"); 
  if( qual[0] != "" ) convQuality_ =StringToEnumValue<reco::Conversion::ConversionQuality>(qual);

  convSelectionCuts_ = ps.getParameter<std::string>("convSelection");
  
  produces<reco::ConversionCollection>("conversions");

  selection_fail = 0;
  algo_fail = 0;
  flag_fail = 0;
  duplicates = 0;
}


void PhotonConversionCandProducer::produce(edm::Event& event, const edm::EventSetup& esetup){
  std::auto_ptr<reco::ConversionCollection> outCollection(new reco::ConversionCollection);

  edm::Handle<reco::ConversionCollection> pConv;
  event.getByLabel(convCollection_,pConv);
  
  edm::Handle<pat::CompositeCandidateCollection> pDiMuons;
  event.getByLabel(diMuonCollection_,pDiMuons);

  edm::Handle<reco::PFCandidateCollection> pPFPhotons;
  event.getByLabel(pfPhotonCollection_,pPFPhotons);

  StringCutObjectSelector<reco::Conversion> *convSelection_ = new StringCutObjectSelector<reco::Conversion>(convSelectionCuts_);

  for(reco::ConversionCollection::const_iterator conv = pConv->begin(); conv != pConv->end(); ++conv){

    if (! ( *convSelection_)(*conv)){
	selection_fail++;	
	continue; // selection string
    }
    if (conv->algo()!= convAlgo_){
	algo_fail++;	
	continue; // select algorithm
    }
    if(convQuality_.size() > 0){
	bool flagsok=true;
	for (std::vector<int>::const_iterator flag = convQuality_.begin(); flag!=convQuality_.end(); ++flag){
	reco::Conversion::ConversionQuality q = (reco::Conversion::ConversionQuality)(*flag);
           if (!conv->quality(q)) {
	      flagsok=false;
	      break;
           }
        }
	if (!flagsok){
	   flag_fail++;
	   continue;
        }
    }

    outCollection->push_back(*conv);
  }

  removeDuplicates(*outCollection);
  event.put(outCollection,"conversions");

  // compatible inner hits ?
  delete convSelection_;
}


/** Put in out collection only those conversion candidates that are not sharing tracks.
    If sharing, keep the one with the best chi2.
 */
void PhotonConversionCandProducer::removeDuplicates(reco::ConversionCollection& c){
  // first sort from high to low chi2 prob
  std::sort(c.begin(),c.end(),ConversionLessByChi2);
  int iter1 = 0;
  // Cycle over all the elements of the collection and compare to all the following, 
  // if two elements have at least one track in common delete the element with the lower chi2
  while(iter1 < (( (int) c.size() ) - 1) ){
     int iter2 = iter1+1;
     while( iter2 < (int) c.size()) if(ConversionEqualByTrack( c[iter1], c[iter2] ) ){
        c.erase( c.begin() + iter2 );
	duplicates++;
        }else{
        iter2++;	// Increment index only if this element is no duplicate. 
			// If it is duplicate check again the same index since the vector rearranged elements index after erasing
        }
     iter1++;
  }

}

void PhotonConversionCandProducer::endJob(){
   std::cout << "Eventi con slection fail: " << selection_fail << std::endl;
   std::cout << "Eventi con algo fail: " << algo_fail << std::endl;
   std::cout << "Eventi con quality fail: " << flag_fail << std::endl;
   std::cout << "Duplicates totali trovati: " << duplicates << std::endl;
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonConversionCandProducer);
