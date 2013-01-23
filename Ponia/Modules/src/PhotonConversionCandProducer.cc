#include <Ponia/Modules/src/PhotonConversionCandProducer.h>
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "Ponia/Modules/interface/VertexReProducer.h"

#include <TMath.h>
#include <boost/foreach.hpp>
#include <vector>

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

bool lt_(std::pair<double,short> a, std::pair<double,short> b) { 
     return a.first < b.first; }

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

  convCollection_          = ps.getParameter<edm::InputTag>("conversions");
  diMuonCollection_        = ps.getParameter<edm::InputTag>("dimuons");
  pfPhotonCollection_      = ps.getParameter<edm::InputTag>("pfphotons");
  thePVs_                  = ps.getParameter<edm::InputTag>("primaryVertexTag");
  wantTkVtxCompatibility_  = ps.getParameter<bool>("wantTkVtxCompatibility");
  sigmaTkVtxComp_          = ps.getParameter<uint32_t>("sigmaTkVtxComp");
  wantCompatibleInnerHits_ = ps.getParameter<bool>("wantCompatibleInnerHits");
  TkMinNumOfDOF_           = ps.getParameter<uint32_t>("TkMinNumOfDOF");

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
  TkVtxC = 0;
  CInnerHits = 0;
}


void PhotonConversionCandProducer::produce(edm::Event& event, const edm::EventSetup& esetup){
  std::auto_ptr<reco::ConversionCollection> outCollection(new reco::ConversionCollection);
   
  edm::Handle<reco::VertexCollection> priVtxs;
  event.getByLabel(thePVs_, priVtxs);
    
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
     
    if (!wantTkVtxCompatibility_ && !wantCompatibleInnerHits_)
       outCollection->push_back(*conv); 
    
    bool flagTkVtxCompatibility  = false;
    bool flagCompatibleInnerHits = false;
	           	
    if (wantTkVtxCompatibility_ && checkTkVtxCompatibility(*conv,*priVtxs.product())) {
          flagTkVtxCompatibility = true; TkVtxC++;}
            
    if (wantCompatibleInnerHits_ && conv->tracks().size()==2) {
         reco::HitPattern hitPatA=conv->tracks().at(0)->hitPattern();
         reco::HitPattern hitPatB=conv->tracks().at(1)->hitPattern();
         if(foundCompatibleInnerHits(hitPatA,hitPatB) && foundCompatibleInnerHits(hitPatB,hitPatA) )
           {
	    flagCompatibleInnerHits = true; CInnerHits++;}
    	 }
    if (!flagTkVtxCompatibility && flagCompatibleInnerHits && conv->tracks().at(0)->ndof() > TkMinNumOfDOF_ && conv->tracks().at(1)->ndof() > TkMinNumOfDOF_ ){
       	outCollection->push_back(*conv);
    }		     
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

bool PhotonConversionCandProducer::
checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection& priVtxs){

  std::vector< std::pair< double, short> > idx[2];
  short ik=-1;
  BOOST_FOREACH(edm::RefToBase<reco::Track> tk, conv.tracks()){
    ik++;
    short count=-1;
    BOOST_FOREACH(const reco::Vertex& vtx,priVtxs){
      count++;
    
      double dz_= tk->dz(vtx.position());
      double dzError_=tk->dzError();
      dzError_=sqrt(dzError_*dzError_+vtx.covariance(2,2));

      if(fabs(dz_)/dzError_ > sigmaTkVtxComp_) continue;
      
      idx[ik].push_back(std::pair<double,short>(fabs(dz_),count));
    }
    if(idx[ik].size()==0) {return false;}
    
    std::stable_sort(idx[ik].begin(),idx[ik].end(),lt_);
  }
  
  if(idx[0][0].second==idx[1][0].second || idx[0][1].second==idx[1][0].second || idx[0][0].second==idx[1][1].second){
    return true;
 } 
  return false;
}

bool PhotonConversionCandProducer::
foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB){
  size_t count=0;
  uint32_t oldSubStr=0;
  for (int i=0; i<hitPatA.numberOfHits() && count<2; i++) {
    uint32_t hitA = hitPatA.getHitPattern(i);
    if (!hitPatA.validHitFilter(hitA) || !hitPatA.trackerHitFilter(hitA)) continue;
    
    if(hitPatA.getSubStructure(hitA)==oldSubStr && hitPatA.getLayer(hitA)==oldSubStr)
      continue;

    if(hitPatB.getTrackerMonoStereo(hitPatA.getSubStructure(hitA),hitPatA.getLayer(hitA)) != 0)
      return true;
    
    oldSubStr=hitPatA.getSubStructure(hitA);
    count++;
  } 
  return false;  
}

void PhotonConversionCandProducer::endJob(){
   std::cout << "Eventi con slection fail: " << selection_fail << std::endl;
   std::cout << "Eventi con algo fail: " << algo_fail << std::endl;
   std::cout << "Eventi con quality fail: " << flag_fail << std::endl;
   std::cout << "Duplicates totali trovati: " << duplicates << std::endl;
   std::cout << "Events with conversion track and vertex compatibility: " << TkVtxC << std::endl;
   std::cout << "Events with compatible inner hits for conversion tracks: " << CInnerHits << std::endl;
}
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonConversionCandProducer);
