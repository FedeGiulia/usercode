// -*- C++ -*-
//
// Package:    DumpConversion
// Class:      DumpConversion
// 
/**\class DumpConversion DumpConversion.cc Ponia/DumpConversion/src/DumpConversion.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alessandro Degano,32 1-C13,+41227678098
//         Created:  Tue Jan  8 15:07:52 CET 2013
// $Id$
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
#include <TMath.h>
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
//
// class declaration
//

class DumpConversion : public edm::EDAnalyzer {
   public:
      explicit DumpConversion(const edm::ParameterSet&);
      ~DumpConversion();

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
      edm::InputTag convCollection_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DumpConversion::DumpConversion(const edm::ParameterSet& conf)

{
   convCollection_ = conf.getParameter<edm::InputTag>("conversions");

}


DumpConversion::~DumpConversion()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DumpConversion::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::ConversionCollection> pConv;
  event.getByLabel(convCollection_,pConv);
  std::cout << "In this event there are " << pConv->size() << " conversions, their vertex chi2 is:" << std::endl;
  for(reco::ConversionCollection::const_iterator conv = pConv->begin(); conv != pConv->end(); ++conv){
	std::cout << TMath::Prob(conv->conversionVertex().chi2(),conv->conversionVertex().ndof()) << " , ";
  }
  std::cout << std::endl;	

}


// ------------ method called once each job just before starting event loop  ------------
void 
DumpConversion::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DumpConversion::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DumpConversion::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DumpConversion::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DumpConversion::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DumpConversion::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DumpConversion::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DumpConversion);
