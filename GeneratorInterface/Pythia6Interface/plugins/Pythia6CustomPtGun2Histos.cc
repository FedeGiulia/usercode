
#include <iostream>

#include "Pythia6CustomPtGun2Histos.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <TH1F.h>
#include <TFile.h>

using namespace edm;
using namespace gen;


/* Generate events with pt that follows the distribution 
   given in user supplied histogram

   Produce alternatively fParticle1Id and fParticle2Id

   Id: Pythia6CustomPtGun2Histos.cc,v 1.1 2011/10/31 10:17:29 argiro Exp
   Author: Stefano Argiro
 */


Pythia6CustomPtGun2Histos::Pythia6CustomPtGun2Histos( const ParameterSet& pset ) :
   Pythia6ParticleGun(pset)
{
   
  using std::string;

   // ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters"); //, defpset ) ;
   fMinEta     = pgun_params.getParameter<double>("MinEta"); // ,-2.2);
   fMaxEta     = pgun_params.getParameter<double>("MaxEta"); // , 2.2);
   fMinPt       = pgun_params.getParameter<double>("MinPt"); // ,  20.);
   fMaxPt       = pgun_params.getParameter<double>("MaxPt"); // , 420.);
   fAddAntiParticle = pgun_params.getParameter<bool>("AddAntiParticle"); //, false) ;  
  
   edm::FileInPath file1(pgun_params.getParameter<string>("HFilename1"));
   edm::FileInPath file2(pgun_params.getParameter<string>("HFilename2"));

   TFile hfile1(file1.fullPath().c_str());
   fUserDist1 = dynamic_cast<TH1F*>(hfile1.Get(pgun_params.getParameter<string>("HName1").c_str()));

   TFile hfile2(file2.fullPath().c_str());
   fUserDist2 = dynamic_cast<TH1F*>(hfile2.Get(pgun_params.getParameter<string>("HName2").c_str()));
   
   if (!fUserDist1) {
     edm::LogError("")<< "Cannot Find histogram " <<
       pgun_params.getParameter<string>("HName1")<< 
       " in file " << pgun_params.getParameter<string>("HFilename1");
       
   }
   if (!fUserDist2) {
     edm::LogError("")<< "Cannot Find histogram " <<
       pgun_params.getParameter<string>("HName2")<< 
       " in file " << pgun_params.getParameter<string>("HFilename2");
       
   }

#ifdef DEBUG
   fTestHisto = new TH1F("test","test",100,0,100);
#endif

   fEvtn=0;


   if (fPartIDs.size() !=2) {
     edm::LogError("CustomGun")<< " Need two particle ids for this module"; 
   }
}

Pythia6CustomPtGun2Histos::~Pythia6CustomPtGun2Histos()
{
  delete fUserDist;
  edm::LogError("ChicGun")<<"Generated Chi's: " << fEvtn;

#ifdef DEBUG
  TFile f("test","recreate");
  fTestHisto->Write();
  delete fTestHisto;
#endif

}

void Pythia6CustomPtGun2Histos::generateEvent()
{
  
 
 
  
   Pythia6Service::InstanceWrapper guard(fPy6Service);	// grab Py6 instance

   // now actualy, start cooking up the event gun 
   //

   // 1st, primary vertex
   //
   HepMC::GenVertex* Vtx = new HepMC::GenVertex( HepMC::FourVector(0.,0.,0.));

   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
     
   int ip=1;
   //   for ( size_t i=0; i<fPartIDs.size(); i++ )
   //   {
   //int particleID = fPartIDs[i]; // this is PDG - need to convert to Py6 !!!
   
   int particleID =0;  
   if (fEvtn%2){
		particleID = fPartIDs[0] ;
		fUserDist = fUserDist1;
	}
   else {
		particleID = fPartIDs[1];
		fUserDist = fUserDist2;
	}
   
   int py6PID = HepPID::translatePDTtoPythia( particleID );
   int dum = 0;
   double pt=0, mom=0, ee=0, the=0, eta=0;
   double mass = pymass_(py6PID);
   
   // fill p(ip,5) (in PYJETS) with mass value right now,
   // because the (hardcoded) mstu(10)=1 will make py1ent
   // pick the mass from there
   pyjets.p[4][ip-1]=mass; 
   
   double phi = (fMaxPhi-fMinPhi)*pyr_(&dum)+fMinPhi;
   
   eta  = (fMaxEta-fMinEta)*pyr_(&dum)+fMinEta;                                                      
   
   the  = 2.*atan(exp(-eta));                                                                          
   // use Von Neuman rejection method to generate a random
   // number according to user-supplied distribution
   
   double hmax = fUserDist->GetMaximum();
   
   while (1){
     pt   = (fMaxPt-fMinPt)*pyr_(&dum)+fMinPt;
     double test = hmax * pyr_(&dum);
     
     int bin = fUserDist->FindBin(pt);
     if  (fUserDist->GetBinContent(bin) > test) break;
   }
#ifdef DEBUG
   fTestHisto->Fill(pt);
#endif
   //
   
   mom = pt/sin(the);
   ee = sqrt(mom*mom+mass*mass);
   
   py1ent_(ip, py6PID, ee, the, phi);
   
   double px     = pyjets.p[0][ip-1]; // pt*cos(phi) ;
   double py     = pyjets.p[1][ip-1]; // pt*sin(phi) ;
   double pz     = pyjets.p[2][ip-1]; // mom*cos(the) ;
   
   HepMC::FourVector p(px,py,pz,ee) ;
   HepMC::GenParticle* Part = 
     new HepMC::GenParticle(p,particleID,1);
   Part->suggest_barcode( ip ) ;
   Vtx->add_particle_out(Part);
   
   if(fAddAntiParticle)
     {
       ip = ip + 1;
       HepMC::GenParticle* APart = addAntiParticle( ip, particleID, ee, eta, phi );
       if ( APart ) Vtx->add_particle_out(APart) ;	    
     }
   ip++;
   //   }
   
   fEvt->add_vertex(Vtx);
   
   // run pythia
   pyexec_();
   

   fEvtn++;
   return;
}

DEFINE_FWK_MODULE(Pythia6CustomPtGun2Histos);
