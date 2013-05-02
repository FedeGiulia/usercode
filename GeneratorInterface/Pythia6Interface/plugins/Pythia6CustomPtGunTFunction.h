#ifndef gen_Pythia6CustomPtGunTFunction_h
#define gen_Pythia6CustomPtGunTFunction_h

#include "Pythia6ParticleGun.h"
#include "TF1.h"

/* Generate events with pt that follows the distribution 
   given in user supplied histogram

   Produce alternatively fParticle1Id and fParticle2Id

   Id: Pythia6CustomPtGunTFunction.h,v 1.1 2011/10/31 10:17:29 argiro Exp
   Author: Stefano Argiro
 */

class TH1F;

namespace gen {

   class Pythia6CustomPtGunTFunction : public Pythia6ParticleGun
   {
   
      public:
      
      Pythia6CustomPtGunTFunction( const edm::ParameterSet& );
      virtual ~Pythia6CustomPtGunTFunction();
      // void produce( edm::Event&, const edm::EventSetup& ) ;
      
      protected:
         void generateEvent() ;
      
      private:
      
         double  fMinEta;
	 double  fMaxEta;
	 double  fMinPt ;
         double  fMaxPt ;
	 bool    fAddAntiParticle;

	 /// Pt distribution to be used to generate particels
//	 TH1F*   fUserDist;
	 TF1*		tfunc;
	 int     fParticle1Id;
	 int     fParticle2Id;  

#define DEBUG 1	   
#ifdef DEBUG	 
	 TH1F*   fTestHisto;
#endif 
	 int     fEvtn;    // event counter
   
   };


}

#endif
