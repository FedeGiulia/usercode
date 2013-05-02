#ifndef gen_Pythia6CustomPtGun2Histos_h
#define gen_Pythia6CustomPtGun2Histos_h

#include "Pythia6ParticleGun.h"

/* Generate events with pt that follows the distribution 
   given in user supplied histogram

   Produce alternatively fParticle1Id and fParticle2Id

   Id: Pythia6CustomPtGun2Histos.h,v 1.1 2011/10/31 10:17:29 argiro Exp
   Author: Stefano Argiro
 */

class TH1F;

namespace gen {

   class Pythia6CustomPtGun2Histos : public Pythia6ParticleGun
   {
   
      public:
      
      Pythia6CustomPtGun2Histos( const edm::ParameterSet& );
      virtual ~Pythia6CustomPtGun2Histos();
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
	 TH1F*   fUserDist;
	 TH1F*   fUserDist1;
	 TH1F*   fUserDist2;

	 int     fParticle1Id;
	 int     fParticle2Id;  
	   
#ifdef DEBUG	 
	 TH1F*   fTestHisto;
#endif 
	 int     fEvtn;    // event counter
   
   };


}

#endif
