#ifndef gen_Pythia6CustomPtGun_h
#define gen_Pythia6CustomPtGun_h

#include "Pythia6ParticleGun.h"

#include "TRandom3.h"
/* Generate events with pt that follows the distribution 
   given in user supplied histogram

   Produce alternatively fParticle1Id and fParticle2Id

   Id: Pythia6CustomPtGun.h,v 1.1 2011/10/31 10:17:29 argiro Exp
   Author: Stefano Argiro
 */

class TH1F;

namespace gen {

   class Pythia6CustomPtGun : public Pythia6ParticleGun
   {
   
      public:
      
      Pythia6CustomPtGun( const edm::ParameterSet& );
      virtual ~Pythia6CustomPtGun();
      // void produce( edm::Event&, const edm::EventSetup& ) ;
      
      protected:
         void generateEvent() ;
      
      private:
      
         double  fMinEta;
	 double  fMaxEta;
	 double  fMinPt ;
         double  fMaxPt ;
	 bool    fAddAntiParticle;
	 double  fRatioPart1;
	 /// Pt distribution to be used to generate particels
	 TH1F*   fUserDist;
	 TRandom3* randomgen;

	 int     fParticle1Id;
	 int     fParticle2Id;  
	   
#ifdef DEBUG	 
	 TH1F*   fTestHisto;
#endif 
	 int     fEvtn;    // event counter
   
   };


}

#endif
