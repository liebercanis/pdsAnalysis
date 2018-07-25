/**
** MG, July 2017 
**/
#ifndef TPMTEVENT_DEFINED
#define TPMTEVENT_DEFINED
#include "TPmtHit.hxx"
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TPmtEvent: public TNamed {
	public:
  enum { TRIGUNKNOWN,TRIG555,TRIG5XX ,TRIG444,TRIG4XX,TRIG111,TRIG1XX,TRIG000,TRIG0XX};
		TPmtEvent();
//		~TPmtEvent();

		void clear();
		void print(int nHitsToPrint=0);
		// data elements
    std::string tag;
    Int_t trigType;
    Int_t    run;
    std::vector<Int_t> rft21;
    std::vector<Int_t> rft22;
    std::vector<Int_t> rft23;
    Int_t    event;
    Int_t    tpcTrig;
    Double_t    pdst;
    Int_t    compSec;
    Long64_t compNano;
    Double_t bclock;  // board 0 clock
    Double_t gapTime; // gap time
    Int_t gapNumber;
    UInt_t  dtime[3];   // caen digitizer time 
    Double_t trigTime[3];  // trigger time for board
    Double_t tPrompt; // one for each board
    Double_t ke;
    Double_t tPromptToRF;
    Double_t tRFave;
    Double_t promptLike;
    std::vector<Double_t> qmax;
    std::vector<Double_t> qsum;
    Int_t  nhits;
		std::vector<TPmtHit>  hit;
		ClassDef(TPmtEvent,8)
};
#endif

