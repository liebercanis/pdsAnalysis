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
		void print(int ihit=0);
		// data elements
    Int_t trigType;
    Int_t    run;
    std::vector<Int_t> rft21;
    std::vector<Int_t> rft22;
    std::vector<Int_t> rft23;
    Int_t    event;
    Int_t    tpcTrig;
    Int_t    pdsTrig;
    UShort_t  gpsYear;
    UShort_t  gpsDay;
    UInt_t    gpsSec;
    UInt_t   gpsNs;
    std::vector<Double_t> qmax;
    std::vector<Double_t> qsum;
    Int_t  nhits;
		std::vector<TPmtHit>  hit;		 
		ClassDef(TPmtEvent,2)
};
#endif

