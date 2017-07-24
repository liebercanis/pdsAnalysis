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
		TPmtEvent();
//		~TPmtEvent();

		void clear();
		void print(int ihit=0);
		// data elements
    Int_t    run;
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
		ClassDef(TPmtEvent,1)
};
#endif

