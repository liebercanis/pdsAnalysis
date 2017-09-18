/**
** MG, July 2017 
**/
#ifndef TPMTHIT_DEFINED
#define TPMTHIT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TVector3.h>
#include <vector>

using namespace std;




// class to store info for the event

class TPmtHit: public TNamed {
	public:
		TPmtHit();
		~TPmtHit();
		void clear();
		void print(int ihit=0);
		// data elements
    Int_t ipmt;
	  Int_t time;
    Int_t tstart;
    Int_t tstop;
    Double_t ratio;
	  Double_t qhit;
	  Double_t qpeak;
	  Double_t qUnhit;
	  Double_t qUnpeak;
    Double_t fwhm; 
    Double_t peakTime;
    Double_t offset;

    Int_t  nsamples;
		std::vector<Double_t>  qsample;		 
		ClassDef(TPmtHit,1)
};
#endif

