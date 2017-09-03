/**
** MG, Sept 1 2017 
**/
#ifndef TPMTSUMMARY_DEFINED
#define TPMTSUMMARY_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>

using namespace std;

// class to store info for the event

class TPmtSummary: public TNamed {
	public:
    enum {NPMT=21};
		TPmtSummary();
		~TPmtSummary();
		void clear();
		void print();
		// data elements
    TString tag;
    Int_t ntrig555;
    Int_t ntrig5xx;
    Int_t ntrig444;
    Int_t ntrig4xx;
    Int_t ntrig111;
    Int_t ntrig1xx;
    Int_t ntrig000;
    Int_t ntrig0xx;
    Double_t qsum[NPMT];
    Double_t eqsum[NPMT];
    Double_t qmax[NPMT];
    Double_t eqmax[NPMT];
    Double_t noise[NPMT];
    Int_t norm[NPMT];
		ClassDef(TPmtSummary,1)
};
#endif

