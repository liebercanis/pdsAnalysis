/**
** MG, Oct 11 2017 
**/
#ifndef TPMTGAINS_DEFINED
#define TPMTGAINS_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>

using namespace std;

// class to store info for the data file labeld by tag 

class TPmtGains: public TNamed {
	public:
    enum {NPMT=21};
		TPmtGains();
		~TPmtGains();
		void clear();
		void print();
    std::string tag;    
    Double_t gain[NPMT];
    Double_t egain[NPMT];
		ClassDef(TPmtGains,1)
};
#endif

