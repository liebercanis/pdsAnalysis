/**
** MG, Jan 12 2018
**/
#ifndef TPMTALIGN_DEFINED
#define TPMTALIGN_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>

using namespace std;

// class to store info for the data file labeld by tag 

class TPmtAlign: public TNamed {
	public:
    enum {NPMT=21};
		TPmtAlign();
		~TPmtAlign();
		void clear();
		void print();
    std::string tag;
    std::vector<Double_t> align0;
    std::vector<Double_t> align1;
    std::vector<Double_t> align2;
		ClassDef(TPmtAlign,1)
};
#endif

