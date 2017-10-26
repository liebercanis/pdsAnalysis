/**
** MG, Sept 1 2017 
**/
#ifndef TPMTSUMMARY_DEFINED
#define TPMTSUMMARY_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>

using namespace std;

// class to store info for the data file labeld by tag 

class TPmtSummary: public TNamed {
	public:
    enum {NPMT=21};
		TPmtSummary();
		~TPmtSummary();
		void clear();
		void print();
    Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
    Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
    Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
    Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}
		// data elements
    std::string tag;
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
    Double_t gain[NPMT];
    Double_t gain_e[NPMT];
    // vectors for times
    std::vector<Int_t>     vtrig;
    std::vector<Int_t>     vrf1;
    std::vector<Int_t>     vrf2;
    std::vector<Int_t>     vrf3;
    std::vector<UInt_t>    vevent;
    std::vector<Long64_t>  ventry;
    std::vector<UInt_t>  vdtime1;   // caen digitizer time 
    std::vector<UInt_t>  vdtime2;   // caen digitizer time 
    std::vector<UInt_t>  vdtime3;   // caen digitizer time 
    std::vector<Int_t>     vcompSec;
    std::vector<Long64_t>  vcompNano;
    //std::vector<UInt_t>    vgpsNs;
    //std::vector<UInt_t>    vgpsSec;
    //std::vector<UShort_t>  vgpsDay;
    //std::vector<UShort_t>  vgpsYear;

		ClassDef(TPmtSummary,2)
};
#endif

