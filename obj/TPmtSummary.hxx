/**
** MG, Sept 1 2017 
**/
#ifndef TPMTSUMMARY_DEFINED
#define TPMTSUMMARY_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TNtuple.h>

using namespace std;

// class to store info for the data file labeld by tag 

class TPmtSummary: public TNamed {
	public:
    const double L=23.2;//m
    const double clight=0.299792458;//m/ns
    const double nmass=939.565;//MeV
    enum {NPMT=21};
		TPmtSummary();
		~TPmtSummary();
		void clear();
		void print();
    Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
    Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
    Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
    Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}
    //TNtuple* ntNeutron= new TNtuple("ntNeutron"," neutrons ","TOF:KE");//ysun

    // neutron spectrum //
    /*void neutron(unsigned iev, Double_t& TOF, Double_t& KE) {
      Double_t offset = 4.0*tZero - L/clight;
      TOF = 4.0*vprompt[iev] - offset;
      Double_t beta = L/TOF/clight;
      Double_t gamma=1.0;
      if(beta>0&&beta<1) gamma = sqrt(1/(1-beta*beta));
      KE=nmass*(gamma-1);
    }

    void fillNeutrons() {
      if(vprompt.size()<1) return;
      Double_t TOF,KE;
      for(unsigned iev=0; iev<vprompt.size(); ++iev) {
        neutron(iev,TOF,KE);
        ntNeutron->Fill(TOF,KE);
      }
    }*///ysun

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
    Double_t tZero;
    Double_t qsum[NPMT];
    Double_t eqsum[NPMT];
    Double_t qmax[NPMT];
    Double_t eqmax[NPMT];
    Double_t noise[NPMT];
    Int_t norm[NPMT];
    Double_t gain[NPMT];
    Double_t gain_e[NPMT];
    Double_t gammapeak;
    // vectors for times
    std::vector<Int_t>     vtrig;
    std::vector<Int_t>     vrf1;
    std::vector<Int_t>     vrf2;
    std::vector<Int_t>     vrf3;
    std::vector<UInt_t>    vevent;
    std::vector<Long64_t>  ventry;
    std::vector<Double_t> vprompt; 
    std::vector<UInt_t>  vdtime1;   // caen digitizer time 
    std::vector<UInt_t>  vdtime2;   // caen digitizer time 
    std::vector<UInt_t>  vdtime3;   // caen digitizer time 
    std::vector<Int_t>     vcompSec;
    std::vector<Long64_t>  vcompNano;
    std::vector<Double_t> tprompt;//ysun
    std::vector<Double_t> tof;//ysun
    std::vector<Double_t> ke;//ysun
    //std::vector<UInt_t>    vgpsNs;
    //std::vector<UInt_t>    vgpsSec;
    //std::vector<UShort_t>  vgpsDay;
    //std::vector<UShort_t>  vgpsYear;

		ClassDef(TPmtSummary,4)//ysun
};
#endif

