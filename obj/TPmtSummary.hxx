/**
** MG, Sept 1 2017 
**/
#ifndef TPMTSUMMARY_DEFINED
#define TPMTSUMMARY_DEFINED
#include <fstream>
#include <iostream>
#include <string>
#include <TNamed.h>
#include <TNtuple.h>

using namespace std;
const double L=23.2;//m
const double clight=0.299792458;//m/ns
const double nmass=939.565;//MeV
const Int_t maxTDiff = 20000000;// max time diff 20 milliseconds

//  compSec compNano RF 1 2 3 digi 1 2 3 tprompt tprompt(torf) tof ke trig nhits beamtrig delta_t
// class to store info for the data file labeled by tag 

class TPmtSummary: public TNamed {
	public:
    enum {NPMT=21};
    enum {NB=3};
		TPmtSummary();
		~TPmtSummary();
		void clear();
    void print(std::ostream &out= std::cout);
    void printFile();
    void printEvent(unsigned it); 
    Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
    Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
    Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
    Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}

    void calcBeamTrig() {
      beamtrig.clear();
      for(UInt_t it=0; it< vrf1.size() ; ++it) {
        Double_t aveRf = Double_t(vrf1[it]+vrf2[it]+vrf3[it])/3.0;
        Int_t b =0;
        if(aveRf > 495 && aveRf < 500) b=1.0;
        beamtrig.push_back(b);
      }
    }

    void checkBeamTrig() {
      //UInt_t last = 1000000000;
      for(UInt_t it=0; it< vrf1.size() ; ++it) {
        bool bTime=false;
        if( abs( int(vdtime1[it])-int(vdtime2[it])) > maxTDiff || 
            abs( int(vdtime1[it])-int(vdtime3[it])) > maxTDiff || 
            abs( int(vdtime2[it])-int(vdtime3[it])) > maxTDiff ) bTime=true;
        //if(bTime && it-last==1) bTime=false;
        //else last=it;
        Int_t check=0;
        bool aZero=false;
        if( vrf1[it] ==0 ||vrf2[it] ==0 || vrf3[it] ==0 ) aZero=true;
        if( vrf1[it] > 480 && vrf1[it] < 500 ) ++check;
        if( vrf2[it] > 480 && vrf2[it] < 500 ) ++check;
        if( vrf3[it] > 480 && vrf3[it] < 500 ) ++check;
        if(bTime&&aZero&&check>0) printEvent(it);
        //if(aZero) printEvent(it);
      }
    }

    void calcDeltaT() {
      Double_t lastRf=0;
      Double_t lastDtime=0;
      deltaT.clear();
      for(UInt_t it=0; it< beamtrig.size() ; ++it) {
        Double_t dtime =Double_t(vdtime1[it]+vdtime2[it]+vdtime3[it])/3.0; 
        if(beamtrig[it]>0) {
          lastRf = Double_t(vrf1[it]+vrf2[it]+vrf3[it])/3.0; 
          lastDtime = dtime;
        }
        Double_t dt = (dtime-lastDtime)*8.0-(459-lastRf)*4.0;
        deltaT.push_back(dt);
      }
    }

    // data elements
    Int_t run;
    Int_t min;
    Int_t seg;
    std::string tag;
    Int_t ntrig555;
    Int_t ntrig5xx;
    Int_t ntrig444;
    Int_t ntrig4xx;
    Int_t ntrig111;
    Int_t ntrig1xx;
    Int_t ntrig000;
    Int_t ntrig0xx;

    Double_t gammapeak;
    Double_t tZero;
    Double_t qsum[NPMT];
    Double_t qrf[NPMT]; // summed charged +/- 100 samples around RF time
    Double_t eqsum[NPMT];
    Double_t qmax[NPMT];
    Double_t eqmax[NPMT];
    Double_t noise[NPMT];
    Int_t norm[NPMT];
    Double_t gain[NPMT];
    Double_t gain_e[NPMT];
    // vectors for times
    std::vector<Int_t>     nhits;  // number of hits in this event
    std::vector<Double_t>  timeToRf; 
    std::vector<Int_t>     vtrig;
    std::vector<Int_t>     vrf1;
    std::vector<Int_t>     vrf2;
    std::vector<Int_t>     vrf3;
    std::vector<UInt_t>    vevent;
    std::vector<Long64_t>  ventry;
    std::vector<Double_t> vprompt1; 
    std::vector<Double_t> vprompt2; 
    std::vector<Double_t> vprompt3; 
    std::vector<UInt_t>  vdtime1;   // caen digitizer time 
    std::vector<UInt_t>  vdtime2;   // caen digitizer time 
    std::vector<UInt_t>  vdtime3;   // caen digitizer time 
    std::vector<Int_t>     vcompSec;
    std::vector<Long64_t>  vcompNano;
    std::vector<Double_t> tprompt;//ysun
    std::vector<Double_t> tof;//ysun
    std::vector<Double_t> ke;//ysun 
    std::vector<Int_t> beamtrig;
    std::vector<Double_t> deltaT;
    //std::vector<UInt_t>    vgpsNs;
    //std::vector<UInt_t>    vgpsSec;
    //std::vector<UShort_t>  vgpsDay;
    //std::vector<UShort_t>  vgpsYear;

		ClassDef(TPmtSummary,8)
};
#endif

