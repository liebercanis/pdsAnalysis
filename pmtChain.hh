/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 12 16:29:00 2017 by ROOT version 5.34/36
// from TTree pmt_tree/photon detector
// found on file: PDSout_06-11-1756_0.root
// ... modified .... M.Gold
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex>
#include <valarray>

#include "TPmtEvent.hxx"
#include "TPmtSummary.hxx"
#include "TPmtGains.hxx"

#include <TROOT.h>
#include <TSystem.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TF1.h>
//   ****back to the V1720 
typedef std::complex<double> Complex;

class pmtChain {
public :
  enum {MAXSAMPLES=2100};
  enum {NB=3,NCPMT=7,NC=NCPMT+1};
  enum {NALLCH=NB*NC};
  enum {NPMT=NCPMT*NB};
  enum {MAXADC=4095};
  enum {THRESHOLDHIGH=1,THRESHOLDLOW=0};
  enum {MAXBUFF=10};
  
   /* yujing 
       1614   2032-2035            B2 has 4 extra events.
       1646   3561-3562            B0 has 2 extra events.
       1714   2058-2059(2054-2055) B2 has 2 extra events.
       1742   2096(2090)           B2 has 1 extra event.    
       1756   347-348(346)         B0 misses 1 event between 347 & 348
       1853   4492-4497            B1 has 6 extra events than B0.
              4499-4500            B1 has 4 extra events than B2.

       >>>  4492 board 1 + 6, board 2 + 2

       1858   1307(1301)           B1 has 1 extra event than B0.
              1310(1303)           B1 has 2 extra events than B2.
       1910   4110(4104)           B2 has 1 extra event
       1914   467(460)             B2 is missing one event
       1914   1110(1104)           B2 has 1 extra event
       *1914   1303                B1 1 extra
       *1914   1303                B2 2 extra
       2020   1066-1067(1059-1060) B1 misses 1 event between 1066 & 1067
              2635-2636(2629-2630) B1 misses 1 event between 2635 & 2636
       2025   467(462)             B1 has 1 extra event than B0.
              471-472(464-465)     B2 has 2 extra events than B0.
       2054   210-212(204-206)     B1 misses 1 event between 210-212.
       2059   852(843)             B2 has 1 extra event
       */
//actually, due to the accumulation of the extra or missing events, the following extra or missing event number should be also shifted.
    enum {NFIX=19};
    enum {MAXEVENT=5000};
    // negative event numbers correspond to missing events
    Int_t skipMin[NFIX]=  {1614,1646,1714,1742,1756,1853,1853,1858, 1858,1910, 1914, 1914, 2020, 2020,2025,2025,2054,2059};
    Int_t skipFirst[NFIX]={2032,3561,2054,2090,-346,4491,4492,1301,-1303,4104, -460, 1104,-1059,-2629, 462, 464,-204, 843};
    Int_t skipLast[NFIX]= {2035,3562,2055,2090,-346,4496,4493,1301,-1303,4104, -460, 1104,-1060,-2630, 462, 465,-206, 843};
    Int_t skipBoard[NFIX]={   2,   0,   2,   2,   0,   1,   2,   1,    2,   2,    2,    2,    1,    1,   1,   2,   1,   2};
  
    enum { MAXRUN=93};
    Int_t misAlignCount[MAXRUN];
  
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree           *outTree; 
   TTree           *tree;
   Int_t           fCurrent; //!current Tree number in a TChain
  


   // Declaration of leaf types
   UInt_t          event_number;
   Int_t           computer_secIntoEpoch;
   Long64_t        computer_nsIntoSec;
   UInt_t          gps_nsIntoSec;
   UInt_t          gps_secIntoDay;
   UShort_t        gps_daysIntoYear;
   UShort_t        gps_Year;
   UShort_t        gps_ctrlFlag;
   UInt_t          digitizer_size[NB];
   UInt_t          digitizer_chMask[NB][NC];
   UInt_t          digitizer_evNum[NB];
   UInt_t          digitizer_time[NB];
   UShort_t        digitizer_waveforms[NB][NC][MAXSAMPLES];
   UInt_t          nDigitizers;
   UInt_t          nChannels;
   UInt_t          nSamples;
   UInt_t          nData;

   // List of branches
   TBranch        *b_event_number;   //!
   TBranch        *b_computer_secIntoEpoch;   //!
   TBranch        *b_computer_nsIntoSec;   //!
   TBranch        *b_gps_nsIntoSec;   //!
   TBranch        *b_gps_secIntoDay;   //!
   TBranch        *b_gps_daysIntoYear;   //!
   TBranch        *b_gps_Year;   //!
   TBranch        *b_gps_ctrlFlag;   //!  
   TBranch        *b_digitizer_size;   //!
   TBranch        *b_digitizer_chMask;   //!
   TBranch        *b_digitizer_evNum;   //!
   TBranch        *b_digitizer_time;   //!
   TBranch        *b_digitizer_waveforms;   //!
   TBranch        *b_nDigitizers;   //!
   TBranch        *b_nChannels;   //!
   TBranch        *b_nSamples;   //!
   TBranch        *b_nData;   //!

   // for output chain
  // Declaration of leaf types
   UInt_t          oevent_number;
   Int_t           ocomputer_secIntoEpoch;
   Long64_t        ocomputer_nsIntoSec;
   UInt_t          ogps_nsIntoSec;
   UInt_t          ogps_secIntoDay;
   UShort_t        ogps_daysIntoYear;
   UShort_t        ogps_Year;
   UShort_t        ogps_ctrlFlag;
   UInt_t          odigitizer_size[NB];
   UInt_t          odigitizer_chMask[NB][NC];
   UInt_t          odigitizer_evNum[NB];
   UInt_t          odigitizer_time[NB];
   UShort_t        odigitizer_waveforms[NB][NC][MAXSAMPLES];
   UInt_t          onDigitizers;
   UInt_t          onChannels;
   UInt_t          onSamples;
   UInt_t          onData;

   // buffering 
   Long64_t        eventNumberBuff[NB];
   UShort_t        waveBuff[NB][NC][MAXSAMPLES];
   UInt_t          sizeBuff[NB];
   UInt_t          chMaskBuff[NB][NC];
   UInt_t          evNumBuff[NB];
   UInt_t          timeBuff[NB];

   

   pmtChain(Int_t maxLoop=0, Long64_t firstEntry=0);
   virtual ~pmtChain();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     InitOutChain();
   UInt_t Loop(UInt_t nToLoop=0, UInt_t firstEntry=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   Int_t eventBuff[NB];
   void getBoard(Long64_t entry,Int_t board);
   void fillBoard(Int_t board);
   
    // file handling
   std::string tag;
   void getTag(std::string fname) { tag = fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")); return;}
   Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
   Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
   Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
   Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}
   std::vector<UInt_t> minList;
   std::vector<UInt_t> segList;
  
   void makeChain();
   std::vector<Int_t> fileTime;
   std::map<int,std::string> timeMap;
   string stag;
   TFile* outFile;

   void copyChain() 
   {
     oevent_number= event_number;
     ocomputer_secIntoEpoch= computer_secIntoEpoch;
     ocomputer_nsIntoSec= computer_nsIntoSec;
     ogps_nsIntoSec= gps_nsIntoSec;
     ogps_secIntoDay= gps_secIntoDay;
     ogps_daysIntoYear= gps_daysIntoYear;
     ogps_Year= gps_Year;
     ogps_ctrlFlag= gps_ctrlFlag;
     for(int ib=0; ib<NB; ++ib) {
       odigitizer_size[ib]=digitizer_size[ib];
       odigitizer_size[ib]=digitizer_size[ib];
       odigitizer_time[ib]=digitizer_time[ib];
       for(int ic=0; ic<NC; ++ic) {
         odigitizer_chMask[ib][ic]=digitizer_chMask[ib][ic];
         for(int is=0; is<MAXSAMPLES; ++is) {
           odigitizer_waveforms[ib][ic][is]=digitizer_waveforms[ib][ic][is];
         }
       }
     }
     onDigitizers=nDigitizers;
     onChannels=nChannels;
     onSamples=nSamples;
     onData=nData;
   }
   
  
   // returns -1 if pmt does not exist 
   // populate 3 boards, each from channel 0-6.  Channel 7 is the RF pulse. 
   // valid pmt are 0 to 20, RF channels are 21,22,23
   void fromPmtNumber(int ipmt, int& ib, int&ic)
   {
     ib=-1; ic=-1;
     if(ipmt<0) return;
     if(ipmt>=NPMT) {
       ib=ipmt-NPMT;
       ic = 7;
     } else {
       ib=(ipmt-ipmt%NCPMT)/NCPMT;
       ic= ipmt%NCPMT;
     }
     return;
   }

   int toPmtNumber(int ib, int ic) 
   {
     int ipmt=-1;
     if(ic<NCPMT) ipmt=ic+NCPMT*ib;
     else ipmt = ib+NPMT; 
     return ipmt;
   }
   // by run addative alignmentes
};
