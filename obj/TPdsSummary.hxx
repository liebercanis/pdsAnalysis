/**
** summary of PDS data 
** MG, Sept 1 2017 
**/
#ifndef TPDSSUMMARY_DEFINED
#define TPDSSUMMARY_DEFINED
#include <iostream>
#include <TNamed.h>
#include <string.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
// local 
#include "TPmtSummary.hxx"
#include "TPmtEvent.hxx"

using namespace std;

// class for making summary
class TPdsSummary: public TNamed {
	public:
    
    enum {MAXSAMPLES=2100};
    enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
    enum {NPMT=NB*NCPMT};
    enum {NALLCH=NB*NC};
    enum {MAXADC=4095};

		TPdsSummary(TString dirName="PDS_beamtime_files");
		~TPdsSummary();
    std::vector<Int_t> findRFTimes(int ipmt, double& step); 
    Int_t triggerInfo();
    void ADCFilter(int iB, int iC);
    void run(Int_t fFirst=0, Int_t maxFiles=0);
    void loop(); 
    void readFile(UInt_t ifile);
    void getTag(std::string fname) { tag = fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")); return;}
    Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
    Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
    Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
    Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}

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

    void printFiles() {
      if(isEmpty) { 
        printf(" file list is empty \n");
        return;
      }
      printf( " have %lu files in list \n",fileList.size());
      for( unsigned ifile =0; ifile < fileList.size() ; ++ifile ) printf(" %i %s \n",ifile,fileList[ifile].c_str());
    }
 
    bool isEmpty; 
    std::string tag;
    TString dirName;
    TString fullDirName;
    std::vector<std::string> fileList;
    TTree *pmt_tree;
    TTree *noBeamTree;
    TTree *lowBeamTree;
    TTree *highBeamTree;
    UInt_t    event;
    Int_t     compSec;
    Long64_t  compNano;
    UInt_t    gpsNs;
    UInt_t    gpsSec;
    UShort_t  gpsDay;
    UShort_t  gpsYear;
    UShort_t digitizer_waveforms[NB][NC][NS];
    
    TFile *summaryFile;
    TTree *summaryTree;
    TPmtSummary* pmtSummary;
    TFile *noBeamFile;
    TFile *lowBeamFile;
    TFile *highBeamFile;
    Int_t badFiles;
    Int_t goodFiles;

		ClassDef(TPdsSummary,2)
};
#endif

