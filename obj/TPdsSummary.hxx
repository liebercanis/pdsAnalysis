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
    void run();
    void loop(); 
    void readFile(TString fileName);
    void getTag(std::string fname) { tag = fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")); return;}
    Int_t getMonth() { return stoi(tag.substr(0,2)) ;}
    Int_t getDay() { return stoi(tag.substr(3,2)) ;}
    Int_t getHour() { return stoi(tag.substr(6,4)) ;}
    Int_t getSegment() { return stoi(tag.substr(11,tag.find(".") -1  - 11));}

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
    UShort_t digitizer_waveforms[NB][NC][NS];
    TFile *summaryFile;
    TTree *summaryTree;
    TPmtSummary* pmtSummary;

		ClassDef(TPdsSummary,1)
};
#endif

