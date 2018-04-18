//Very simple functions to help plot PMT waveforms
#ifndef DISPLAY_DEFINED
#define DISPLAY_DEFINED

#include <unistd.h>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <TNamed.h>
#include <TObjArray.h>
#include <TChainElement.h>
#include "TDirectory.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

const UShort_t bogus = -999;
const int cWidth=800, cHeight=700;
const int ADCrange = 4095;

using namespace std;

class display : public TNamed {
  public :
    enum{NSAMPLES = 2100, NBOARDS = 3, NPMTS = 8};
    enum {MAXADC=4095};
    display();
    virtual ~display();
    void TestRange(TH1F *histogram, int &startBin, int &stopBin);
    void SetFilter(bool forBaseline = true, bool forDisplay = true);
    void OpenFile(const char *infile, const char *treename = "pmt_tree");
    void AddChain(const char *infiles, const char *treename = "pmt_tree");
    void AddTag(TString tag="07-31-1518_0",const char *treename = "pmt_tree");
    void ClearChain();
    void ClearTree(); 
    void SetupHistos(float offsetstepADC);
    void ADCfilter(int iB, int iC);
    std::vector<Int_t> getRFTimes(int ib);
    float GetBaseline(int iB, int iC, int baselinestart=10, int baselinewide=50);
    float GetBaseline(TH1F *h1, int baselinestart=10, int baselinewide=50);
    float GetSigma(TH1F *h1, int baselinestart=10, int baselinewide=50);
    bool FillHistos(int EvNum, float offsetstepADC = 50.);
    bool SumEvents(bool trim = false, float offsetstepADC = 50.);
    bool DrawEvent(int EvNum, int mask0 = 255, int mask1 = 255, int mask2 = 255, bool showSum = true, float offsetstepADC = 50.);
    float EstimateTS(int EvNum, float threshold);
    int GetT0Bin(TH1F *histogram, float threshold, int startBin, int stopBin);
    int GetMinimum(TH1F *histogram, float &minVal, int startBin, int stopBin, bool StepForward=true);
    void print(TString format=".jpeg")  {if(c1) c1->Print(format.Data()) ;}
    void SetTagList();
    TH1F *histo[NBOARDS][NPMTS];
    TH1F *histoDraw[NBOARDS][NPMTS];
    TH1F *hsum ;
    TH1F *hsumDraw ;
    TH1F *hmaster ;
    std::vector<std::string> tagList;
  private:
    UShort_t waveforms[NBOARDS][NPMTS][NSAMPLES];
    TString tag; // event tag
    TFile *tf;
    TTree *pmt_tree;
    TChain *chain;
    bool histosetup;
    bool prefilter, postfilter;
    int Trims[NBOARDS][NPMTS];
    int Mask0;
    int Mask1;
    int Mask2;
    TCanvas *c1;

    ClassDef(display,1);
   
};
#endif
