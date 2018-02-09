//Very simple functions to help plot PMT waveforms

#include <unistd.h>
#include <cstdio>
#include <iostream>

#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"

const int NSAMPLES = 2100, NBOARDS = 3, NPMTS = 8;
const UShort_t bogus = -999;
UShort_t waveforms[NBOARDS][NPMTS][NSAMPLES];
const int cWidth=800, cHeight=700;
int ADCrange = 4095;
TString tag; // event tag
TFile *tf = 0;
TTree *pmt_tree = 0;
TChain *chain = 0;

bool histosetup = false;
TH1F *histo[NBOARDS][NPMTS];
TH1F *histoDraw[NBOARDS][NPMTS];
TH1F *hsum = 0;
TH1F *hsumDraw = 0;
TH1F *hmaster = 0;

bool prefilter = true, postfilter = true;

TCanvas *c1 = 0;

int Mask0 = 255;
int Mask1 = 255;
int Mask2 = 255;

int Trims[NBOARDS][NPMTS];

void TestRange(TH1F *histogram, int &startBin, int &stopBin) {
    
    if (stopBin < startBin) {
        int temp = stopBin;
        stopBin = startBin;
        startBin = temp;
    }
    
    int NBinsX = histogram->GetNbinsX();
    
    if (startBin < 1) { startBin = 1; }
    if (stopBin > NBinsX) { stopBin = NBinsX; }
    
}

void SetFilter(bool forBaseline = true, bool forDisplay = true) {
    
    prefilter = forBaseline;
    postfilter = forDisplay;
    
}

void OpenFile(const char *infile, const char *treename = "pmt_tree") {
    
    tf = new TFile(infile, "READ");
    pmt_tree = (TTree *)tf->Get(treename);
    pmt_tree->SetBranchAddress("digitizer_waveforms", &waveforms);

}

void AddChain(const char *infiles, const char *treename = "pmt_tree") {
    
    if (!chain) { chain = new TChain(treename); }
    int numfiles = chain->Add(infiles);
    
    chain->SetBranchAddress("digitizer_waveforms", &waveforms);
    
    if (numfiles == 0) { cout << "No files added" << endl; }
    cout << " chain has " << chain->GetEntries() << " entries " << endl;
    //chain->ls();
}

void ClearChain() {
    
    delete chain;
    chain = 0;
}

void ClearTree() {
    
    delete pmt_tree;
    pmt_tree = 0;
}

void SetupHistos(float offsetstepADC) {
    
    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            
            if (!histo[iB][iC]) {
                histo[iB][iC] = new TH1F(Form("B%iC%i_0",iB, iC), Form("B%i C%i",iB, iC), NSAMPLES, 0., 1.*NSAMPLES);
                histo[iB][iC]->SetStats(0);
                
                histoDraw[iB][iC] = new TH1F(Form("B%iC%i",iB, iC), Form("B%i C%i",iB, iC), NSAMPLES, 0., 1.*NSAMPLES);
                histoDraw[iB][iC]->SetStats(0);
                
                if (iC ==7) {
                    histoDraw[iB][iC]->SetLineColor(kGreen);
                }
                
            }
        }
    }
    
    hmaster = new TH1F("hmaster", "ADC Plot Event -1;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
    hmaster->GetYaxis()->SetRangeUser(-4.*offsetstepADC, 24*offsetstepADC);
    hmaster->GetYaxis()->SetTitleOffset(1.30);
    hmaster->SetStats(0);
    
    hsum = new TH1F("hsum", "PMT Sum;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
    hsumDraw = new TH1F("hsumDraw", "PMT Sum;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
    
    hsumDraw->SetLineColor(kRed);
    
    histosetup = true;
}

void ADCfilter(int iB, int iC) {
    
    for (int iS = 0; iS<NSAMPLES; ++iS) {
        if (waveforms[iB][iC][iS] > ADCrange) {
            if (iS > 0) { waveforms[iB][iC][iS] = waveforms[iB][iC][iS-1];}
            else {
                int iS2 = 0;
                while (waveforms[iB][iC][iS2] > ADCrange) {
                    waveforms[iB][iC][0] = waveforms[iB][iC][iS2+1];
                    ++iS2;
                }
            }
        }
    }
}

float GetBaseline(int iB, int iC, int baselinestart=10, int baselinewide=50) {
    float baseline = 0.;
    for (int iS = baselinestart; iS<(baselinestart+baselinewide); ++iS) {
        baseline += waveforms[iB][iC][iS]/(1.*baselinewide);
    }
    
    return baseline;
    
}

float GetBaseline(TH1F *h1, int baselinestart=10, int baselinewide=50) {
    
    int baselinestop = baselinestart+baselinewide;
    TestRange(h1, baselinestart, baselinestop);
    
    float baseline = 0.;
    for (int iS = baselinestart; iS<baselinestop; ++iS) {
        baseline += h1->GetBinContent(iS)/(1.*baselinewide);
    }

    return baseline;
}

float GetSigma(TH1F *h1, int baselinestart=10, int baselinewide=50) {
    
    int baselinestop = baselinestart+baselinewide;
    TestRange(h1, baselinestart, baselinestop);
    
    float tot = 0., tot2 = 0.;
    for (int iS = baselinestart; iS<baselinestop; ++iS) {
        tot += h1->GetBinContent(iS)/(1.*baselinewide);
        tot2 += h1->GetBinContent(iS)*h1->GetBinContent(iS)/(1.*baselinewide);
    }
    
    return sqrt(tot2-tot*tot);
}

bool FillHistos(int EvNum, float offsetstepADC = 50.) {
    
    if (!histosetup) { SetupHistos(offsetstepADC); }
    
    int NumEvents = 0;
    
    // default to a chain
    if (chain) { NumEvents = chain->GetEntries(); }
    else if (pmt_tree) { NumEvents = pmt_tree->GetEntries(); }
    else { cout << "No TChain or TTree" << endl; return false; }
    
    if (EvNum >= NumEvents) { cout << "Entry " << EvNum << " > " << NumEvents << " entries" << endl; return false; }
    
    // default to a chain
    if (chain) { chain->GetEntry(EvNum); }
    else if (pmt_tree) { pmt_tree->GetEntry(EvNum); }
    
    // now fill histograms
    float baseline = 0.;
    float offset = 0.;
    bool needS0 = false;

    float sumAll=0;

    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            
            if (prefilter) { ADCfilter(iB, iC); }
            baseline = GetBaseline(iB, iC);
            if (postfilter) { ADCfilter(iB, iC); }
            
            histo[iB][iC]->Reset();
            histoDraw[iB][iC]->Reset();
            float sum=0;
            for (int iS = 0; iS<NSAMPLES; ++iS) {
                histo[iB][iC]->Fill(iS+0.5, (1.*waveforms[iB][iC][iS]-baseline) );
                if (iC == 7) {
                    histoDraw[iB][iC]->Fill(iS+0.5, ((1.*waveforms[iB][iC][iS]-baseline)*offsetstepADC/(1.*ADCrange+1.)+offset) );
                }
                else {
                    histoDraw[iB][iC]->Fill(iS+0.5, (1.*waveforms[iB][iC][iS]-baseline+offset) );
                  sum-= float(waveforms[iB][iC][iS]-baseline);
                  sumAll-= float(waveforms[iB][iC][iS]-baseline);
                }

            }
            if(iC!=7) printf(" ib %i ic %i sum %f \n",iB,iC,sum);
            
            offset += offsetstepADC;
        }
    }
    
    hmaster->SetTitle(Form("ADC Plot File %s Event %i;Sample (4 ns per sample);ADC Counts",tag.Data(),EvNum));
    printf(" %s %f \n",hmaster->GetTitle(),sumAll);
    
    return true;

}

bool SumEvents(bool trim = false, float offsetstepADC = 50.) {
    
    hsum->Reset();
    hsumDraw->Reset();
    
    int sample;
    
    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            if (iC == 7) { continue; }
            
            for (int iS = 1; iS<=NSAMPLES; ++iS) {
                
                if (trim) {
                    sample = iS + Trims[iB][iC];
                    if (sample < 1) { sample = 1; }
                    else if (sample > NSAMPLES ) { sample = NSAMPLES; }
                }
                else { sample = iS; }
                hsum->Fill(iS-0.5, histo[iB][iC]->GetBinContent(sample) );
                hsumDraw->Fill(iS-0.5, histo[iB][iC]->GetBinContent(sample) );
            }
        }
    }
    
    for (int iS = 1; iS<=NSAMPLES; ++iS) { hsumDraw->Fill(iS-0.5, -offsetstepADC); }

    return true;
}

bool DrawEvent(int EvNum, int mask0 = 255, int mask1 = 255, int mask2 = 255, bool showSum = true, float offsetstepADC = 50.) {
    
  bool status = FillHistos(EvNum, offsetstepADC);
  status &= SumEvents(false, offsetstepADC);

  if (!status) { return false; }
  int mask = 0;

  bool redraw = false;
  if ( (mask0 != Mask0) || (mask1 != Mask1) || (mask2 != Mask2) ) { redraw = true; }
  Mask0 = mask0; Mask1 = mask1; Mask2 = mask2;

  if (!c1) { c1 = new TCanvas("c1", "", cWidth, cHeight); }

  hmaster->Draw();

  for (int iB = 0; iB<NBOARDS; ++iB) {

    if (iB == 0) { mask = Mask0; }
    else if (iB == 1) { mask = Mask1; }
    else if (iB == 2) { mask = Mask2; }

    for (int iC = 0; iC<NPMTS; ++iC) {

      if (mask & (1 << iC)) {
        histoDraw[iB][iC]->Draw("hist same");
      }
    }
  }
    
  if (showSum) { hsumDraw->Draw("hist same"); }
    
  c1->Modified(); c1->Update();
    
  return true;
    
}

float EstimateTS(int EvNum, float threshold) {
    
  float peak=0, total=0, F90=0;
  int pulse0 = 0;
  int iS;

  bool status = FillHistos(EvNum);
  status &= SumEvents(false);

  if (!status) { return -1.; }

  for (iS = 400; iS < 550; ++iS) {
    if (hsum->GetBinContent(iS) < -threshold) {
      pulse0 = iS;
      break;
    }
  }
    
  if (iS == 550) { return -1; }

  for (int iS = pulse0-5; iS < pulse0 + 23; ++iS) {

    peak += hsum->GetBinContent(iS);
    total += hsum->GetBinContent(iS);
  }
  for (int iS = pulse0 + 23; iS < pulse0+1200; ++iS) total += hsum->GetBinContent(iS);
  

  if(total>0) F90 = peak / total;

  printf("EstimateTS event %i %f \n",EvNum,F90); 

  return F90;

}

int GetT0Bin(TH1F *histogram, float threshold, int startBin, int stopBin) {
    

    TestRange(histogram, startBin, stopBin);
    
    for (int iS = startBin; iS <= stopBin; ++iS) {
        if (histogram->GetBinContent(iS) < -threshold) { return iS; }
    }
    
    return -1;
    
}

int GetMinimum(TH1F *histogram, float &minVal, int startBin, int stopBin, bool StepForward=true) {
    
    int iStart = startBin, iStop = stopBin, step = 1;
    TestRange(histogram, startBin, stopBin);
    
    if (!StepForward) { iStop = startBin; iStart = stopBin; step = -1; }
    
    minVal = ADCrange;
    float val = ADCrange;
    int minBin = -1;
    
    for (int iS = iStart; iS <= iStop; iS += step) {
        val = histogram->GetBinContent(iS);
        if (val < minVal) {
            minVal = val;
            minBin = iS;
        }
    }
    
    return minBin;
    
}


void vis(int iEvent=0, TString theTag= "07-22-1408_0") 
{
  tag=theTag;
   TString inputFileName = TString("pdsData/PDSout_")+tag+TString(".root");
   AddChain(inputFileName);
  // open ouput file and make some histograms
  TString outputFileName = TString("visualizer-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  DrawEvent(iEvent);
  printf(" DrawEvent(i) \n");
}

