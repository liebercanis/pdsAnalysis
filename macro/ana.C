#include <vector>
#include "TPmtEvent.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
  TH1D *hOcc;
  TH1D *hNoise;
  TH1D *hBase;
  TH1D* hSamples[NPMT];
  TH1D* hPeaks[NPMT];
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hNHits[NPMT];
  TH1D* hQMax[NPMT];
  TH1D* hCounts[NPMT];
  TH1D* hBaseline[NPMT];

  
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



void getFileHistograms(TFile *infile)
{
   // get histograms from file
  hOcc = (TH1D*)infile->Get("occupancy");
  hNoise = (TH1D*)infile->Get("noise");
  hBase= (TH1D*)infile->Get("base");
 
  TString hname;
  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("Samples_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hSamples[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("Peaks_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hPeaks[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("Counts_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hCounts[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("Baseline_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hBaseline[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("hitQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hHitQ[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("QMax_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hQMax[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("NHits_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hNHits[ipmt] = (TH1D*)infile->Get(hname);
      hname.Form("FFTPmt%i",ipmt);
      hFFT[ipmt] = (TH1D*)infile->Get(hname);
    }
  }

}

void canPlots(TString tag)
{
    TString canname;
    enum {NCAN=7};
    TCanvas *can1[NCAN];
    TCanvas *can2[NCAN];
    TCanvas *can3[NCAN];
    TCanvas *can4[NCAN];
    TCanvas *can5[NCAN];

    int ican=-1;
    int ip=0;
    for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
      if(ipmt%3==0) {
        ip=0;
        ++ican;
        canname.Form("FFT-set%i-run-%s",ican,tag.Data());
        can1[ican] = new TCanvas(canname,canname);
        can1[ican]->Divide(1,3);
        canname.Form("counts-set%i-run-%s",ican,tag.Data());
        can2[ican] = new TCanvas(canname,canname);
        can2[ican]->Divide(1,3);
        canname.Form("samples-set%i-run-%s",ican,tag.Data());
        can3[ican] = new TCanvas(canname,canname);
        can3[ican]->Divide(1,3);
        canname.Form("hitCharge-set%i-run-%s",ican,tag.Data());
        can4[ican] = new TCanvas(canname,canname);
        can4[ican]->Divide(1,3);
        canname.Form("qMax-set%i-run-%s",ican,tag.Data());
        can5[ican] = new TCanvas(canname,canname);
        can5[ican]->Divide(1,3);
      }
      can1[ican]->cd(ip+1); hFFT[ipmt]->Draw();
      can4[ican]->cd(ip+1); gPad->SetLogy(); hHitQ[ipmt]->Draw();
      can5[ican]->cd(ip+1); gPad->SetLogy(); hQMax[ipmt]->Draw();
      can3[ican]->cd(ip+1); 
      hPeaks[ipmt]->Draw();
      hSamples[ipmt]->Draw("sames");
      can2[ican]->cd(ip+1);  gPad->SetLogy(); hCounts[ipmt]->Draw();
      ++ip;
    }

    for(int ican=0; ican<NCAN; ++ican) {
      can1[ican]->Print(".pdf");
      can2[ican]->Print(".pdf");
      can3[ican]->Print(".pdf");
      can4[ican]->Print(".pdf");
      can5[ican]->Print(".pdf");
    }
}

void ana(TString tag= "07-21-1740_0")
{
  TString inputFileName = TString("../pdsOutput/pmtAna_")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *pmtTree=NULL;

  // tree has to be in file
  Long64_t aSize=0;
  pmtTree = (TTree*) infile->Get("pmtTree");
  if(pmtTree) aSize=pmtTree->GetEntriesFast();
  printf(" pmtTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("ana-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  
  getFileHistograms(infile);
  canPlots(tag);

  
  TPmtEvent *ev = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&ev);

  for(unsigned entry =0; entry < 200; ++entry ) {
    pmtTree->GetEntry(entry);
    if(entry%100==0) printf("\t entry %i hits %i \n",entry,ev->nhits);
  }

  // end of ana 
  outfile->Write();
}
