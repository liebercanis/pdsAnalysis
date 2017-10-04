#include <vector>
#include "TPmtEvent.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
  TH1D *hTimeToRF;
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

  TH1D* hQinHit; 
  TH1D* hQHitOn[NPMT];
  TH1D* hQHitOff[NPMT];
  TH1D* hQHitNoBeam[NPMT];
  TH1D* hQHitLength[NPMT];
  TNtuple *ntHit;
  //TH1D* hQHitTime[NPMT];
  TPmtEvent *ev;
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
      can2[ican]->cd(ip+1);  gPad->SetLogy(); if(hCounts[ipmt]) hCounts[ipmt]->Draw();
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


void newCanPlots(TString tag)
{
    TString canname;
    enum {NCAN=7};
    TCanvas *can1[NCAN];
    TCanvas *can2[NCAN];

    int ican=-1;
    int ip=0;
    for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
      if(ipmt%3==0) {
        ip=0;
        ++ican;
        canname.Form("Qhit-set%i-run-%s",ican,tag.Data());
        can1[ican] = new TCanvas(canname,canname);
        can1[ican]->Divide(1,3);
        canname.Form("QHitNoBeam-set%i-run-%s",ican,tag.Data());
        can2[ican] = new TCanvas(canname,canname);
        can2[ican]->Divide(1,3);

      }
      can1[ican]->cd(ip+1); gPad->SetLogy(); 
      hQHitOn[ipmt]->SetNormFactor();
      hQHitOff[ipmt]->SetNormFactor();
      hQHitOn[ipmt]->Draw();
      hQHitOff[ipmt]->SetLineColor(kRed);
      hQHitOff[ipmt]->Draw("sames");
      can2[ican]->cd(ip+1); gPad->SetLogy(); hQHitNoBeam[ipmt]->Draw();
      ++ip;
    }

    for(int ican=0; ican<NCAN; ++ican) {
      can1[ican]->Print(".pdf");
      can2[ican]->Print(".pdf");
    }
}



void ana(TString tag="07-31-1555_0")
    //"07-26-0914_0") // "07-31-1555_0") //
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
  hTimeToRF = new TH1D("TimeToRF"," hit time to RF time",MAXSAMPLES/2,0,MAXSAMPLES+1);
  ntHit = new TNtuple("ntAnaHit", " ana hits ","ipmt:rft:peakt:length:ratio:q:qp");
  hQinHit = new TH1D("QinHit"," ADC counts of bins in hit ",2000,0,20);
  outfile->ls();
    
  TString hname;
  TString htitle;

  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      
      hname.Form("QhitOff_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit between RF board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitOff[ipmt] = new TH1D(hname,htitle,200,0,200);
      hQHitOff[ipmt]->SetXTitle(" ADC counts ");
      hQHitOff[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

      hname.Form("QhitNoBeam_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit between RF board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitNoBeam[ipmt] = new TH1D(hname,htitle,200,0,200);
      hQHitNoBeam[ipmt]->SetXTitle(" ADC counts ");
      hQHitNoBeam[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

   
      hname.Form("QhitOn_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit at RF board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitOn[ipmt] = new TH1D(hname,htitle,200,0,200);
      hQHitOn[ipmt]->SetXTitle(" ADC counts ");
      hQHitOn[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

      hname.Form("QhitLength_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit/length board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitLength[ipmt] = new TH1D(hname,htitle,100,0,100);
      hQHitLength[ipmt]->SetXTitle(" qhit/length (ADC counts) ");
      hQHitLength[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

    }
  }
  
  getFileHistograms(infile);
  //canPlots(tag);

  
  ev = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&ev);

  std::vector<TPmtHit>  hit;
  for(unsigned entry =0; entry < aSize; ++entry ) {
    pmtTree->GetEntry(entry);
    hit.clear();
    hit = ev->hit;
    if(entry%1000==0) printf("...entry %i hits %lu trig type %i (%lu,%lu,%lu) \n",entry,hit.size(), ev->trigType,
        ev->rft21.size(),ev->rft23.size(),ev->rft23.size());
    for(int ihit =0; ihit < hit.size(); ++ihit) {
      TPmtHit* phit = &(ev->hit[ihit]);
      for(int is=0; is<phit->nsamples; ++is) hQinHit->Fill(phit->qsample[is]);
      Int_t timeToRF = phit->timeToRF; 
      if(ev->trigType==TPmtEvent::TRIG111||ev->trigType==TPmtEvent::TRIG444||ev->trigType==TPmtEvent::TRIG555) 
        hTimeToRF->Fill(timeToRF);
      //if(timeToRF<MAXSAMPLES) printf(" trig %i hit %i pmt %i time %i \n",ev->trigType, ihit, phit->ipmt,timeToRF);
      int length = TMath::Abs(phit->tstop-phit->tstart)+1;
      //printf(" \t %i %i ipmt %i length %i qhit %f \n",entry,ihit,phit->ipmt,length,phit->qhit);
      // cut on time since RF pulse
      if(timeToRF<100) {
        hQHitOn[phit->ipmt]->Fill(phit->qhit);
        hQHitLength[phit->ipmt]->Fill(phit->qhit/double(length));
      } else if(timeToRF>500) {// off RF pulse 
        hQHitOff[phit->ipmt]->Fill(phit->qhit);  
      }

      if(ev->trigType==TPmtEvent::TRIG000)  hQHitNoBeam[phit->ipmt]->Fill(phit->qhit);  

      
      ntHit->Fill(phit->ipmt,phit->timeToRF, phit->peakTime, length, phit->ratio, phit->qhit, phit->qpeak);
    }
  }

  // end of ana 
  newCanPlots(tag);
  outfile->Write();

  for(int ipmt=0; ipmt<NPMT; ++ ipmt) printf(" mean[%i]=%f ; \n",ipmt, hQHitNoBeam[ipmt]->GetMean());
}
