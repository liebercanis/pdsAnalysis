#include <vector>
#include "TPmtEvent.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
double off;
double spacer;
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

  TH1D* hQHitCut[NPMT];
  TH1D* hQHitLength[NPMT];
  //TH1D* hQHitTime[NPMT];
  
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
      hname.Form("Samples_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hSamples[ipmt] = (TH1D*)infile->Get(hname);
      //if(hSamples[ipmt]) cout << " got " << hSamples[ipmt]->GetName() << endl;
      //else cout << " cannot find hSamples pmt " << ipmt << endl;
      if(ipmt<0||ipmt>=NPMT) continue;
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

void setHistStyle(TH1D* h, int color=kBlack, bool stat=false)
{
  if(!h) return;
  //h->SetLineWidth(0.5);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(1);
  h->SetStats(stat);
  h->SetTitleOffset(1,"Y");
  h->GetYaxis()->SetLabelSize(0.02);
  h->SetAxisRange(-4*spacer,spacer*NB*NC+4*spacer,"Y");
}



TH1D *histOffset(TH1D* h, float off)
{
  if(!h) { cout << " Oy, no hist " << endl; return NULL;} 
  TString hname;
  hname.Form("Plot-%s",h->GetName());
  TH1D* h1 = (TH1D*) h->Clone(hname);
  // baseline 
  double base=0;
  for(int ibin=1; ibin<= h->GetNbinsX(); ++ibin) base+= h->GetBinContent(ibin);
  base/=double(h->GetNbinsX());

  for(int ibin=1; ibin<= h->GetNbinsX(); ++ibin) {
    h1->SetBinContent(ibin+1,h->GetBinContent(ibin)+off-base);  // note that histo bins start with 1
  }
  return h1;
}

void plotSummary(int iEvent=0, TString tag= "07-22-1408_0")
{
  spacer=10.;
 
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
  TString outputFileName = TString("eventPlot-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());



  
  getFileHistograms(infile);
  
  TPmtEvent *ev = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&ev);

  // 8*3 colors
  int ci1 = TColor::GetColor("#73bac6");
  int ci2 = TColor::GetColor("#ef617b");
  int ci3 = TColor::GetColor("#ef8210");
  int colors[NC*NB]={1,2,5,3,4,6,7,ci1,1,2,5,3,4,6,7,ci2,1,2,5,3,4,6,7,ci3};

  TString eventName;
  eventName.Form("file%s-event-%d",tag.Data(),iEvent);
  TCanvas *canA = new TCanvas(eventName,eventName);

 /* 
  TPad *pad[NB*NC];
  TString padName;
   for(int ib=0; ib<NB; ++ib) {    
    for(int ic=0; ic<NC; ++ic) {
      padName.Form("bd%i-ch%i",ib,ic);
      pad[ib*ic] = new TPad(padName,padName,0.02,0.02,0.48,0.83,kWhite);
      pad[ib*ic]->Draw();
    }
   }
*/

  //canA->Range(0,0,25,18);
  //canA->SetFillColor(17);
  
 
  TString histTitle;
  TString opt;
  bool first;
  off=spacer+spacer*double(NB*NC);

  for(int ib=0; ib<NB; ++ib) {    
    for(int ic=NC-1; ic>=0; --ic) {
      if(first) {
        opt=TString("");
        first=false;
      } else opt=TString("SAME");
      off -= spacer;
      int ipmt = toPmtNumber(ib,ic);

      TH1D *hist = hSamples[ipmt];
      if(!hist) {
        printf(" no hist board %i ch %i \n",ib,ic); 
        continue;
      }
      else printf("%s offset %f \n",hist->GetName(),off);
      cout << hist->GetName() << endl;
      if(ic==NC-1) setHistStyle(hist,kRed);
      else setHistStyle(hist,kBlack);
      TH1D* h1 = histOffset(hist,off);
      histTitle.Form("%s-%i",tag.Data(),iEvent); 
      h1->SetTitle(histTitle);
      h1->SetXTitle("4ns samples");
      h1->SetYTitle(" ADC couts (arb offset)  ");
      h1->Draw(opt);
    }
  }
  canA->Print(".pdf");
}

