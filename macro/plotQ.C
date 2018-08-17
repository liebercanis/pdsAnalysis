#include <vector>
#include "TPmtEvent.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};

TH1D* hQHit[NPMT];
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
 
  TString hname;
  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      hname.Form("QhitCut_b%u_ch%u_pmt%u",ib,ic,ipmt);
      hQHit[ipmt] = (TH1D*)infile->Get(hname);
      if(hQHit[ipmt]) cout << " got " << hQHit[ipmt]->GetName() << endl;
      else cout << " cannot find hSamples pmt " << ipmt << endl;
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
}



void plotQ(int iEvent=0, TString tag= "07-31-1555_0")
{
  TString inputFileName = TString("ana-")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  getFileHistograms(infile);

  // loop to find peak below 200 
  Int_t down = 0;
  Int_t izero = hQHit[0]->FindBin(5);
  Int_t istart = hQHit[0]->FindBin(100);
  Int_t ilast =  hQHit[0]->GetNbinsX();
  printf("starting bin is %i \n",istart);
 
  double pmtGain[NPMT];
  for(int ipmt=0; ipmt<NPMT; ++ipmt ) {
    pmtGain[ipmt]=0;
    // find first peak
    Int_t ibin=0;
    Int_t pbin=0;
    double peak=0;
    for(ibin=izero; ibin <  ilast; ++ibin) {
      double ni = hQHit[ipmt]->GetBinContent(ibin);
      if(ni>peak) {
        peak=ni; pbin=ibin;
      } 
    }
    // find first dip
    Int_t iend = pbin;
    double nmin=peak;
    int up=0;
    for(ibin=pbin; ibin <  ilast; ++ibin) {
      double ni = hQHit[ipmt]->GetBinContent(ibin);
      if(ni<nmin) {
        nmin=ni; iend=ibin; up=0;
      } else ++up;
      if(up>2) break;
    }
    printf(" pmt %i peak ADC %f peak num  %f \n min ADC %f nmin %f  \n",ipmt,hQHit[ipmt]->GetBinLowEdge(pbin),peak,hQHit[ipmt]->GetBinLowEdge(iend),nmin );

    //printf(" \t pmt  %i istart %i  qstart %f \n",ipmt,istart,hQHit[ipmt]->GetBinContent(istart));
    double nmax=0; 
    int mbin=iend;
    down=0;
    for(ibin=istart; ibin>iend ; --ibin) {
      double ni = hQHit[ipmt]->GetBinContent(ibin);
      if(ni>nmax) {
        nmax=ni; mbin=ibin; down=0; 
      } else ++down;
      if(down>5) break;
    }
    pmtGain[ipmt]=hQHit[ipmt]->GetBinLowEdge(mbin);
    printf(" \t pmt %i max bin %i max ADC %f nmax %f down %i \n",ipmt,mbin,hQHit[ipmt]->GetBinLowEdge(mbin),nmax,down);
  }

  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" pmt %i gain %0.1f \n",ipmt,pmtGain[ipmt]);
}
