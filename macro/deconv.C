// test polya distributions
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TKey.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TVirtualFFT.h"

using namespace TMath;
TString htag;
std::vector<TH1F*> hlist;
std::vector<TH1F*> hdlist;
std::vector<TH1F*> hflist;
std::vector<double> ddigi;
std::vector<double> fdigi;
TH1F *hResponse140;
TH1F *hResponse;
/**** weiner filter *****/
void WFilter(TVirtualFFT *fFFTD, TVirtualFFT *fInverseFFTD, TVirtualFFT *fFFTR, TH1F *hFData) 
{
    int nFFTSize = hFData->GetNbinsX();
    for (int i = 1; i<nFFTSize; ++i) {
       double rl, im;
       fFFTD->GetPointComplex(i,rl,im);
       std::complex<double>  c(rl,im);
       double rrl, rim;
       fFFTR->GetPointComplex(i,rrl,rim);
       std::complex<double>  rc(rrl,rim);
       double fresp = std::abs(rc);
       if(fresp!=0) c /= fresp;
       fInverseFFTD->SetPoint(i,c.real(), c.imag());
    }
    fInverseFFTD->Transform();
    // Set the samples into the calibrated pulse digit.
    for (int ibin =0; ibin<nFFTSize; ++ibin) {
      double v=fInverseFFTD->GetPointReal(ibin)/double(nFFTSize);
      if(isnan(v)) { printf(" i %i v NAN \n",ibin); break; }
<<<<<<< HEAD
      //printf(" i %i v %f \n",ibin, v);
=======
      printf(" i %i v %f \n",ibin, v);
>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1
      hFData->SetBinContent(ibin+1,v);
    }
}
/*******/

/**** from cooper *****/
float Chi2Overlap(TH1F *data, TH1F *response, int startBin, int overlapRange) {
    
    int nbinsD = data->GetNbinsX();
    int nbinsR = response->GetNbinsX();
    
    if (overlapRange > nbinsR) { overlapRange = nbinsR; }
    
    if (startBin > nbinsD) { return 0.; }
    
    int endBin = startBin + overlapRange;
    if (endBin > nbinsD) { endBin = nbinsD; }
    
    float valD, valR, valDR=0., valRR=0.;
    
    for (int i = startBin; i<endBin; ++i) {
        
        valD = data->GetBinContent(i);
        valR = response->GetBinContent(i-startBin+1);
        
        valDR += valD*valR;
        valRR += valR*valR;
        
    }
    
    return valDR/(valRR+0.0000001);
}
/**** from cooper *****/

void AddWave(TH1F *data, TH1F *response, float scale, int startBin) {

    int nbinsD = data->GetNbinsX();
    int nbinsR = response->GetNbinsX();
    
    if (startBin > nbinsD) { return; }
    
    int endBin = startBin + nbinsR;
    if (endBin > nbinsD) { endBin = nbinsD; }
    
    float valD, valR;
    
    for (int i = startBin; i<endBin; ++i) {
        
        valD = data->GetBinContent(i);
        valR = response->GetBinContent(i-startBin+1);
        
        data->SetBinContent(i, (valD+scale*valR));

    }
    
    //data->SetBinContent(startBin, scale);
    
}

/**** from cooper *****/
void Deconvolve(TH1F *data, TH1F *response, TH1F *decon1, float threshold=4, int offset=5, int overlapRange=140) 
{

  TH1F *data2 =  dynamic_cast<TH1F *>(data->Clone("data2"));
  TH1F *decon2 = dynamic_cast<TH1F *>(data->Clone("decon2"));

    
  int nbinsD = data->GetNbinsX();
  data2->SetBins(nbinsD, 0., (float)nbinsD);

  decon1->Reset();
  decon1->SetBins(nbinsD, 0., (float)nbinsD);
  decon2->Reset();
  decon2->SetBins(nbinsD, 0., (float)nbinsD);
  for (int i=1;i<=nbinsD;++i){ data2->SetBinContent(i, data->GetBinContent(i)); }

  float over;
  float maxval = 0.;

  for (int i=1;i<=nbinsD;++i) {
    over = Chi2Overlap(data2, response, i, overlapRange);
    decon2->SetBinContent(i, over);


    // real bump and is rising
    if ( (over > threshold) && (over > maxval) ) { maxval = over; }

    // real bump in maxval and now falling; previous bin is peak
    if ( (maxval > threshold) && (over < maxval) ) {

      AddWave(data2, response, -maxval, i-1);

      decon1->SetBinContent(i-1, -maxval);
      maxval = 0;

      // want to jump ahead to get "down off the bump"
      i += offset;
    }
  }
    
}




/*****************************************************
   read all histograms 
*******************************************************/

// neils filter
std::vector<Double_t> MovingAverageFilter(std::vector<Double_t> signal,Int_t aveN)
{
  Int_t N = aveN;
  if(aveN%2==0) ++N; 
  std::vector<Double_t> filter;
  Int_t N2 = std::floor(N/2);
  for(int i = N2; i < Int_t(signal.size())-N2; i++){
    Double_t sum = 0;
    for(int j = i-N2; j <= i+N2; j++){
      sum += signal[j];
    }
    sum /= N;
    filter.push_back(sum);
  }
  
  for(int i = 0; i < N2 ; i++){
    std::vector<Double_t>::iterator it;
    it = filter.begin();
    filter.insert(it,0.);
    filter.push_back(0);
  }
  return filter;
}


void reading(TDirectory *fdir) 
{
  TObject *obj=NULL;
  TKey* key;
  TIter nextkey(fdir->GetListOfKeys());
  htag = TString("QhitNoBeam");
  //htag = TString("QhitOff");
  while ( (key = (TKey*)nextkey()) )  {
    fdir->cd();
    obj = key->ReadObj();
    if(obj->IsA()->InheritsFrom("TH1F")) { //case of TH1 or TProfile
      TH1F* h1 = dynamic_cast<TH1F*>(obj);
      if(TString(h1->GetName()).Contains("response")) hResponse140=h1;
      if(TString(h1->GetName()).Contains("ev")) hlist.push_back(h1);
    } 
  }
}

void deconv(TString tag="led-pulse-events-cooper")
{
  enum {NPMT=21};
  TString inputFileName = tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
<<<<<<< HEAD

  
  //infile->ls();
  reading(infile);
  for(UInt_t ih=0; ih<hlist.size(); ++ih) printf(" %i %s \n",ih,hlist[ih]->GetName());
  hResponse140->Smooth(5);

  TString outputFileName = tag+TString("-deconv.root");
  printf(" opening file %s \n",outputFileName.Data()); 
  TFile *outfile = new TFile(outputFileName,"recreate");

  
=======
  //infile->ls();
  reading(infile);
  for(UInt_t ih=0; ih<hlist.size(); ++ih) printf(" %i %s \n",ih,hlist[ih]->GetName());

  

>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1
  // for deconvolution, make clones
  for(UInt_t ih=0; ih<hlist.size(); ++ih) hdlist.push_back(dynamic_cast<TH1F *>(hlist[ih]->Clone(Form("decEv%ih",ih))));
   // for filtesr deconvolution, make clones
  for(UInt_t ih=0; ih<hlist.size(); ++ih) hflist.push_back(dynamic_cast<TH1F *>(hlist[ih]->Clone(Form("fdecEv%ih",ih))));
<<<<<<< HEAD
  hResponse = dynamic_cast<TH1F *>(hlist[0]->Clone("fResponse"));
  hResponse->Clear();
  hResponse->SetName("fReponse");
  hResponse->SetTitle(" full response function");
  hResponse->Sumw2(false);
  double hnorm =0;
  for(int ibin=0; ibin<=hResponse->GetNbinsX(); ++ibin) {
    hResponse->SetBinContent(ibin,0);
    hResponse->SetBinError(ibin,0);
  }
  for(int ibin=1; ibin<=hResponse->GetNbinsX(); ++ibin) {
    if(ibin<hResponse140->GetNbinsX()) hResponse->SetBinContent(ibin-2, hResponse140->GetBinContent(ibin));
    //if(ibin<60) hResponse->SetBinContent(ibin-2, hResponse140->GetBinContent(ibin));
  }
  printf(" \n\n\t response integral %f \n\n",abs(hResponse->Integral()));
  hnorm = abs(hResponse->Integral());
  for(int ibin=1; ibin<=hResponse->GetNbinsX(); ++ibin)  hResponse->SetBinContent(ibin, hResponse->GetBinContent(ibin)/hnorm);
  printf(" \n\n\t response integral %f \n\n",abs(hResponse->Integral()));

  outfile->Add(hResponse140);
  outfile->Add(hResponse);
  //hResponse140->Print("all");
  //hResponse->Print("all");
=======
  hResponse = dynamic_cast<TH1F *>(hlist[0]->Clone("fullResponse"));
  hResponse->Clear();
  double hnorm =0;
  for(int ibin=1; ibin<=hResponse->GetNbinsX(); ++ibin) {
    if(ibin<hResponse140->GetNbinsX()) hResponse->SetBinContent(ibin, hResponse140->GetBinContent(ibin));
    else hResponse->SetBinContent(ibin,0);
    hResponse->SetBinError(ibin,0); 
    hnorm += hResponse->GetBinContent(ibin);
  }
  for(int ibin=1; ibin<=hResponse->GetNbinsX(); ++ibin)  hResponse->SetBinContent(ibin, hResponse->GetBinContent(ibin)/hnorm);

>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1
 
  TString canTitle;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  TCanvas* canr;
  printf(" \n *********** pmt response  %s ************ \n",tag.Data());
  canTitle.Form("%s-response",tag.Data());
  canr = new TCanvas(canTitle,canTitle);
  hResponse->SetLineColor(kBlack);
  hResponse->Draw();
  //canr->Print(".pdf");

  // histogram for response function
  for(UInt_t ih=0; ih<hdlist.size(); ++ih){ 
    Deconvolve(hlist[ih], hResponse, hdlist[ih]);
    TCanvas* candecon;
    printf(" \n *********** pmt response  %s ************ \n",tag.Data());
    canTitle.Form("%s-response-ev-%i",tag.Data(),ih);
    candecon= new TCanvas(canTitle,canTitle);
    hlist[ih]->SetLineColor(kBlack);
    hlist[ih]->SetMarkerStyle(3);
    hlist[ih]->SetMarkerColor(kBlack);
    hlist[ih]->SetMarkerSize(0.5);
    hlist[ih]->GetXaxis()->SetRangeUser(400,600);	
    hlist[ih]->Draw();
<<<<<<< HEAD


    for(int ibin=1; ibin<=hdlist[ih]->GetNbinsX(); ++ibin)  hdlist[ih]->SetBinError(ibin, sqrt( hdlist[ih]->GetBinContent(ibin)));
=======
>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1
    
    hdlist[ih]->SetLineColor(kRed);
    hdlist[ih]->SetMarkerStyle(4);
    hdlist[ih]->SetMarkerSize(0.5);
    hdlist[ih]->SetMarkerColor(kRed);
    hdlist[ih]->GetXaxis()->SetRangeUser(400,600);	
    hdlist[ih]->Draw("same");
    candecon->Print(".pdf");
  }


  
  // smoothed
  int MAXSAMPLES = hResponse->GetNbinsX();
  for(int ibin=1; ibin<=MAXSAMPLES ; ++ibin) ddigi.push_back( double(hResponse->GetBinContent(ibin)));
  fdigi = MovingAverageFilter(ddigi,1); // parameter is top hat window which should be odd

  // makei fft plots 
  int ientry=1;
  int ipmt=0;
  TString hname,htitle;
  hname.Form("raw_pmt%i_ev%i",ipmt,Int_t(ientry));
  htitle.Form("raw samples pmt %i ev %i",ipmt,Int_t(ientry));
  TH1D* hDigi = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
  hname.Form("filtered_pmt%i_ev%i",ipmt,Int_t(ientry));
  htitle.Form("filtered digi samples pmt %i ev %i",ipmt,Int_t(ientry));
  TH1D* hFDigi = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
  for(unsigned idigi=0; idigi<ddigi.size(); ++idigi) {
    hDigi->SetBinContent(int(idigi+1),ddigi[idigi]);
    hFDigi->SetBinContent(int(idigi+1),fdigi[idigi]);
  }

  TCanvas* canFilt;
  canTitle.Form("%s-filtered-response",tag.Data());
  canFilt = new TCanvas(canTitle,canTitle);
  hDigi->SetLineColor(kBlack);
  hDigi->Draw();
  hFDigi->SetLineColor(kRed);
  hFDigi->Draw("same");

  //canr->Print(".pdf");


  // FFT initialize 
  Int_t nFFTSize = int(fdigi.size());
  TVirtualFFT *fFFTR = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
  //TVirtualFFT *fInverseFFTR = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");
  for(int is =0; is<nFFTSize; ++is) {
    fFFTR->SetPoint(is,fdigi[is]);
  }

  fFFTR->Transform();
  TH1D* hfft = new TH1D("FFTResponse","FFT response",nFFTSize/2,0,nFFTSize/2);

  // fill samples FFT histogram && elec response in time domain
  printf(" created %s %s  size = %i \n",hfft->GetName(),hfft->GetTitle(), nFFTSize );
  // skip first bin which is pedestal
  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFTR->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im);
    hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
  } 
 
  gStyle->SetOptStat(0);
  TCanvas* canfft;
  printf(" \n *********** pmt response fft  %s ************ \n",tag.Data());
  canTitle.Form("%s-response-FFT",tag.Data());
  canfft = new TCanvas(canTitle,canTitle);
  hfft->SetLineColor(kBlack);
  hfft->Draw();
  //canr->Print(".pdf");
  //
  // FFT of signals loop over all signal histograms
 

  TH1D* hfftSig[MAXSAMPLES];
  for(UInt_t ih=0; ih<hflist.size(); ++ih) {
     TVirtualFFT *fFFTD = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
     TVirtualFFT *fInverseFFTD = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");

    // fill ddigi from histogram
    ddigi.clear();
    for(int ibin=1; ibin<=MAXSAMPLES ; ++ibin) ddigi.push_back( double(hlist[ih]->GetBinContent(ibin)));
     for(int is =0; is<nFFTSize; ++is) {
       fFFTD->SetPoint(is,ddigi[is]);
     }
     fFFTD->Transform();
     hfftSig[ih] = new TH1D(Form("FFTPmt%i",ih),Form("FFT PMT %i",ih),nFFTSize/2,0,nFFTSize/2);
     // skip first bin which is pedestal
     for (int i = 1; i<nFFTSize/2; ++i) {
       double rl, im;
       fFFTD->GetPointComplex(i,rl,im);
       std::complex<double>  c(rl,im);
       hfftSig[ih]->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
     } 

     printf(" fft size fftD %i fftR %i fftInv %i %i \n",*fFFTD->GetN(),  *fFFTR->GetN(), *fInverseFFTD->GetN(),nFFTSize); 
     WFilter(fFFTD,fInverseFFTD,fFFTR,hflist[ih]);
<<<<<<< HEAD
=======
     delete fFFTD;
     delete fInverseFFTD;
>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1

  }

  
  // histogram for response function
  for(UInt_t ih=0; ih<hdlist.size(); ++ih){ 
    printf(" \n *********** pmt fft filter   %s ************ \n",tag.Data());
    canTitle.Form("%s-filter-ev-%i",tag.Data(),ih);
    TCanvas* candecon= new TCanvas(canTitle,canTitle);
    hlist[ih]->SetLineColor(kBlack);
    hlist[ih]->SetMarkerStyle(3);
    hlist[ih]->SetMarkerColor(kBlack);
    hlist[ih]->SetMarkerSize(0.5);
    hlist[ih]->GetXaxis()->SetRangeUser(400,600);	
    hlist[ih]->Draw();
    
    hflist[ih]->SetLineColor(kRed);
<<<<<<< HEAD
    hflist[ih]->SetMarkerStyle(21);
=======
    hflist[ih]->SetMarkerStyle(4);
>>>>>>> 8da87108be6515f302dd9761483619c348f6c5c1
    hflist[ih]->SetMarkerSize(0.5);
    hflist[ih]->SetMarkerColor(kRed);
    hflist[ih]->GetXaxis()->SetRangeUser(400,600);	
    hflist[ih]->Draw("same");
    candecon->Print(".pdf");
  }


}

