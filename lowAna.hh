//////////////////////////////////////////////////////////
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
#include <complex>
#include <valarray>
#include <string>

#include "TPmtEvent.hxx"
#include "TPmtGains.hxx"
#include "TPmtSummary.hxx"
#include "readClock.hxx"
#include <TROOT.h>
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

class lowAna {
public :
  enum {MAXSAMPLES=2100};
  enum {NB=3,NCPMT=7,NC=NCPMT+1};
  enum {NALLCH=NB*NC};
  enum {NPMT=NCPMT*NB};
  enum {MAXADC=4095};
  enum {THRESHOLDHIGH=1,THRESHOLDLOW=0};
  //peak finding
  enum {minLength=3,maxHalfLength=100};
  const double L=23.2;//m
  const double clight=0.299792458;//m/ns
  const double nmass=939.565;//MeV
  const Double_t GAMMAPEAK=-628.089;//ysun
  const UInt_t gate = 4500000; // PDS gate 
  const float toMicro = 1.0E-3;
 
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
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

   lowAna(Int_t maxLoop=0,Int_t firstEntry=0);
   virtual ~lowAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   UInt_t Loop(UInt_t nToLoop=0, UInt_t firstEntry=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // Top Hat filter algorithm
   // parameter is top hat window which should be odd 
   std::vector<Double_t> MovingAverageFilter(std::vector<Double_t> signal,Int_t aveN);
   // nearest RF time to hit peak time
   void getTimeToRF(UInt_t board); 
   
   TTree* pmtTree;
   TPmtEvent* pmtEvent;
   TTree* summaryTree;
   TPmtSummary *pmtSummary;
   TFile* outFile;
   TPmtGains *pmtGains;   // new gains from this run
   TPmtGains *goodGains; // the ones used
   TFile *gainFile;
   readClock *rclock;
   
   // get trigger type 
   Int_t triggerInfo();
   // trigger type counters
   
   std::vector<Int_t> findMaxPeak(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
   std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
   Int_t findHits(Int_t ipmt, Double_t sum,  std::vector<Int_t> peakTime, std::vector<Double_t> ddigi, std::vector<Double_t> ddigiUn, Int_t type); 
   Int_t readGainConstants(TString fileName="gainConstants.txt"); // returns number of gains read
   double getBaseline(int ipmt ) { return hBase->GetBinContent(ipmt+1); }
   std::vector<Int_t> findRFTimes(int ipmt,double& digiMin);
   void ADCFilter(int iB, int iC);
   void qualitySummary();
   Double_t getPromptTime();
   Double_t getPromptTimeToRF();
   Double_t getTimeToRF();


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

   void integralHist(TH1F *hist, TH1F* ihist) {
     int nbins = hist->GetNbinsX();
     Float_t val=0;
     for(int i=0; i<nbins; ++i) {
       val += hist->GetBinContent(i);
       ihist->SetBinContent(i,val);
     }
   }

   Float_t getTime(TH1F *ihist,int& maxbin, Float_t& delta) {
     int nbins = ihist->GetNbinsX();
     maxbin=0;
     delta=0;
     for(int i=1; i<nbins; ++i) {
       Float_t d = ihist->GetBinContent(i) - ihist->GetBinContent(i-1);
       if(d>delta) {
         delta=d;
         maxbin=i;
       }
     }
     Float_t time = ihist->GetXaxis()->GetBinLowEdge(maxbin);
     return time;
   }

   // by run addative alignmentes
   Long64_t addAlign1, addAlign2;
   Double_t baselineNominal[NPMT];
   Double_t gain[NPMT]; //Q = ADC/gain
   Double_t egain[NPMT]; //Q = ADC/gain
   TNtuple *ntTrig;
   //FFT
   Int_t nFFTSize;
   TH1D* FFTFilter(Int_t pmt); 
   /// The fft class to take the fourier transform.
   TVirtualFFT *fFFT;
   /// The fft class to take the inverse fourier transform.
   TVirtualFFT *fInverseFFT;

   std::vector<Int_t> rftime21;
   std::vector<Int_t> rftime22;
   std::vector<Int_t> rftime23;

   // histogram pointers
   TH1D* hSamples[NALLCH];  // include RF
   TH1D* hSamplesSum;

   TH1D* hTPrompt;
   TH1D* hTPromptEvent;
   TDirectory *promptDir;
   TDirectory *histDir;
   TDirectory *trigTimeDir;

   TH1D* hSamplesPDS[NPMT];  // include RF
   TH1D* hSamplesPDSSum;

   TH1D* hPeaks[NPMT];
   TH1D* hFFT[NPMT];
   TH1D* hHitQ[NPMT];
   TH1D* hNHits[NPMT];
   TH1D* hQMax[NPMT];
   TH1D* hRawQ[NPMT];
   
   TH1D* hCounts[NPMT];
   TH1D* hBaseline[NPMT];
   TH1D* hOcc;
   TH1D* hNoise;
   TH1D* hBase; 

  TH1F* hQTime[NB];
  TH1F*  hQTimeIntegral[NB];
};
