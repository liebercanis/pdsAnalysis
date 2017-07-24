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

#include "TPmtEvent.hxx"

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

typedef std::complex<double> Complex;


class pmtAna {
public :
  enum {MAXSAMPLES=4200};
  enum {NB=3,NC=16,NS=MAXSAMPLES};
  enum {NPMT=NB*NC/2};

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
   UShort_t        digitizer_waveforms[NB][NC][NS];
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

   pmtAna(TString tag="07-12-1900_0",Int_t maxLoop=0);
   virtual ~pmtAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   UInt_t Loop(UInt_t nToLoop=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TTree* pmtTree;
   TPmtEvent* pmtEvent;
   
   std::vector<Int_t> findMaxPeak(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
   std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
   Int_t findHits(Int_t ipmt, Double_t sum,  std::vector<Int_t> peakTime, std::vector<Double_t> ddigi); 

   double getBaseline(int ipmt ) { return hBase->GetBinContent(ipmt+1); } 

   // returns -1 if pmt does not exist 
   int toPmtNumber(int ib, int ic) 
   {
     int ipmt;
     if(ic%2!=0) ipmt=-1;
     else ipmt = ib*NC/2+ ic/2;
     return ipmt;
   }

   // valid pmt are 0 to NPMT-1
   void fromPmtNumber(int ipmt, int& ib, int&ic)
   {
     ib=-1; ic=-1;
     if(ipmt>=NPMT&&ipmt<0) return;
     int jpmt = 2*ipmt;
     ib=(jpmt-jpmt%NC)/NC;
     ic= jpmt - NC*ib;
   }

   Double_t baselineNominal[NPMT];
   TNtuple *ntPmt;
   TNtuple *ntDigi;
   TNtuple *ntHit;
   //FFT
   Int_t nFFTSize;
   TH1D* FFTFilter(Int_t pmt); 
   /// The fft class to take the fourier transform.
   TVirtualFFT *fFFT;
   /// The fft class to take the inverse fourier transform.
   TVirtualFFT *fInverseFFT;


   // histogram pointers
   TH1D* hSamples[NPMT];
   TH1D* hPeaks[NPMT];
   TH1D* hFFT[NPMT];
   TH1D* hHitQ[NPMT];
   TH1D* hQMax[NPMT];

   TH1D* hCounts[NPMT];
   TH1D* hBaseline[NPMT];
   TH1D* hOcc;
   TH1D* hNoise;
   TH1D* hBase;   
};
