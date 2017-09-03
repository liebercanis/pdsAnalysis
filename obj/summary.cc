////////////////////////////////////////////////////////
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

TTree *pmt_tree;
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
enum {NALLCH=NB*NC};
enum {MAXADC=4095};
UShort_t digitizer_waveforms[NB][NC][NS];
TTree *summaryTree;
TPmtSummary* pmtSummary;

// store info for summary

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

std::vector<Int_t> findRFTimes(int ipmt, double& step) 
{
  std::vector<Int_t> rftimes;
  int ib; int ic;
  fromPmtNumber(ipmt,ib,ic);

  // find baseline
  std::vector<UShort_t> udigi; 
  for (UInt_t is=0; is<NS; ++is) udigi.push_back(digitizer_waveforms[ib][ic][is]);
  std::sort(udigi.begin(), udigi.end());
  UShort_t baseline = udigi[0.5*double(NS)];
  
  // looking for negative values.  
  UShort_t digiMin=MAXADC;
  for (UInt_t is=0; is<NS; ++is) {
    digitizer_waveforms[ib][ic][is]=TMath::Min( baseline , digitizer_waveforms[ib][ic][is]);
    if(digitizer_waveforms[ib][ic][is]<digiMin) digiMin=digitizer_waveforms[ib][ic][is];
  }
  
  step = double(digiMin) - double(baseline);
  // return if step down is too small
  if(step>-500) return rftimes;
  // pick off start of rising edge
  bool isRF=false;
  for (UInt_t is=0; is<NS; ++is){
    double digi = double(digitizer_waveforms[ib][ic][is]) - double(baseline);
    if(digi<0.75*step&&!isRF) {
      rftimes.push_back(is);
      isRF=true;
    } else if(digi>0.75*step) {
      isRF=false;
    }
  }
  return rftimes;
}


// trigger information
Int_t triggerInfo()
{
  Int_t type = TPmtEvent::TRIGUNKNOWN; // unknosn 
  // RF channels 
  double s1,s2,s3;

  std::vector<Int_t> rftime21 = findRFTimes(21,s1);
  std::vector<Int_t> rftime22 = findRFTimes(22,s2);
  std::vector<Int_t>  rftime23 = findRFTimes(23,s3);
  //UInt_t totalTimes = rftime21.size()+rftime21.size()+rftime21.size();
  double t1 = 0; if(rftime21.size()>0) t1 = rftime21[0];
  double t2 = 0; if(rftime22.size()>0) t2 = rftime22[0];
  double t3 = 0; if(rftime23.size()>0) t3 = rftime23[0];

  Int_t r1 = Int_t(rftime21.size());
  Int_t r2 = Int_t(rftime22.size());
  Int_t r3 = Int_t(rftime23.size());

  // determine trigger type and count types 
  
  if(r1==0&&r2==0&r3==0) { // zero 
    type = TPmtEvent::TRIG000;
    ++pmtSummary->ntrig000;
  } else if (r1==0||r2==0||r3==0) { 
    type = TPmtEvent::TRIG0XX;
    ++pmtSummary->ntrig5xx; 
  } else if (r1==5&&r2==5&&r3==5) { //five 
    type = TPmtEvent::TRIG555; 
    ++pmtSummary->ntrig555;
  } else if ( r1==5 || r2==5 || r3==5 ) {
   type = TPmtEvent::TRIG5XX; 
    ++pmtSummary->ntrig5xx;
  } else if (r1==4&&r2==4&&r3==4) {  //four
    type = TPmtEvent::TRIG444; 
    ++pmtSummary->ntrig444;
  } else if ( r1==4 || r2==4 || r3==4 ) {
   type = TPmtEvent::TRIG4XX; 
    ++pmtSummary->ntrig4xx;
  }
  return type;
 }


void ADCFilter(int iB, int iC) 
{
  for (int is = 0; is<NS; ++is) {
    if (digitizer_waveforms[iB][iC][is] > MAXADC) {
      if (is > 0) { digitizer_waveforms[iB][iC][is] = digitizer_waveforms[iB][iC][is-1];}
      else {
        int is2 = 0;
        while (digitizer_waveforms[iB][iC][is2] > MAXADC) {
          digitizer_waveforms[iB][iC][0] = digitizer_waveforms[iB][iC][is2+1];
          ++is2;
        }
      }
    }
  }
}


void Loop() 
{
 
  Long64_t entries = pmt_tree->GetEntries(); 
  cout << " pmt_tree has entries = " << entries << endl;
  pmtSummary->clear();

  std::vector<Double_t> sdigi;  // source
  std::vector<Double_t> ddigi;  // baseline subtracted

  // collect info for summary 
  std::vector<Int_t> vtrig;
  std::vector<Int_t> vpmt;
  std::vector<Int_t> vtmax;
  std::vector<Double_t> vqmax;
  std::vector<Double_t> vsum;
  std::vector<Double_t> vnoise;

  // loop over events 
  for(Long64_t ientry=0; ientry<entries; ++ientry) {
    pmt_tree->GetEntry(ientry);
    if(ientry%1000==0) printf(" ..... %lld \n",ientry);
    Int_t trigType = triggerInfo();
    for(UInt_t ib=0; ib<NB; ++ib) {
      for(UInt_t ic=0; ic<NC; ++ic) {
          // filter waveforms for stuck bits
          ADCFilter(ib,ic);
          // get pmt number
          int ipmt = toPmtNumber(ib,ic);
          if(ipmt<0||ipmt>=NPMT) continue;
              
          // make a vector of samples for sorting.
          sdigi.clear();
          ddigi.clear();
          
          // Find the sample median and it's "sigma".
          for (UInt_t is=0; is<NS; ++is) {
            sdigi.push_back(double(digitizer_waveforms[ib][ic][is]));
          }

          std::sort(sdigi.begin(), sdigi.end());
          double baselineMedian = sdigi[0.5*double(NS)];
          double baselineSigma = sdigi[0.16*double(NS)];
          baselineSigma = std::abs(baselineSigma-baselineMedian);
          double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);
          
          // loop over digitizations
          Int_t tmax=0;
          double qmax=0;
          double sum=0; 
           for(UInt_t is=0 ; is<NS; ++is) {
            double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])-baselineMedian);
            if(digi>qmax) {
              qmax=digi;
              tmax=Int_t(is)+1;
            }
            ddigi.push_back(digi);
            if(digi>3.0*noise) { 
              //if(is>450&&is<470) 
              sum+=digi;
            }
          } // loop over digitizations
          if(ientry%1000==0) printf(" \t ipmt %i trig %i tmax %i qmax %f sum %f noise %f \n",ipmt,trigType,tmax,qmax,sum,noise);
          vtrig.push_back(trigType);
          vpmt.push_back(ipmt);
          vtmax.push_back(tmax);
          vqmax.push_back(qmax);
          vsum.push_back(sum);
          vnoise.push_back(noise);
      }
    }
  }
  // fill rest of summary
  // ugly utility arrays used in summary
  Double_t x[NPMT], y[NPMT], z[NPMT],y2[NPMT],z2[NPMT],ex[NPMT], ey[NPMT], ez[NPMT];
  Double_t norm[NPMT], aveNoise[NPMT]; 

  // zero utility arrays
  for(Int_t j=0; j<NPMT; ++j) {
    x[j]=Double_t(j); ex[j]=0;  
    y[j]=0; z[j]=0; y2[j]=0; z2[j]=0; ey[j]=0; ez[j]=0;
    norm[j]=0; aveNoise[j]=0;
  }
  
  printf(" \n \n \t calculating averages with %lu entries \n",vpmt.size());
 
  for (UInt_t k=0 ;k<vpmt.size(); k++) {
    int ipmt = vpmt[k];
    // cosmic cut
    bool cut = vtmax[k]>440&&vtmax[k]<480&&vsum[k]<5000;
    //if(k%100==0) printf(" \t ipmt %i norm %.0f trig %i tmax %i qmax %f sum %f noise %f \n",ipmt,norm[ipmt],vpmt[k],vtrig[k],vtmax[k],vqmax[k],vsum[k],vnoise[k]);
    if(cut) { 
      y[ipmt]+=vqmax[k];
      y2[ipmt]+=pow(vqmax[k],2.);
      z[ipmt]+=vsum[k];
      z2[ipmt]+=pow(vsum[k],2.);
      aveNoise[ipmt]+=vnoise[k];
      norm[ipmt]+=1.0;
    }
  }
  for(Int_t j=0; j<NPMT; ++j) {
    y[j]/= norm[j]; z[j]/=norm[j]; y2[j]/=norm[j]; z2[j]/=norm[j]; aveNoise[j]/=norm[j];
  }

  for(Int_t j=0; j<NPMT; ++j) {
    ey[j]= sqrt( (y2[j]-pow(y[j],2.))/norm[j]);
    ez[j]= sqrt( (z2[j]-pow(z[j],2.))/norm[j]);
  }

  //for(Int_t j=0; j<NPMT; ++j) 
  //  printf(" ipmt %i norm %i qmax %.2f +/- %.2f sum %.2f +/- %.2f \n",j,int(norm[j]),y[j],ey[j],z[j],ez[j]);

  // store averages in summary
  for(Int_t j=0; j<NPMT; ++j) {
    pmtSummary->norm[j]=norm[j];
    pmtSummary->qmax[j]=y[j];
    pmtSummary->eqmax[j]=ey[j];
    pmtSummary->qsum[j]=z[j];
    pmtSummary->eqsum[j]=ez[j];
    pmtSummary->noise[j]=aveNoise[j];
  }

  pmtSummary->print();
  summaryTree->Fill();
}

void readFile(TString fileName) 
{
  cout << " reading file " << fileName << endl;
  TFile*  fin = new TFile(fileName, "READ");
  pmt_tree = (TTree *)fin->Get("pmt_tree");
  if(!pmt_tree) {
    printf(" cannot find pmt_tree in file %s\n",fileName.Data());
    fin->Close();
    return;
  }

  pmt_tree->SetBranchAddress("digitizer_waveforms", &digitizer_waveforms);
  
  Loop();

  fin->Close();
  
}

void summary(TString dirName="PDS_beamtime_files")
{
  // get list of files
  std::vector<std::string> fileList;
  TString fullDirName= TString("2017/")+dirName;
  void *dirp = gSystem->OpenDirectory(fullDirName);
  if (!dirp) return;
  char *direntry;
  Long_t id, size,flags,modtime;
  //loop on all entries of this directory
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) { 
    //cout << direntry << endl;
    string fname = string(direntry);
    if ( strstr(fname.c_str(), "PDSout" )==NULL ) continue;
    if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
    //string tag= fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_"));
    fileList.push_back(fname);
  }

  // structure for holding pmt info 
  //open output file
  TString summaryFileName = TString("pdsOutput/pdsSummary_")+dirName+ TString(".root");
  TFile* summaryFile = new TFile(summaryFileName,"recreate");
  summaryFile->cd();
  printf(" opening summary file %s \n",summaryFileName.Data());
  summaryTree = new TTree("summaryTree","summaryTree");
  pmtSummary  = new TPmtSummary();
  summaryTree->Branch("pmtSummary",&pmtSummary);


  printf(" total of files in %s is %lu \n list of tags: \n",dirName.Data(),fileList.size());

  // loop over files
  for( unsigned ifile =0; ifile < fileList.size() ; ++ifile ) {
    printf(" %i %s \n",ifile,fileList[ifile].c_str());
    TString fullName = fullDirName+TString("/")+TString(fileList[ifile].c_str());
    readFile(fullName);
  }

  summaryFile->Write();

}

