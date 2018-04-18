#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TPmtSummary.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
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


void summary(TString tag= "PDS_beamtime_files")
{
  //TString inputFileName = TString("../pdsOutput/pdsSummary_")+tag+TString(".root");
  //printf(" opening file %s \n",inputFileName.Data()); 
  //TFile *infile = new TFile(inputFileName);
  TChain *sumTree= new TChain("summaryTree");
  //sumTree->Add("../pdsOutput/pmtSummary_07-31-1555_0.root");
  //sumTree->Add("../pdsOutput/pmtSummary_07-31-2133_0.root");
  //sumTree->Add("../pdsOutput/pmtSummary_07-31-1853_0.root");
  sumTree->Add("../pdsOutput/pmtSummary_07-31-1600_0.root");
  Long64_t aSize=0;
  if(sumTree) aSize=sumTree->GetEntries();
  else  printf(" no summaryTree  \n");
  printf(" summaryTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("summary-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  TNtuple *ntime = new TNtuple("ntime","time","ev:sec:ns:rf1:rf2:rf3:drf1:dt1:ddt1:dt2:dt3:tof1:tof2:tof3");
  TH1F* hRF12 = new TH1F("hRF12"," rf1 -rf2 (ns)",100,-50,50);
  TH1F* hRF13 = new TH1F("hRF13"," rf1 -rf3 (ns)",100,-50,50);

  TPmtSummary *pmtSum = new TPmtSummary();
  sumTree->SetBranchAddress("pmtSummary",&pmtSum);
  outfile->cd();
  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
  TString name, title;

  unsigned late=0;
  printf(" sizeof compSec %lu sizeof compNano %lu \n",sizeof(pmtSum->vcompSec[0]),sizeof(pmtSum->vcompNano[0]));
  Int_t ayear = 60*60*24*365;
  Int_t aday = 60*60*24;

  std::vector<Long64_t> vnave;   // nano averages 
  std::vector<Long64_t> vtave;   // nano averages 
  

  printf(" year = %i s day %i s \n",ayear,aday);
  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t min = pmtSum->getMin();
    Int_t seg = pmtSum->getSegment();
    Int_t time = 60*min+60*60*24*day+60*60*24*30*month;
    
    if(entry%100==0) printf("...entry %u tag %s  month %i day %i min %i seg %i time %i (s) vector size %i \n",
        entry,tag.c_str(),month,day,min,seg,time,int(pmtSum->vtrig.size()) );

    // average to track the drift 
    vnave.clear();
    vtave.clear();
    unsigned nstep = pmtSum->vcompSec.size()/5;
    printf(" size of sec %zu nstep %u  \n",pmtSum->vcompSec.size(),nstep);
    
    Long64_t na=0;
    Long64_t norm=0;
    for(unsigned i=0; i<pmtSum->vcompSec.size(); ++i ){
      na +=pmtSum->vcompNano[i];
      ++norm;
      if(i%nstep==0&&i>0) {
        na /= norm;
        vnave.push_back(na);
        vtave.push_back(i);
        na=0;
        norm=0;
      } 
    }


    //for(unsigned i=0; i<vnave.size(); ++i ) printf(" %u time %lld ave %E (ms) \n",i,vtave[i],double(vnave[i])*1.0E-6);


    float adcToNs = 4.0;   // 250 MHz digi
    float trigToNs = 8.0;  // 150 MHz trigger clock

    Long64_t jstep =0;
    float drf1=0;
    float ddt1=0;
    for(unsigned i=0; i<pmtSum->vcompSec.size(); ++i ){
      if(i%nstep==0&&i>0&&jstep<vnave.size()-2) ++jstep;
      float tns  = float(pmtSum->vcompNano[i]);
      float tave = float(vnave[jstep])- float(vnave[0]);
      float tsec= pmtSum->vcompSec[i]-pmtSum->vcompSec[0];
      float tsecInNs = tsec*1.0E9;
      float tcms = tns - tave;
      float dt1 = float(pmtSum->vdtime1[i]-pmtSum->vdtime1[0])*trigToNs; // convert to ns
      if(i>0) ddt1 = float(pmtSum->vdtime1[i]-pmtSum->vdtime1[i-1])*trigToNs; // convert to ns
      float dt2 = float(pmtSum->vdtime2[i]-pmtSum->vdtime2[0])*trigToNs; // convert to ns
      float dt3 = float(pmtSum->vdtime3[i]-pmtSum->vdtime3[0])*trigToNs; // convert to ns
      float rf1 = float(pmtSum->vrf1[i])*adcToNs; // convert to ns
      if(i>0) drf1 = float(pmtSum->vrf1[i]- pmtSum->vrf1[i-1])*adcToNs; 
      float rf2 = float(pmtSum->vrf2[i])*adcToNs; // convert to ns
      float rf3 = float(pmtSum->vrf3[i])*adcToNs; // convert to ns
      float tof1 = float(pmtSum->vdtime1[i])*trigToNs -  float(pmtSum->vrf1[i])*adcToNs;
      float tof2 = float(pmtSum->vdtime2[i])*trigToNs -  float(pmtSum->vrf2[i])*adcToNs;
      float tof3 = float(pmtSum->vdtime3[i])*trigToNs -  float(pmtSum->vrf3[i])*adcToNs;

      if(rf1>0) {
        hRF12->Fill(rf1-rf2);
        hRF13->Fill(rf1-rf3);
      }
       
      // convert differnces to microsec
      drf1 /=1.0E3;
      ddt1 /=1.0E3;
      if(i%1000==0) printf(" %i %f tns %E (ns) dt %E rf %E  \n",i,tsec,tns,dt1,rf1);
      ntime->Fill(float(i),tsec,tns,rf1,rf2,rf3,drf1,dt1,ddt1,dt2,dt3,tof1,tof2,tof3);
    }
  }
  outfile->Write();
}
