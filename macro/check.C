#include <vector>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMultiGraph.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TNtuple.h>
//#include "TPmtEvent.hxx"

enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};

enum {NTRK=0,NCLU,NSHO,UNKN,ALL};
const UInt_t gate = 4500000; // PDS gate 
const UInt_t pulse = 180;  // time between micro pulses
const float toMicro = 1.0E-3;
const float toMilli= 1.0E-6;
const ULong_t timeZero   = 1501535896380432212;  // starting from run 0 event 0
//const ULong_t timeZero = 1501543704222496000;
//1501543650000000000; 
const ULong_t timeZeroSec = 1501543650; 
const ULong_t NANO = 8;
const float jump = 195.0E6;


std::vector<ULong_t> pdsNano;
std::vector<ULong_t> pdsSec;
std::vector<float> pdsEnergy;
std::vector<float> pdsMicro;


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

void check()
{

  TString fileTag=TString("lowAnaNoAlign-0-0");

  printf(" timeZero %lu %0.9E \n",timeZero,float(timeZero));
  ULong_t STEP = pow(2.0,31);
  ULong_t HALF =pow (2.0,30); 
  
  
  // open PDS file
  //TString fileTag=TString("lowAnaNoAlign-0-0");
  TString inputPDSFileName = TString("../pdsOutput/")+fileTag+TString(".root");
  printf(" opening file %s \n",inputPDSFileName.Data()); 
  TFile *inpdsfile = new TFile(inputPDSFileName);
  TTree *pmtTree=NULL;

  // tree has to be in file
  Long64_t bSize=0;
  pmtTree = (TTree*) inpdsfile->Get("pmtTree");
  if(pmtTree) bSize=pmtTree->GetEntriesFast();
  printf(" pmtTree with %i entries \n",int(bSize));
  TPmtEvent *pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&pmtEvent);
  std::vector<TPmtHit> hit;



  // open ouput file and make some histograms
  TString outputFileName = TString("check_")+ fileTag+ TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  outfile->cd();
  TH1F* hQinHit = new TH1F("QinHit"," ADC counts of bins in hit ",2000,0,20);
  TNtuple* ntPmt = new TNtuple("ntPmt","pmt by board","ev:run:pdst:dt0:dt1:dt2:nrf:delta0:delta1:delta2:off:rt0:rt1:rt2");
  TNtuple* ntRaw = new TNtuple("ntRaw","board time diffs","ev:run:pdst:dpdst:d0:d1:d2:delta0:delta1:delta2");
  Int_t nQTimeBins = Int_t( gate/1000);
  TH1F* hQTime = new TH1F("QTime","summed charge versus time(micro-sec)",nQTimeBins,0,float(gate)*toMicro);
  TNtuple *ntMatch = new TNtuple("ntMatch"," TPC PDS matching ","itpc:ipds:tpcsec:pdssec:tpcnano:pdsnano:tb0:tb1:tb2:diff:diff0:qsum");


  // look at corresponding PDS data
  std::vector<float> qbt0;
  std::vector<float> qbt1;
  std::vector<float> qbt2;
  std::vector<float> qb0;
  std::vector<float> qb1;
  std::vector<float> qb2;


  float qboard[NB];
  ULong_t dtime[NB]; 
  ULong_t dlast[NB];
  ULong_t offset[NB];
  ULong_t dstime[NB];
  ULong_t dtnano[NB];
  ULong_t dstartNano[NB];
  std::vector<ULong_t> timeb0;
  std::vector<ULong_t> timeb1;
  std::vector<ULong_t> timeb2;

  for(int ib=0; ib<NB; ++ib) {
    dtime[ib]=0; dlast[ib]=0; offset[ib]=0; dstime[ib]=0; dtnano[ib]=0; dstartNano[ib]=0;
  }

  qb0.clear();qb1.clear();qb2.clear();
  qbt0.clear(); qbt1.clear(); qbt2.clear();
  Int_t nrf=0;

  pmtTree->GetEntry(0);
  ULong_t pdsCompSecf = pmtEvent->compSec;
  ULong_t pdsCompNanof = pmtEvent->compNano; 
  ULong_t pdsCompTimef = pdsCompSecf*1000000000 + pdsCompNanof;
  pmtTree->GetEntry(bSize-1);
  ULong_t pdsCompSecl = pmtEvent->compSec;
  ULong_t pdsCompNanol = pmtEvent->compNano; 
  ULong_t pdsCompTimel = pdsCompSecl*1000000000 + pdsCompNanol;

  
  Int_t thisRun = -1;

  unsigned firstEntry = 0;
  float ebt0=0;
  float ebt1=0;
  float ebt2=0;

  double rawLast[4]={0,0,0,0};
  double rdiff[4]={0,0,0,0};


  for(unsigned entry = 224221 ; entry < bSize; ++entry ) {
    pmtTree->GetEntry(entry);
 
    // look for hits in time with TPC event
    Int_t run = pmtEvent->run;
    ULong_t pdsCompSec = pmtEvent->compSec;
    ULong_t pdsCompNano = pmtEvent->compNano; 
    ULong_t pdsCompTime = pdsCompSec*1000000000 + pdsCompNano ;
    if( entry==0) printf(" .... entry %u run %i pdsCompTime %lu %E\n",entry,run,pdsCompTime,double(pdsCompTime));
    if(pdsCompTime < timeZero) continue;
    pdsCompTime -= timeZero;
    if(run>44) break;
    if( entry%5000 == 0) printf(" .... %s entry %u run %i pdsCompTime %lu %E\n",pmtEvent->tag.c_str(),entry,run,pdsCompTime,double(pdsCompTime));

    // raw dtimes
    for(int ib=0; ib<NB; ++ib) {
      rdiff[ib]= double(pmtEvent->dtime[ib]) - rawLast[ib];
      rawLast[ib]=double(pmtEvent->dtime[ib]);
    }
    rdiff[3] = double(pdsCompTime) - rawLast[3];
    rawLast[3]=double(pdsCompTime);
    ntRaw->Fill( float(entry),float(run),float(pdsCompTime),rdiff[3],rawLast[0],rawLast[1],rawLast[2],rdiff[0],rdiff[1],rdiff[2]);
    
    // new run 
    if(thisRun ==0 ) {
      thisRun=run;
      for(int ib=0; ib<NB; ++ib) {
        offset[ib]=0;
        dtime[ib] = ULong_t(pmtEvent->dtime[ib]);
        dstime[ib]= ULong_t(pmtEvent->dtime[ib]);
        dstartNano[ib]=NANO*dstime[ib]; 
      }
      printf("\t entry %u run %i pdsTime %lu setting board dstartNano %lu %lu %lu  \n",entry,run,pdsCompTime,
          dstartNano[0],dstartNano[1],dstartNano[2]); 
    }


    // take out steps 
    for(int ib=0; ib<NB; ++ib) {
      dtime[ib] = ULong_t(pmtEvent->dtime[ib]);
      if( dtime[ib]-dlast[ib]>HALF && dlast[ib]>0 ) offset[ib] += STEP;
      dstime[ib]=dtime[ib]+offset[ib];
    }

    float dtime01 = float(dstime[0]) - float(dstime[1]);
    if(dtime01>5000.) printf(" BIG>>>>> entry %u pds %lu dt0 %lu dt1 %lu dt2 %lu \n",entry,pdsCompTime,dstime[0],dstime[1],dstime[2]);


    Double_t tRFave = pmtEvent->tRFave;
    if(tRFave>1) ++nrf; // count RF pulses

   // convert to ns
    for(int ib=0;ib<NB; ++ib) dtnano[ib]=dstime[ib]*NANO;
    // these times in micro-seconds
    if(dtnano[0] < dstartNano[0]) printf(" WARNING : %lu %lu \n",dtnano[0] , dstartNano[0]);
   
    //printf(" entry %i dtime[0] %ul dtime[1] %ul dtnano[0] %ul dtnano[1] %ul diff %0.9E   \n",
    //  entry,dtime[0],dtime[1],dtnano[0]-dstartNano[0],dtnano[1]-dstartNano[1],fbt0-fbt1); 
    // check board alignments
    if(rdiff[0]>0&&rdiff[1]>0&&rdiff[2]>0) {
      if(rdiff[3] > jump && ( rdiff[0]*8<jump || rdiff[1]*8<jump || rdiff[2]*8<jump)) 
        printf(" MISALIGN event %u run %i %f (%f ,%f , %f) \n",entry,run,rdiff[3],rdiff[0]*8,rdiff[1]*8,rdiff[2]*8); 
    }
    //

    printf(" CHECK event %u run %i pdstime= %lu diff_pds = %.0f board digi = (%lu , %lu , %lu) board diffs = (%.0f ,%.0f , %.0f) (ns) \n",
        entry,run,pdsCompTime,rdiff[3],
        pmtEvent->dtime[0],pmtEvent->dtime[1],pmtEvent->dtime[2],
        rdiff[0]*8,rdiff[1]*8,rdiff[2]*8); 

    // save last digi times
    // take out steps 
    for(int ib=0; ib<NB; ++ib) dlast[ib]= Long_t(pmtEvent->dtime[ib]);

  }

  printf(" number of PDS triggers is %lu \n",qbt0.size());
  outfile->Write();
}
