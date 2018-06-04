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
ULong_t timeZero   = 1501535896380432212;  // starting from run 0 event 0
//const ULong_t timeZero = 1501543704222496000;
//1501543650000000000; 
ULong_t timeZeroSec = 1501543650; 
const ULong_t NANO = 8;
const float jump = 190.0E6;


std::vector<ULong_t> pdsNano;
std::vector<ULong_t> pdsSec;
std::vector<float> pdsEnergy;
std::vector<float> pdsMicro;
std::vector<ULong_t> bsum1(32);
std::vector<ULong_t> bsum2(32);
std::vector<ULong_t> bsum3(32);

enum {NRUNS=90,NCHECKS=7};
int checkSum[NRUNS][NCHECKS];

void bitSum(UInt_t word, std::vector<ULong_t> &bsum) 
{
  ULong_t lword = ULong_t(word);
  std::bitset<32> foo(lword);
  for(unsigned i=0; i<bsum.size(); ++i) if(foo.test(i)) ++bsum[i];
}


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

void check(int icheck=1)
{
  TString fileTag;
  if(icheck==0) fileTag=TString("lowAna-pmtChain-fix4-0-0");
  else fileTag=TString("lowAnaNoAlign-0-0");
  printf(" calling checkNo with icheck %i fileTag %s \n",icheck,fileTag.Data());

  ULong_t STEP = pow(2.0,31);
  ULong_t HALF =pow (2.0,30); 
  
  
  // open PDS file
  TString inputPDSFileName = TString("../pdsOutput/")+fileTag+TString(".root");
  printf(" check is opening file %s \n",inputPDSFileName.Data()); 
  TFile *inpdsfile = new TFile(inputPDSFileName);
  TTree *pmtTree=NULL;

  // tree has to be in file
  Long64_t bSize=0;
  pmtTree = (TTree*) inpdsfile->Get("pmtTree");
  if(!pmtTree) return;
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

  TH1D* hBits0 = new TH1D("Bits0"," board 0 bits",32,0,32);
  TH1D* hBits1 = new TH1D("Bits1"," board 1 bits",32,0,32);
  TH1D* hBits2 = new TH1D("Bits2"," board 2 bits",32,0,32);


  TH1F* hQinHit = new TH1F("QinHit"," ADC counts of bins in hit ",2000,0,20);
  TNtuple* ntPmt = new TNtuple("ntPmt","pmt by board","ev:run:pdst:dt0:dt1:dt2:nrf:delta0:delta1:delta2:off:rt0:rt1:rt2");
  TNtuple* ntRaw = new TNtuple("ntRaw","board time diffs","ev:run:pdst:dpdst:d0:d1:d2:delta0:delta1:delta2");
  Int_t nQTimeBins = Int_t( gate/1000);
  TH1F* hQTime = new TH1F("QTime","summed charge versus time(micro-sec)",nQTimeBins,0,float(gate)*toMicro);
  TNtuple *ntMatch = new TNtuple("ntMatch"," TPC PDS matching ","itpc:ipds:tpcsec:pdssec:tpcnano:pdsnano:tb0:tb1:tb2:diff:diff0:qsum");
  TNtuple *ntJump = new TNtuple("ntJump"," holdoffs ","run:ev:njump:bjump0:bjump1:bjump2:ncheck:pdst:dpdst:d0:d1:d2:rft0:rft1:rft2");
  TNtuple *ntClock = new TNtuple("ntClock"," clocks ","run:ev:pdst:dpdst:tb0:tb1:tb2:d0:d1:d2:rft0:rft1:rft2");

  typedef struct {
    Int_t  run;
    Int_t  event;
    Int_t  compSec;
    Int_t rf0;
    Int_t rf1;
    Int_t rf2;
    Long64_t compNano;
    UInt_t  dt0;
    UInt_t  dt1;
    UInt_t  dt2;   // caen digitizer time 
    Double_t tPrompt; // one for each board
    Double_t tPromptToRF;
  } CLOCK;

  static CLOCK clock;
  /*
    C : a character string terminated by the 0 character
    B : an 8 bit signed integer (Char_t)
    b : an 8 bit unsigned integer (UChar_t)
    S : a 16 bit signed integer (Short_t)
    s : a 16 bit unsigned integer (UShort_t)
    I : a 32 bit signed integer (Int_t)
    i : a 32 bit unsigned integer (UInt_t)
    F : a 32 bit floating point (Float_t)
    D : a 64 bit floating point (Double_t)
    L : a 64 bit signed integer (Long64_t)
    l : a 64 bit unsigned integer (ULong64_t)
    O : [the letter o, not a zero] a boolean (Bool_t)
    */


  TTree *tClock = new TTree("TClk"," PDS clocks ");
  tClock->Branch("clk",&clock,"run/I:event/I:sec/I:rf0/I:rf1/I:rf2/I:nano/l:dt0/i:dt1/i:dt2/i:tprompt/D:tPromptToRF/D");


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

  for(unsigned i=0; i<bsum1.size(); ++i) { bsum1[i]=0;bsum2[i]=0;bsum3[i]=0;}

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

  unsigned njump =0;
  unsigned ncheck =0;
  unsigned bjump[NB]={0,0,0};
  Int_t rft[NB]={0,0,0};

  unsigned ifirst = 0; //5000 * 44;
  unsigned ilast = bSize; //5000 * 45;
  TString mess;
  for(int irun=0; irun<NRUNS; ++irun) for(itype=0; itype<NCHECKS; ++itype) checkSum[irun][itype]=0;

  for(unsigned entry = ifirst ; entry < ilast; ++entry ) {
    pmtTree->GetEntry(entry);
    mess.Clear();
 
    // look for hits in time with TPC event
    Int_t run = pmtEvent->run;
    ULong_t pdsCompSec = pmtEvent->compSec;
    ULong_t pdsCompNano = pmtEvent->compNano; 
    ULong_t pdsCompTime = pdsCompSec*1000000000 + pdsCompNano ;
    if( entry==0) printf(" .... entry %u run %i pdsCompTime %lu %E\n",entry,run,pdsCompTime,double(pdsCompTime));
    if(entry==ifirst) {
      timeZero = pdsCompTime;
      printf(" timeZero %lu %0.9E \n",timeZero,float(timeZero));
    }

    // gather RF clocks
    for(int ib=0; ib< NB; ++ib) rft[ib]=0;
    if( pmtEvent->rft21.size()>0) rft[0]=pmtEvent->rft21[0];
    if( pmtEvent->rft22.size()>0) rft[1]=pmtEvent->rft22[0];
    if( pmtEvent->rft23.size()>0) rft[2]=pmtEvent->rft23[0];

    // fill tree
    clock.run      = pmtEvent->run;
    clock.event    = pmtEvent->event;
    clock.compSec  = pmtEvent->compSec;
    clock.rf0      = rft[0];
    clock.rf1      = rft[1];
    clock.rf2      = rft[2];
    clock.compNano = pmtEvent->compNano;
    clock.dt0      = pmtEvent->dtime[0];
    clock.dt1      = pmtEvent->dtime[1];
    clock.dt2      = pmtEvent->dtime[2];
    clock.tPrompt  = pmtEvent->tPrompt; // one for each board
    clock.tPromptToRF = pmtEvent->tPromptToRF;
    tClock->Fill();

    pdsCompTime -= timeZero;
    //if(run<44) continue;
    //if(run>44) break;
    if( entry%1 == 0) printf(" XXXXX %s entry %u run %i pdsCompTime %lu %E %E %E \n",pmtEvent->tag.c_str(),entry,run,pdsCompTime,double(pdsCompTime),
        clock.tPrompt,clock.tPromptToRF);

    bitSum(pmtEvent->dtime[0],bsum1);
    bitSum(pmtEvent->dtime[1],bsum2);
    bitSum(pmtEvent->dtime[2],bsum3);

    // raw dtimes
    for(int ib=0; ib<NB; ++ib) {
      rdiff[ib]= double(pmtEvent->dtime[ib]) - rawLast[ib];
      rawLast[ib]=double(pmtEvent->dtime[ib]);
    }
    rdiff[3] = double(pdsCompTime) - rawLast[3];
    rawLast[3]=double(pdsCompTime);
    ntRaw->Fill( float(entry),float(run),float(pdsCompTime),rdiff[3],rawLast[0],rawLast[1],rawLast[2],rdiff[0],rdiff[1],rdiff[2]);

    if(rdiff[0]<0) bitSum(pmtEvent->dtime[0],bsum1);
    if(rdiff[1]<0) bitSum(pmtEvent->dtime[1],bsum2);
    if(rdiff[2]<0) bitSum(pmtEvent->dtime[2],bsum3);

    
    /// new run 
    if(thisRun ==0 ) {
      thisRun=run;
      for(int ib=0; ib<NB; ++ib) {
        offset[ib]=0;
        dtime[ib] = ULong_t(pmtEvent->dtime[ib]);
        dstime[ib]= ULong_t(pmtEvent->dtime[ib]);
        dstartNano[ib]=NANO*dstime[ib]; 
      }
      printf("***** NEW RUN entry %u run %i pdsTime %lu setting board dstartNano %lu %lu %lu  \n",entry,run,pdsCompTime,
          dstartNano[0],dstartNano[1],dstartNano[2]); 
      printf("\t run %2i MISSA %4i MISSB %4i MISSC %4i MISSD %4i MISSE %4i comp jumps %4i board jumps %4i \n",
        run,checkSum[run][0],checkSum[run][1],checkSum[run][2],checkSum[run][3],checkSum[run][4],checkSum[run][5],checkSum[run][6] );

    }


    // take out steps 
    for(int ib=0; ib<NB; ++ib) {
      dtime[ib] = ULong_t(pmtEvent->dtime[ib]);
      if( dtime[ib]-dlast[ib]>HALF && dlast[ib]>0 ) offset[ib] += STEP;
      dstime[ib]=dtime[ib]+offset[ib];
    }

    float dtime01 = float(dstime[0]) - float(dstime[1]);
    //if(dtime01>5000.) printf(" BIG>>>>> entry %u pds %lu dt0 %lu dt1 %lu dt2 %lu \n",entry,pdsCompTime,dstime[0],dstime[1],dstime[2]);


    Double_t tRFave = pmtEvent->tRFave;
    if(tRFave>1) ++nrf; // count RF pulses

   // convert to ns
    for(int ib=0;ib<NB; ++ib) dtnano[ib]=dstime[ib]*NANO;
    // these times in micro-seconds
    if(dtnano[0] < dstartNano[0]) printf(" WARNING : %lu %lu \n",dtnano[0] , dstartNano[0]);
   
    //printf(" entry %i dtime[0] %ul dtime[1] %ul dtnano[0] %ul dtnano[1] %ul diff %0.9E   \n",
    //  entry,dtime[0],dtime[1],dtnano[0]-dstartNano[0],dtnano[1]-dstartNano[1],fbt0-fbt1); 
    // check board alignments
 
    int njumps=0;
    int nzeros=0;
    for(int ib=0; ib<NB; ++ib) {
      if(rdiff[ib]*8>jump) ++njumps;
      if(rdiff[ib]<0) ++nzeros;
    }

    /*
       MISS_____
       A == computer time jump but missing at least 2 board jumps
       B == >2 board jumps, no computer time jump
       C == at least 1 board jump but less than 2
       D == missing RF ~500
       E == missing RF ~600
       */

    // rf checks
    for(unsigned ib=0; ib<NB; ++ib) {
      if( float(rdiff[ib]*8) >1.0E8 && float(rdiff[ib]*8) < 3.0E8 ) ++bjump[ib]; 
    }

      int nrf500=0;
    int nrf600=0;
    for(int ib=0; ib<NB; ++ib) {
      if(rft[ib]>450&&rft[ib]<550) ++nrf500;
      if(rft[ib]>550&&rft[ib]<650) ++nrf600;
    }

    if(rdiff[3] > jump && njumps<2&&nzeros==0) {
      mess = mess + TString("A");
      ++checkSum[run][0];
    }
    if( rdiff[3]<jump && njumps>=2){
      mess = mess + TString("B");
      ++checkSum[run][1];
    }
    if( njumps>0&&njumps<2) {
      mess = mess + TString("C");
      ++checkSum[run][2];
    }
    if(nrf500>0&&nrf500<3) {
      ++checkSum[run][3];
      mess = mess + TString("D");
    }
    if(nrf600>0&&nrf600<3) {
      ++checkSum[run][4];
      mess = mess + TString("E");
    }
    if(mess.Sizeof()>1) mess = TString("MISS")+mess;
    mess.Resize(9);

    if(njumps>=2)  ++checkSum[run][6];      


    if(rdiff[3] > jump) {
      ncheck=0;
      ++checkSum[run][5];      
      mess = TString("J ")+mess;
    } else 
      mess = TString("  ")+mess;

    printf(" %s %3u ev %3u run %3i pdst (%10lu,%10lu)  dpdt = %10.0f rft(%3u,%3u,%3u) digi = (%10u,%10u,%10u) diffs = (%10.0f,%10.0f,%10.0f) (ns) \n",
        mess.Data(),++ncheck,// bjump[0],bjump[1],bjump[2],
        entry,run,pdsCompSec,pdsCompNano,rdiff[3], rft[0],rft[1],rft[2],
        pmtEvent->dtime[0],pmtEvent->dtime[1],pmtEvent->dtime[2],
        rdiff[0]*8,rdiff[1]*8,rdiff[2]*8); //,  clock.tPrompt, clock.tPromptToRF);

    ntJump->Fill(float(run),float(entry),float(njump),float(bjump[0]), float(bjump[1]), float(bjump[2]), float(ncheck),
        float(pdsCompTime),float(rdiff[3]),
        float(rdiff[0]*8),float(rdiff[1]*8),float(rdiff[2]*8),
        float(rft[0]),float(rft[1]),float(rft[2]) );
    ntClock->Fill(float(run),float(entry),float(pdsCompTime),float(rdiff[3]),float(dtnano[0]),float(dtnano[1]),float(dtnano[2]),
        float(pmtEvent->dtime[0]), float(pmtEvent->dtime[1]), float(pmtEvent->dtime[2]), 
        float(rft[0]),float(rft[1]),float(rft[2]) );
    // if(entry==76607) { printf(" BLAH %u %lu %.0f \n",entry,pmtEvent->dtime[0],float(pmtEvent->dtime[0]))}
    
    //= new TNtuple("ntJump"," holdoffs ","run:ev:njump:ncheck:pdst:dpdst:d0:d1:d2");

  // save last digi times
  // take out steps 
  for(int ib=0; ib<NB; ++ib) dlast[ib]= Long_t(pmtEvent->dtime[ib]);

 
  }

  for(unsigned i=0; i<bsum1.size(); ++i)  hBits0->SetBinContent(i+1,bsum1[i]);
  for(unsigned i=0; i<bsum2.size(); ++i)  hBits1->SetBinContent(i+1,bsum2[i]);
  for(unsigned i=0; i<bsum3.size(); ++i)  hBits2->SetBinContent(i+1,bsum3[i]);
  

  printf(" number of PDS triggers is %lld tclock size %lld summary of errors : \n",ntClock->GetEntries(),tClock->GetEntries());

  for(int irun=0; irun<NRUNS; ++irun) 
    printf("\t run %2i MISSA %4i MISSB %4i MISSC %4i MISSD %4i MISSE %4i comp jumps %4i board jumps %4i \n",
        irun,checkSum[irun][0],checkSum[irun][1],checkSum[irun][2],checkSum[irun][3],checkSum[irun][4],checkSum[irun][5],checkSum[irun][6] );
  
  outfile->Write();
}
