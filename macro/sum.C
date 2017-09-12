#include <vector>
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


void sum(TString tag= "PDS_beamtime_files_max1000")
{
  TString inputFileName = TString("../pdsOutput/pdsSummary_")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *sumTree=NULL;

  // tree has to be in file
  Long64_t aSize=0;
  sumTree = (TTree*) infile->Get("summaryTree");
  if(sumTree) aSize=sumTree->GetEntriesFast();
  printf(" smmaryTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("sum-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  TPmtSummary *pmtSum = new TPmtSummary();
  sumTree->SetBranchAddress("pmtSummary",&pmtSum);

   for(unsigned entry =0; entry < aSize; ++entry ) {
     sumTree->GetEntry(entry);
     string tag = pmtSum->tag;
     Int_t month = pmtSum->getMonth();
     Int_t day = pmtSum->getDay();
     Int_t hour = pmtSum->getHour();
     Int_t seg = pmtSum->getSegment();
     if(entry%1==0) printf("...entry %i tag %s  month %i day %i hour %i seg %i \n",entry,tag.c_str(),month,day,hour,seg);
    }

  outfile->Write();
}
