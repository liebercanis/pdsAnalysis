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
std::string tag;
void getTag(std::string fname) { tag = fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")); return;}
Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}


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


void sum()
{
  TString fullDirName;
  std::vector<std::string> fileList;
  
  // get list of files
  printf("Low intensity runs range from 0731_1518 to 0731_2130. They are all in one day.\n");
  TString sumTag("low-intensity");
  // get list of files
  fullDirName= TString("/data1/gold/pdsOutput/");
  void *dirp = gSystem->OpenDirectory(fullDirName);
  cout << " pdsOut full directory name is " << fullDirName << endl;
  if (!dirp) {
    cout << " returning with NULL directory pointer ! \n"; 
    return;
  }
  fileList.clear();
  char* direntry;
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) { 
    //cout << direntry << endl;
    string fname = string(direntry);
    if ( strstr(fname.c_str(), "pmtSummary" )==NULL ) continue;
    if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
    getTag(fname); 
    Int_t month = getMonth();
    Int_t day =  getDay();
    Int_t min =  getMin();
    Int_t segment =  getSegment();
    //printf(" %i %i %i %i \n",month,day,min,segment);
    if( month==7 && day == 31 && min>1517 && min < 2131 ) fileList.push_back(fname);
  }

  printf(" total of files in %s is %lu  \n",fullDirName.Data(),fileList.size());
  if(fileList.size()<1) return;
  for(unsigned ifile = 0; ifile < fileList.size(); ++ ifile) printf("\t %i %s \n",ifile,fileList[ifile].c_str());

  //TString inputFileName = TString("../pdsOutput/pdsSummary_")+tag+TString(".root");
  //printf(" opening file %s \n",inputFileName.Data()); 
  //TFile *infile = new TFile(inputFileName);
  TChain *sumTree= new TChain("summaryTree");
  TPmtSummary *pmtSum = new TPmtSummary();
  sumTree->SetBranchAddress("pmtSummary",&pmtSum);
  
  for(unsigned ifile = 0; ifile < fileList.size(); ++ ifile) {
    TString addName = fullDirName + TString(fileList[ifile].c_str());
    sumTree->Add(addName);
  }

  //sumTree = (TTree*) infile->Get("summaryTree");
  Long64_t aSize=0;
  if(sumTree) aSize=sumTree->GetEntries();
  else  printf(" no summaryTree  \n");
  printf(" summaryTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("sum-")+sumTag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  TNtuple *ntTrigSum = new TNtuple("ntTrigSum"," trigger sum ","time:n555:n5xx:n444:n4xx:n111:n1xx:n000:n0xx");
  TNtuple *ntTime = new TNtuple("ntTime","timing","run:trig:type:tzero:rf1:rf2:rf3:tp1:tp2:tp3:dt1:dt2:dt3");

  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
 for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    std::string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t min = pmtSum->getMin();
    Int_t seg = pmtSum->getSegment();
    double time = double(min)+double(60*24*day+60*24*30*(month-7)) - 60*24*20;

    ntTrigSum->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);
    
    if(entry%100==0) printf("...entry %u tag %s  month %i day %i min %i seg %i time %0.f \n",entry,tag.c_str(),month,day,min,seg,time);

    for(unsigned itr = 0; itr< pmtSum->vtrig.size(); ++ itr) ntTime->Fill(entry,itr,pmtSum->vtrig[itr],pmtSum->tZero[0],
        pmtSum->vrf1[itr],pmtSum->vrf2[itr],pmtSum->vrf3[itr],
        pmtSum->vprompt1[itr],pmtSum->vprompt2[itr],pmtSum->vprompt3[itr],
        pmtSum->vdtime1[itr],pmtSum->vdtime2[itr],pmtSum->vdtime3[itr]);
 }
} 
