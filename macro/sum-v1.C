#include <vector>
#include <map>
#include "TTree.h"
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
  std::map<int,std::string> timeMap;
  
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
    Int_t time = min +60*24*day+60*24*30*(month-7)- 60*24*20;
    
    if( month==7 && day == 31 && min>1517 && min < 2131 ) {
      timeMap.insert ( std::pair<int,std::string>(time,fname) );
    }
  }

  printf(" total of files in %s is %lu  \n",fullDirName.Data(),timeMap.size());
  if(timeMap.size()<1) return;

  //TString inputFileName = TString("../pdsOutput/pdsSummary_")+tag+TString(".root");
  //printf(" opening file %s \n",inputFileName.Data()); 
  //TFile *infile = new TFile(inputFileName);
  TChain *sumTree= new TChain("summaryTree");
  TPmtSummary *pmtSum = new TPmtSummary();
  sumTree->SetBranchAddress("pmtSummary",&pmtSum);
  

  std::cout << "timeMap contains:\n";
  std::map<int,std::string>::iterator iter;
  for (iter=timeMap.begin(); iter!=timeMap.end(); ++iter) {
    TString addName = fullDirName + TString(iter->second.c_str());
    cout << "adding to tree " << iter->first << " => " << iter->second << "  file " << addName << endl;
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
  TNtuple *ntTime = new TNtuple("ntTime","timing","run:trig:type:tzero:rf1:rf2:rf3:dt1:dt2:dt3:ddt1:ddt2:ddt3");
  TTree *trSortTime = new TTree("SortTime","sort timing");
  Double_t sortTime1; 
  Double_t sortTime2;
  Double_t sortTime3;
  Double_t sortDtime1; 
  Double_t sortDtime2;
  Double_t sortDtime3;
  Double_t sortOff1; 
  Double_t sortOff2;
  Double_t sortOff3;

  Double_t sortTrig;
  trSortTime->Branch("trig",&sortTrig,"trig/D");
  trSortTime->Branch("t1",&sortTime1,"t1/D");
  trSortTime->Branch("t2",&sortTime2,"t2/D");
  trSortTime->Branch("t3",&sortTime3,"t3/D");
  trSortTime->Branch("dt1",&sortDtime1,"dt1/D");
  trSortTime->Branch("dt2",&sortDtime2,"dt2/D");
  trSortTime->Branch("dt3",&sortDtime3,"dt3/D");
  trSortTime->Branch("o1",&sortOff1,"o1/D");
  trSortTime->Branch("o2",&sortOff2,"o2/D");
  trSortTime->Branch("o3",&sortOff3,"o3/D");

  //trSortTime->Print();
  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
  UInt_t skip[3]={0,0,0};
  UInt_t sskip[3]={0,0,0};
  UInt_t nskip[3]={0,0,0};

  std::vector<Double_t> strig;
  std::vector<Double_t> st1;
  std::vector<Double_t> st2;
  std::vector<Double_t> st3;
  std::vector<Double_t> sdt1;
  std::vector<Double_t> sdt2;
  std::vector<Double_t> sdt3;
  std::vector<Double_t> ot1;
  std::vector<Double_t> ot2;
  std::vector<Double_t> ot3;


  Double_t offset1 =0;
  Double_t offset2 =0;
  Double_t offset3 =0;

  Int_t treeNumber=0;
  UInt_t trigsize=0;
  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    if(sumTree->GetTreeNumber()!=treeNumber) {
      treeNumber=sumTree->GetTreeNumber();
      printf("opening tree number %i \n",treeNumber);
      for(int iskip=0; iskip<3; ++iskip) skip[iskip]=0;
    }
    std::string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t min = pmtSum->getMin();
    Int_t seg = pmtSum->getSegment();
    double time = double(min)+double(60*24*day+60*24*30*(month-7)) - 60*24*20;

    ntTrigSum->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);

    if(entry%100==0) printf("...entry %u tag %s  month %i day %i min %i seg %i time %0.f \n",entry,tag.c_str(),month,day,min,seg,time);

    for(unsigned itr = 0; itr< pmtSum->vtrig.size(); ++ itr) {
      // check to see if we have reached overflow.
      if(itr>0&& pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr] > 1E6) offset1 += pmtSum->vdtime1[itr-1]; 
      if(itr>0&& pmtSum->vdtime2[itr-1]-pmtSum->vdtime2[itr] > 1E6) offset2 += pmtSum->vdtime2[itr-1]; 
      if(itr>0&& pmtSum->vdtime3[itr-1]-pmtSum->vdtime3[itr] > 1E6) offset3 += pmtSum->vdtime3[itr-1]; 
      Double_t dt1 = pmtSum->vdtime1[itr]+skip[0]+offset1;
      Double_t dt2 = pmtSum->vdtime2[itr]+skip[1]+offset2;
      Double_t dt3 = pmtSum->vdtime3[itr]+skip[2]+offset3;
      // save starting shifts
      for(int iskip=0; iskip<3; ++iskip) sskip[iskip] =  skip[iskip];

      strig.push_back(Double_t(itr + treeNumber*(pmtSum->vtrig.size())));
      st1.push_back(Double_t(pmtSum->vdtime1[itr])+offset1);
      st2.push_back(Double_t(pmtSum->vdtime2[itr])+offset2);
      st3.push_back(Double_t(pmtSum->vdtime3[itr])+offset3);
      if(dt1==0||dt2==0||dt3==0) cout << "\t file " << treeNumber << " entry " << itr << " " << dt1 << ", " << dt2<<  ", " << dt3 << endl;
      if(dt1!=dt2||dt1!=dt3||dt2!=dt3) {
        if(dt1<dt2) {skip[0]+= dt2-dt1; ++nskip[0]; dt1 = Double_t(pmtSum->vdtime1[itr]+skip[0])+offset1;}
        if(dt1<dt3) {skip[0]+= dt3-dt1; ++nskip[0]; dt1 = Double_t(pmtSum->vdtime1[itr]+skip[0])+offset1;}
        if(dt2<dt1) {skip[1]+= dt1-dt2; ++nskip[1]; dt2 = Double_t(pmtSum->vdtime2[itr]+skip[1])+offset2;}
        if(dt2<dt3) {skip[1]+= dt3-dt2; ++nskip[1]; dt2 = Double_t(pmtSum->vdtime2[itr]+skip[1])+offset2;}
        if(dt3<dt1) {skip[2]+= dt1-dt3; ++nskip[2]; dt3 = Double_t(pmtSum->vdtime3[itr]+skip[2])+offset3;}
        if(dt3<dt2) {skip[2]+= dt2-dt3; ++nskip[2]; dt3 = Double_t(pmtSum->vdtime3[itr]+skip[2])+offset3;}
      }
      sdt1.push_back(dt1);
      sdt2.push_back(dt2);
      sdt3.push_back(dt3);

      // convert to microseconds
      ot1.push_back(Double_t(skip[0]-sskip[0])*0.008);
      ot2.push_back(Double_t(skip[1]-sskip[1])*0.008);
      ot3.push_back(Double_t(skip[2]-sskip[2])*0.008);

       
      if(dt1!=dt2||dt1!=dt3||dt2!=dt3) 
        printf(" ERROR!!! %u skip number %u,%u,%u :  %u %u %u \n",itr,skip[0],skip[1],skip[2],dt1,dt2,dt3); 
      unsigned gtrig = itr + entry*trigsize;
      ntTime->Fill(entry,gtrig,pmtSum->vtrig[itr],pmtSum->tZero[0],
          pmtSum->vrf1[itr],pmtSum->vrf2[itr],pmtSum->vrf3[itr],
          //pmtSum->vprompt1[itr],pmtSum->vprompt2[itr],pmtSum->vprompt3[itr],
          pmtSum->vdtime1[itr],pmtSum->vdtime2[itr],pmtSum->vdtime3[itr],
          dt1,dt2,dt3);
    }
    trigsize= pmtSum->vtrig.size();
    printf("  file %i skips %u %u %u skip number %u,%u,%u  \n",treeNumber,nskip[0],nskip[1],nskip[2],skip[0],skip[1],skip[2]);
  }

  //std::sort(st1.begin(), st1.end());
  //std::sort(st2.begin(), st2.end());
  //std::sort(st3.begin(), st3.end());
  for(unsigned it =0; it< st1.size(); ++it) {
    sortTrig = Double_t(it);
    sortTime1 = st1[it];
    sortTime2 = st2[it];
    sortTime3 = st3[it];
    sortDtime1 = sdt1[it]; 
    sortDtime2 = sdt2[it]; 
    sortDtime3 = sdt3[it]; 
    sortOff1 = ot1[it];
    sortOff2 = ot2[it];
    sortOff3 = ot3[it];
    if(st1[it]==0||st2[it]==0||st3[it]==0) cout << sortTrig << "  "<< st1[it] << "  " << st2[it] << "  " << st3[it] << endl;
    trSortTime->Fill();
  }

  TString name,title;
  TGraph *gtr1 = new TGraph(strig.size(),&strig[0],&st1[0]);
  name.Form("time1");
  title.Form("time PDS board 1 ");
  gtr1->SetNameTitle(name,title);
  gtr1->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  gtr1->GetHistogram()->GetYaxis()->SetTitle(" PDS time 1 ");
  gtr1->SetMarkerStyle(20); gtr1->SetMarkerSize(0.5);gtr1 ->SetMarkerColor(kBlack);


  TGraph *gtr2 = new TGraph(strig.size(),&strig[0],&st2[0]);
  name.Form("time2");
  title.Form("time PDS board 2 ");
  gtr2->SetNameTitle(name,title);
  gtr2->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  gtr2->GetHistogram()->GetYaxis()->SetTitle(" PDS time 2 ");
  gtr2->SetMarkerStyle(20); gtr2->SetMarkerSize(0.5);gtr2 ->SetMarkerColor(kBlue);

  TGraph *gtr3 = new TGraph(strig.size(),&strig[0],&st3[0]);
  name.Form("time3");
  title.Form("time PDS board 3 ");
  gtr3->SetNameTitle(name,title);
  gtr3->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  gtr3->GetHistogram()->GetYaxis()->SetTitle(" PDS time 3 ");
  gtr3->SetMarkerStyle(20); gtr3->SetMarkerSize(0.5);gtr3 ->SetMarkerColor(kRed);

  TCanvas *can1 = new TCanvas("t1-trig","t1-trig");
  gtr1->Draw();
  gtr2->Draw("same");
  gtr3->Draw("same");

  printf(" write file %s \n",outfile->GetName());
  outfile->Write();
} 
