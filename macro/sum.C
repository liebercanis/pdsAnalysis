#include <vector>
#include <map>
#include <bitset> 
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
const Double_t STEP=34359738.368; // this is 2^32/125 step in micro-seconds
// returns -1 if pmt does not exist 
// populate 3 boards, each from channel 0-6.  Channel 7 is the RF pulse. 
// valid pmt are 0 to 20, RF channels are 21,22,23
std::string tag;
void getTag(std::string fname) { tag = fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")); return;}
Int_t getMonth() { return atoi(tag.substr(0,2).c_str()) ;}
Int_t getDay() { return atoi(tag.substr(3,2).c_str()) ;}
Int_t getMin() { return atoi(tag.substr(6,4).c_str()) ;}
Int_t getSegment() { return atoi(tag.substr(11,tag.find(".") -1  - 11).c_str());}

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


void sum()
{
  TString fullDirName;
  std::map<int,std::string> timeMap;
  std::vector<ULong_t> bsum1(32);
  std::vector<ULong_t> bsum2(32);
  std::vector<ULong_t> bsum3(32);
  for(unsigned i=0; i<bsum1.size(); ++i) { bsum1[i]=0;bsum2[i]=0;bsum3[i]=0;}
  
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
    if ( strstr(fname.c_str(), "pmtSummary2_" )==NULL ) continue;
    if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
    //cout << fname << endl;
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
  TTree *trSortTime = new TTree("sortTime","sort timing");
  TH1D* hBits1 = new TH1D("Bits1"," board 1 bits",32,0,32);
  TH1D* hBits2 = new TH1D("Bits2"," board 2 bits",32,0,32);
  TH1D* hBits3 = new TH1D("Bits3"," board 3 bits",32,0,32);
  Double_t sortTime1; 
  Double_t sortTime2;
  Double_t sortTime3;
  Double_t sortDtime1; 
  Double_t sortDtime2;
  Double_t sortDtime3;
  Double_t sortOff1; 
  Double_t sortOff2;
  Double_t sortOff3;
  Double_t sortDiff1; 
  Double_t sortDiff2;
  Double_t sortDiff3;

  Double_t sortDiffRf1; 
  Double_t sortDiffRf2;
  Double_t sortDiffRf3;


  Double_t sortRun;
  Double_t sortTrig;
  Double_t sortType;
  Double_t sortTzero;
  Double_t sortRf1;
  Double_t sortRf2;
  Double_t sortRf3;
  
  Double_t sortOrf1;
  Double_t sortOrf2;
  Double_t sortOrf3;

  Double_t sortDrf1;
  Double_t sortDrf2;
  Double_t sortDrf3;



  Double_t sortCtime;
  trSortTime->Branch("time",&sortCtime,"time/D");
  trSortTime->Branch("run",&sortRun,"run/D");
  trSortTime->Branch("trig",&sortTrig,"trig/D");
  trSortTime->Branch("type",&sortType,"type/D");
  trSortTime->Branch("tzero",&sortTzero,"tzero/D");
  trSortTime->Branch("rf1",&sortRf1,"rf1/D");
  trSortTime->Branch("rf2",&sortRf2,"rf2/D");
  trSortTime->Branch("rf3",&sortRf3,"rf3/D");

  trSortTime->Branch("drf1",&sortDrf1,"drf1/D");
  trSortTime->Branch("drf2",&sortDrf2,"drf2/D");
  trSortTime->Branch("drf3",&sortDrf3,"drf3/D");
  trSortTime->Branch("orf1",&sortOrf1,"orf1/D");
  trSortTime->Branch("orf2",&sortOrf2,"orf2/D");
  trSortTime->Branch("orf3",&sortOrf3,"orf3/D");


  trSortTime->Branch("t1",&sortTime1,"t1/D");
  trSortTime->Branch("t2",&sortTime2,"t2/D");
  trSortTime->Branch("t3",&sortTime3,"t3/D");
  trSortTime->Branch("dt1",&sortDtime1,"dt1/D");
  trSortTime->Branch("dt2",&sortDtime2,"dt2/D");
  trSortTime->Branch("dt3",&sortDtime3,"dt3/D");
  trSortTime->Branch("o1",&sortOff1,"off1/D");
  trSortTime->Branch("o2",&sortOff2,"off2/D");
  trSortTime->Branch("o3",&sortOff3,"off3/D");

  trSortTime->Branch("diff1",&sortDiff1,"diff1/D");
  trSortTime->Branch("diff2",&sortDiff2,"diff2/D");
  trSortTime->Branch("diff3",&sortDiff3,"diff3/D");
  trSortTime->Branch("diffRf1",&sortDiffRf1,"diffRf1/D");
  trSortTime->Branch("diffRf2",&sortDiffRf2,"diffRf2/D");
  trSortTime->Branch("diffRf3",&sortDiffRf3,"diffRf3/D");




  //trSortTime->Print();
  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
  
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

  std::vector<Double_t> rf1;
  std::vector<Double_t> rf2;
  std::vector<Double_t> rf3;

  std::vector<Double_t> orf1;
  std::vector<Double_t> orf2;
  std::vector<Double_t> orf3;

  std::vector<Double_t> drf1;
  std::vector<Double_t> drf2;
  std::vector<Double_t> drf3;

  std::vector<Double_t> srun;
  std::vector<Double_t> stype;
  std::vector<Double_t> stzero;
  std::vector<Double_t> stime;
  std::vector<UInt_t> vdt1;
  std::vector<UInt_t> vdt2;
  std::vector<Uint_t> vdt3;
 


  Double_t offset[3]= {0,0,0};
  Double_t compTime =0;

  Int_t treeNumber=0;
  UInt_t trigsize=0;
  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    if(sumTree->GetTreeNumber()!=treeNumber) {
      treeNumber=sumTree->GetTreeNumber();
      printf("opening tree number %i \n",treeNumber);
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
      // calculate computer time
      compTime += pmtSum->vnano[itr]*1.0E-3; // nano to micro
      //if(itr<100) cout << itr << "  " <<  compTime << "  " <<  pmtSum->vnano[itr] << endl;
      stime.push_back(compTime);

      // check to see if we have reached overflow. convert everything to microseconds
      //if(itr>0&& pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr]) 
      //cout << " offset " << ( pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr])/125 << endl; 
      if(itr>0&& pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr] > 1E6) offset[0] += STEP;
      if(itr>0&& pmtSum->vdtime2[itr-1]-pmtSum->vdtime2[itr] > 1E6) offset[1] += STEP;
      if(itr>0&& pmtSum->vdtime3[itr-1]-pmtSum->vdtime3[itr] > 1E6) offset[2] += STEP;
      strig.push_back(Double_t(itr + treeNumber*(pmtSum->vtrig.size())));
      // bit sum
      st1.push_back(Double_t(pmtSum->vdtime1[itr])/125.+offset[0]);
      st2.push_back(Double_t(pmtSum->vdtime2[itr])/125.+offset[1]);
      st3.push_back(Double_t(pmtSum->vdtime3[itr])/125.+offset[2]);

      vdt1.push_back(pmtSum->vdtime1[itr]);
      vdt2.push_back(pmtSum->vdtime2[itr]);
      vdt3.push_back(pmtSum->vdtime3[itr]);

      rf1.push_back(Double_t(pmtSum->vrf1[itr])/250.+Double_t(pmtSum->vdtime1[itr])/125.+offset[0]);
      rf2.push_back(Double_t(pmtSum->vrf2[itr])/250.+Double_t(pmtSum->vdtime2[itr])/125.+offset[1]);
      rf3.push_back(Double_t(pmtSum->vrf3[itr])/250.+Double_t(pmtSum->vdtime3[itr])/125.+offset[2]);
      srun.push_back(treeNumber);
      stype.push_back(pmtSum->vtrig[itr]);
      stzero.push_back(pmtSum->tZero[0]);
    }
  }
  
  /*
  printf(" bit sum 1 \n");
  for(unsigned i=0; i<bsum1.size(); ++i) cout << i << " " << bsum1[i] << " ; ";
  cout << endl;
  */

  // no calculate offsets
  Double_t shift[3]={0,0,0};
  Double_t dshift[3]={0,0,0};
  Double_t rfshift[3]={0,0,0};
  Double_t drfshift[3]={0,0,0};

  UInt_t nshift[3]={0,0,0};

  Double_t dt1,dt2,dt3;
  ot1.push_back(0);
  ot2.push_back(0);
  ot3.push_back(0);

  for(unsigned it =0; it< st1.size(); ++it) {

    dt1 = st1[it] + shift[0];
    dt2 = st2[it] + shift[1];
    dt3 = st3[it] + shift[2];

     // zero incremental shifts
    for(int iskip=0; iskip<3; ++iskip) dshift[iskip] = 0; 

    // align everything with board 1
    if(dt1!=dt2||dt1!=dt3||dt2!=dt3) {
      //if(dt1<dt2) {dshift[0]+= dt2-dt1; ++nshift[0]; dt1 += dt2-dt1;}
      //if(dt1<dt3) {dshift[0]+= dt3-dt1; ++nshift[0]; dt1 += dt3-dt1;}
      if(dt1<dt2) {dshift[1]-= dt2-dt1; ++nshift[0]; dt2 -= dt2-dt1;}
      if(dt1<dt3) {dshift[2]-= dt3-dt1; ++nshift[0]; dt3 -= dt3-dt1;}
      if(dt2<dt1) {dshift[1]+= dt1-dt2; ++nshift[1]; dt2 += dt1-dt2;}
      if(dt2<dt3) {dshift[1]+= dt3-dt2; ++nshift[1]; dt2 += dt3-dt2;}
      if(dt3<dt1) {dshift[2]+= dt1-dt3; ++nshift[2]; dt3 += dt1-dt3;}
      if(dt3<dt2) {dshift[2]+= dt2-dt3; ++nshift[2]; dt3 += dt2-dt3;}
    }
   
    /*
    Double_t ddt1 = st1[it];
    Double_t ddt2=  st2[it];
    Double_t ddt3 = st3[it];
    if(ddt1==0||ddt2==0||ddt3==0) cout << "\t file " << treeNumber << " entry " << it << " " << ddt1 << ", " << ddt2<<  ", " << ddt3 << endl;
    if(ddt1!=ddt2||ddt1!=ddt3||ddt2!=ddt3) {
      if(ddt1<ddt2) ddshift[0]+= ddt2-ddt1;
      if(ddt1<ddt3) ddshift[0]+= ddt3-ddt1; 
      if(ddt2<ddt1) ddshift[1]+= ddt1-ddt2;
      if(ddt2<ddt3) ddshift[1]+= ddt3-ddt2;
      if(ddt3<ddt1) ddshift[2]+= ddt1-ddt3;
      if(ddt3<ddt2) ddshift[2]+= ddt2-ddt3; 
    }
    */

 
    // save incremental shifts 
    // .. convert to microseconds by divide by 125 MHz
    ot1.push_back(Double_t(dshift[0]));
    ot2.push_back(Double_t(dshift[1]));
    ot3.push_back(Double_t(dshift[2]));

    if(dt1!=dt2||dt1!=dt3||dt2!=dt3) 
      printf(" ERROR!!! %u skip number %E,%E,%E :  %E %E %E \n",it,shift[0],shift[1],shift[2],dt1,dt2,dt3);

    // increment shifts
    for(int iskip=0; iskip<3; ++iskip) shift[iskip]+=dshift[iskip];

    sdt1.push_back(dt1);
    sdt2.push_back(dt2);
    sdt3.push_back(dt3);

    // do the same for RF pulses align everything with board 1
    //     // zero incremental shifts
    for(int iskip=0; iskip<3; ++iskip) drfshift[iskip] = 0;


    Double_t crf1 = rf1[it] + rfshift[0];
    Double_t crf2 = rf2[it] + rfshift[1];
    Double_t crf3 = rf3[it] + rfshift[2];

    if(crf1!=crf2||crf1!=crf3||crf2!=crf3) {
      //if(crf1<crf2) {dshift[0]+= crf2-crf1; ++nshift[0]; crf1 += crf2-crf1;}
      //if(crf1<crf3) {dshift[0]+= crf3-crf1; ++nshift[0]; crf1 += crf3-crf1;}
      if(crf1<crf2) {drfshift[1]-= crf2-crf1; crf2 -= crf2-crf1;}
      if(crf1<crf3) {drfshift[2]-= crf3-crf1; crf3 -= crf3-crf1;}
      if(crf2<crf1) {drfshift[1]+= crf1-crf2; crf2 += crf1-crf2;}
      if(crf2<crf3) {drfshift[1]+= crf3-crf2; crf2 += crf3-crf2;}
      if(crf3<crf1) {drfshift[2]+= crf1-crf3; crf3 += crf1-crf3;}
      if(crf3<crf2) {drfshift[2]+= crf2-crf3; crf3 += crf2-crf3;}
    }
    // increment shifts
    for(int iskip=0; iskip<3; ++iskip) rfshift[iskip]+=drfshift[iskip];

    drf1.push_back(crf1);
    drf2.push_back(crf2);
    drf3.push_back(crf3);
    orf1.push_back(Double_t(drfshift[0]));
    orf2.push_back(Double_t(drfshift[1]));
    orf3.push_back(Double_t(drfshift[2]));

    if(crf1!=crf2||crf1!=crf3||crf2!=crf3) 
      printf(" ERROR!!! %u skip number %E,%E,%E :  %E %E %E \n",it,drfshift[0],drfshift[1],drfshift[2],crf1,crf2,crf3);

  }

  // fill tree
  for(unsigned it =0; it< st1.size(); ++it) {
    sortDiff1=0;
    sortDiff2=0;
    sortDiff3=0;
    if(it>0) {
      sortDiff1=(st1[it]-st1[it-1]);
      sortDiff2=(st2[it]-st2[it-1]);
      sortDiff3=(st3[it]-st3[it-1]);
      if(st1[it-1]>st1[it]) bitSum(vdt1[it],bsum1);
      if(st2[it-1]>st2[it]) bitSum(vdt2[it],bsum2);
      if(st3[it-1]>st3[it]) bitSum(vdt3[it],bsum3);
    }
    sortDiffRf1=0;
    sortDiffRf2=0;
    sortDiffRf3=0;
    if(it>0) {
      sortDiffRf1=(rf1[it]-rf1[it-1]);
      sortDiffRf2=(rf2[it]-rf2[it-1]);
      sortDiffRf3=(rf3[it]-rf3[it-1]);
    }

    sortTrig = Double_t(it);
    sortCtime=stime[it];
    sortRun=srun[it];
    sortType=stype[it];
    sortTzero=stzero[it];
    sortRf1=rf1[it];
    sortRf2=rf2[it];
    sortRf3=rf3[it];
    sortTime1 = st1[it];
    sortTime2 = st2[it];
    sortTime3 = st3[it];
    sortDtime1 = sdt1[it]; 
    sortDtime2 = sdt2[it]; 
    sortDtime3 = sdt3[it]; 
    sortOff1 = ot1[it];
    sortOff2 = ot2[it];
    sortOff3 = ot3[it];

    sortDrf1=drf1[it];
    sortDrf2=drf2[it];
    sortDrf3=drf3[it];
    sortOrf1 = orf1[it];
    sortOrf2 = orf2[it];
    sortOrf3 = orf3[it];
    
    if(st1[it]==0||st2[it]==0||st3[it]==0) cout << sortTrig << "  "<< st1[it] << "  " << st2[it] << "  " << st3[it] << endl;
    trSortTime->Fill();
  }

  for(unsigned i=0; i<bsum1.size(); ++i)  hBits1->SetBinContent(i+1,bsum1[i]);
  for(unsigned i=0; i<bsum2.size(); ++i)  hBits2->SetBinContent(i+1,bsum2[i]);
  for(unsigned i=0; i<bsum3.size(); ++i)  hBits3->SetBinContent(i+1,bsum3[i]);

  outfile->Write();
  printf(" write file %s \n",outfile->GetName());
  return;

  TString name,title;

  TGraph *gtr2 = new TGraph(strig.size(),&strig[0],&ot2[0]);
  name.Form("time2");
  title.Form("time offset PDS board 2 ");
  gtr2->SetNameTitle(name,title);
  gtr2->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  gtr2->GetHistogram()->GetYaxis()->SetTitle("  offset mu-sec 2");
  gtr2->SetMarkerStyle(20); gtr2->SetMarkerSize(0.5); gtr2 ->SetMarkerColor(kBlue);

  TGraph *gtr3 = new TGraph(strig.size(),&strig[0],&ot3[0]);
  name.Form("time3");
  title.Form("time offset  PDS board 3 ");
  gtr3->SetNameTitle(name,title);
  gtr3->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  gtr3->GetHistogram()->GetYaxis()->SetTitle(" offset mu-sec 3");
  gtr3->SetMarkerStyle(20); gtr3->SetMarkerSize(0.5);gtr3 ->SetMarkerColor(kRed);


  TGraph *grf2 = new TGraph(strig.size(),&strig[0],&orf2[0]);
  name.Form("rf-offset2");
  title.Form("RF offset  PDS board 2 ");
  grf2->SetNameTitle(name,title);
  grf2->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  grf2->GetHistogram()->GetYaxis()->SetTitle(" RF offset mu-sec 2");
  grf2->SetMarkerStyle(20); grf2->SetMarkerSize(0.5);grf2 ->SetMarkerColor(kBlue);

  TGraph *grf3 = new TGraph(strig.size(),&strig[0],&orf3[0]);
  name.Form("rf-offset3");
  title.Form("RF offset  PDS board 3 ");
  grf3->SetNameTitle(name,title);
  grf3->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  grf3->GetHistogram()->GetYaxis()->SetTitle(" RF offset mu-sec 3");
  grf3->SetMarkerStyle(20); grf3->SetMarkerSize(0.5);grf3 ->SetMarkerColor(kRed);



  TCanvas *can2 = new TCanvas("offset2-trig","offset2-trig"); gtr2->Draw();
  TCanvas *can3 = new TCanvas("offset3-trig","offset3-trig"); gtr3->Draw();

  TCanvas *canrf2 = new TCanvas("rf-offset2-trig","rf-offset2-trig"); grf2->Draw();
  TCanvas *canrf3 = new TCanvas("rf-offset3-trig","rf-offset3-trig"); grf3->Draw();

} 
