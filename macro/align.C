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
#include "TPmtAlign.hxx"
  
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


void align()
{
  TString fullDirName;
  Double_t STEP= (pow(2.0,32) - pow (2.0,31))/125.0; // step in microseconds
  std::map<int,std::string> timeMap;
  std::vector<ULong_t> bsum1(32);
  std::vector<ULong_t> bsum2(32);
  std::vector<ULong_t> bsum3(32);
  std::vector<Int_t> fileTime;
  for(unsigned i=0; i<bsum1.size(); ++i) { bsum1[i]=0;bsum2[i]=0;bsum3[i]=0;}
  
  // get list of files
  printf("Low intensity runs range from 0731_1518 to 0731_2130. They are all in one day. step %E \n",STEP);
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
      fileTime.push_back(time);
      timeMap.insert ( std::pair<int,std::string>(time,fname) );
    }
  }

  printf(" total of files in %s is %lu  \n",fullDirName.Data(),timeMap.size()-1);
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
  // printout some events


  printf(" size %u \n",timeMap.size());
  /*
  unsigned mapCount=0;
  for (iter=timeMap.begin(); iter!=timeMap.end(); ++iter) {
    cout << " map entry " << mapCount++ << " => " << iter->second << endl;
  }
  */

  // open ouput file for alignment constants
  TString alignFileName = TString("align-")+sumTag+TString(".root");
  TFile *alignFile = new TFile(alignFileName,"recreate");
  printf(" opening output file %s \n",alignFileName.Data());
  TTree *alignTree = new TTree("alignTree","alignTree");
  TPmtAlign* pmtAlign  = new TPmtAlign();
  alignTree->Branch("pmtAlign",&pmtAlign);
  
  // open ouput file and make some histograms
  TString outputFileName = TString("sum-")+sumTag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  TNtuple *ntTrigSum = new TNtuple("ntTrigSum"," trigger sum ","time:n555:n5xx:n444:n4xx:n111:n1xx:n000:n0xx");
  TTree *trDiff = new TTree("trDiff"," diff ");
  Double_t trdRun, trdEvent,trdTrig, trdD1,trdD2,trdD3,trdS1,trdS2,trdS3;
  trDiff->Branch("run",&trdRun,"run/D");
  trDiff->Branch("ev",&trdEvent,"ev/D");
  trDiff->Branch("trig",&trdTrig,"trig/D");
  trDiff->Branch("d1",&trdD1,"d1/D");
  trDiff->Branch("d2",&trdD2,"d2/D");
  trDiff->Branch("d3",&trdD3,"d3/D");
  trDiff->Branch("s1",&trdS1,"s1/D");
  trDiff->Branch("s2",&trdS2,"s2/D");
  trDiff->Branch("s3",&trdS3,"s3/D");



  TTree *trSortTime = new TTree("sortTime","sort timing");
  TH1D* hBits1 = new TH1D("Bits1"," board 1 bits",32,0,32);
  TH1D* hBits2 = new TH1D("Bits2"," board 2 bits",32,0,32);
  TH1D* hBits3 = new TH1D("Bits3"," board 3 bits",32,0,32);
  TH1D* hDiff1 = new TH1D("Diff1"," board time differences ",200,-.50E8,.50E8);
  TH1D* hDiff2 = new TH1D("Diff2"," board time differences ",200,-.50E8,.50E8);
  TH1D* hDiff3 = new TH1D("Diff3"," board time differences ",200,-.50E8,.50E8);
  hDiff1->GetXaxis()->SetTitle(" diff in sequential board digitization time (micro-sec) ");
  hDiff2->GetXaxis()->SetTitle(" diff in sequentialboard digitization time (micro-sec) ");
  hDiff3->GetXaxis()->SetTitle(" diff in sequentialboard digitization time (micro-sec) ");

  TH1D* hADiff1 = new TH1D("ADiff1"," board time differences ",200,-.50E8,.50E8);
  TH1D* hADiff2 = new TH1D("ADiff2"," board time differences ",200,-.50E8,.50E8);
  TH1D* hADiff3 = new TH1D("ADiff3"," board time differences ",200,-.50E8,.50E8);
  hADiff1->GetXaxis()->SetTitle(" diff in sequential board digitization time (micro-sec) ");
  hADiff2->GetXaxis()->SetTitle(" diff in sequentialboard digitization time (micro-sec) ");
  hADiff3->GetXaxis()->SetTitle(" diff in sequentialboard digitization time (micro-sec) ");


  hDiff1->SetLineColor(kBlue);
  hDiff2->SetLineColor(kBlack);
  hDiff3->SetLineColor(kRed);
  hADiff1->SetLineColor(kBlue);
  hADiff2->SetLineColor(kBlack);
  hADiff3->SetLineColor(kRed);

  

  Double_t sortTime1; 
  Double_t sortTime2;
  Double_t sortTime3;
  Double_t sortShiftTime; 
  Double_t sortShiftRFTime; 
  Double_t sortShiftTimeDiff;
  Double_t sortOff1; 
  Double_t sortOff2;
  Double_t sortOff3;
  Double_t sortDiff1; 
  Double_t sortDiff2;
  Double_t sortDiff3;

  Double_t sortDiffRf1; 
  Double_t sortDiffRf2;
  Double_t sortDiffRf3;

  Double_t sortRTime1;
  Double_t sortRTime2;
  Double_t sortRTime3;

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


  Double_t sortCtime;

  Double_t sortRunShift1;
  Double_t sortRunShift2;
  Double_t sortRunShift3;

  trSortTime->Branch("time",&sortCtime,"time/D");
  trSortTime->Branch("run",&sortRun,"run/D");
  trSortTime->Branch("trig",&sortTrig,"trig/D");
  trSortTime->Branch("type",&sortType,"type/D");
  trSortTime->Branch("tzero",&sortTzero,"tzero/D");

  trSortTime->Branch("runshift1",&sortRunShift1,",runshift1/D"); // vectors runShift1,2,3
  trSortTime->Branch("runshift2",&sortRunShift2,",runshift2/D");
  trSortTime->Branch("runshift3",&sortRunShift3,",runshift3/D");

  trSortTime->Branch("rf1",&sortRf1,"rf1/D");
  trSortTime->Branch("rf2",&sortRf2,"rf2/D");
  trSortTime->Branch("rf3",&sortRf3,"rf3/D");

  trSortTime->Branch("orf1",&sortOrf1,"orf1/D");
  trSortTime->Branch("orf2",&sortOrf2,"orf2/D");
  trSortTime->Branch("orf3",&sortOrf3,"orf3/D");

  trSortTime->Branch("tr1",&sortRTime1,"tr1/D");
  trSortTime->Branch("tr2",&sortRTime2,"tr2/D");
  trSortTime->Branch("tr3",&sortRTime3,"tr3/D");

  trSortTime->Branch("t1",&sortTime1,"t1/D");  // is st1
  trSortTime->Branch("t2",&sortTime2,"t2/D");
  trSortTime->Branch("t3",&sortTime3,"t3/D");
  trSortTime->Branch("stime",&sortShiftTime,"stime/D"); // adjusted time
  trSortTime->Branch("difftime",&sortShiftTimeDiff,"difftime/D"); // adjusted time
  
  trSortTime->Branch("rftime",&sortShiftRFTime,"rftime/D"); // adusted RF 
  trSortTime->Branch("o1",&sortOff1,"off1/D");  // is ot1 incremental shifts
  trSortTime->Branch("o2",&sortOff2,"off2/D");
  trSortTime->Branch("o3",&sortOff3,"off3/D");

  trSortTime->Branch("diff1",&sortDiff1,"diff1/D");      // st1 seqential differences 
  trSortTime->Branch("diff2",&sortDiff2,"diff2/D");
  trSortTime->Branch("diff3",&sortDiff3,"diff3/D");
  trSortTime->Branch("diffRf1",&sortDiffRf1,"diffRf1/D");// rf1 sequential differences
  trSortTime->Branch("diffRf2",&sortDiffRf2,"diffRf2/D");
  trSortTime->Branch("diffRf3",&sortDiffRf3,"diffRf3/D");




  //trSortTime->Print();
  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
  
  std::vector<Double_t> strig;

  // raw times
  std::vector<UInt_t> vdt1;
  std::vector<UInt_t> vdt2;
  std::vector<UInt_t> vdt3;
 
  // continuous times
  std::vector<Double_t> st1;
  std::vector<Double_t> st2;
  std::vector<Double_t> st3;

  // shifted times
  std::vector<Double_t> sdt1;
  std::vector<Double_t> sdt2;
  std::vector<Double_t> sdt3;

  // start of run offsets 
  std::vector<ULong64_t> runShift1;
  std::vector<ULong64_t> runShift2;
  std::vector<ULong64_t> runShift3;

  // offsets 
  std::vector<Double_t> ot1;
  std::vector<Double_t> ot2;
  std::vector<Double_t> ot3;

  // shifted RF times
  std::vector<Double_t> rf1;
  std::vector<Double_t> rf2;
  std::vector<Double_t> rf3;

  //RF offsets
  std::vector<Double_t> orf1;
  std::vector<Double_t> orf2;
  std::vector<Double_t> orf3;

  // RF differences 
  std::vector<Double_t> drf1;
  std::vector<Double_t> drf2;
  std::vector<Double_t> drf3;

  std::vector<Double_t> srun;
  std::vector<Double_t> stype;
  std::vector<Double_t> stzero;
  std::vector<Double_t> stime;

  std::vector<unsigned> stree;
  std::vector<unsigned> sevent;

 

  Double_t offset[3]= {0,0,0};
  Double_t lastTime[3]={0,0,0};
  Double_t compTime =0;

  Int_t treeNumber=0;
  UInt_t trigsize=0;
  Double_t check[3]={0,0,0};

  pmtAlign->clear();
  
  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    if(sumTree->GetTreeNumber()!=treeNumber) {
      treeNumber=sumTree->GetTreeNumber();
      printf("\t opening tree number %i entry %u vtrig size %zu \n",treeNumber,entry, pmtSum->vtrig.size());
    }
    std::string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t min = pmtSum->getMin();
    Int_t seg = pmtSum->getSegment();
    double time = double(min)+double(60*24*day+60*24*30*(month-7)) - 60*24*20;
    
    ntTrigSum->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);

    //if(entry%100==0) printf("...entry %u tag %s  month %i day %i min %i seg %i time %0.f \n",entry,tag.c_str(),month,day,min,seg,time);
    pmtAlign->clear();
    pmtAlign->tag=tag;

    for(unsigned itr = 0; itr< pmtSum->vtrig.size(); ++ itr) {
      // calculate computer time
      compTime += pmtSum->vnano[itr]*1.0E-3; // nano to micro
      //if(itr%1000==0)  cout << " compTime " << itr << "  " <<  compTime << "  " <<  pmtSum->vnano[itr] << endl;
      stime.push_back(compTime);

      // align the boards at start of run to board 1
      if(itr==0) {
        Double_t delta1 = (Double_t(pmtSum->vdtime1[0])- lastTime[0])/125.0;
        Double_t delta2 = (Double_t(pmtSum->vdtime2[0])- lastTime[1])/125.0;
        Double_t delta3 = (Double_t(pmtSum->vdtime3[0])- lastTime[2])/125.0;
        offset[1] = offset[0] + delta1 - delta2;
        offset[2] = offset[0] + delta1 - delta3;
      }
          
      // check to see if we have reached overflow. convert everything to microseconds
      //if(itr>0&& pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr]) 
      //cout << " offset " << ( pmtSum->vdtime1[itr-1]-pmtSum->vdtime1[itr])/125 << endl;
      for(int icheck=0; icheck<3; ++icheck) check[icheck]=0.0;
      if(itr>0) {
        check[0]= Double_t(pmtSum->vdtime1[itr-1])- Double_t(pmtSum->vdtime1[itr]);
        check[1]= Double_t(pmtSum->vdtime2[itr-1])- Double_t(pmtSum->vdtime2[itr]);
        check[2]= Double_t(pmtSum->vdtime3[itr-1])- Double_t(pmtSum->vdtime3[itr]);
      }
      for(int icheck=0; icheck<3; ++icheck) if(check[icheck]> 1.0E9) {
        //cout << " BIG STEP  " << icheck << "  " << itr 
        //<< " check " << check[icheck] << " , " <<  pmtSum->vdtime1[itr-1] << "  " << pmtSum->vdtime1[itr] << endl; 
        offset[icheck] +=STEP;  
      }


      strig.push_back(Double_t(itr + treeNumber*(pmtSum->vtrig.size())));
      // bit sum
      stree.push_back(entry);
      sevent.push_back(itr);

      st1.push_back(Double_t(pmtSum->vdtime1[itr])/125.+offset[0]);
      st2.push_back(Double_t(pmtSum->vdtime2[itr])/125.+offset[1]);
      st3.push_back(Double_t(pmtSum->vdtime3[itr])/125.+offset[2]);

      runShift1.push_back(ULong64_t(offset[0]));
      runShift2.push_back(ULong64_t(offset[1]));
      runShift3.push_back(ULong64_t(offset[2]));

      pmtAlign->align0.push_back(offset[0]);
      pmtAlign->align1.push_back(offset[1]);
      pmtAlign->align2.push_back(offset[2]);

  
      unsigned ilast = st1.size() -1;
      //if(ilast>0) if( itr==0 )                    printf(" ------- %i %u %u %E %E %E %E \n",treeNumber,itr,ilast,st1[ilast],st1[ilast-1],
        //  Double_t(pmtSum->vdtime1[itr]),lastTime[0]);
      
      if(ilast>0) if( st1[ilast]-st1[ilast-1]<0 ) printf(" WOOFWOOF %i %u %u %E %E %E %E \n",treeNumber,itr,ilast,st1[ilast],st1[ilast-1],
          Double_t(pmtSum->vdtime1[itr]),lastTime[0]);

      if(itr>0) {
        hDiff1->Fill( (Double_t(pmtSum->vdtime1[itr]) - lastTime[0])/125.0 );
        hDiff2->Fill( (Double_t(pmtSum->vdtime2[itr]) - lastTime[1])/125.0 );
        hDiff3->Fill( (Double_t(pmtSum->vdtime3[itr]) - lastTime[2])/125.0 );
      }
      
      if(itr>0) {
         trdRun=  Double_t(treeNumber);
         trdEvent=Double_t(entry);
         trdTrig= Double_t(itr);
         trdD1=   Double_t(pmtSum->vdtime1[itr] - pmtSum->vdtime1[itr-1]);
         trdD2=   Double_t(pmtSum->vdtime2[itr] - pmtSum->vdtime2[itr-1]);
         trdD3=   Double_t(pmtSum->vdtime3[itr] - pmtSum->vdtime3[itr-1]);
         trdS1= st1[ilast] - st1[ilast-1];
         trdS2= st2[ilast] - st2[ilast-1];
         trdS3= st3[ilast] - st3[ilast-1];
         trDiff->Fill();
       }
      
      if(itr>0) {
        vdt1.push_back(pmtSum->vdtime1[itr]);
        vdt2.push_back(pmtSum->vdtime2[itr]);
        vdt3.push_back(pmtSum->vdtime3[itr]);
        rf1.push_back(Double_t(pmtSum->vrf1[itr])/250.+Double_t(pmtSum->vdtime1[itr])/125.+offset[0]);
        rf2.push_back(Double_t(pmtSum->vrf2[itr])/250.+Double_t(pmtSum->vdtime2[itr])/125.+offset[1]);
        rf3.push_back(Double_t(pmtSum->vrf3[itr])/250.+Double_t(pmtSum->vdtime3[itr])/125.+offset[2]);
      } else {
        vdt1.push_back(0);
        vdt2.push_back(0);
        vdt3.push_back(0);
        rf1.push_back(0);
        rf2.push_back(0);
        rf3.push_back(0);
      }


      srun.push_back(treeNumber);
      stype.push_back(pmtSum->vtrig[itr]);
      stzero.push_back(pmtSum->tZero[0]);

      // store the last times 
      lastTime[0]=Double_t(pmtSum->vdtime1[itr]);
      lastTime[1]=Double_t(pmtSum->vdtime2[itr]);
      lastTime[2]=Double_t(pmtSum->vdtime3[itr]);
    }
    alignTree->Fill();
  }
  for(unsigned it =1; it< st1.size(); ++it) if((st1[it]-st1[it-1])<0) printf(" BLAHBLAH %u %E %E \n",it,st1[it],st1[it-1]);
 
  /*
  printf(" bit sum 1 \n");
  for(unsigned i=0; i<bsum1.size(); ++i) cout << i << " " << bsum1[i] << " ; ";
  cout << endl;
  */

  // now calculate offsets
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
    if(it>0) {
      hADiff1->Fill(Double_t(st1[it]-st1[it-1]));
      hADiff2->Fill(Double_t(st2[it]-st2[it-1]));
      hADiff3->Fill(Double_t(st3[it]-st3[it-1]));
    }

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
    
    
    if(it%5000!=0) {
      ot1.push_back(Double_t(dshift[0]));
      ot2.push_back(Double_t(dshift[1]));
      ot3.push_back(Double_t(dshift[2]));
    } else {
      ot1.push_back(0); ot2.push_back(0); ot3.push_back(0);
    }

    if(dt1!=dt2||dt1!=dt3||dt2!=dt3) 
      printf(" ERROR!!! %u skip number %E,%E,%E :  %E %E %E \n",it,shift[0],shift[1],shift[2],dt1,dt2,dt3);

    // increment shifts
    for(int iskip=0; iskip<3; ++iskip) shift[iskip]+=dshift[iskip];

     if(it%5000!=0) {
       sdt1.push_back(dt1);
       sdt2.push_back(dt2);
       sdt3.push_back(dt3);
     } else {
      sdt1.push_back(0); sdt2.push_back(0); sdt3.push_back(0);
     }

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

     if(it%5000!=0) {
       drf1.push_back(crf1);
       drf2.push_back(crf2);
       drf3.push_back(crf3);
       orf1.push_back(Double_t(drfshift[0]));
       orf2.push_back(Double_t(drfshift[1]));
       orf3.push_back(Double_t(drfshift[2]));
     } else {
      drf1.push_back(0); drf2.push_back(0); drf3.push_back(0);
      orf1.push_back(0); orf2.push_back(0); orf3.push_back(0);
     }

    if(crf1!=crf2||crf1!=crf3||crf2!=crf3) 
      printf(" ERROR!!! %u skip number %E,%E,%E :  %E %E %E \n",it,drfshift[0],drfshift[1],drfshift[2],crf1,crf2,crf3);

  }

  // fill tree
  for(unsigned it =0; it< st1.size(); ++it) {
    bitSum(vdt1[it],bsum1);
    bitSum(vdt2[it],bsum2);
    bitSum(vdt3[it],bsum3);

    sortRTime1= Double_t(vdt1[it]);
    sortRTime2= Double_t(vdt2[it]);
    sortRTime3= Double_t(vdt3[it]);

    sortRunShift1=Double_t(runShift1[it]);
    sortRunShift2=Double_t(runShift2[it]);
    sortRunShift3=Double_t(runShift3[it]);
    
    sortDiff1=0;
    sortDiff2=0;
    sortDiff3=0;
    bool skip = false;
    if(it>0) skip = it%5000==0;
    if(it>0) {
      Double_t sdiff = st1[it]-st1[it-1];
      if(sdiff>1.0E5&&it<1) printf(" skip %u %E \n",it,sdiff);
    }
    if(it>0&&!skip) {
      sortDiff1=(st1[it]-st1[it-1]);
      if(sortDiff1<0) printf(" BLAHWOOF  %u %E %E %E \n",it,sortDiff1,st1[it],st1[it-1]);
      
      sortDiff2=(st2[it]-st2[it-1]);
      sortDiff3=(st3[it]-st3[it-1]);
      //if(st1[it-1]>st1[it]) bitSum(vdt1[it],bsum1);
      //if(st2[it-1]>st2[it]) bitSum(vdt2[it],bsum2);
      //if(st3[it-1]>st3[it]) bitSum(vdt3[it],bsum3);
    }
    sortDiffRf1=0;
    sortDiffRf2=0;
    sortDiffRf3=0;
    if(it>0) skip = rf1[it]==0||rf1[it-1]==0;
    if(it>0&&!skip) {
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
    sortShiftTime = sdt1[it]; 
    sortShiftRFTime = drf1[it]; 
    
    if(it>0) sortShiftTimeDiff = sdt1[it] - sdt1[it-1];

    sortOff1 = ot1[it];
    sortOff2 = ot2[it];
    sortOff3 = ot3[it];

    sortOrf1 = orf1[it];
    sortOrf2 = orf2[it];
    sortOrf3 = orf3[it];
    
    if(st1[it]==0||st2[it]==0||st3[it]==0) cout << sortTrig << "  "<< st1[it] << "  " << st2[it] << "  " << st3[it] << endl;
    trSortTime->Fill();

    // printout some events
    
    //for(unsigned ij = 0; ij<fileTime.size() ; ++ij) cout << " file " << ij << " time "  
    //  << fileTime[ij] << " is " << timeMap[fileTime[ij]]<< endl;
    //printf(" size %u \n",timeMap.size());
    //for(ifile = 0; ifile<timeMap.size(); ++ifile) cout << " file " << ifile << " is " << timeMap[ifile] << endl;
    if(sortDiff3>1E6) {
      int utree = int(stree[it]);
      std::string sfile = timeMap[fileTime[utree]];
      unsigned ilocal = it%5000;
      printf(" LARGESHIFT3  event %u,%u tree %u  %s %E \n",ilocal,sevent[it],utree,sfile.c_str(),sortDiff3);
    }
     
  }

  for(unsigned i=0; i<bsum1.size(); ++i)  hBits1->SetBinContent(i+1,bsum1[i]);
  for(unsigned i=0; i<bsum2.size(); ++i)  hBits2->SetBinContent(i+1,bsum2[i]);
  for(unsigned i=0; i<bsum3.size(); ++i)  hBits3->SetBinContent(i+1,bsum3[i]);

  alignFile->Write();
  printf(" write file %s \n",alignFile->GetName());
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


  TGraph *grf1 = new TGraph(strig.size(),&strig[0],&orf1[0]);
  name.Form("rf-1");
  title.Form("RF  PDS board 1 ");
  grf1->SetNameTitle(name,title);
  grf1->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  grf1->GetHistogram()->GetYaxis()->SetTitle(" RF  time ");
  grf1->SetMarkerStyle(20); grf1->SetMarkerSize(0.5);grf1 ->SetMarkerColor(kBlue);


  TGraph *grf2 = new TGraph(strig.size(),&strig[0],&orf2[0]);
  name.Form("rf-2");
  title.Form("RF  PDS board 2 ");
  grf2->SetNameTitle(name,title);
  grf2->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  grf2->GetHistogram()->GetYaxis()->SetTitle(" RF  time ");
  grf2->SetMarkerStyle(20); grf2->SetMarkerSize(0.5);grf2 ->SetMarkerColor(kBlue);

  TGraph *grf3 = new TGraph(strig.size(),&strig[0],&orf3[0]);
  name.Form("rf-3");
  title.Form("RF  PDS board 3 ");
  grf3->SetNameTitle(name,title);
  grf3->GetHistogram()->GetXaxis()->SetTitle(" trigger number ");
  grf3->GetHistogram()->GetYaxis()->SetTitle(" RF  time ");
  grf3->SetMarkerStyle(20); grf3->SetMarkerSize(0.5);grf3->SetMarkerColor(kBlue);
  



 // TCanvas *can2 = new TCanvas("offset2-trig","offset2-trig"); gtr2->Draw();
 // TCanvas *can3 = new TCanvas("offset3-trig","offset3-trig"); gtr3->Draw();\

  TCanvas *canrf1 = new TCanvas("rf1-trig","rf1-trig"); grf1->Draw();
  TCanvas *canrf2 = new TCanvas("rf2-trig","rf2-trig"); grf2->Draw();
  TCanvas *canrf3 = new TCanvas("rf3-trig","rf3-trig"); grf3->Draw();

} 
