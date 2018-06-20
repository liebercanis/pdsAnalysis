// clock info 
typedef struct {
  Int_t  run;
  Int_t  event;
  Int_t  compSec;
  Int_t rf0;
  Int_t rf1;
  Int_t rf2;
  Long64_t compNano;
  UInt_t  dt0,dt1,dt2;   // caen digitizer time 
  Double_t tPrompt; // one for each board
  Double_t tPromptToRF;
} CLOCK;

static CLOCK clk;

TTree *tree;

typedef struct {
  Long64_t ip;
  Long64_t ev;
  Long64_t dt;
  double bt;
  double dbt;
  double jt;
  double rf;
} TIMEMAP;

static TIMEMAP tmap0;
static TIMEMAP tmap1;
static TIMEMAP tmap2;


enum {NB=3};
const double jump = 195.0E6;
const double jumpMax = 201.0E6;

std::vector<ULong_t> bsum2(32);
std::vector<ULong_t> bsumAll2(32);
std::vector<ULong_t> bsumAll0(32);

void bitSum(UInt_t word, std::vector<ULong_t> &bsum) 
{
  ULong_t lword = ULong_t(word);
  std::bitset<32> foo(lword);
  for(unsigned i=0; i<bsum.size(); ++i) if(foo.test(i)) ++bsum[i];
}

void rclock(Int_t ifile=0, Long64_t max=0, Long64_t first=0)
{
  TString fileInputName;
  if(ifile==0) fileInputName=TString("check_lowAnaNoAlign-0-0.root");
  else {
    cout << " invalid "  << ifile << endl;
    return;
  }

  // open the file
  TFile *finput = new TFile(fileInputName);
  if (!finput || !finput->IsOpen()) {
    cout << " could not find file " << fileInputName << endl;
    return;
  }
  cout << " opened file " << fileInputName << endl;
  finput->GetObject("TClk",tree);
  Long64_t fsize = tree->GetEntriesFast();
  cout << " size of TClk tree " << fsize << endl;

  tree->SetBranchAddress("clk",&clk);

  for(unsigned i=0; i<bsum2.size(); ++i) bsum2[i]=0;

  double MAXVAL = pow(2.0,32);
  double STEP = pow(2.0,31);
  double HALF = pow(2.0,30);

  printf(" time cuts step %.0f half %.0f %.0f  %X %X %X \n",STEP, HALF, MAXVAL, unsigned(STEP) , unsigned(HALF), unsigned(MAXVAL) );
  double btime[NB];
  double blast[NB];
  double pdsLast;
  double boff[NB];
  double bslope[NB];
  double dtime[NB];
  double rft[NB];
  double dlast[NB];
  double delta[NB];
  int nskips[NB];
  double bstart[NB];
  double jumpTime[NB];
  double rfTime[NB];
  int njump[NB];
  int nRF[NB];
  int jcount[NB];
  int rfcount[NB];
  int rfcountEvent[NB];
  int isRF[NB];
  for(int ib=0; ib<NB; ++ib) { 
    btime[ib]=0;  blast[ib]=0;  boff[ib]=0; dtime[ib]=0; dlast[ib]=0;  delta[ib]=0; nskips[ib]=0; bstart[ib]=0; bslope[ib]=0;
    rft[ib]=0;
    jumpTime[ib]=0;
    rfTime[ib]=0;
    njump[ib]=0;
    nRF[ib]=0;
    jcount[ib]=0;
    rfcount[ib]=0;
    rfcountEvent[ib]=0;
    isRF[ib]=0;
  }

  bslope[1]=1.0+3.30E-6; // 
  bslope[2]=1.0+1.6E-6; // 
  

  TString outFileName;
  outFileName.Form("clocks-%i.root",ifile);

 
  TFile *fout = new TFile(outFileName.Data(),"RECREATE");
  fout->cd();

  TH1D* hBits2 = new TH1D("Bits2"," board 2 bits",32,0,32);
  TH1D* hBitsAll2 = new TH1D("BitsAll2"," board 2 bits",32,0,32);
  TH1D* hBitsAll0 = new TH1D("BitsAll0"," board 2 bits",32,0,32);
  TH1D* hDiff0 = new TH1D("Diff0"," time since last jump board 0 ns",800,0.,80000.);
  TH1D* hDiff1 = new TH1D("Diff1"," time since last jump board 1 ns",800,0.,80000.);
  TH1D* hDiff2 = new TH1D("Diff2"," time since last jump board 2 ns",800,0.,80000.);


  TTree *bClock = new TTree("bclk"," PDS clocks ");
  boardClock *bclk = new boardClock;
  bClock->Branch("c",&bclk);
  //,"run/I:event/I:sec/I:rf0/I:rf1/I:rf2/I:nrf0/I:nrf1/I:nrf2/I:nj0/I:nj1/I:nj2/I:nano/l:dt0/i:dt1/i:dt2/i:tprompt/D:tPromptToRF/D:bt0/D:bt1/D:bt2/D:db0/D:db1/D:db2/D:rftime0/D:rftime1/D:rftime2/D:jtime0/D:jtime1/D:jtime2/D");
  bClock->Print("all");

  //TNtuple *ntSort  = new TNtuple("ntsort"," sorted clocks","ev:bt0:bt1:bt2");
  TNtuple *ntStep  = new TNtuple("ntstep"," steps ","ev:steps:s0:s1:s2:delta0:delta1:delta2:bt0:bt1:bt2");
  TNtuple *ntJump  = new TNtuple("ntjump"," jumps ","ev:bt0:delta0:delta1:delta2:n0:n1:n2:djt0:djt1:djt2");

  TTree *ntMap0 = new TTree("ntmap0"," time map 0 ");
  ntMap0->Branch("map0",&tmap0,"ip/L:ev/L:dt/L:bt/D:dbt/D:jt/D:rf/D");
  TTree *ntMap1 = new TTree("ntmap1"," time map 1 ");
  ntMap1->Branch("map1",&tmap1,"ip/L:ev/L:dt/L:bt/D:dbt/D:jt/D:rf/D");
  TTree *ntMap2 = new TTree("ntmap2"," time map 12 ");
  ntMap2->Branch("map2",&tmap2,"ip/L:ev/L:dt/L:bt/D:dbt/D:jt/D:rf/D");

  ntMap0->Print("all");

  TNtuple *ntEvCount  = new TNtuple("evCount","   event count ","ev:rf0:rf1:rf2:j0:j1:j2");
  TNtuple *ntRunCount  = new TNtuple("runCount"," run count  ","run:rf0:rf1:rf2:j0:j1:j2");
  TNtuple *ntRunTime = new TNtuple("runTime"," run times ","run:start:end:dstart:dend");
  TNtuple *ntRFevent = new TNtuple("RFevent"," RF event numbers ","ev0:ev1:ev2:n0:n1:n2");
  TNtuple *ntGap = new TNtuple("gap"," gap-rF correlation ","ev:bt0:n0:n1:n2:r0:r1:r2:jt0:rft0:jt1:rft1:jt2:rft2");

  std::vector<double> c0;
  std::vector<double> c1;
  std::vector<double> c2;


  std::vector<double> rt0;
  std::vector<double> rt1;
  std::vector<double> rt2;
  std::vector<double> jtime0;
  std::vector<double> jtime1;
  std::vector<double> jtime2;
  std::vector<double> rftime0;
  std::vector<double> rftime1;
  std::vector<double> rftime2;


  std::vector<double> jcountRun0;
  std::vector<double> jcountRun1;
  std::vector<double> jcountRun2;
  std::vector<double> rfcountRun0;
  std::vector<double> rfcountRun1;
  std::vector<double> rfcountRun2;
  std::vector<double> runCount;
  int runNumber=0;

  std::vector<int> isrf0;

  std::map<double,Long64_t> timeMap0; 
  std::map<double,Long64_t> timeMap1;
  std::map<double,Long64_t> timeMap2;

  Long64_t nentries = max;
  if(nentries==0) nentries = fsize;

  std::bitset<3> steps;

  std::vector<Long64_t> rfev0;
  std::vector<Long64_t> rfev1;
  std::vector<Long64_t> rfev2;

  std::vector<Long64_t> rfevCount0;
  std::vector<Long64_t> rfevCount1;
  std::vector<Long64_t> rfevCount2;


  double offSet =0;

  // want number of events after the gap
  Int_t ngapCount=0;
  std::vector<Int_t> ngapTrig;
  pdsLast=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    tree->GetEntry(jentry);   
    bclk->clear();
    // load digi times 
    //
    dtime[0]=clk.dt0; dtime[1]=clk.dt1; dtime[2]=clk.dt2;
    rft[0]=clk.rf0; rft[1]=clk.rf1; rft[2]=clk.rf2;

    double pdst = static_cast<double>(clk.compSec)*1000000000 + static_cast<double>(clk.compNano);
    double dpdst=0;
    if(pdsLast>0) dpdst = pdst - pdsLast;
    if( dpdst*1e-6 > 155 ) { 
      ngapTrig.push_back(ngapCount);
      ngapCount=0;
    } else {
      ++ngapCount;
    }
    pdsLast = pdst;
  }

  printf(" pds gap count is ngapTrig size %lu \n",ngapTrig.size());




  Int_t ngap=0;
  pdsLast=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //if(jentry>20000) break;
    tree->GetEntry(jentry);   
    bclk->clear();
    // load digi times 
    //
    dtime[0]=clk.dt0; dtime[1]=clk.dt1; dtime[2]=clk.dt2;
    rft[0]=clk.rf0; rft[1]=clk.rf1; rft[2]=clk.rf2;

    double pdst = static_cast<double>(clk.compSec)*1000000000 + static_cast<double>(clk.compNano);
    double dpdst=0;
    if(pdsLast>0) dpdst = pdst - pdsLast;
    if( dpdst*1e-6 > 155 ) { //&& dpdst*1e-6 < 260) {
      bclk->ngap = ++ngap;
      if(unsigned(ngap) < ngapTrig.size()) bclk->ngapTrig= ngapTrig[ngap];
    } else {
      bclk->ngap=0;
    }

    if( dpdst*1e-9 > 17) offSet += dpdst/8;  // units of board digitizer clock 

    for(int ib=0; ib<NB; ++ib) {
      isRF[ib]=0;
      ++nRF[ib];
      if(rft[ib]>494&&rft[ib]<500) {
        isRF[ib]=1;
        nRF[ib]=0;
        ++rfcount[ib];
        ++rfcountEvent[ib];
      }
    }

    // lists of rf events by board
    if(isRF[0]==1) rfev0.push_back(jentry);
    if(isRF[1]==1) rfev1.push_back(jentry);
    if(isRF[2]==1) rfev2.push_back(jentry);

    if(isRF[0]==1) {
      rfevCount0.push_back(rfcountEvent[0]);
      rfevCount1.push_back(rfcountEvent[1]);
      rfevCount2.push_back(rfcountEvent[2]);
    }


    isrf0.push_back( isRF[0] );
   
    // take out steps
    steps.reset();
    for(int ib=0; ib<NB; ++ib) {
      delta[ib] =  dtime[ib] - dlast[ib];
      if( -delta[ib] > HALF && dlast[ib]>0 ) {
        boff[ib] += STEP;
        steps.set(size_t(ib));
        //cout <<" **** " << ib << "  " <<  steps << endl; 
      } else
        ++nskips[ib];
      btime[ib]=dtime[ib]+boff[ib]+offSet;
    }


    // if(jentry>5000) break;
    // start of run
    /*
    if(jentry%5000==0) {
      bstart[1]= - btime[1] + btime[0];
      bstart[2]= - btime[2] + btime[0];
    }
    for(int ib=0; ib<NB; ++ib) {
      btime[ib] += bstart[ib];
      blast[ib] += bstart[ib];
    }
    */

    // slope 
    //btime[1] /= bslope[1];
    //btime[2] /= bslope[2];


    //for(int ib=0; ib<NB; ++ib) btime[ib]=dtime[ib]+boff[ib]+bstart[ib];

    if( steps.test(2) ) bitSum(dtime[2],bsum2);
    bitSum(dtime[2],bsumAll2);
    bitSum(dtime[0],bsumAll0);


    c0.push_back(btime[0]);
    c1.push_back(btime[1]);
    c2.push_back(btime[2]);
    rt0.push_back(dtime[0]);
    rt1.push_back(dtime[1]);
    rt2.push_back(dtime[2]);
    timeMap0.insert ( std::pair<double,Long64_t>(btime[0],jentry) );
    timeMap1.insert ( std::pair<double,Long64_t>(btime[1],jentry) );
    timeMap2.insert ( std::pair<double,Long64_t>(btime[2],jentry) );

    
    if(steps.to_ulong()>0) {
      ntStep->Fill(float(jentry),float(steps.to_ulong()),float(nskips[0]),float(nskips[1]),float(nskips[2]),
          delta[0],delta[1],delta[2],dtime[0],dtime[1],dtime[2],btime[0],btime[1],btime[2]);
      //printf(" DLAST (%10.0f %10.0f %10.0f) \n",dlast[0],dlast[1],dlast[2]);
    }

    //reset steps if bit is set 
    for(size_t istep=0; istep<3; ++ istep) if( steps.test(istep) ) nskips[istep]=0;


    // save jump time, RF time
    for(int ib=0; ib<NB; ++ib) {
      if(delta[ib]*8>jump && delta[ib]*8<jumpMax )   jumpTime[ib]=btime[ib];
      if(isRF[ib]==1) rfTime[ib]=btime[ib];
    }

    rftime0.push_back( rfTime[0] );
    rftime1.push_back( rfTime[1] );
    rftime2.push_back( rfTime[2] );

    
    for(int ib=0; ib<NB; ++ib) {
      if(delta[ib]*8<jump && delta[ib]*8<jumpMax  ) {
        ++njump[ib];
      } else  { 
        njump[ib]=0;
        ++jcount[ib];
      }
    }

    
    ntJump->Fill(jentry,btime[0],btime[0]-blast[0],btime[1]-blast[1],btime[2]-blast[2],njump[0],njump[1],njump[2],btime[0]-jumpTime[0],btime[1]-jumpTime[1],btime[2]-jumpTime[2]);
    ntGap->Fill(jentry,btime[0]*8,nRF[0],nRF[1],nRF[2],rft[0]*8,rft[1]*8,rft[2]*8,
        jumpTime[0]*8,rfTime[0]*8,jumpTime[1]*8,rfTime[1]*8,jumpTime[2]*8,rfTime[2]*8);


    hDiff0->Fill(btime[0]-jumpTime[0]);
    hDiff1->Fill(btime[1]-jumpTime[1]);
    hDiff2->Fill(btime[2]-jumpTime[2]);

      // save jump times
    
    jtime0.push_back(jumpTime[0]);
    jtime1.push_back(jumpTime[1]);
    jtime2.push_back(jumpTime[2]);


    //if( (isRF[0]&!isRF[1])||(isRF[1]&!isRF[0]) ) printf(" MISSINGRF (%.0f,%.0f,%.0f) entry %lli run %i  ev %lli \n",rft[0],rft[1],rft[2],jentry,runNumber,jentry%5000);

    if(jentry%5000==0) {
      jcountRun0.push_back(jcount[0]);
      rfcountRun0.push_back(rfcount[0]);
      jcountRun1.push_back(jcount[1]);
      rfcountRun1.push_back(rfcount[1]);
      jcountRun2.push_back(jcount[2]);
      rfcountRun2.push_back(rfcount[2]);
      runCount.push_back(double(runNumber));
        //if(jcount[0]!=jcount[1]||jcount[2]!=rfcount[0])
        printf(" run %i %10lld  jcount (%6i %6i %6i )  rfcount  (%6i %6i %6i ) total rfcount (%6i %6i %6i )  \n",runNumber,jentry,jcount[0],jcount[1],jcount[2],
            rfcount[0],rfcount[1],rfcount[2],rfcountEvent[0],rfcountEvent[1],rfcountEvent[2]);
      ntRunCount->Fill( float(runNumber), float(rfcount[0]), float(rfcount[1]), float(rfcount[2]), float(jcount[0]), float(jcount[1]), float(jcount[2]) ); 
      // time in seconds
      ntRunTime->Fill(double(runNumber),btime[0]*8.0E-9,blast[0]*8.0E-9,dtime[0]*8.0E-9,dlast[0]*8.0E-9);
      ++runNumber;
    }

     ntEvCount->Fill( float(jentry),float(rfcount[0]), float(rfcount[1]), float(rfcount[2]), float(jcount[0]), float(jcount[1]), float(jcount[2]) );


    // fill tree
    bclk->run      = clk.run;
    bclk->event    = jentry;
    bclk->compSec  = clk.compSec;
    bclk->rf0      = clk.rf0;
    bclk->rf1      = clk.rf1;
    bclk->rf2      = clk.rf2;
    bclk->nrf0     = nRF[0]; 
    bclk->nrf1     = nRF[1]; 
    bclk->nrf2     = nRF[2];
    bclk->nj0      = njump[0]; 
    bclk->nj1      = njump[1]; 
    bclk->nj2      = njump[2];
    bclk->pdst =  pdst;
    bclk->dpdst = dpdst;
    bclk->dt0      = clk.dt0;
    bclk->dt1      = clk.dt1;
    bclk->dt2      = clk.dt2;
    bclk->tPrompt  = clk.tPrompt; // one for each board
    bclk->tPromptToRF = clk.tPromptToRF;
    bclk->bt0=btime[0];
    bclk->bt1=btime[1];
    bclk->bt2=btime[2];
    bclk->dbt0=btime[0]-blast[0];
    bclk->dbt1=btime[1]-blast[1];
    bclk->dbt2=btime[2]-blast[2];
    bclk->rftime0=rfTime[0];
    bclk->rftime1=rfTime[1];
    bclk->rftime2=rfTime[2];
    bclk->jtime0=jumpTime[0];
    bclk->jtime1=jumpTime[1];
    bclk->jtime2=jumpTime[2];
    bClock->Fill();

    //if( jentry>15300&&jentry<15400) {
    //  printf(" XXXXXXXXXXXXXXX  %llu %u %.0f \n",jentry,clk.dt0,btime[0]);
    //}


    // save last 
    pdsLast = pdst;
    for(int ib=0; ib<NB; ++ib) {
      dlast[ib] = dtime[ib];
      blast[ib] = btime[ib];
    }
  } // loop over entries

  printf(" run %i %10lld  jcount (%6i %6i %6i )  rfcount  (%6i %6i %6i ) total rfcount (%6i %6i %6i )  \n",runNumber,nentries,jcount[0],jcount[1],jcount[2],
            rfcount[0],rfcount[1],rfcount[2],rfcountEvent[0],rfcountEvent[1],rfcountEvent[2]);

  // time maps
  std::vector<double> jlist0;
  std::vector<double> jlist1;
  std::vector<double> jlist2;
  
  std::vector<int> nlist0;
  std::vector<int> nlist1;
  std::vector<int> nlist2;

  
  std::map<double,Long64_t>::iterator iter;

  Long64_t place=0;
  double last0=-1;
  for (iter=timeMap0.begin(); iter!=timeMap0.end(); ++iter) {
    if(last0==-1) last0=iter->first;
    tmap0.ip=place;
    tmap0.ev=iter->second;
    tmap0.dt = rt0[place];
    tmap0.bt=iter->first;
    tmap0.dbt=iter->first-last0;
    tmap0.jt=jtime0[unsigned(iter->second)];
    tmap0.rf=rftime0[unsigned(iter->second)];
    //printf(" %llu %.0f %.0f \n",place,iter->first,c0[place]); the same
    ntMap0->Fill();
    if( (iter->first-last0)*8<jump && (iter->first-last0)*8<jumpMax  ) {
      jlist0.push_back(8*iter->first); 
      nlist0.push_back( iter->second );
    }
    last0=iter->first;
    ++place;
    //if(place<100) printf(" %.0f %.0f %.0f %.0f \n",float(place),float(iter->second),float(iter->first),float(rt0[unsigned(iter->second)]));
  }


  place=0;
  double last1=-1;
  for (iter=timeMap1.begin(); iter!=timeMap1.end(); ++iter) {
    if(last1==-1) last0=iter->first;
    tmap1.ip=place;
    tmap1.ev=iter->second;
    tmap1.dt = rt1[place];
    tmap1.bt=iter->first;
    tmap1.dbt=iter->first-last1;
    tmap1.jt=jtime1[unsigned(iter->second)];
    tmap1.rf=rftime1[unsigned(iter->second)];
    ntMap1->Fill();

    if( (iter->first-last1)*8<jump && (iter->first-last1)*8<jumpMax  ){
      jlist1.push_back(8*iter->first); 
      nlist1.push_back( iter->second );
    }
    last1=iter->first;
    ++place;
    //if(place<100) printf(" %.0f %.0f %.0f %.0f \n",float(place),float(iter->second),float(iter->first),float(rt0[unsigned(iter->second)]));
  }

  place=0;
  double last2=-1;
  
  for (iter=timeMap2.begin(); iter!=timeMap2.end(); ++iter) {
    if(last1==-1) last0=iter->first;
    tmap2.ip=place;
    tmap2.ev=iter->second;
    tmap2.dt = rt2[place];
    tmap2.bt=iter->first;
    tmap2.dbt=iter->first-last2;
    tmap2.jt=jtime2[unsigned(iter->second)];
    tmap2.rf=rftime2[unsigned(iter->second)];
    ntMap2->Fill();

    if( (iter->first-last2)*8<jump && (iter->first-last2)*8<jumpMax  ) {
      jlist2.push_back(8*iter->first); 
      nlist2.push_back( iter->second );
    }
    last2=iter->first;
    ++place;
    //if(place<100) printf(" %.0f %.0f %.0f %.0f \n",float(place),float(iter->second),float(iter->first),float(rt0[unsigned(iter->second)]));
  }


  unsigned minList = jlist0.size();
  if(jlist1.size()<minList) minList=jlist1.size();
  if(jlist2.size()<minList) minList=jlist2.size();

  std::vector<double> jlistDiff01;
  std::vector<double> jlistDiff02;
  std::vector<int> nlistDiff01;
  std::vector<int> nlistDiff02;

  for(unsigned il=0; il<minList; ++ il) {
    jlistDiff01.push_back( jlist0[il] - jlist1[il] );
    jlistDiff02.push_back( jlist0[il] - jlist2[il] );
    nlistDiff01.push_back( nlist0[il] - nlist1[il] );
    nlistDiff02.push_back( nlist0[il] - nlist2[il] );

  }

  unsigned maxList = rfev0.size();
  if(rfev1.size()>maxList) maxList=rfev1.size();
  if(rfev2.size()<maxList) maxList=rfev2.size();

  for(unsigned il=0; il<maxList; ++ il) {
    unsigned index0 = il;
    if(index0>rfev0.size()-1) index0 = rfev0.size()-1;
    unsigned index1 = il;
    if(index1>rfev1.size()-1) index1 = rfev1.size()-1;
    unsigned index2 = il;
    if(index2>rfev2.size()-1) index2 = rfev2.size()-1;
    //printf(" %u  %lld %lld %lld \n",il,rfev0[index0],rfev1[index1],rfev2[index2]);
    ntRFevent->Fill( float(rfev0[index0]) , float(rfev1[index1]) , float(rfev2[index2]),float(rfevCount0[index0]),float(rfevCount1[index1]),float(rfevCount2[index2]));
  }




  //for(unsigned ir=0; ir<runCount.size(); ++ir)  ntRunCount->Fill( runCount[ir], rfcountRun0[ir], rfcountRun1[ir], rfcountRun2[ir], jcountRun0[ir], jcountRun1[ir], jcountRun2[ir] ); 

  /*
  TMultiGraph* gmultiTime = new TMultiGraph();

  TGraph *grf0 =   new TGraph(runCount.size(),&runCount[0],&rfcountRun0[0]);
  TGraph *grf1 =   new TGraph(runCount.size(),&runCount[0],&rfcountRun1[0]);
  TGraph *grf2 =   new TGraph(runCount.size(),&runCount[0],&rfcountRun2[0]);
  grf0->SetMarkerColor(kBlack); grf0->SetMarkerStyle(20); grf0->SetMarkerSize(0.5); gmultiTime->Add(grf0);
  grf1->SetMarkerColor(kBlue); grf1->SetMarkerStyle(20); grf1->SetMarkerSize(0.5); gmultiTime->Add(grf1);
  grf2->SetMarkerColor(kRed); grf2->SetMarkerStyle(20); grf2->SetMarkerSize(0.5); gmultiTime->Add(grf2);


  TString multiTitle;
  multiTitle.Form(" rf time count ; run ; count");
  gmultiTime->SetTitle(multiTitle.Data());
  TCanvas *cmultiTime = new TCanvas("rf-count","rf-count");
  gmultiTime->Draw("ap");
  cmultiTime->Print(".gif");


  TGraph *gjump0 =   new TGraph(runCount.size(),&runCount[0],&jcountRun0[0]);
  TGraph *gjump1 =   new TGraph(runCount.size(),&runCount[0],&jcountRun1[0]);
  TGraph *gjump2 =   new TGraph(runCount.size(),&runCount[0],&jcountRun2[0]);
  */


  /*
  TGraph *gr01 = new TGraph(minList,&jlist0[0],&jlistDiff01[0]);
  TGraph *gr02 = new TGraph(minList,&jlist0[0],&jlistDiff02[0]);

  TCanvas *can1 = new TCanvas("jump01","jump01");
  gr01->Draw("ap");
  TCanvas *can2 = new TCanvas("jump02","jump02");
  gr02->Draw("ap");


  TGraph *gn01 = new TGraph(minList,&nlist0[0],&nlistDiff01[0]);
  TGraph *gn02 = new TGraph(minList,&nlist0[0],&nlistDiff02[0]);
  gn01->SetTitle("jump diff board 0,1");
  gn02->SetTitle("jump diff board 0,2");
  gn01->GetHistogram()->GetXaxis()->SetTitle(" event ");
  gn01->GetHistogram()->GetYaxis()->SetTitle(" jump event number diff");
  gn02->GetHistogram()->GetXaxis()->SetTitle(" event ");
  gn02->GetHistogram()->GetYaxis()->SetTitle(" jump event number diff");




  TCanvas *can3 = new TCanvas("nev-jump01","nev-jump01");
  gn01->Draw("ap");
  TCanvas *can4 = new TCanvas("nev-jump02","nev-jump02");
  gn02->Draw("ap");


  fout->Append(can1);
  fout->Append(can2);
  fout->Append(can3);
  fout->Append(can4);

 */

  cout << "\t\t size of c0 " << c0.size() << " size of c1 " << c1.size() << " size of c2 " << c2.size()<< endl;
  cout << "\t\t size of jist0 " << jlist0.size() << " size of jlist1 " << jlist1.size() << " size of jlist2 " << jlist2.size()<< endl;
  cout << "\t\t size of timeMap0 " << timeMap0.size() << " size of timeMap1 " << timeMap1.size() << " size of timeMap2 " << timeMap2.size()<< endl;


  for(unsigned i=0; i<bsum2.size(); ++i)  hBits2->SetBinContent(i+1,bsum2[i]);
  for(unsigned i=0; i<bsumAll2.size(); ++i)  hBitsAll2->SetBinContent(i+1,bsumAll2[i]);
  for(unsigned i=0; i<bsumAll0.size(); ++i)  hBitsAll0->SetBinContent(i+1,bsumAll0[i]);

  fout->Write();

  char buff[250];
  std::string fileName=std::string("map0.txt");
  cout << " make file " << fileName << endl;
  std::ofstream ofs;
  ofs.open (fileName.c_str(), std::ofstream::out );
  place=0;
  Long64_t lastEv=-1;
  ofs << "# time ordering of board 0 " << endl;
  for (iter=timeMap0.begin(); iter!=timeMap0.end(); ++iter) {
    //if( iter->second<lastEv) {
      sprintf(buff," place %lli sorted time %.0f digitime  %.0f event %lli  \n",place,c0[unsigned(iter->second)], rt0[unsigned(iter->second)],iter->second);
      ofs << buff;
    //}
    lastEv=iter->second;
    ++place;
  }
  ofs.close();

  fileName=std::string("map1.txt");
  cout << " make file " << fileName << endl;
  ofs.open (fileName.c_str(), std::ofstream::out );
  place=0;
  lastEv=-1;


  ofs << "# time ordering of board 1 " << endl ;
  for (iter=timeMap1.begin(); iter!=timeMap1.end(); ++iter) {
    if(  iter->second<lastEv ) {
      sprintf(buff," place %lli sorted time %.0f digitime  %.0f event %lli  \n",place,c1[unsigned(iter->second)], rt1[unsigned(iter->second)],iter->second);
      ofs << buff;
    }
    lastEv=iter->second;
    ++place;
  }
  
  cout << " boardClock has  " << bClock->GetEntriesFast() << endl;
  ofs.close();
}
