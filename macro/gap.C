
enum {NB=3};
const double jump = 194.0E6;
const double jumpMax = 202.0E6;
const unsigned long gapTimeDrift = 2056; // ns
std::vector<int> tpcRun;
std::vector<int> tpcEvent;
std::vector<unsigned long long> tpcTime;
Int_t nSamples;
TVirtualFFT *fFFT;
TVirtualFFT *fInverseFFT;





void readTPCTimes()
{
  tpcRun.clear();
  tpcEvent.clear();
  tpcTime.clear();

  TString fileName("summary_TPC.txt");
  ifstream in;
  in.open(fileName);
  Int_t events=0;
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return;
  }
  printf(" readTPCTimes from file %s \n",fileName.Data());
  string line;
  int run;
  int event;
  unsigned long long time;
  while (in.good()) {
    if(in.eof()) break;
    // look for comment or blank line
    char c=in.peek();
    if(c=='#') {  // ASCII # 
      getline(in,line); // throw away line
      //cout << " peek "  << c << " " << line <<endl;
      continue;
    }
    in >> run >> event >> time;
    getline(in,line); // throw away to end of line
    tpcRun.push_back(run);
    tpcEvent.push_back(event);
    tpcTime.push_back(time);
    unsigned iev = tpcTime.size()-1;
    //if(iev%5000==0) printf(" %10u run %10i event %10i tpcTime %12llu \n",iev,tpcRun[iev],tpcEvent[iev],tpcTime[iev]);
  } 
  in.close();
}


void gap() 
{


  TString fileInputName=TString("clocks.root");

  // open the file
  TFile *finput = new TFile(fileInputName);
  if (!finput || !finput->IsOpen()) {
    cout << " could not find file " << fileInputName << endl;
    return;
  }
  cout << " opened file " << fileInputName << endl;

  TTree *tree;
  finput->GetObject("bclk",tree);
  
  unsigned long long fsize = tree->GetEntriesFast();
  cout << " size of TClk tree " << fsize << endl;

  boardClock *bclk = new boardClock;
  tree->SetBranchAddress("c",&bclk);

  readTPCTimes();
  printf(" have read %lu TPC times \n",tpcTime.size());

  
   // data elements
  std::vector<unsigned long long> glist0; 
  std::vector<unsigned long long> glist1; 
  std::vector<unsigned long long> glist2; 
  std::vector<UInt_t> nj0; 
  std::vector<UInt_t> nj1; 
  std::vector<UInt_t> nj2;
  std::vector<UInt_t> rft0;
  std::vector<unsigned long long> pdstime;
  std::vector<unsigned long long> gtime0;
  std::vector<unsigned long long> mtime;
  std::vector<unsigned long long> gaptime0;
  std::vector<unsigned long long> gorder0;
  std::vector<unsigned long long> pulsetime0;
  std::vector<unsigned long long> rftime0;
  std::vector<unsigned long long> rflist0;
  std::vector<unsigned long long> rftime1;
  std::vector<unsigned long long> rflist1;
  std::vector<unsigned long long> rftime2;
  std::vector<unsigned long long> rflist2;

 std::map<unsigned long long,unsigned int> timeMap0; 



  TFile *fout = new TFile("gap.root","RECREATE");

  
  TTree *tClock = new TTree("tclk"," TPC clock ");
  tpcClock *tclk = new tpcClock;
  tClock->Branch("c",&tclk);


  TH1D *hDiff = new TH1D("hdiff"," time diff tpc - gap ",400,-200,200);
  hDiff->GetXaxis()->SetTitle(" time difference milli-sec");

  TH1D *hDiffRF = new TH1D("hdiffRF"," time diff tpc - RF ",400,-200,200);
  hDiffRF->GetXaxis()->SetTitle(" time difference milli-sec");

  TNtuple *nt = new TNtuple("nt"," tpc-gap match ","entry:run:itpc:igap:irf:tpctime:gpdstime:rfpdstime:gtime:rftime:dtpc:dgtime:tmod:gmod");
  TNtuple *ntGap = new TNtuple("ntgap"," gaps ","ngap:gtime:diff:nj");
  TNtuple *ntGapDiff = new TNtuple("ntgdiff"," gap diffs ","ngap:diff:gsum");
  TNtuple *ntTpcDiff = new TNtuple("nttpcdiff"," tpc diffs ","ngap:diff:gsum");
  TNtuple *ntPdsTime = new TNtuple("pdsTime"," pds times ","run:event:pdst:dpdst:bt0:dbt0");
  TNtuple *ntMtime = new TNtuple("mtime"," tpc computer ","run:time:tdiff:mod");
  TNtuple *ntGorder = new TNtuple("gorder"," pds board 0  ","ev:bt0:tdiff:mod");


  bool isGap[NB];
  bool isRF[NB];

  tree->GetEntry(0);

  unsigned long long btime0 = static_cast<unsigned long long>(bclk->bt0);
  unsigned long long pdstime0 = static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano);
  unsigned long long mtime0 = tpcTime[0];
  cout << " tpc time zero " <<  mtime0 << " pdstime0  " << pdstime0 << " board 0 time zero " << btime0 << endl;


  for(unsigned iev=0; iev<tpcTime.size(); ++iev ) mtime.push_back( tpcTime[iev] -mtime0);

  for(unsigned iev=1; iev<tpcTime.size(); ++iev ) ntMtime->Fill(double(tpcRun[iev]),double(iev), double(mtime[iev])*1.0e-6, double(mtime[iev]-mtime[iev-1])*1.0e-6,  double( mtime[iev]%200000000 )*1.0e-6 );

  printf(" TPC time zero %10E last time %10E elapsed time %10E \n", double(mtime0)*1.0e-6, double(tpcTime[tpcTime.size()-1] )*1.0e-6, double(mtime[tpcTime.size()-1] )*1.0e-6);


  /*/ do FFT
  nSamples=mtime.size();
  // initialize fft
  fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");

  vector<double> signal;
  for(unsigned iev=0; iev<nSamples; ++iev ) signal.push_back(  double(mtime[iev])*1.0e-9 );

  for(int is =0; is<nSamples; ++is) fFFT->SetPoint(is, signal[is]);
  fFFT->Transform();
  TH1D* hfft = new TH1D("fftTPC","FFT TPC times",nSamples/2,0,nSamples/2);
  hfft->SetXTitle(" frequency HZ");
  hfft->SetYTitle(" power ");

  //std::vector<std::complex<double> > VectorComplex;
  //std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<nSamples; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im); //.real or .imag accessors
    hfft->SetBinContent(i,hfft->GetBinContent(i)+std::abs(c));
    //VectorComplex.push_back(c);
    //realVec.push_back( VectorComplex[i].real());
    //imVec.push_back(VectorComplex[i].imag() );
  }
*/



  unsigned long long totalDrift=0;
  unsigned long long pdsLast=0;
  unsigned long long bt0Last=0;
  // loop
  for (Long64_t jentry=0; jentry<fsize ;jentry++) {
    tree->GetEntry(jentry);
    for(int ib=0; ib<NB; ++ib) { isGap[ib]=false; isRF[ib]=false;}

    unsigned pdsEv = unsigned(bclk->event);

    unsigned long long lpdst =  static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano) -pdstime0;
    unsigned long long lbt0  =  static_cast<unsigned long long>(bclk->bt0)*8 - btime0*8 - totalDrift;
    ntPdsTime->Fill( double(bclk->run),double(pdsEv),double(lpdst)*1E-6, double(bclk->dpdst)*1E-6, double(lbt0)*1E-6 , double(lbt0)*1E-6 -double(bt0Last)*1E-6   );

    double bdelta0 = bclk->dbt0*8;
    double bdelta1 = bclk->dbt1*8;
    double bdelta2 = bclk->dbt2*8;

    if( bdelta0 >jump && bdelta0<jumpMax  ) {
      isGap[0]=true;
      glist0.push_back(pdsEv);
      gtime0.push_back( lbt0);
      nj0.push_back(bclk->nj0);
      rft0.push_back(bclk->rf0);
      pdstime.push_back( static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano) -pdstime0);
    }
    if( bdelta1 >jump && bdelta1 <jumpMax  ) {
      isGap[1]=true;
      glist1.push_back(pdsEv);
      nj1.push_back(bclk->nj1);
    }

    if( bdelta2>jump && bdelta2<jumpMax  ) {
      isGap[2]=true;
      glist2.push_back(pdsEv);
      nj2.push_back(bclk->nj2);
    }

    // RF list

    if( bclk->rf0>494 && bclk->rf0<500  ) isGap[0]=true;
    if( bclk->rf1>494 && bclk->rf1<500  ) isGap[1]=true;
    if( bclk->rf2>494 && bclk->rf2<500  ) isGap[2]=true;

    if(isGap[0] ) {
      rflist0.push_back(pdsEv);
      rftime0.push_back(lbt0 );
    }
    if(isGap[1] ) {
      rflist1.push_back(pdsEv);
      rftime1.push_back( lbt0);
    }
    if(isGap[2] ) {
      rflist2.push_back(pdsEv);
      rftime2.push_back( lbt0);
    }

    if(isGap[0]) totalDrift += gapTimeDrift;
    pdsLast = lpdst;

    bt0Last = lbt0;


  }

  printf(" gap counts : %lu %lu %lu \n",glist0.size(),glist1.size(),glist2.size() ) ;
  printf(" RF  counts : %lu %lu %lu \n",rflist0.size(),rflist1.size(),rflist2.size() ) ;
  if(glist0.size()==0 || rflist0.size()==0 ) return; 

  

  for(unsigned jev=0; jev<gtime0.size() ; ++jev)  gaptime0.push_back( gtime0[jev]-gtime0[0]);
  for(unsigned jev=0; jev<rftime0.size() ; ++jev) pulsetime0.push_back( rftime0[jev]-rftime0[0]);

  for(unsigned jev=0; jev<gtime0.size() ; ++jev) timeMap0.insert ( std::pair<unsigned long long,unsigned int>(gaptime0[jev],unsigned(jev)) );

  // fill time ordered gaps
  std::map<unsigned long long, unsigned int>::iterator iter;
  for (iter=timeMap0.begin(); iter!=timeMap0.end(); ++iter) gorder0.push_back(iter->first);

   printf(" PDS time zero %10E elapsed time %10E board time zero %10E elapsed time %10E \n", pdstime0*1.0e-6 , double(pdstime[pdstime.size()-1] )*1.0e-6, 
        btime0*1.0e-6 , double(gtime0[gtime0.size()-1] )*1.0e-6);


   /* fill running clock */

  double tdiff=0;
  double rtime=0;
  tclk->clear();
  for(unsigned iev=1; iev<gaptime0.size() ; ++iev) {
    tdiff = ( double(gaptime0[iev]) - double(gaptime0[iev-1]))*1E-6;
    if(tdiff>194&&tdiff<202){ 
      //rtime += tdiff;
      rtime += 200.0;
      //printf(" %u %E %E tdiff %E seconds \n",iev,double(gaptime0[iev])*1E-6,double(gaptime0[iev-1])*1E-6,tdiff*10-3);
    } else 
      rtime += tdiff;
    tclk->event.push_back(iev);
    tclk->time.push_back(double(gaptime0[iev])*1E-6);
    tclk->rtime.push_back(rtime);
    tclk->tdiff.push_back(tdiff);
  }


  tClock->Fill();
  Long64_t ntpcentries = tClock->GetEntries();
  printf(" tpcClock filled number of runs %lld  \n",ntpcentries);
  fout->Write();





  double gapSum=0;
  double gapDiff=0;
  unsigned igapDiff=0;
  unsigned long long glast=gtime0[0];
  for(unsigned jev=0; jev<gtime0.size() ; ++jev) {
    if(jev>0) glast=gtime0[jev-1];
    ntGap->Fill( jev, double(gtime0[jev])*1E-6 , double(gtime0[jev])*1E-6 - double(glast)*1E-6,double(nj0[jev]));
    if( double(gtime0[jev]-glast)*1E-6 > 1E12) printf(" gap %10u time %llu previous time %llu diff %f \n", jev, gtime0[jev], glast ,  double(gtime0[jev])*1E-6 - double(glast)*1E-6);
    gapDiff = double(gtime0[jev]-glast)*1E-6;
    // sequential gap 
    if(gapDiff>180&&gapDiff<220) { 
      gapSum += gapDiff;
      ntGapDiff->Fill(++igapDiff,gapDiff,gapSum);
    }
  
  }

  for(unsigned iev=0; iev<mtime.size() ; ++iev) if(iev%5000==0) printf(" %10u run %10i event %10i tpcTime %12llu mtime %12llu \n",iev,tpcRun[iev],tpcEvent[iev],tpcTime[iev],mtime[iev]);
  for(unsigned iev=0; iev<gorder0.size() ; ++iev)  if(iev%5000==0) printf(" %10u gorder %10llu \n",iev,gorder0[iev]);


  for(unsigned iev=1; iev<gorder0.size(); ++iev ) ntGorder->Fill(double(iev), double(gorder0[iev])*1.0e-6, double(gorder0[iev]-gorder0[iev-1])*1.0e-6, double( gorder0[iev]%200000000 )*1.0e-6 );
  return;


  // output table 
  std::string tableFile=std::string("matchSummary.txt");
  std::ofstream ofs;
  ofs.open (tableFile.c_str(), std::ofstream::out );
  cout << " writing table " << tableFile << endl;
  char buff[500];
  TDatime stamp;
  sprintf(buff,"# table created %s times in ms \n# entry:run:itpc:igap:irf:tpctime:gpdstime:rfpdstime:gtime:rftime:dtpc:dgtime:drftime:tpctime-gtime\n",stamp.AsString());
  ofs << buff;


  unsigned jmin=0;
  unsigned rfmin=0;
  unsigned maxGtime = 80000;
  unsigned maxRFtime = 80000;
  unsigned long long mlast=0;
  unsigned jminMinus=0;
  unsigned rfminMinus=0;
  unsigned nbad=0;

  vector<double> tmicro;
  vector<double> gmicro;
  vector<double> tmicroErr;
  vector<double> gmicroErr;

  vector<double> gapnum;
  vector<double> gapnumErr;
  vector<double> gdelta;

  printf(" looping over %lu TPC events \n",mtime.size());


  double tpcSum=0;
  double tpcDiff=0;
  unsigned itpcDiff=0;
  for(unsigned iev=1; iev<mtime.size() ; ++iev) {
    tpcDiff = double(mtime[iev])*1E-6 -  double(mtime[iev-1])*1E-6;
    if(tpcDiff>180&&tpcDiff<220) { 
      tpcSum += tpcDiff;
      ntTpcDiff->Fill(++itpcDiff,tpcDiff,tpcSum);
    }
  }

  for(unsigned iev=0; iev<mtime.size() ; ++iev) {
  //for(unsigned iev=0; iev<100 ; ++iev) {
    if(iev>0) mlast = mtime[iev-1];
    tmicro.push_back( double(mtime[iev])*1E-6 ); 
    gapnum.push_back(double(iev));
    gapnumErr.push_back(0);
    /* find nearest gap */
    double mindiff = 1000000000000;
    double diff = mindiff;
    for(unsigned jev=0; jev< maxGtime ; ++jev) {
      diff = abs( double(mtime[iev]) - double(gorder0[jev]));
      if(diff<mindiff) { mindiff=diff; jmin=jev; }
      if(jmin>0) jminMinus = jmin -1;
    }

    double gdiff = mindiff;

    //if(iev%1000==0) printf(" GAP %u  %u mtime %llu gtime %llu diff %lli  tpcdiff %f\n",iev,jmin,mtime[iev],gorder0[jmin], mtime[iev] - gorder0[jmin],  double(mtime[iev])*1E-6 -  double(mlast)*1E-6);
    hDiff->Fill(  double(mtime[iev])*1E-6 - double(gorder0[jmin])*1E-6 );
    /*
    if( abs( Long64_t(gorder0[jmin]) - Long64_t(mtime[iev]) )  > 100000000 &&jmin>0) {
       printf(" \t BAD %u %u  %u mtime %llu gtime (%llu ,%llu ,%llu) diff %lli tpcdiff %f \n",
           ++nbad,iev,jmin,mtime[iev],gorder0[jmin-1],gorder0[jmin],gorder0[jmin+1], Long64_t(mtime[iev]) - Long64_t(gorder0[jmin]), double(mtime[iev])*1E-6 -  double(mlast)*1E-6);
    }
    */
    /* find nearest RF */
    mindiff = 1000000000000;
    diff = mindiff;
    for(unsigned kev=0; kev< maxRFtime ; ++kev) {
      diff = abs( double(mtime[iev]) - double(pulsetime0[kev]));
      if(diff<mindiff) { mindiff=diff; rfmin=kev; }
      if(rfmin>0) rfminMinus = rfmin -1; 
    }

    double rfdiff = mindiff;
    //if(iev%1000==0) printf(" RF %u  %u mtime %llu rftime %llu diff %lli  tpcdiff %f  tpctime %10E pdstime %10E \n",
    //    iev,jmin,mtime[iev],pulsetime0[rfmin], mtime[iev] - pulsetime0[rfmin], double(mtime[iev])*1E-6 -  double(mlast)*1E-6,   double(mtime[iev])*1E-6 , double(pdstime[jmin])*1E-6 );
    hDiffRF->Fill(  double(mtime[iev])*1E-6 - double(pulsetime0[rfmin])*1E-6 );
    nt->Fill(iev,tpcRun[iev],tpcEvent[iev],double(glist0[jmin]),double(rflist0[rfmin]), 
        double(mtime[iev])*1E-6,  
        double(pdstime[jmin])*1E-6, 
        double(pdstime[rfmin])*1E-6, 
        double(gorder0[jmin])*1E-6,
        double(pulsetime0[rfmin])*1E-6,
        double(mtime[iev])*1E-6 -  double(mlast)*1E-6, 
        double(gorder0[jmin])*1E-6 - double(gorder0[jminMinus])*1E-6,
        double(mtime[iev]%200000000 )*1.0e-6,
        double(gorder0[jmin]%200000000 )*1.0e-6);
        gmicro.push_back( (gorder0[jmin])*1E-6 ); 
        tmicroErr.push_back(10.0);
        gmicroErr.push_back(1.0);
        gdelta.push_back( double(mtime[iev])*1E-6 - double(gorder0[jmin])*1E-6 );


    sprintf(buff," %9u %6u %6u %9llu %9llu %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f \n ",
        iev,tpcRun[iev],tpcEvent[iev],glist0[jmin],rflist0[rfmin], 
        double(mtime[iev])*1E-6,  
        double(pdstime[jmin])*1E-6, 
        double(pdstime[rfmin])*1E-6, 
        double(gorder0[jmin])*1E-6,
        double(pulsetime0[rfmin])*1E-6,
        double(mtime[iev])*1E-6 -  double(mlast)*1E-6, 
        double(gorder0[jmin])*1E-6 - double(gorder0[jminMinus])*1E-6,
        double(pulsetime0[rfmin])*1E-6 - double(pulsetime0[rfminMinus])*1E-6,
        double(mtime[iev])*1E-6-double(gorder0[jmin])*1E-6);
    if(iev%1000==0) cout << buff;
    ofs<< buff;
  }

  // graph
  TMultiGraph* gmultiTime = new TMultiGraph();

  TGraphErrors *gtpc  =   new TGraphErrors(500,&gapnum[7000],&tmicro[7000],&gapnumErr[7000],&tmicroErr[7000]);
  TGraphErrors *gpds  =   new TGraphErrors(500,&gapnum[7000],&gmicro[7000],&gapnumErr[7000],&gmicroErr[7000]);
  TGraph *grdelta = new TGraph(gapnum.size(),&gapnum[0],&gdelta[0]);
  grdelta->SetMarkerColor(kBlue);   grdelta->SetMarkerSize(.5); grdelta->SetMarkerStyle(4); 
  grdelta->SetTitle("tpctime - gtime ");


  gtpc->SetMarkerColor(kBlack); gtpc->SetMarkerStyle(6); gtpc->SetMarkerSize(0.5); gmultiTime->Add(gtpc);
  gpds->SetMarkerColor(kBlue);  gpds->SetMarkerStyle(5); gpds->SetMarkerSize(1.0); gmultiTime->Add(gpds);


  TString multiTitle;
  multiTitle.Form(" gaps ; event; time in milli-seconds ");
  gmultiTime->SetTitle(multiTitle.Data());
  TCanvas *cmultiTime = new TCanvas("gap-time","gap-time");
  gmultiTime->Draw("ap");
  cmultiTime->Print(".gif");

  TCanvas *cdelta = new TCanvas("gap-delta","gap-delta");
  grdelta->Draw("ap");
  cdelta->Print(".gif");
  
  fout->Add(cmultiTime);
  fout->Add(cdelta);


  fout->Write();
  ofs.close();

}
