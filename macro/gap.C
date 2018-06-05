
enum {NB=3};
const double jump = 194.0E6;
const double jumpMax = 202.0E6;
std::vector<int> tpcRun;
std::vector<int> tpcEvent;
std::vector<unsigned long long> tpcTime;


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
  std::vector<unsigned long long> pdstime;
  std::vector<unsigned long long> gtime0;
  std::vector<unsigned long long> mtime;
  std::vector<unsigned long long> gaptime0;
  std::vector<unsigned long long> pulsetime0;
  std::vector<unsigned long long> rftime0;
  std::vector<unsigned long long> rflist0;
  std::vector<unsigned long long> rftime1;
  std::vector<unsigned long long> rflist1;
  std::vector<unsigned long long> rftime2;
  std::vector<unsigned long long> rflist2;





  TFile *fout = new TFile("gap.root","RECREATE");

  TH1D *hDiff = new TH1D("hdiff"," time diff tpc - gap ",400,-200,200);
  hDiff->GetXaxis()->SetTitle(" time difference milli-sec");

  TH1D *hDiffRF = new TH1D("hdiffRF"," time diff tpc - RF ",400,-200,200);
  hDiffRF->GetXaxis()->SetTitle(" time difference milli-sec");

  TNtuple *nt = new TNtuple("nt"," tpc-gap match ","entry:run:itpc:igap:irf:tpctime:gpdstime:rfpdstime:gtime:rftime:dtpc:dgtime:drftime");
  TNtuple *ntGap = new TNtuple("ntgap"," gaps ","ngap:gtime:diff:nj");
  TNtuple *ntPdsTime = new TNtuple("pdsTime"," pds times ","ipds:pdst:bt0");


  bool isGap[NB];
  bool isRF[NB];

  tree->GetEntry(0);

  unsigned long long btime0 = static_cast<unsigned long long>(bclk->bt0);
  unsigned long long pdstime0 = static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano);
  unsigned long long mtime0 = tpcTime[0];
  cout << " tpc time zero " <<  mtime0 << " pdstime0  " << pdstime0 << " board 0 time zero " << btime0 << endl;


  for(unsigned iev=0; iev<tpcTime.size(); ++iev ) mtime.push_back( tpcTime[iev] -mtime0);

  printf(" TPC time zero %10E last time %10E elapsed time %10E \n", double(mtime0)*1.0e-6, double(tpcTime[tpcTime.size()-1] )*1.0e-6, double(mtime[tpcTime.size()-1] )*1.0e-6);



  // loop
  for (Long64_t jentry=0; jentry<fsize ;jentry++) {
    tree->GetEntry(jentry);
    for(int ib=0; ib<NB; ++ib) { isGap[ib]=false; isRF[ib]=false;}

     unsigned long long lpdst =  static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano) -pdstime0;
     unsigned long long lbt0  =  static_cast<unsigned long long>(bclk->bt0)*8 - btime0 ;
     ntPdsTime->Fill( double(jentry),double(lpdst)*1E-6, double(lbt0)*1E-6 );



    if( bclk->bdelta0*8>jump && bclk->bdelta0*8<jumpMax  ) {
      isGap[0]=true;
      glist0.push_back(jentry);
      gtime0.push_back( static_cast<unsigned long long>(bclk->bt0)*8 - btime0 );
      nj0.push_back(bclk->nj0);
      pdstime.push_back( static_cast<unsigned long long>(bclk->compSec)*1000000000 + static_cast<unsigned long long>(bclk->compNano) -pdstime0);
    }
    if( bclk->bdelta1*8>jump && bclk->bdelta1*8<jumpMax  ) {
      isGap[1]=true;
      glist1.push_back(jentry);
      nj1.push_back(bclk->nj1);
    }

    if( bclk->bdelta2*8>jump && bclk->bdelta2*8<jumpMax  ) {
      isGap[2]=true;
      glist2.push_back(jentry);
      nj2.push_back(bclk->nj2);
    }

    // RF list

    if( bclk->rf0>494 && bclk->rf0<500  ) isGap[0]=true;
    if( bclk->rf1>494 && bclk->rf1<500  ) isGap[1]=true;
    if( bclk->rf2>494 && bclk->rf2<500  ) isGap[2]=true;

    if(isGap[0] ) {
      rflist0.push_back(jentry);
      rftime0.push_back( static_cast<unsigned long long>(bclk->bt0)*8 - btime0 );
    }
    if(isGap[1] ) {
      rflist1.push_back(jentry);
      rftime1.push_back( static_cast<unsigned long long>(bclk->bt0)*8 - btime0 );
    }
    if(isGap[2] ) {
      rflist2.push_back(jentry);
      rftime2.push_back( static_cast<unsigned long long>(bclk->bt0)*8 - btime0 );
    }


  }

  printf(" gap counts : %lu %lu %lu \n",glist0.size(),glist1.size(),glist2.size() ) ;
  printf(" RF  counts : %lu %lu %lu \n",rflist0.size(),rflist1.size(),rflist2.size() ) ;
  for(unsigned jev=0; jev<gtime0.size() ; ++jev) gaptime0.push_back( gtime0[jev]-gtime0[0]);
  for(unsigned jev=0; jev<rftime0.size() ; ++jev) pulsetime0.push_back( rftime0[jev]-rftime0[0]);

   printf(" PDS time zero %10E elapsed time %10E board time zero %10E elapsed time %10E \n", pdstime0*1.0e-6 , double(pdstime[pdstime.size()-1] )*1.0e-6, 
        btime0*1.0e-6 , double(gtime0[gtime0.size()-1] )*1.0e-6);


  unsigned long long glast=gtime0[0];
  for(unsigned jev=0; jev<gtime0.size() ; ++jev) {
    if(jev>0) glast=gtime0[jev-1];
    ntGap->Fill( jev, double(gtime0[jev])*1E-6 , double(gtime0[jev])*1E-6 - double(glast)*1E-6,double(nj0[jev]));
    if( double(gtime0[jev]-glast)*1E-6 > 1E12) printf(" gap %10u time %llu previous time %llu diff %f \n", jev, gtime0[jev], glast ,  double(gtime0[jev])*1E-6 - double(glast)*1E-6);
  }

  for(unsigned iev=0; iev<mtime.size() ; ++iev) if(iev%5000==0) printf(" %10u run %10i event %10i tpcTime %12llu mtime %12llu \n",iev,tpcRun[iev],tpcEvent[iev],tpcTime[iev],mtime[iev]);
  for(unsigned iev=0; iev<gaptime0.size() ; ++iev)  if(iev%5000==0) printf(" %10u gtime %10llu \n",iev,gaptime0[iev]);

  // output table 
  std::string tableFile=std::string("matchSummary.txt");
  std::ofstream ofs;
  ofs.open (tableFile.c_str(), std::ofstream::out );
  cout << " writing table " << tableFile << endl;
  char buff[500];
  TDatime stamp;
  sprintf(buff,"# table created %s \n# entry:run:itpc:igap:irf:tpctime:gpdstime:rfpdstime:gtime:rftime:dtpc:dgtime:drftime\n",stamp.AsString());
  ofs << buff;


  unsigned jmin=0;
  unsigned rfmin=0;
  unsigned maxGtime = 80000;
  unsigned maxRFtime = 80000;
  unsigned long long mlast=0;
  unsigned jminMinus=0;
  unsigned rfminMinus=0;
  unsigned nbad=0;
  printf(" looping over %lu TPC events \n",mtime.size());
  for(unsigned iev=0; iev<mtime.size() ; ++iev) {
  //for(unsigned iev=0; iev<100 ; ++iev) {
    /* find nearest gap */
    if(iev>0) mlast = mtime[iev-1];
    double mindiff = 1000000000000;
    double diff = mindiff;
    for(unsigned jev=0; jev< maxGtime ; ++jev) {
      diff = abs( double(mtime[iev]) - double(gaptime0[jev]));
      if(diff<mindiff) { 
        mindiff=diff; 
        jmin=jev;       
        if(rfmin>0) jminMinus=jmin-1;
      }
    }
    //if(iev%1000==0) printf(" GAP %u  %u mtime %llu gtime %llu diff %lli  tpcdiff %f\n",iev,jmin,mtime[iev],gaptime0[jmin], mtime[iev] - gaptime0[jmin],  double(mtime[iev])*1E-6 -  double(mlast)*1E-6);
    hDiff->Fill(  double(mtime[iev])*1E-6 - double(gaptime0[jmin])*1E-6 );
    /*
    if( abs( Long64_t(gaptime0[jmin]) - Long64_t(mtime[iev]) )  > 100000000 &&jmin>0) {
       printf(" \t BAD %u %u  %u mtime %llu gtime (%llu ,%llu ,%llu) diff %lli tpcdiff %f \n",
           ++nbad,iev,jmin,mtime[iev],gaptime0[jmin-1],gaptime0[jmin],gaptime0[jmin+1], Long64_t(mtime[iev]) - Long64_t(gaptime0[jmin]), double(mtime[iev])*1E-6 -  double(mlast)*1E-6);
    }
    */
    /* find nearest RF */
    mindiff = 1000000000000;
    diff = mindiff;
    for(unsigned jev=0; jev< maxRFtime ; ++jev) {
      diff = abs( double(mtime[iev]) - double(pulsetime0[jev]));
      if(diff<mindiff) { mindiff=diff; rfmin=jev; }
      if(rfmin>0) rfminMinus = rfmin -1; 
    }
    //if(iev%1000==0) printf(" RF %u  %u mtime %llu rftime %llu diff %lli  tpcdiff %f  tpctime %10E pdstime %10E \n",
    //    iev,jmin,mtime[iev],pulsetime0[rfmin], mtime[iev] - pulsetime0[rfmin], double(mtime[iev])*1E-6 -  double(mlast)*1E-6,   double(mtime[iev])*1E-6 , double(pdstime[jmin])*1E-6 );
    hDiffRF->Fill(  double(mtime[iev])*1E-6 - double(pulsetime0[rfmin])*1E-6 );
    nt->Fill(iev,tpcRun[iev],tpcEvent[iev],jmin, rfmin, 
        double(mtime[iev])*1E-6,  
        double(pdstime[jmin])*1E-6, 
        double(pdstime[rfmin])*1E-6, 
        double(gaptime0[jmin])*1E-6,
        double(pulsetime0[rfmin])*1E-6,
        double(mtime[iev])*1E-6 -  double(mlast)*1E-6, 
        double(gtime0[jmin])*1E-6 - double(gtime0[jminMinus])*1E-6,
        double(pulsetime0[rfmin])*1E-6 - double(pulsetime0[rfminMinus])*1E-6);

    sprintf(buff," %9u %6u %6u %6u %6u %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f \n ",
        iev,tpcRun[iev],tpcEvent[iev],jmin, rfmin, 
        double(mtime[iev])*1E-6,  
        double(pdstime[jmin])*1E-6, 
        double(pdstime[rfmin])*1E-6, 
        double(gaptime0[jmin])*1E-6,
        double(pulsetime0[rfmin])*1E-6,
        double(mtime[iev])*1E-6 -  double(mlast)*1E-6, 
        double(gtime0[jmin])*1E-6 - double(gtime0[jminMinus])*1E-6,
        double(pulsetime0[rfmin])*1E-6 - double(pulsetime0[rfminMinus])*1E-6);
    if(iev%1000==0) cout << buff;
    ofs<< buff;
  }

  fout->Write();
  ofs.close();

}
