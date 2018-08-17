#include "readClock.hxx"

void readClock::Init() 
{
  TString fileInputName=TString("clocks-fix5.root");
  // open the file
  finput = new TFile(fileInputName);
  if (!finput || !finput->IsOpen()) {
    printf(" readClock:: could not find file %s\n",fileInputName.Data());
    return;
  }
  printf(" readClock:: opened file %s\n",fileInputName.Data());
  finput->GetObject("bclk",bTree);
  fsize = bTree->GetEntriesFast();
  printf(" size of bclk tree %lld \n",fsize);
  bclk = new boardClock;
  bTree->SetBranchAddress("c",&bclk);
  // event map
  finput->GetObject("ntEvMap",evMap);
  printf(" size of ntEvMap %lld \n",evMap->GetEntries());
  evMap->SetBranchAddress("entryOld",&fentryOld);
  evMap->SetBranchAddress("entryNew",&fentryNew);
}

void readClock::getClock(Long64_t jentry, double& pdst, double& time, double& jtime, int& jnumber )
{    
  time=0;
  jtime=0;
  jnumber=0;
   evMap->GetEntry(jentry);
   if(fentryNew<0) {
     //printf(" readClock returning 0 for entry %lld\n",jentry);
     return;
   }
   if(jentry!=fentryOld) printf(" readClock ERROR entry %lld %lld\n",jentry,fentryOld);
   bTree->GetEntry(fentryNew);
   pdst = bclk->pdst;
   time  = bclk->bt0*8;
   jtime = bclk->jtime0*8;
   jnumber = bclk->nj0;
}
void readClock::readEntry(Long64_t jentry ) 
{
  if(jentry>fsize) return;
  double t0= double(timeZero);
  double pdst,time,jtime;
  int number;
  getClock(jentry,pdst,time,jtime,number);
  printf(" %lld %20.4f %20.4f %i \n", jentry,(pdst-t0)*1e-9,(jtime-t0)*1e-9, number);
}

ClassImp(readClock)
