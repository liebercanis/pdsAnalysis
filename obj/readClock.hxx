#ifndef READCLOCK_DEFINED
#define READCLOCK_DEFINED
#include <TObject.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TFile.h>
#include "boardClock.hxx"
static const ULong_t timeZero = 1501535896380432212;
class readClock:  public TObject {
  public :
    readClock() { 
      Init(); 
    }
    virtual ~readClock() { 
      finput->Close();
    }
    void Init();
    void getClock(Long64_t jentry, double& pdst, double& time, double& jtime, int& jnumber );
    void readEntry(Long64_t jentry = 0); 

    boardClock* bclk;
    TNtuple *evMap;
    TTree *bTree;
    Long64_t fsize;
    TFile *finput;
    Float_t fentryOld;
    Float_t fentryNew;

  ClassDef(readClock,1)

 };
#endif
