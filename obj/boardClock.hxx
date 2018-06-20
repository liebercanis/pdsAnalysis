#ifndef BOARDCLOCK_DEFINED
#define BOARDCLOCK_DEFINED
#include <TObject.h>
class boardClock:  public TObject {
  public :

    boardClock() { printf(" boardClock instance \n"); }
    virtual ~boardClock() {;} 

    Int_t  run;
    Int_t  event;
    Int_t  compSec;
    Int_t  ngap;
    Int_t  ngapTrig;
    Int_t rf0;
    Int_t rf1;
    Int_t rf2;
    Int_t nrf0;
    Int_t nrf1;
    Int_t nrf2;
    Int_t nj0;
    Int_t nj1;
    Int_t nj2;
    Long64_t compNano;
    UInt_t dt0;
    UInt_t dt1;
    UInt_t dt2;   // caen digitizer time 
    Double_t tPrompt; // one for each board
    Double_t tPromptToRF;
    Double_t pdst;
    Double_t dpdst;
    Double_t bt0;
    Double_t bt1;
    Double_t bt2;
    Double_t dbt0;
    Double_t dbt1;
    Double_t dbt2;
    Double_t rftime0;
    Double_t rftime1;
    Double_t rftime2;
    Double_t jtime0;
    Double_t jtime1;
    Double_t jtime2;
   

    void clear() {
      // fill tree
      run      =0; 
      event    =0; 
      compSec  =0; 
      ngap     =0;
      ngapTrig =0;
      rf0      =0; 
      rf1      =0; 
      rf2      =0; 
      nrf0     =0; 
      nrf1     =0; 
      nrf2     =0; 
      nj0      =0;  
      nj1      =0;  
      nj2      =0; 
      compNano =0; 
      dt0      =0; 
      dt1      =0; 
      dt2      =0; 
      tPrompt  =0;
      tPromptToRF   =0;
      pdst     =0;
      dpdst    =0;
      bt0           =0;
      bt1           =0;
      bt2           =0;
      dbt0          =0;
      dbt1          =0;
      dbt2          =0;
      rftime0       =0;
      rftime1       =0;
      rftime2       =0;
      jtime0        =0;
      jtime1        =0;
      jtime2        =0;
    }

	ClassDef(boardClock,6)

 };
#endif

