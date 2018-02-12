#include "TPmtSummary.hxx"
ClassImp(TPmtSummary)

TPmtSummary::TPmtSummary(): TNamed("TPmtSummary","TPmtSummary")
{
  clear();
}


TPmtSummary::~TPmtSummary(){}

void TPmtSummary::clear()
{
  tZero=0;
  tag.clear(); // std string
  vtrig.clear();
  vevent.clear();
  ventry.clear();
  vcompSec.clear();
  vcompNano.clear();
  vrf1.clear(); 
  vrf2.clear(); 
  vrf3.clear(); 
  vprompt1.clear();
  vprompt2.clear();
  vprompt3.clear();
  vdtime1.clear(); 
  vdtime2.clear(); 
  vdtime3.clear(); 
  tprompt.clear();//ysun
  tof.clear();//ysun
  ke.clear();//ysun
  nhits.clear();  // number of hits in this event
  beamtrig.clear();  
  deltaT.clear(); 
  timeToRf.clear(); 
  

  run=0; min=0; seg=0; gammapeak=0;

  ntrig555=0; ntrig5xx=0; ntrig444=0; ntrig4xx=0; ntrig111=0;ntrig1xx=0; ntrig000=0; ntrig0xx=0;
  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    qsum[ipmt]=0; eqsum[ipmt]=0;
    qmax[ipmt]=0; eqmax[ipmt]=0;
    norm[ipmt]=0; noise[ipmt]=0;
    gain[ipmt]=0; gain_e[ipmt]=0;
    qrf[ipmt]=0;
  }
  
 }

void TPmtSummary::print()
{
   printf(" \n SSSSSSSSSSS summary %s: total events %lu SSSSSSSSSSSSS ",tag.c_str(),vtrig.size());
   printf(" n555 %i  n5xx= %i n444= %i n4xx= %i n111= %i n1xx= %i n000= %i n0xx=  %i   \n",
        ntrig555,ntrig5xx,ntrig444,ntrig4xx,ntrig111,ntrig1xx,ntrig000,ntrig1xx);

    Int_t nrfTrig = ntrig555+ntrig444+ntrig111;
    printf(" PDS triggers %i  RF  triggers %i \n",ntrig000,nrfTrig);
    //printf(" gammapeak %f  tZero %f \n",gammapeak,tZero); 
    printf(" gammapeak %f  \n",gammapeak); 

   
    //printf(" event compSec compNano RF 1 2 3 digi 1 2 3 tprompt tprompt(torf) tof ke trig nhits beamtrig delta_t \n ");
    printf(" row event entry trig compSec compNano RF 1 2 3 dtime 1 2 3 tprompt timeToRf tof ke trig nhits beamtrig delta_t \n ");

    for(unsigned it=0; it<vevent.size() ; ++it ) {
      printf(" %4i  %4i  %4i  %4i  %lld  %5i %5i %5i %9u %9u %9u %10.3f %10.3f %10.3f %10.3f %2d %9d %2d %10.3f  \n ", 
          (int) it,  vevent[it], (int) ventry[it],  vcompSec[it], vcompNano[it],
          vrf1[it], vrf2[it], vrf3[it],  vdtime1[it],vdtime2[it], vdtime3[it],
          tprompt[it], timeToRf[it], tof[it] , ke[it], 
          vtrig[it], nhits[it], beamtrig[it],deltaT[it] );  
    }
}

void TPmtSummary::printEvent(unsigned it) 
{
     
  printf(" %4i  %4i  %4i  %4i  %lld  %5i %5i %5i %9u %9u %9u %10.3f %10.3f %10.3f %10.3f %2d %9d %2d %10.3f  \n ",tag.c_str(), 
      (int) it,  vevent[it], (int) ventry[it],  vcompSec[it], vcompNano[it],
      vrf1[it], vrf2[it], vrf3[it],  vdtime1[it],vdtime2[it], vdtime3[it],
      tprompt[it], timeToRf[it], tof[it] , ke[it], 
      vtrig[it], nhits[it], beamtrig[it],deltaT[it] );  
}


