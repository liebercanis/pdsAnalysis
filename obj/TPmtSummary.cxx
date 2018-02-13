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

void TPmtSummary::printFile()
{
  std::string fileName=std::string("pmtSummary-")+tag+std::string(".txt");
  std::ofstream ofs;
  ofs.open (fileName, std::ofstream::out );
  print(ofs);
  ofs.close();
}

void TPmtSummary::print(std::ostream &out)
{
   char buff[250];
   sprintf(buff," \n \t\t PDS run summary %s: total events %lu ",tag.c_str(),vtrig.size()); out << buff;
   sprintf(buff," n555 %i  n5xx= %i n444= %i n4xx= %i n111= %i n1xx= %i n000= %i n0xx=  %i   \n",
        ntrig555,ntrig5xx,ntrig444,ntrig4xx,ntrig111,ntrig1xx,ntrig000,ntrig1xx); out << buff;

    Int_t nrfTrig = ntrig555+ntrig444+ntrig111;
    sprintf(buff," PDS triggers %i  RF  triggers %i \n",ntrig000,nrfTrig); out << buff;
    //printf(" gammapeak %f  tZero %f \n",gammapeak,tZero); 
    sprintf(buff," gammapeak %f  \n",gammapeak); out << buff; 

   
    //printf(" event compSec compNano RF 1 2 3 digi 1 2 3 tprompt tprompt(torf) tof ke trig nhits beamtrig delta_t \n ");
    sprintf(buff," row event entry trig compSec compNano RF 1 2 3 dtime 1 2 3 tprompt timeToRf tof ke trig nhits beamtrig delta_t \n "); out <<  buff;

    for(unsigned it=0; it<vevent.size() ; ++it ) {
 
      sprintf(buff," %4i  %4i  %4i  %4i  %lld  %5i %5i %5i %9u %9u %9u %10.3f %10.3f %10.3f %10.3f %2d %9d %2d %10.3f  \n ", 
          (int) it,  vevent[it], (int) ventry[it],  vcompSec[it], vcompNano[it],
          vrf1[it], vrf2[it], vrf3[it],  vdtime1[it],vdtime2[it], vdtime3[it],
          tprompt[it], timeToRf[it], tof[it] , ke[it], 
          vtrig[it], nhits[it], isBeamTrig(it) ,deltaT[it] );  
      out << buff;
    }
}

void TPmtSummary::printEvent(unsigned it) 
{
     
  printf(" %4i  %4i  %4i  %4i  %lld  %5i %5i %5i %9u %9u %9u %10.3f %10.3f %10.3f %10.3f %2d %9d %2d %10.3f  \n ",tag.c_str(), 
      (int) it,  vevent[it], (int) ventry[it],  vcompSec[it], vcompNano[it],
      vrf1[it], vrf2[it], vrf3[it],  vdtime1[it],vdtime2[it], vdtime3[it],
      tprompt[it], timeToRf[it], tof[it] , ke[it], 
      vtrig[it], nhits[it], isBeamTrig(it),deltaT[it] );  
}


