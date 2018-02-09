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
  vsec.clear();
  vnano.clear();
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

    printf(" PDS triggers %i \n",ntrig000);
    Int_t nrfTrig = ntrig555+ntrig444;
    printf(" RF  triggers %i \n",nrfTrig);
   /*
    printf(" \n PMT averages \n");
    for(Int_t j=0; j<NPMT; ++j){ 
      printf(" ipmt %i norm %i qmax %.2f +/- %.2f sum %.2f +/- %.2f gain %.2f +/- %.2f \n",j,int(norm[j]),qmax[j],eqmax[j],qsum[j],eqsum[j],gain[j],gain_e[j]);
    }
    */

    // print out time info
    /*for(unsigned it=0; it<vtrig.size() ; ++it ) {
      printf(" event (%i,%i) trig %i sec %i ns ",(int) vevent[it], (int) ventry[it],vtrig[it], (int) vcompSec[it]);
      cout << vcompNs[it] << endl;
    }
    */
}
