#include "TPmtHit.hxx"
ClassImp(TPmtHit)

TPmtHit::TPmtHit(): TNamed("TPmtHit","TPmtHit")
{
  clear();
}


TPmtHit::~TPmtHit(){}

void TPmtHit::clear()
{
  ipmt=0;
  timeToRF=0;
  tstart=0;
  tstop=0;
  ratio=0;
  qhit=0;
  qpeak=0;
  qUnhit=0;
  qUnpeak=0;
  fwhm=0; 
  peakTime=0;
  offset=0;
  nsamples=0;
  qsample.clear();	
 }

void TPmtHit::print(int ihit)
{
  printf(" \n \n TPmtHit %i  timeToRF %i qhit %f qmax %f \n",ihit,timeToRF,qhit,qpeak);
  printf("      nsamples %i start %i stop %i \n",nsamples,tstart,tstop);
  for(int i=0; i<nsamples ; ++i) printf(" %i %f ; ",i,qsample[i]);
  printf(" end of TPmtHit \n\n");
}

