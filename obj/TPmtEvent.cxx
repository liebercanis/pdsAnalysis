#include "TPmtEvent.hxx"
ClassImp(TPmtEvent)

TPmtEvent::TPmtEvent(): TNamed("TPmtEvent","TPmtEvent")
{
  clear();
}

//TPmtEvent::~TPmtEvent(){}

void TPmtEvent::clear()
{
  trigType = TRIGUNKNOWN; 
  run=0;
  event=0;
  tpcTrig=0;
  pdsTrig=0;
  gpsYear=0;
  gpsDay=0;
  gpsSec=0;
  gpsNs=0;
  nhits=0;
  hit.clear();	 
  qsum.clear();
  qmax.clear();
 }

void TPmtEvent::print(int ipmt)
{
  printf(" \n  TPmtEvent ipmt %i  run %i nHits \n",ipmt,run,nhits);
  for(Int_t i=0; i<nhits ; ++i) hit[i].print();
  printf(" end of TPmtEvent \n");
}

