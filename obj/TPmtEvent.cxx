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
  tag.clear();
  run=0;
  rft21.clear();
  rft22.clear();
  rft23.clear();
  event=0;
  tpcTrig=0;
  pdsTrig=0;
  compSec=0;
  compNano=0;
  tRFave=0;
  promptLike=0;
  nhits=0;
  hit.clear();	 
  qsum.clear();
  qmax.clear();
  tPrompt=0;   
  for(int i=0; i<3; ++i) {
    dtime[i]=0;
  }
}

void TPmtEvent::print(int nHitsToPrint) 
{
  int r21=0; if(rft21.size()>0) r21=rft21[0];
  int r22=0; if(rft22.size()>0) r22=rft22[0];
  int r23=0; if(rft23.size()>0) r23=rft23[0];
  printf(" \n  TPmtEvent %s run %i sec %i nano %i rf21 %i rf22 %i rf23 %i nhits %i \n",tag.c_str(),run,compSec,int(compNano),r21,r22,r23,nhits);
  for(Int_t i=0; i<nHitsToPrint; ++i) hit[i].print();
  printf(" end of TPmtEvent \n");
}

