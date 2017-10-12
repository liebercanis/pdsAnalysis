#include "TPmtGains.hxx"
ClassImp(TPmtGains)

TPmtGains::TPmtGains(): TNamed("TPmtGains","TPmtGains")
{
  clear();
}


TPmtGains::~TPmtGains(){}

void TPmtGains::clear()
{
  tag.clear(); // std string
  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    gain[ipmt]=0; egain[ipmt]=0;
  }
  
 }

void TPmtGains::print()
{
   printf(" \n\t GGGGGGGGGG PMT Gains GGGGGGGGGGGGGG  %s \n",tag.c_str());

      for(Int_t j=0; j<NPMT; ++j){ 
      printf(" gain[%i]=%f ; egain[%i]= %f ;\n",j,gain[j],j,egain[j]);
    }
}
