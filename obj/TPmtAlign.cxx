#include "TPmtAlign.hxx"
ClassImp(TPmtAlign)

TPmtAlign::TPmtAlign(): TNamed("TPmtAlign","TPmtAlign")
{
  clear();
}


TPmtAlign::~TPmtAlign(){}

void TPmtAlign::clear()
{
  tag.clear();
  align0.clear();
  align1.clear();
  align2.clear();
}

void TPmtAlign::print()
{
  printf(" \n\t AAAAAAAAAAA PMT Align AAAAAAAAAA run  %s \n",tag.c_str());
}
