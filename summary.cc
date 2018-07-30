////////////////////////////////////////////////////////
#include <string.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
// local 
#include "TPmtSummary.hxx"


void summary(TString dirName="PDS_beamtime_files")
{
 

  // get list of files
  std::vector<std::string> taglist;
  TString fullDirName= TString("2017/")+dirName;
  void *dirp = gSystem->OpenDirectory(fullDirName);
  if (!dirp) return;
  char *direntry;
  Long_t id, size,flags,modtime;
  //loop on all entries of this directory
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) { 
    //cout << direntry << endl;
    string fname = string(direntry);
    if ( strstr(fname.c_str(), "PDSout" )==NULL ) continue;
    if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
    string tag= fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_"));
    taglist.push_back(tag);
  }

  printf(" total of files in %s is %lu \n list of tags: \n",dirName.Data(),taglist.size());
  for( unsigned itag =0; itag < taglist.size() ; ++itag ) {
    printf(" %i %s \n",itag,taglist[itag].c_str());
  }

  
  TString summaryFileName = TString("pdsOutput/pdsSummary_")+dirName+ TString(".root");
  summaryFile = new TFile(summaryFileName,"recreate");
  summaryFile->cd();
  printf(" opening summary file %s \n",summaryFileName.Data());
  TTree *summaryTree = new TTree("summaryTree","summaryTree");
  pmtSummary  = new TPmtSummary();
  summaryTree->Branch("pmtSummary",&pmtSummary);


}
