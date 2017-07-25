#include <vector>
#include "TPmtEvent.hxx"

void ana(TString tag= "07-12-1900_0")
{
  TString inputFileName = TString("../pmtAna_")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *pmtTree=NULL;

  // tree has to be in file
  Long64_t aSize=0;
  pmtTree = (TTree*) infile->Get("pmtTree");
  if(pmtTree) aSize=pmtTree->GetEntriesFast();
  printf(" pmtTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("ana-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  TH1D *hNHits = new TH1D("NHits"," pmt hits per event ",200,0,200);
 
  
  TPmtEvent *ev = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&ev);

  for(unsigned entry =0; entry < 200; ++entry ) {
    pmtTree->GetEntry(entry);
    if(entry%100==0) printf("\t entry %i hits %i \n",entry,ev->nhits);
    hNHits->Fill(ev->nhits);
  }

  // end of ana 
  outfile->Write();
}
