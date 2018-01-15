#include "TPmtAlign.hxx"
void readAlign(TString fileName="align-low-intensity.root") 
{
  cout << " reading file " << fileName << endl;
  TFile*  fin = new TFile(fileName, "READ");
  if(fin->IsZombie()) {
    printf(" cannot read file %s\n",fileName.Data());
    return;
  }
  TTree* atree = (TTree *)fin->Get("alignTree");
  if(!atree) {
    printf(" cannot find alignTree in file %s\n",fileName.Data());
    fin->Close();
    return;
  }
  
  // tree to look at alignments
  TFile*  fout = new TFile("readAlign.root", "recreate"); 
  fout->cd();
  TTree *tAlign = new TTree("tAlign"," PDS alignments ");
  Double_t a0,a1,a2;
  Double_t trig=0;
  Int_t run,event;
  tAlign->Branch("run",&run,"run/I");
  tAlign->Branch("event",&event,"event/I");
  tAlign->Branch("trig",&trig,"trig/D");
  tAlign->Branch("a0",&a0,"a0/D");
  tAlign->Branch("a1",&a1,"a1/D");
  tAlign->Branch("a2",&a2,"a2/D");
  
  TPmtAlign *pmtAlign = new TPmtAlign();
  atree->SetBranchAddress("pmtAlign",&pmtAlign);
  ULong64_t nEntry = atree->GetEntries();
  printf(" have %d alignments runs \n",int(nEntry));
  for(ULong64_t entry=0; entry< nEntry; ++entry){
    atree->GetEntry(entry);
    cout << pmtAlign->tag << " size " << pmtAlign->align0.size() << endl;
    for(unsigned i=0; i < pmtAlign->align0.size(); ++i) {
      run=int(entry);
      event=int(i);
      a0=pmtAlign->align0[i];
      a1=pmtAlign->align1[i];
      a2=pmtAlign->align2[i];
      tAlign->Fill();
      ++trig;
    }
  }
  cout << " tAlign entries = " << tAlign->GetEntries() << endl;
}


