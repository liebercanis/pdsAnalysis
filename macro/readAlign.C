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
  TH1D* hShift1 = new TH1D("Shift1"," board 1 alignment by run ",58,1,58);
  hShift1->SetXTitle(" run number ");
  hShift1->SetYTitle(" digitization time additive shift ");
  TH1D* hShift2 = new TH1D("Shift2"," board 2 alignment by run ",58,1,58);
  hShift2->SetXTitle(" run number ");
  hShift2->SetYTitle(" digitization time additive shift ");
  
  printf("run & start 0 & start 1 & start 2 & start0-start1 & start0-start2 \n"); 

  for(ULong64_t entry=0; entry< nEntry; ++entry){
    atree->GetEntry(entry);
    cout << " run "  << entry << " tag " << pmtAlign->tag << endl;
    for(unsigned i=0; i < pmtAlign->align0.size(); ++i) {
      run=int(entry);
      event=int(i);
      trig=pmtAlign->trig[i];
      a0=pmtAlign->align0[i];
      a1=pmtAlign->align1[i];
      a2=pmtAlign->align2[i];
      Long64_t rf0 = Long64_t(pmtAlign->rf0[i]*250.0);
      Long64_t rf1 = Long64_t(pmtAlign->rf1[i]*250.0);
      Long64_t rf2 = Long64_t(pmtAlign->rf2[i]*250.0);
      if(run==13) if(abs(rf0-rf1)>100||abs(rf0-rf2)>100||abs(rf1-rf2)>100) printf(" %u & %lli & %lli & %lli \\\\ \n",i,rf0,rf1,rf2);
      tAlign->Fill();
    }
    // for now all are the same, so take first entry
    // convert shifts to units of digitization time
    Long64_t start0 = Long64_t(pmtAlign->start0);
    Long64_t start1 = Long64_t(pmtAlign->start1);
    Long64_t start2 = Long64_t(pmtAlign->start2);

      
    hShift1->SetBinContent(entry+1, Double_t(start0-start1));
    hShift2->SetBinContent(entry+1, Double_t(start0-start2));
  }
  cout << " tAlign entries = " << tAlign->GetEntries() << endl;
}


