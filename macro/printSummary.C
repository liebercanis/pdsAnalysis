void printSummary() 
{
  TFile *fin = new TFile("../pdsOutput/lowAna-50-0.root","READONLY");
  TTree *tree;
  fin->GetObject("summaryTree",tree);
  if(!tree) { printf(" summary tree not found \n"); return 0;}
  Long64_t nruns = tree->GetEntries();
  printf(" have %lld runs \n",nruns);
  TPmtSummary *psum = new TPmtSummary();
  tree->SetBranchAddress("pmtSummary",&psum);
  for(Long64_t entry=0; entry<nruns ; ++entry) {
    tree->GetEntry(entry);
    psum->print();
    psum->printFile();
  }
}
