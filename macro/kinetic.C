
void kinetic()
{
  double L=23.2;//m
  double clight=0.299792458;//m/ns
  double nmass=939.565;//MeV

  cout << " L " << L << " clight "  << clight << " nmass " << nmass;

  TString fileName ; fileName.Form("../pdsOutput/lowAna-pmtChain-fix5.root");

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(_file0->IsZombie()) {
    cout << " not finding file " << fileName << endl;
    return;
  }

  TTree* tree;
  _file0->GetObject("summaryTree",tree);
  cout << " size of " << tree->GetEntries() << endl;

  TH1F* hgamma = new TH1F("gamma"," timeToRf ",250,-700,-200);
  TH1F* hkinetic = new TH1F("kinetic"," kinetic ",1000,0,1000);
  TH1F* hkinetica = new TH1F("kinetica"," kinetic ",1000,0,1000);
  TH1F* hkineticb = new TH1F("kineticb"," kinetic ",1000,0,1000);

  TCanvas *cgamma = new TCanvas("gamma","gamma");
  tree->Draw("timeToRf>>gamma","timeToRf>-700&&timeToRf<-200");

  double gpeak = hgamma->GetBinLowEdge(hgamma->GetMaximumBin());
  cout << " peak bin " <<  gpeak << endl;

  double gpeak1 = -610;
  double gpeak2 = -628.089;

  Long64_t nruns = tree->GetEntries();
  printf(" have %lld runs \n",nruns);
  TPmtSummary *psum = new TPmtSummary();
  tree->SetBranchAddress("pmtSummary",&psum);
  for(Long64_t entry=0; entry<nruns ; ++entry) {
    tree->GetEntry(entry);
    for(unsigned j=0; j< psum->timeToRf.size(); ++j) {
      double tof=psum->timeToRf[j]-gpeak+L/clight;
      double tof1=psum->timeToRf[j]-gpeak1+L/clight;
      double tof2=psum->timeToRf[j]-gpeak2+L/clight;
      if(tof<0) continue;
      double beta = L/tof/clight;
      double betaa = L/tof1/clight;
      double betab = L/tof2/clight;
      double gamma2 = 1./(1.-beta*beta);
      double gamma2a = 1./(1.-betaa*betaa);
      double gamma2b = 1./(1.-betab*betab);
      if(gamma2>0) hkinetic->Fill( nmass*(sqrt(gamma2)-1.0));
      if(gamma2a>0) hkinetica->Fill( nmass*(sqrt(gamma2a)-1.0));
      if(gamma2b>0) hkineticb->Fill( nmass*(sqrt(gamma2b)-1.0));
    }
  }


  TCanvas *ckinetic = new TCanvas("kinetic","kinetic");
  ckinetic->SetLogy();
  hkinetic->SetLineWidth(2);
  hkinetic->Draw();
  //hkinetica->SetLineColor(kRed);
  hkineticb->SetLineColor(kRed);
  hkineticb->SetLineWidth(2);

  //hkinetica->Draw("same");
  hkineticb->Draw("same");
  hkinetic->Draw("same");

}

///pmtSummary->timeToRf
