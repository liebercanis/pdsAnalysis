void plotSummary() 
{
  TFile *fin = new TFile("../pdsOutput/lowAna-0-0.root","READONLY");
  TTree *tree;
  fin->GetObject("summaryTree",tree);
  if(!tree) { printf(" summary tree not found \n"); return 0;}

  // open ouput file
  TString outFileName = TString("pdsSummary-lowAna.root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  TNtuple *ntTime = new TNtuple("ntTime"," time ntuple ","trig:ev:run:comps:compmu:dt0:dt1:dt2");


  Long64_t nruns = tree->GetEntries();
  printf(" have %lld runs \n",nruns);
  TPmtSummary *pmtSum = new TPmtSummary();
  tree->SetBranchAddress("pmtSummary",&pmtSum);

  Double_t STEP= pow (2.0,31); // step in microseconds
  Double_t offset[3];
  Double_t check[3];
  Double_t compTimeZero=0;
 
  for(Long64_t entry=0; entry<nruns ; ++entry) {

    // start of run
    tree->GetEntry(entry);
    Double_t compTime=0;
    Double_t dtime0Zero=Double_t(pmtSum->vdtime1[0]);
    Double_t dtime1Zero=Double_t(pmtSum->vdtime2[0]);
    Double_t dtime2Zero=Double_t(pmtSum->vdtime3[0]);
    for(int icheck=0; icheck<3; ++icheck) {
      offset[icheck]=0.0;
      check[icheck]=0.0;
    }
    
    // loop over events in tis run
    for(unsigned itr = 0; itr< pmtSum->vtrig.size(); ++ itr) {
      // look for steps
      if(itr>0) {
        check[0]= Double_t(pmtSum->vdtime1[itr-1])- Double_t(pmtSum->vdtime1[itr]);
        check[1]= Double_t(pmtSum->vdtime2[itr-1])- Double_t(pmtSum->vdtime2[itr]);
        check[2]= Double_t(pmtSum->vdtime3[itr-1])- Double_t(pmtSum->vdtime3[itr]);
      }
      for(int icheck=0; icheck<3; ++icheck) if(check[icheck]> 1.0E9) offset[icheck] +=STEP;  
        
      compTime= Double_t(pmtSum->vcompSec[itr])-Double_t(pmtSum->vcompSec[0]) + Double_t(pmtSum->vcompNano[itr])*1E9;
      Double_t compsec = Double_t(pmtSum->vcompSec[itr])-Double_t(pmtSum->vcompSec[0]);
      Double_t compmicro = compsec*1.0E6+Double_t(pmtSum->vcompNano[itr])*1.0E-3;
      //if(itr>0) compTime = Double_t(pmtSum->vcompNano[itr])/1000.0 - Double_t(pmtSum->vcompNano[itr-1])/1000.0 ;//+  Double_t(pmtSum->vcompSec[itr])*1.0E6; 
      Double_t dt0 = Double_t(pmtSum->vdtime1[itr])-dtime0Zero+offset[0];
      Double_t dt1 = Double_t(pmtSum->vdtime2[itr])-dtime1Zero+offset[1];
      Double_t dt2 = Double_t(pmtSum->vdtime3[itr])-dtime2Zero+offset[2];
      ntTime->Fill(Double_t(pmtSum->beamtrig[itr]),Double_t(itr*entry),Double_t(entry),compsec,compmicro,dt0/125.0,dt1/125.0,dt2/125.0);
    }

  }

  outFile->Write();
}
