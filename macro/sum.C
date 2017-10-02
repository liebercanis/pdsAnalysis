#include <vector>
#include "TPmtSummary.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
   
// returns -1 if pmt does not exist 
// populate 3 boards, each from channel 0-6.  Channel 7 is the RF pulse. 
// valid pmt are 0 to 20, RF channels are 21,22,23
void fromPmtNumber(int ipmt, int& ib, int&ic)
{
  ib=-1; ic=-1;
  if(ipmt<0) return;
  if(ipmt>=NPMT) {
    ib=ipmt-NPMT;
    ic = 7;
  } else {
    ib=(ipmt-ipmt%NCPMT)/NCPMT;
    ic= ipmt%NCPMT;
  }
  return;
}

int toPmtNumber(int ib, int ic) 
{
  int ipmt=-1;
  if(ic<NCPMT) ipmt=ic+NCPMT*ib;
  else ipmt = ib+NPMT; 
  return ipmt;
}


void sum(TString tag= "PDS_beamtime_files")
{ 
  //TString inputFileName = TString("../pdsOutput/pdsSummary_")+tag+TString(".root");
  //printf(" opening file %s \n",inputFileName.Data()); 
  //TFile *infile = new TFile(inputFileName);
  TChain *sumTree= new TChain("summaryTree");
  sumTree->Add("../pdsOutput/pdsSummary_files-0-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-1000-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-1500-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-2000-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-2500-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-3000-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-3500-500_PDS_beamtime_files.root");
  sumTree->Add("../pdsOutput/pdsSummary_files-500-500_PDS_beamtime_files.root");

  //sumTree = (TTree*) infile->Get("summaryTree");
  Long64_t aSize=0;
  if(sumTree) aSize=sumTree->GetEntries();
  else  printf(" no summaryTree  \n");
  printf(" summaryTree with %i entries \n",int(aSize));

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("sum-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());

  TH2F *hSumPeak = new TH2F("SumPeak"," sum versus peak ",1000,0,4000,150,0,150);
  hSumPeak->GetXaxis()->SetTitle(" charge sum ");
  hSumPeak->GetYaxis()->SetTitle(" peak charge ");

  TNtuple *ntTrig = new TNtuple("ntTrig"," triggers ","time:n555:n5xx:n444:n4xx:n111:n1xx:n000:n0xx");


  TPmtSummary *pmtSum = new TPmtSummary();
  sumTree->SetBranchAddress("pmtSummary",&pmtSum);
 
  outfile->cd();
  int icolor[NPMT]={1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,1,2,3};
  TString name, title;

  TGraphErrors *grsum[NPMT];
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    grsum[ipmt] = new TGraphErrors(aSize);
    name.Form("qmax_pmt%i",ipmt);
    title.Form("qmax_pmt%i",ipmt);
    grsum[ipmt]->SetNameTitle(name,title);
    grsum[ipmt]->SetMarkerStyle(20+ipmt%3);
    grsum[ipmt]->SetMarkerColor(icolor[ipmt]);
    grsum[ipmt]->SetLineColor(icolor[ipmt]);
    grsum[ipmt]->SetMarkerSize(0.5);
  }
  TGraphErrors *grpeak[NPMT];
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    grpeak[ipmt] = new TGraphErrors(aSize);
    name.Form("qsum_pmt%i",ipmt);
    title.Form("qsum_pmt%i",ipmt);
    grpeak[ipmt]->SetNameTitle(name,title);
    grpeak[ipmt]->SetMarkerStyle(20+ipmt%3);
    grpeak[ipmt]->SetMarkerColor(icolor[ipmt]);
    grpeak[ipmt]->SetLineColor(icolor[ipmt]);
    grpeak[ipmt]->SetMarkerSize(0.5);
  }
  

  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t hour = pmtSum->getHour();
    Int_t seg = pmtSum->getSegment();
    double time = double(hour)+double(24*day+24*30*(month-7));
    ntTrig->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);
    
    if(entry%100==0) printf("...entry %i tag %s  month %i day %i hour %i seg %i time %0.f \n",entry,tag.c_str(),month,day,hour,seg,time);
    double etime=0;
    for(int ipmt=0; ipmt<NPMT; ++ipmt) {
      hSumPeak->Fill(pmtSum->qsum[ipmt],pmtSum->qmax[ipmt]);
      grsum[ipmt]->SetPoint(entry,time,pmtSum->qsum[ipmt]); 
      grsum[ipmt]->SetPointError(entry,etime,pmtSum->eqsum[ipmt]);  
      grpeak[ipmt]->SetPoint(entry,time,pmtSum->qmax[ipmt]); 
      grpeak[ipmt]->SetPointError(entry,etime,pmtSum->eqmax[ipmt]);  

    }
  }
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {outfile->Append(grsum[ipmt]);outfile->Append(grpeak[ipmt]);} 
  TString can1Name; can1Name.Form("%s-%s",tag.Data(),"sum");
  cout << " making " << can1Name << endl;
  TCanvas *c1 = new TCanvas(can1Name,can1Name);
  grsum[0]->GetHistogram()->GetXaxis()->SetTitle(" time in hours ");
  grsum[0]->SetTitle("sum charge");
  grsum[0]->Draw("AP");
  for(int ipmt=1; ipmt<NPMT; ++ipmt) grsum[ipmt]->Draw("PSAME");

  TString can2Name; can2Name.Form("%s-%s",tag.Data(),"peak");
  cout << " making " << can2Name << endl;
  TCanvas *c2 = new TCanvas(can2Name,can2Name);
  grpeak[0]->GetHistogram()->GetXaxis()->SetTitle(" time in hours ");
  grpeak[0]->SetTitle("peak charge");
  grpeak[0]->Draw("AP");
  for(int ipmt=1; ipmt<NPMT; ++ipmt) grpeak[ipmt]->Draw("PSAME");


  outfile->Write();
}
