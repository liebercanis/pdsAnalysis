#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraphErrors.h"
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
  Double_t dnorm = double(aSize);
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
    name.Form("qsum_pmt%i",ipmt);
    title.Form("qsum_pmt%i",ipmt);
    grsum[ipmt]->SetNameTitle(name,title);
    grsum[ipmt]->SetMarkerStyle(20+ipmt%3);
    grsum[ipmt]->SetMarkerColor(icolor[ipmt]);
    grsum[ipmt]->SetLineColor(icolor[ipmt]);
    grsum[ipmt]->SetMarkerSize(0.5);
  }
  TGraphErrors *grpeak[NPMT];
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    grpeak[ipmt] = new TGraphErrors(aSize);
    name.Form("qmax_pmt%i",ipmt);
    title.Form("qmax_pmt%i",ipmt);
    grpeak[ipmt]->SetNameTitle(name,title);
    grpeak[ipmt]->SetMarkerStyle(20+ipmt%3);
    grpeak[ipmt]->SetMarkerColor(icolor[ipmt]);
    grpeak[ipmt]->SetLineColor(icolor[ipmt]);
    grpeak[ipmt]->SetMarkerSize(0.5);
  }
  
  Double_t qmaxAve[NPMT];
  Double_t qsumAve[NPMT];
  Double_t qcooper[NPMT];

  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    qmaxAve[ipmt]=0; qsumAve[ipmt]=0;
  }

  unsigned norm=0;
  Int_t amonth=7, aday=31, ahour=1555;
  double atime = double(ahour)+double(24*aday+24*30*(amonth-7));
  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t hour = pmtSum->getHour();
    Int_t seg = pmtSum->getSegment();
    double time = double(hour)+double(24*day+24*30*(month-7));

    ntTrig->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);
    
    if(entry%100==0) printf("...entry %u tag %s  month %i day %i hour %i seg %i time %0.f \n",entry,tag.c_str(),month,day,hour,seg,time);
    double etime=0;
    for(int ipmt=0; ipmt<NPMT; ++ipmt) {
      if( TMath::IsNaN(pmtSum->qsum[ipmt])  ) continue;
      if( TMath::IsNaN(pmtSum->qmax[ipmt])  ) continue;
      ++norm;
      qmaxAve[ipmt] += pmtSum->qmax[ipmt];
      qsumAve[ipmt] += pmtSum->qsum[ipmt];
      //if(time==atime) printf("  asum[%i]=%.0f ; \n",ipmt,qsumAve[ipmt]);
      hSumPeak->Fill(pmtSum->qsum[ipmt],pmtSum->qmax[ipmt]);
      grsum[ipmt]->SetPoint(entry,time,pmtSum->qsum[ipmt]); 
      grsum[ipmt]->SetPointError(entry,etime,pmtSum->eqsum[ipmt]);  
      grpeak[ipmt]->SetPoint(entry,time,pmtSum->qmax[ipmt]); 
      grpeak[ipmt]->SetPointError(entry,etime,pmtSum->eqmax[ipmt]);  
    }
  }
  //for(int ipmt=0; ipmt<NPMT; ++ipmt) {outfile->Append(grsum[ipmt]);outfile->Append(grpeak[ipmt]);} 
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

  printf(" averages over %u runs \n",norm);
  TGraph* gSumAve = new TGraph(NPMT);
  TGraph* gMaxAve = new TGraph(NPMT);
  for(int ipmt=1; ipmt<NPMT; ++ipmt) {
    qmaxAve[ipmt] /= qmaxAve[0];
    qsumAve[ipmt] /= qsumAve[0];
    gSumAve->SetPoint(ipmt,double(ipmt),qsumAve[ipmt]);
    gMaxAve->SetPoint(ipmt,double(ipmt),qmaxAve[ipmt]);
  }
  qmaxAve[0]=1.0;
  qsumAve[0]=1.0;
  gSumAve->SetPoint(0,double(0),qsumAve[0]);
  gMaxAve->SetPoint(0,double(0),qmaxAve[0]);
  

 for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" gsumAve[%i]= %.3f ;\n",ipmt,qsumAve[ipmt]);
  
 qcooper[0]= 1.000000 ;
 qcooper[1]=  1.028461 ;
 qcooper[2]=  0.883047 ;
 qcooper[3]=  0.899605 ;
 qcooper[4]=  1.144030 ;
 qcooper[5]=  1.021451 ;
 qcooper[6]=  0.813883 ;
 qcooper[7]=  0.972044 ;
 qcooper[8]=  0.830097 ;
 qcooper[9] = 0.924417 ;
 qcooper[10] = 0.818511 ;
 qcooper[11] = 0.867445 ;
 qcooper[12] = 0.560126 ;
 qcooper[13] = 0.578796 ;
 qcooper[14] = 0.708039 ;
 qcooper[15] = 0.778511 ;
 qcooper[16] = 0.794114 ;
 qcooper[17] = 0.680632 ;
 qcooper[18] = 1.106116 ;
 qcooper[19] = 1.103592 ;
 qcooper[20] = 0.776326 ;


  Double_t mean[NPMT];

  mean[0]=3.840878 ; 
  mean[1]=2.986571 ; 
  mean[2]=3.334651 ; 
  mean[3]=3.234994 ; 
  mean[4]=3.124421 ; 
  mean[5]=3.406811 ; 
  mean[6]=2.704257 ; 
  mean[7]=3.564039 ; 
  mean[8]=3.372040 ; 
  mean[9]=3.937634 ; 
  mean[10]=3.797565 ; 
  mean[11]=4.234440 ; 
  mean[12]=3.897901 ; 
  mean[13]=3.996275 ; 
  mean[14]=3.471692 ; 
  mean[15]=4.167746 ; 
  mean[16]=3.513075 ; 
  mean[17]=3.427285 ; 
  mean[18]=3.561007 ; 
  mean[19]=3.272667 ; 
  mean[20]=3.924667 ; 

  for(int ipmt=1; ipmt<NPMT; ++ipmt) mean[ipmt] /= mean[0];
  mean[0]=1;

  Double_t asum[NPMT];

  asum[0]=2365500 ; 
  asum[1]=2264227 ; 
  asum[2]=2575016 ; 
  asum[3]=2614957 ; 
  asum[4]=1588514 ; 
  asum[5]=2114388 ; 
  asum[6]=2283447 ; 
  asum[7]=2416328 ; 
  asum[8]=1359574 ; 
  asum[9]=1435025 ; 
  asum[10]=1477013 ; 
  asum[11]=1754151 ; 
  asum[12]=1422394 ; 
  asum[13]=1936981 ; 
  asum[14]=1812810 ; 
  asum[15]=2248298 ; 
  asum[16]=2182574 ; 
  asum[17]=2035239 ; 
  asum[18]=2892739 ; 
  asum[19]=2859542 ; 
  asum[20]=2015205 ; 

 for(int ipmt=1; ipmt<NPMT; ++ipmt) asum[ipmt] /= asum[0];
  asum[0]=1;

  double fitOff[NPMT]; 
   // fit to QhitOff with ZERO gain
  // fits QhitOff
 fitOff[0]=  5.930 ;
 fitOff[1]=  7.187 ;
 fitOff[2]=  4.572 ;
 fitOff[3]=  8.728 ;
 fitOff[4]=  3.866 ;
 fitOff[5]=  8.443 ;
 fitOff[6]=  3.891 ;
 fitOff[7]=  9.556 ;
 fitOff[8]=  6.472 ;
 fitOff[9]=  8.277 ;
 fitOff[10]=  4.880 ;
 fitOff[11]=  4.970 ;
 fitOff[12]=  3.384 ;
 fitOff[13]=  5.789 ;
 fitOff[14]=  3.758 ;
 fitOff[15]=  6.147 ;
 fitOff[16]=  5.535 ;
 fitOff[17]=  8.462 ;
 fitOff[18]=  5.826 ;
 fitOff[19]=  9.195 ;
 fitOff[20]=  5.694 ;

 double fitNoBeam[NPMT]; 
 // fit to QhitNoBeam with trigger000 and ZERO gain
// fits QhitNoBeam
 fitNoBeam[0]=  5.969 ;
 fitNoBeam[1]=  7.237 ;
 fitNoBeam[2]=  4.597 ;
 fitNoBeam[3]=  8.709 ;
 fitNoBeam[4]=  3.839 ;
 fitNoBeam[5]=  8.466 ;
 fitNoBeam[6]=  3.809 ;
 fitNoBeam[7]=  9.602 ;
 fitNoBeam[8]=  6.518 ;
 fitNoBeam[9]=  8.196 ;
 fitNoBeam[10]=  4.912 ;
 fitNoBeam[11]=  4.949 ;
 fitNoBeam[12]=  3.369 ;
 fitNoBeam[13]=  5.771 ;
 fitNoBeam[14]=  3.753 ;
 fitNoBeam[15]=  6.122 ;
 fitNoBeam[16]=  5.651 ;
 fitNoBeam[17]=  8.491 ;
 fitNoBeam[18]=  5.802 ;
 fitNoBeam[19]=  9.217 ;
 fitNoBeam[20]=  5.652 ;



  for(int ipmt=1; ipmt<NPMT; ++ipmt) fitOff[ipmt] /= fitOff[0];
  fitOff[0]=1;
  
  
  for(int ipmt=1; ipmt<NPMT; ++ipmt) fitNoBeam[ipmt] /= fitNoBeam[0];
  fitNoBeam[0]=1;



 TGraph* gcooper = new TGraph(NPMT);
 for(int ipmt=0; ipmt<NPMT; ++ipmt) {
   gcooper->SetPoint(ipmt,double(ipmt),qcooper[ipmt]);
 }
 
 
  TGraph* gasum = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gasum->SetPoint(ipmt,double(ipmt),asum[ipmt]);
  }
 

  TGraph* gfitOff = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gfitOff->SetPoint(ipmt,double(ipmt),fitOff[ipmt]);
  }
 
  TGraph* gfitNoBeam = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gfitNoBeam->SetPoint(ipmt,double(ipmt),fitNoBeam[ipmt]);
  }
 

  TString can3Name; can3Name.Form("%s-%s",tag.Data(),"relative-gain");
  //gSumAve->Print();
  cout << " making " << can3Name << endl;
  TCanvas *c3 = new TCanvas(can3Name,can3Name);
  gSumAve->GetHistogram()->GetXaxis()->SetTitle(" sum (blue squ) max (black tri) rc (red tri) fit off (gr circ)  fit nobeam (purp star)   : ipmt ");
  gSumAve->GetHistogram()->GetYaxis()->SetTitle(" relative gain ");
  gSumAve->SetMarkerSize(1);
  gSumAve->SetMarkerColor(kBlue);
  gSumAve->SetMarkerStyle(21);
  gMaxAve->SetMarkerSize(1);
  gMaxAve->SetMarkerStyle(22);
  gSumAve->SetTitle("gain");
  gSumAve->GetHistogram()->SetAxisRange(0,2.0,"Y");
  gSumAve->GetHistogram()->SetTitle("relative gain");
  

  gcooper->SetMarkerSize(1);
  gcooper->SetMarkerColor(kRed);
  gcooper->SetMarkerStyle(23);

  gfitOff->SetMarkerSize(1.3);
  gfitOff->SetMarkerColor(kGreen);
  gfitOff->SetMarkerStyle(20);

  gfitNoBeam->SetMarkerSize(1.5);
  gfitNoBeam->SetMarkerColor(6);
  gfitNoBeam->SetMarkerStyle(29);

  
  
  gSumAve->Draw("AP");
  gMaxAve->Draw("PSAME");
  gcooper->Draw("PSAME");
  gfitOff->Draw("PSAME");
  gfitNoBeam->Draw("PSAME");
  
}
