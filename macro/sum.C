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

  unsigned qnorm[NPMT];
  for(int ipmt=0; ipmt<NPMT; ++ipmt) qnorm[ipmt]=0;
  unsigned late=0;
  Int_t amonth=7, aday=31, amin=1555;
  double atime = double(amin)+double(24*aday+24*30*(amonth-7))-60*24*20;

  for(unsigned entry =0; entry < aSize; ++entry ) {
    sumTree->GetEntry(entry);
    string tag = pmtSum->tag;
    Int_t month = pmtSum->getMonth();
    Int_t day = pmtSum->getDay();
    Int_t min = pmtSum->getMin();
    Int_t seg = pmtSum->getSegment();
    double time = double(min)+double(60*24*day+60*24*30*(month-7)) - 60*24*20;

    ntTrig->Fill(time,pmtSum->ntrig555,pmtSum->ntrig5xx,pmtSum->ntrig444,pmtSum->ntrig4xx,pmtSum->ntrig111,pmtSum->ntrig1xx,pmtSum->ntrig000,pmtSum->ntrig0xx);
    
    if(entry%100==0) printf("...entry %u tag %s  month %i day %i min %i seg %i time %0.f \n",entry,tag.c_str(),month,day,min,seg,time);
    double etime=0;
    ++late;
    for(int ipmt=0; ipmt<NPMT; ++ipmt) {
      if( TMath::IsNaN(pmtSum->qsum[ipmt])  ) continue;
      if( TMath::IsNaN(pmtSum->qmax[ipmt])  ) continue;
      if(time>22600) { // average only late runs
        ++qnorm[ipmt];
        qmaxAve[ipmt] += pmtSum->qmax[ipmt];
        qsumAve[ipmt] += pmtSum->qsum[ipmt];
      }
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
  grsum[0]->GetHistogram()->GetXaxis()->SetTitle(" time in mins ");
  grsum[0]->SetTitle("sum charge");
  grsum[0]->Draw("AP");
  for(int ipmt=1; ipmt<NPMT; ++ipmt) grsum[ipmt]->Draw("PSAME");

  TString can2Name; can2Name.Form("%s-%s",tag.Data(),"peak");
  cout << " making " << can2Name << endl;
  TCanvas *c2 = new TCanvas(can2Name,can2Name);
  grpeak[0]->GetHistogram()->GetXaxis()->SetTitle(" time in mins ");
  grpeak[0]->SetTitle("peak charge");
  grpeak[0]->Draw("AP");
  for(int ipmt=1; ipmt<NPMT; ++ipmt) grpeak[ipmt]->Draw("PSAME");


  outfile->Write();

  printf(" averages over %u runs of %lld total \n",late,aSize);
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


  Double_t gain[NPMT];
  //PMT averages thrushold 1 side 0
 gain[0]=14.234858; 
 gain[1]=17.412120 ;
 gain[2]=14.131400 ;
 gain[3]=15.326554 ;
 gain[4]=20.247572 ;
 gain[5]=18.639789 ;
 gain[6]=13.448876 ;
 gain[7]=16.289118 ;
 gain[8]=14.963304 ;
 gain[9]=14.884147 ;
 gain[10]=14.343309 ;
 gain[11]=13.432883 ;
 gain[12]=10.920633 ;
 gain[13]=10.440021 ;
 gain[14]=11.908392 ;
 gain[15]=13.387776 ;
 gain[16]=12.343357 ;
 gain[17]=11.165986 ;
 gain[18]=17.630539 ;
 gain[19]=17.289387 ;
 gain[20]=14.914158 ;

 
  for(int ipmt=1; ipmt<NPMT; ++ipmt) gain[ipmt] /= gain[0];
  gain[0]=1;

   double fitOff[NPMT]; 
   // fit to QhitOff with ZERO gain
  // fits QhitOff
  //
  //
// fits QhitOff
 fitOff[0]=  46.216 ;
 fitOff[1]=  46.741 ;
 fitOff[2]=  43.052 ;
 fitOff[3]=  40.486 ;
 fitOff[4]=  45.727 ;
 fitOff[5]=  44.279 ;
 fitOff[6]=  37.363 ;
 fitOff[7]=  43.323 ;
 fitOff[8]=  39.452 ;
 fitOff[9]=  36.738 ;
 fitOff[10]=  36.328 ;
 fitOff[11]=  37.600 ;
 fitOff[12]=  29.199 ;
 fitOff[13]=  28.911 ;
 fitOff[14]=  28.590 ;
 fitOff[15]=  40.030 ;
 fitOff[16]=  30.553 ;
 fitOff[17]=  30.273 ;
 fitOff[18]=  50.722 ;
 fitOff[19]=  50.168 ;
 fitOff[20]=  36.464 ;
 
 double fitNoBeam[NPMT]; 
 // fit to QhitNoBeam with trigger000 and ZERO gain
// fits QhitNoBeam
 
// fits QhitNoBeam
 fitNoBeam[0]=  45.699 ;
 fitNoBeam[1]=  46.325 ;
 fitNoBeam[2]=  42.689 ;
 fitNoBeam[3]=  39.519 ;
 fitNoBeam[4]=  45.633 ;
 fitNoBeam[5]=  44.245 ;
 fitNoBeam[6]=  36.687 ;
 fitNoBeam[7]=  42.372 ;
 fitNoBeam[8]=  38.983 ;
 fitNoBeam[9]=  36.003 ;
 fitNoBeam[10]=  36.028 ;
 fitNoBeam[11]=  36.485 ;
 fitNoBeam[12]=  28.265 ;
 fitNoBeam[13]=  28.211 ;
 fitNoBeam[14]=  28.074 ;
 fitNoBeam[15]=  39.351 ;
 fitNoBeam[16]=  30.212 ;
 fitNoBeam[17]=  29.326 ;
 fitNoBeam[18]=  50.523 ;
 fitNoBeam[19]=  49.844 ;
 fitNoBeam[20]=  35.698 ;

  for(int ipmt=1; ipmt<NPMT; ++ipmt) fitOff[ipmt] /= fitOff[0];
  fitOff[0]=1;
  
  
  for(int ipmt=1; ipmt<NPMT; ++ipmt) fitNoBeam[ipmt] /= fitNoBeam[0];
  fitNoBeam[0]=1;


 Double_t fitRaw[NPMT];
 Double_t fitRawE[NPMT];

 fitRaw[0]=  16.585 ; fitRawE[0]= 0.067 ;
 fitRaw[1]=  19.137 ; fitRawE[1]= 0.066 ;
 fitRaw[2]=  16.014 ; fitRawE[2]= 0.071 ;
 fitRaw[3]=  17.370 ; fitRawE[3]= 0.080 ;
 fitRaw[4]=  21.603 ; fitRawE[4]= 0.093 ;
 fitRaw[5]=  20.474 ; fitRawE[5]= 0.102 ;
 fitRaw[6]=  14.843 ; fitRawE[6]= 0.051 ;
 fitRaw[7]=  19.574 ; fitRawE[7]= 0.121 ;
 fitRaw[8]=  17.317 ; fitRawE[8]= 0.109 ;
 fitRaw[9]=  16.498 ; fitRawE[9]= 0.084 ;
 fitRaw[10]=  15.760 ; fitRawE[10]= 0.077 ;
 fitRaw[11]=  15.610 ; fitRawE[11]= 0.075 ;
 fitRaw[12]=  12.375 ; fitRawE[12]= 0.075 ;
 fitRaw[13]=  12.508 ; fitRawE[13]= 0.077 ;
 fitRaw[14]=  13.323 ; fitRawE[14]= 0.074 ;
 fitRaw[15]=  15.080 ; fitRawE[15]= 0.056 ;
 fitRaw[16]=  13.809 ; fitRawE[16]= 0.076 ;
 fitRaw[17]=  13.593 ; fitRawE[17]= 0.146 ;
 fitRaw[18]=  19.926 ; fitRawE[18]= 0.068 ;
 fitRaw[19]=  19.017 ; fitRawE[19]= 0.062 ;
 fitRaw[20]=  16.804 ; fitRawE[20]= 0.088 ;

  for(int ipmt=1; ipmt<NPMT; ++ipmt) fitRaw[ipmt] /= fitRaw[0];
  fitRaw[0]=1;



  cout << endl;
 for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" gain[%i]= %.3f ;\n",ipmt,gain[ipmt]);
  cout << endl;
 for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" fitRaw[%i]= %.3f ;\n",ipmt,fitRaw[ipmt]);
  cout << endl;
 
  
  
 TGraph* gcooper = new TGraph(NPMT);
 for(int ipmt=0; ipmt<NPMT; ++ipmt) {
   gcooper->SetPoint(ipmt,double(ipmt),qcooper[ipmt]);
 }
 
 
  TGraph* ggain = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    ggain->SetPoint(ipmt,double(ipmt),gain[ipmt]);
  }
 

  TGraph* gfitOff = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gfitOff->SetPoint(ipmt,double(ipmt),fitOff[ipmt]);
  }
 
  TGraph* gfitNoBeam = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gfitNoBeam->SetPoint(ipmt,double(ipmt),fitNoBeam[ipmt]);
  }

  TGraph* gfitRaw = new TGraph(NPMT);
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gfitRaw->SetPoint(ipmt,double(ipmt),fitRaw[ipmt]);
  }
 
 

  TString can3Name; can3Name.Form("%s-%s",tag.Data(),"relative-gain");
  //gSumAve->Print();
  cout << " making " << can3Name << endl;
  TCanvas *c3 = new TCanvas(can3Name,can3Name);
  gcooper->GetHistogram()->GetXaxis()->SetTitle(" pmt number ");
  gcooper->GetHistogram()->GetYaxis()->SetTitle(" relative gain ");
  gSumAve->SetMarkerSize(1);
  gSumAve->SetMarkerColor(kBlue);
  gSumAve->SetMarkerStyle(21);
  gMaxAve->SetMarkerSize(1);
  gMaxAve->SetMarkerStyle(22);
  gcooper->SetTitle("gain");
  gcooper->GetHistogram()->SetAxisRange(.5,1.5,"Y");
  gcooper->GetHistogram()->SetTitle("relative gain");
  

  gcooper->SetMarkerSize(1);
  gcooper->SetMarkerColor(kRed);
  gcooper->SetMarkerStyle(23);

  ggain->SetMarkerSize(1.3);
  ggain->SetMarkerColor(kBlack);
  ggain->SetMarkerStyle(4);

  gfitNoBeam->SetMarkerSize(1.5);
  gfitNoBeam->SetMarkerColor(6);
  gfitNoBeam->SetMarkerStyle(28);

  gfitNoBeam->SetMarkerSize(1.5);
  gfitNoBeam->SetMarkerColor(kGreen);
  gfitNoBeam->SetMarkerStyle(29);
  
  gfitRaw->SetMarkerSize(1.5);
  gfitRaw->SetMarkerColor(kGreen);
  gfitRaw->SetMarkerStyle(29);

  
  //gSumAve->Draw("AP");
  //gMaxAve->Draw("PSAME");
  gcooper->Draw("AP");
  //gfitOff->Draw("PSAME");
  gfitRaw->Draw("PSAME");
  ggain->Draw("PSAME");
    
}
