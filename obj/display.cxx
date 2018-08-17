//Very simple functions to help plot PMT waveforms
#include "display.hxx"
ClassImp(display)

display::display(): TNamed("TPmtHit","TPmtHit")
{
  histosetup = false;
  prefilter = true, postfilter = true;
  chain=NULL;
  Mask0=255;
  Mask1=255;
  Mask2=255;
  c1=NULL;
  tagList.clear();
    
  hsum = new TH1F("hsum", "PMT Sum;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
  hsumDraw = new TH1F("hsumDraw", "PMT Sum;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
    
  hsumDraw->SetLineColor(kRed);
  
  for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
          histo[iB][iC]=NULL;
          histoDraw[iB][iC]=NULL;
        }
  }

  SetupHistos(50.);
  
}
display::~display()
{
  if (!chain) return;
  delete chain->GetCurrentFile();
}

void display::TestRange(TH1F *histogram, int &startBin, int &stopBin) {
    
    if (stopBin < startBin) {
        int temp = stopBin;
        stopBin = startBin;
        startBin = temp;
    }
    
    int NBinsX = histogram->GetNbinsX();
    
    if (startBin < 1) { startBin = 1; }
    if (stopBin > NBinsX) { stopBin = NBinsX; }
    
}

void display::SetFilter(bool forBaseline, bool forDisplay) {
    
    prefilter = forBaseline;
    postfilter = forDisplay;
    
}

void display::OpenFile(const char *infile, const char *treename ) {
    
    tf = new TFile(infile, "READ");
    pmt_tree = (TTree *)tf->Get(treename);
    pmt_tree->SetBranchAddress("digitizer_waveforms", &waveforms);

}
void display::AddTag(TString tag, const char *treename) 
{   
  cout << "AddTag tag "<< tag << endl;
  if (!chain) { chain = new TChain(treename); } 
  TString fileName = TString("pdsData/PDSout_") + tag + TString(".root");
  int numfiles = chain->Add(fileName);
  chain->SetBranchAddress("digitizer_waveforms", &waveforms);
  if (numfiles == 0) { cout << "No files added" << endl; }
  cout << " chain has " << chain->GetEntries() << " entries " << endl;
}

void display::SetTagList(){
  TObjArray *fileElements=chain->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=NULL;
  while (( chEl=(TChainElement*)next() )) {
    string fname = string(chEl->GetTitle());
    tagList.push_back(fname.substr( fname.find("_")+1, fname.find(".") -1  - fname.find("_")));
  }
  for(unsigned it=0; it< tagList.size() ; ++ it) cout << it << "  " << tagList[it] << endl;
}

void display::AddChain(const char *infiles, const char *treename ) {
    if (!chain) { chain = new TChain(treename); }
    int numfiles = chain->Add(infiles);
    chain->SetBranchAddress("digitizer_waveforms", &waveforms);
    if (numfiles == 0) { cout << "No files added" << endl; }
    cout << " chain has " << chain->GetEntries() << " entries " << endl;
    //chain->ls();
}

void display::ClearChain() {
    
    delete chain;
    chain = 0;
}

void display::ClearTree() {
    
    delete pmt_tree;
    pmt_tree = 0;
}

void display::SetupHistos(float offsetstepADC) {
    
    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            if (!histo[iB][iC]) {
                histo[iB][iC] = new TH1F(Form("B%iC%i_0",iB, iC), Form("B%i C%i",iB, iC), NSAMPLES, 0., 1.*NSAMPLES);
                histo[iB][iC]->SetStats(0);
                
                histoDraw[iB][iC] = new TH1F(Form("B%iC%i",iB, iC), Form("B%i C%i",iB, iC), NSAMPLES, 0., 1.*NSAMPLES);
                histoDraw[iB][iC]->SetStats(0);
                
                if (iC ==7) {
                    histoDraw[iB][iC]->SetLineColor(kGreen);
                }
                
            }
        }
    }
    
    hmaster = new TH1F("hmaster", "ADC Plot Event -1;Sample (4 ns per sample);ADC Counts", NSAMPLES, 0., 1.*NSAMPLES);
    hmaster->GetYaxis()->SetRangeUser(-4.*offsetstepADC, 24*offsetstepADC);
    hmaster->GetYaxis()->SetTitleOffset(1.30);
    hmaster->SetStats(0);
    
    histosetup = true;
    gDirectory->ls();    
}

void display::ADCfilter(int iB, int iC) {
    
    for (int iS = 0; iS<NSAMPLES; ++iS) {
        if (waveforms[iB][iC][iS] > ADCrange) {
            if (iS > 0) { waveforms[iB][iC][iS] = waveforms[iB][iC][iS-1];}
            else {
                int iS2 = 0;
                while (waveforms[iB][iC][iS2] > ADCrange) {
                    waveforms[iB][iC][0] = waveforms[iB][iC][iS2+1];
                    ++iS2;
                }
            }
        }
    }
}

float display::GetBaseline(int iB, int iC, int baselinestart, int baselinewide) {
    float baseline = 0.;
    for (int iS = baselinestart; iS<(baselinestart+baselinewide); ++iS) {
        baseline += waveforms[iB][iC][iS]/(1.*baselinewide);
    }
    
    return baseline;
    
}

float display::GetBaseline(TH1F *h1, int baselinestart, int baselinewide) {
    
    int baselinestop = baselinestart+baselinewide;
    TestRange(h1, baselinestart, baselinestop);
    
    float baseline = 0.;
    for (int iS = baselinestart; iS<baselinestop; ++iS) {
        baseline += h1->GetBinContent(iS)/(1.*baselinewide);
    }

    return baseline;
}

float display::GetSigma(TH1F *h1, int baselinestart, int baselinewide) {
    
    int baselinestop = baselinestart+baselinewide;
    TestRange(h1, baselinestart, baselinestop);
    
    float tot = 0., tot2 = 0.;
    for (int iS = baselinestart; iS<baselinestop; ++iS) {
        tot += h1->GetBinContent(iS)/(1.*baselinewide);
        tot2 += h1->GetBinContent(iS)*h1->GetBinContent(iS)/(1.*baselinewide);
    }
    
    return sqrt(tot2-tot*tot);
}

bool display::FillHistos(int EvNum, float offsetstepADC ){ 

    printf(" fill histos for %i \n",EvNum);
    
    if (!histosetup) { SetupHistos(offsetstepADC); }
    
    int NumEvents = 0;
    
    // default to a chain
    if (chain) { NumEvents = chain->GetEntries(); }
    else if (pmt_tree) { NumEvents = pmt_tree->GetEntries(); }
    else { cout << "No TChain or TTree" << endl; return false; }
    
    if (EvNum >= NumEvents) { cout << "Entry " << EvNum << " > " << NumEvents << " entries" << endl; return false; }
    
    // default to a chain
    if (chain) { chain->GetEntry(EvNum); }
    else if (pmt_tree) { pmt_tree->GetEntry(EvNum); }
    
    // now fill histograms
    float baseline = 0.;
    float offset = 0.;
    bool needS0 = false;

    float sumAll=0;

    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            
            if (prefilter) { ADCfilter(iB, iC); }
            baseline = GetBaseline(iB, iC);
            if (postfilter) { ADCfilter(iB, iC); }
            
            histo[iB][iC]->Reset();
            histoDraw[iB][iC]->Reset();
            float sum=0;
            for (int iS = 0; iS<NSAMPLES; ++iS) {
                histo[iB][iC]->Fill(iS+0.5, (1.*waveforms[iB][iC][iS]-baseline) );
                if (iC == 7) {
                    histoDraw[iB][iC]->Fill(iS+0.5, ((1.*waveforms[iB][iC][iS]-baseline)*offsetstepADC/(1.*ADCrange+1.)+offset) );
                }
                else {
                    histoDraw[iB][iC]->Fill(iS+0.5, (1.*waveforms[iB][iC][iS]-baseline+offset) );
                  sum-= float(waveforms[iB][iC][iS]-baseline);
                  sumAll-= float(waveforms[iB][iC][iS]-baseline);
                }

            }
            if(iC!=7) printf(" ib %i ic %i sum %f \n",iB,iC,sum);
            
            offset += offsetstepADC;
        }
    }
    UInt_t evFile  = EvNum%5000;  
    UInt_t iFile = EvNum/5000;
    hmaster->SetTitle(Form("ADC Plot File %s Event %i;Sample (4 ns per sample);ADC Counts",tagList[iFile].c_str(),evFile));
    printf(" %s %f \n",hmaster->GetTitle(),sumAll);
    std::vector<Int_t> rftimes0 =getRFTimes(0); 
    std::vector<Int_t> rftimes1 =getRFTimes(1); 
    std::vector<Int_t> rftimes2 =getRFTimes(2); 
    Int_t rf0=-1;
    Int_t rf1=-1;
    Int_t rf2=-1;
    if(rftimes0.size()>0) rf0 = rftimes0[0];
    if(rftimes1.size()>0) rf1 = rftimes1[0];
    if(rftimes2.size()>0) rf2 = rftimes2[0];
    //printf(" sizes of rftimes %u %u %u \n",rftimes0.size(),rftimes1.size(),rftimes2.size());
    printf(" rftimes %i %i %i  \n",rf0,rf1,rf2);
    
    return true;

}

std::vector<Int_t> display::getRFTimes(int ib) 
{
  printf(" getting RF times board %i \n",ib);
  std::vector<Int_t> rftimes;
  rftimes.clear();
  int ic=7;

  // find baseline
  std::vector<UShort_t> udigi; 
  for (UInt_t is=0; is<NSAMPLES; ++is) udigi.push_back(waveforms[ib][ic][is]);
  std::sort(udigi.begin(), udigi.end());
  UShort_t baseline = udigi[0.5*double(NSAMPLES)];

  // looking for negative values.  
  UShort_t digiMin=MAXADC;
  for (UInt_t is=0; is<NSAMPLES; ++is) {
    waveforms[ib][ic][is]=TMath::Min( baseline,waveforms[ib][ic][is]);
    if(waveforms[ib][ic][is]<digiMin) digiMin=waveforms[ib][ic][is];
  }

  double step = double(digiMin) - double(baseline);
  // return if step down is too small
  if(step>-500) {
    printf(" step too small %f %u \n",step,rftimes.size());
    return rftimes;
  }
  // pick off start of rising edge
  bool isRF=false;
  for (UInt_t is=0; is<NSAMPLES; ++is){
    double digi = double(waveforms[ib][ic][is]) - double(baseline);
    if(digi<0.75*step&&!isRF) {
      rftimes.push_back(is);
      isRF=true;
    } else if(digi>0.75*step) {
      isRF=false;
    }
  }
  printf(" board %i has %u rftimes\n",ib,rftimes.size());
  for(unsigned it =0; it< rftimes.size(); ++it ) printf(" board %i %i %i",ib,it,rftimes[it]);
  return rftimes;
}

bool display::SumEvents(bool trim , float offsetstepADC) {
    
    hsum->Reset();
    hsumDraw->Reset();
    
    int sample;
    
    for (int iB = 0; iB<NBOARDS; ++iB) {
        for (int iC = 0; iC<NPMTS; ++iC) {
            if (iC == 7) { continue; }
            
            for (int iS = 1; iS<=NSAMPLES; ++iS) {
                
                if (trim) {
                    sample = iS + Trims[iB][iC];
                    if (sample < 1) { sample = 1; }
                    else if (sample > NSAMPLES ) { sample = NSAMPLES; }
                }
                else { sample = iS; }
                hsum->Fill(iS-0.5, histo[iB][iC]->GetBinContent(sample) );
                hsumDraw->Fill(iS-0.5, histo[iB][iC]->GetBinContent(sample) );
            }
        }
    }
    
    for (int iS = 1; iS<=NSAMPLES; ++iS) { hsumDraw->Fill(iS-0.5, -offsetstepADC); }

    return true;
}

bool display::DrawEvent(int EvNum, int mask0 , int mask1 , int mask2 , bool showSum, float offsetstepADC ) 
{
    
  bool status = FillHistos(EvNum, offsetstepADC);
  status &= SumEvents(false, offsetstepADC);

  if (!status) { return false; }
  int mask = 0;

  bool redraw = false;
  if ( (mask0 != Mask0) || (mask1 != Mask1) || (mask2 != Mask2) ) { redraw = true; }
  Mask0 = mask0; Mask1 = mask1; Mask2 = mask2;

  if (!c1) { c1 = new TCanvas("c1", "", cWidth, cHeight); }
  UInt_t evFile  = EvNum%5000;  
  UInt_t iFile = EvNum/5000;
  c1->SetName(Form("Run%sEvent%u",tagList[iFile].c_str(),evFile));
  c1->SetTitle(Form("ADC Plot File %s Event %i;Sample (4 ns per sample);ADC Counts",tagList[iFile].c_str(),evFile));
  hmaster->Draw();

  for (int iB = 0; iB<NBOARDS; ++iB) {

    if (iB == 0) { mask = Mask0; }
    else if (iB == 1) { mask = Mask1; }
    else if (iB == 2) { mask = Mask2; }

    for (int iC = 0; iC<NPMTS; ++iC) {

      if (mask & (1 << iC)) {
        histoDraw[iB][iC]->Draw("hist same");
      }
    }
  }
    
  if (showSum) { hsumDraw->Draw("hist same"); }
    
  //c1->Modified(); c1->Update();
    
  return true;
    
}

float display::EstimateTS(int EvNum, float threshold) {
    
  float peak=0, total=0, F90=0;
  int pulse0 = 0;
  int iS;

  bool status = FillHistos(EvNum);
  status &= SumEvents(false);

  if (!status) { return -1.; }

  for (iS = 400; iS < 550; ++iS) {
    if (hsum->GetBinContent(iS) < -threshold) {
      pulse0 = iS;
      break;
    }
  }
    
  if (iS == 550) { return -1; }

  for (int iS = pulse0-5; iS < pulse0 + 23; ++iS) {

    peak += hsum->GetBinContent(iS);
    total += hsum->GetBinContent(iS);
  }
  for (int iS = pulse0 + 23; iS < pulse0+1200; ++iS) total += hsum->GetBinContent(iS);
  

  if(total>0) F90 = peak / total;

  printf("EstimateTS event %i %f \n",EvNum,F90); 

  return F90;

}

int display::GetT0Bin(TH1F *histogram, float threshold, int startBin, int stopBin) {
    

    TestRange(histogram, startBin, stopBin);
    
    for (int iS = startBin; iS <= stopBin; ++iS) {
        if (histogram->GetBinContent(iS) < -threshold) { return iS; }
    }
    
    return -1;
    
}

int display::GetMinimum(TH1F *histogram, float &minVal, int startBin, int stopBin, bool StepForward) {
    
    int iStart = startBin, iStop = stopBin, step = 1;
    TestRange(histogram, startBin, stopBin);
    
    if (!StepForward) { iStop = startBin; iStart = stopBin; step = -1; }
    
    minVal = ADCrange;
    float val = ADCrange;
    int minBin = -1;
    
    for (int iS = iStart; iS <= iStop; iS += step) {
        val = histogram->GetBinContent(iS);
        if (val < minVal) {
            minVal = val;
            minBin = iS;
        }
    }
    
    return minBin;
    
}
