#include "lowAna.hh"
#include <TF1.h>
#include <TPaveStats.h>
#include <TText.h>

lowAna::lowAna(Int_t maxLoop, Int_t firstEntry)
{
  if(readGainConstants()==0) {
    printf(" cannot read gain constants file so abort \n");
    return;
  }
 
  fChain=NULL;
  TString fileName("pdsOutput/pmtChainLow.root");
  printf(" looking for file %s\n",fileName.Data());
  TFile *f = new TFile(fileName,"readonly");
  if(f->IsZombie()) {
    printf(" couldnt open file %s so abort.\n",fileName.Data());
    return;
  }
  TTree *tree;
  f->GetObject("pdsTree",tree);
  tree->ls();
  Init(tree);
  if(!fChain) return;
  Long64_t nentries = fChain->GetEntries();

  printf(" pdsTree has %lld entries \n",nentries);
  
  // initicalize fft 
  nFFTSize = int(MAXSAMPLES);
  fFFT = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");

 /* 
  TString gainFileName = TString("pdsOutput/pmtGains_")+tag+ TString(".root");
  gainFile = new TFile(gainFileName,"recreate");
  gainFile->cd();
  printf(" opening summary file %s \n",summaryFileName.Data());
  TTree *gainsTree = new TTree("gainsTree","gainsTree");
  pmtGains  = new TPmtGains();
  gainsTree->Branch("pmtGains",&pmtGains);
  pmtGains->tag=tag;
  */
  
  // open ouput file and make some histograms
  TString outputFileName;
  outputFileName.Form("pdsOutput/lowAna-%i-%i.root",maxLoop,firstEntry);
  TString(".root");
  outFile = new TFile(outputFileName,"recreate");
  //promptDir = outFile->mkdir("promptDir");
  outFile->cd();
  printf(" opening output file %s \n",outputFileName.Data());
  summaryTree = new TTree("summaryTree","summaryTree");
  pmtSummary  = new TPmtSummary();
  summaryTree->Branch("pmtSummary",&pmtSummary);
  

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);

  

  //pmtTree->ls();
  //
  //ntuples
  ntTrig = new TNtuple("ntTrig"," trigger ","run:event:r1:r2:r3:t1:t2:t3");

  promptDir =  outFile->mkdir("promptDir");
  histDir =  outFile->mkdir("histDir");
  histDir->cd();
  // histos 
  hOcc =  new TH1D("occupancy","occupancy by pmt",NPMT,0,NPMT);
  hOcc->SetXTitle(" pmt number ");
  hOcc->SetYTitle(" hits per event ");
  hNoise = new TH1D("noise","baseline subtracted noise by pmt",NPMT,0,NPMT);
  hNoise->SetXTitle(" pmt number ");
  hBase = new TH1D("base","baseline by pmt",NPMT,0,NPMT);
  hBase->SetXTitle(" pmt number ");
  hBase->Sumw2();


  TString hname;
  TString htitle;


  hTPrompt = new TH1D("TPrompt"," peak of charge weighted pulse times",MAXSAMPLES+500,-500,MAXSAMPLES);
  hTPrompt->SetXTitle(" prompt peak (sample time) ");

  hTPromptEvent = new TH1D("TPromptEvent"," peak of charge weighted pulse times, single event",MAXSAMPLES,0,MAXSAMPLES);
  hTPromptEvent->SetXTitle(" prompt peak (sample time) ");


  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0) continue;
      hname.Form("Samples_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Samples board%u channel%u pmt%u",ib,ic,ipmt);
      hSamples[ipmt] = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
      hSamples[ipmt]->SetXTitle(" sample number ");
    }
  }
  hSamplesSum = new TH1D("SampleSum"," samples summed over PMTs",MAXSAMPLES,0,MAXSAMPLES);
  hSamplesSum->SetXTitle(" sample number ");

  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("SamplesPDS_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Samples board%u channel%u pmt%u",ib,ic,ipmt);
      hSamplesPDS[ipmt] = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
      hSamplesPDS[ipmt]->SetXTitle(" sample number ");
    }
  }
  hSamplesPDSSum = new TH1D("SampleSumPDS"," samples summed over PMTs",MAXSAMPLES,0,MAXSAMPLES);
  hSamplesPDSSum->SetXTitle(" sample number ");


  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("Peaks_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Peaks board%u channel%u pmt%u",ib,ic,ipmt);
      hPeaks[ipmt] = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
      hPeaks[ipmt]->SetXTitle(" sample number ");
      hPeaks[ipmt]->SetLineColor(kRed);
      hPeaks[ipmt]->SetMarkerColor(kRed);
      hPeaks[ipmt]->SetFillColor(kRed);
      hPeaks[ipmt]->SetFillStyle(3002);

      hname.Form("Sum_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" counts board%u channel%u pmt %u",ib,ic,ipmt);
      hCounts[ipmt] = new TH1D(hname,htitle,500,0,5000);
      hCounts[ipmt]->SetXTitle(" baseline subtracted summed ADC counts ");

      hname.Form("Baseline_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" Baseline board%u channel%u pmt %u",ib,ic,ipmt);
      hBaseline[ipmt] = new TH1D(hname,htitle,100,-50,50);
      hBaseline[ipmt]->SetXTitle(" baseline fluctuation (ADC counts) ");

      hname.Form("HitQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" hit charge board%u channel%u pmt %u",ib,ic,ipmt);
      hHitQ[ipmt] = new TH1D(hname,htitle,60,0,60);
      hHitQ[ipmt]->SetXTitle(" hit charge in pulse time  (ADC counts) ");

      hname.Form("RawQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      //htitle.Form(" Raw Charge board%u channel%u pmt%u",ib,ic,ipmt);
      htitle.Form("");
      hRawQ[ipmt] = new TH1D(hname,htitle,60,0,60);
      hRawQ[ipmt]->SetXTitle(Form(" Raw Charge [b%u c%u pmt%u] (ADC)",ib,ic,ipmt));
      hRawQ[ipmt]->SetYTitle(" Entries / [ADC]");


      hname.Form("QMax_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" max ADC c board%u channel%u pmt %u",ib,ic,ipmt);
      hQMax[ipmt] = new TH1D(hname,htitle,50,0,100);
      hQMax[ipmt]->SetXTitle(" q max (ADC counts) ");

      hname.Form("NHits_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" Hits per event c board%u channel%u pmt %u",ib,ic,ipmt);
      hNHits[ipmt] = new TH1D(hname,htitle,20,0,20);
      hNHits[ipmt]->SetXTitle(" number of hits per event ");
    }
  }
  //gDirectory->ls();

  /***  loop over entries zero = all ***/
  UInt_t nLoop = Loop(maxLoop,firstEntry);

  //gainsTree->Fill();
  //gainFile->Write();
  pmtSummary->print();
  summaryTree->Fill();

  outFile->Write();
  printf(" wrote output file %s \n",outFile->GetName());
  //printf(" wrote gains file %s \n",gainFile->GetName());

}


lowAna::~lowAna()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

UInt_t lowAna::Loop(UInt_t nToLoop,UInt_t firstEntry)
{
  if (fChain == 0) return 0;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes=0;
  std::vector<Double_t> sdigi;  // source
  std::vector<Double_t> ddigi;  // baseline subtracted
  std::vector<Double_t> fdigi;  // baseline subtracted and filterd

  // no gain applied
  std::vector<Double_t> sdigiUn;  // source
  std::vector<Double_t> ddigiUn;  // baseline subtracted

  UInt_t nloop=nentries;
  if(nToLoop!=0) nloop = nToLoop;
  printf(" entries %lld looping %d first %d \n",nentries,nloop,firstEntry);
  std::vector<Double_t> vpromptLike;
  Int_t currentRun = -1;
  // loop over entries
  for (Long64_t jentry=firstEntry; jentry<nloop+firstEntry; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) { printf(" load tree returns %lld\n",ientry); break;}
    nbytes += fChain->GetEntry(jentry);
    if(jentry%1000==0) printf(" \t.... entry %lld nbytes %lld pmtTree entries %lld \n",jentry,nbytes,pmtTree->GetEntries());
    // is this a new run?
    Int_t thisRun = Int_t(gps_ctrlFlag);
    std::string tag=std::string("07-31-")+to_string(int(gps_nsIntoSec))+std::string("-")+to_string(int(gps_secIntoDay)); 
    if( thisRun!= currentRun) {
      if(currentRun>=0) {
        qualitySummary();
        pmtSummary->print();
        summaryTree->Fill();
      }
      currentRun=thisRun;
      pmtSummary->clear();
      pmtSummary->tag = tag;
      pmtSummary->run=currentRun;
      pmtSummary->min=int(gps_nsIntoSec);
      pmtSummary->seg=int(gps_secIntoDay);
    }

    // clear the event
    pmtEvent->clear();
   
    // save event info
    // construct the tag 
    pmtEvent->tag = tag; 
    pmtEvent->run= currentRun;
    pmtEvent->event=event_number;
    pmtEvent->compSec=computer_secIntoEpoch;
    pmtEvent->compNano=computer_nsIntoSec;
    // trigger type
    pmtEvent->trigType = triggerInfo();
    if(event_number%5000==0) printf(" run %i event %i tag %s \n",pmtEvent->run,pmtEvent->event,(pmtEvent->tag).c_str());
    //pmtEvent.tpcTrig;
    //pmtEvent.pdsTrig;
    pmtEvent->rft21=rftime21;
    pmtEvent->rft22=rftime22;
    pmtEvent->rft23=rftime23;
     // summary info
    pmtSummary->vsec.push_back(computer_secIntoEpoch);
    pmtSummary->vnano.push_back(computer_nsIntoSec);
    pmtSummary->vdtime1.push_back(digitizer_time[0]);
    pmtSummary->vdtime2.push_back(digitizer_time[1]);
    pmtSummary->vdtime3.push_back(digitizer_time[2]);
    pmtSummary->vtrig.push_back(pmtEvent->trigType);
    pmtSummary->vevent.push_back(event_number);
    pmtSummary->ventry.push_back(jentry);    
    pmtSummary->vcompSec.push_back(computer_secIntoEpoch);
    pmtSummary->vcompNano.push_back(computer_nsIntoSec);
    
  
    UInt_t rftime[3];
    rftime[0]=0; if(rftime21.size()>0) rftime[0]=UInt_t(rftime21[0]);
    rftime[1]=0; if(rftime22.size()>0) rftime[1]=UInt_t(rftime22[0]);
    rftime[2]=0; if(rftime23.size()>0) rftime[2]=UInt_t(rftime23[0]);


    for(UInt_t ib=0; ib<NB; ++ib) {
      UInt_t time = digitizer_time[ib];
      //printf(" board %u time %u \n",ib,time);
      for(UInt_t ic=0; ic<NC; ++ic) {

        // filter waveforms for stuck bits
        ADCFilter(ib,ic);
        // get pmt number
        int ipmt = toPmtNumber(ib,ic);
        if(ipmt<0||ipmt>=NPMT) continue;

        // make a vector of samples for sorting.
        sdigi.clear();
        ddigi.clear();
        fdigi.clear();
        sdigiUn.clear();
        ddigiUn.clear();

        double sum=0;
        double sumUn=0;

        // Find the sample median and it's "sigma".
        for (UInt_t is=0; is<MAXSAMPLES; ++is) {
          sdigi.push_back(double(digitizer_waveforms[ib][ic][is])/gain[ipmt]);
          sdigiUn.push_back(double(digitizer_waveforms[ib][ic][is]));
        }

        std::sort(sdigi.begin(), sdigi.end());
        double baselineMedian = sdigi[0.5*double(MAXSAMPLES)];
        double baselineSigma = sdigi[0.16*double(MAXSAMPLES)];
        baselineSigma = std::abs(baselineSigma-baselineMedian);

        std::sort(sdigiUn.begin(), sdigiUn.end());
        double baselineMedianUn = sdigiUn[0.5*double(MAXSAMPLES)];
        double baselineSigmaUn = sdigiUn[0.16*double(MAXSAMPLES)];
        baselineSigmaUn = std::abs(baselineSigmaUn-baselineMedianUn);


        //baselineSigma = std::abs(baselineSigma-baselineMedian);
        //noise = sdigi[0.68*sdigi.size()];/
        hBase->SetBinContent(ipmt+1,hBase->GetBinContent(ipmt+1)+baselineMedian);
        hBase->SetBinError(ipmt+1,hBase->GetBinError(ipmt+1)+baselineSigma);
        double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);
        hNoise->SetBinContent(ipmt+1,hNoise->GetBinContent(ipmt+1)+noise);
        if(ientry==0) baselineNominal[ipmt]= baselineMedian;
        else hBaseline[ipmt]->Fill(baselineMedian-baselineNominal[ipmt]);

        //if(ientry==0) hFFT[ipmt]=FFTFilter(ipmt);

        UInt_t tmax=0;
        double qmax=0;
        Double_t qrf=0;
        UInt_t tmaxUn=0;
        double qmaxUn=0;
        UInt_t rfLow=0;
        if( Int_t(rftime[ipmt])-100 >= 0) rfLow = rftime[ipmt];
        for(UInt_t is=0 ; is<MAXSAMPLES; ++is) {
          double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])/gain[ipmt]-baselineMedian);
          if(digi>qmax) {
            qmax=digi;
            tmax=is+1;
          } 
          // RF window sum 
          if(is>rfLow&&is<rftime[ipmt]+100) qrf += digi; 
          // witout gain
          double digiUn = -1.0*(double(digitizer_waveforms[ib][ic][is])-baselineMedianUn);
          if(digiUn>qmaxUn) {
            qmaxUn=digiUn;
            tmaxUn=is+1;
          }

          ddigi.push_back(digi);
          ddigiUn.push_back(digiUn);
          if(pmtEvent->trigType == TPmtEvent::TRIG000) hSamplesPDS[ipmt]->SetBinContent(int(is+1),hSamplesPDS[ipmt]->GetBinContent(int(is+1))+digi);
          else hSamples[ipmt]->SetBinContent(int(is+1),hSamples[ipmt]->GetBinContent(int(is+1))+digi);
          // here I am not worrying about the difference between noise and gain-corrected noise.  just using gain-corrected noise
          if(digi>3.0*noise) { 
            //if(is>450&&is<470) 
            sum+=digi;
            sumUn+=digiUn;
          }
        } // loop over digitizations

        // filtered 
        fdigi = MovingAverageFilter(ddigi,5); // parameter is top hat window which should be odd

        // make some single event sample plots 
        if(jentry<2) {
          histDir->cd();
          TString hname,htitle;
          hname.Form("Digi_pmt%i_ev%i",ipmt,Int_t(ientry));
          htitle.Form("digi samples pmt %i ev %i",ipmt,Int_t(ientry));
          TH1D* hDigi = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
          hname.Form("FDigi_pmt%i_ev%i",ipmt,Int_t(ientry));
          htitle.Form("filtered digi samples pmt %i ev %i",ipmt,Int_t(ientry));
          TH1D* hFDigi = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
          for(unsigned idigi=0; idigi<ddigi.size(); ++idigi) {
            hDigi->SetBinContent(int(idigi+1),ddigi[idigi]);
            hFDigi->SetBinContent(int(idigi+1),fdigi[idigi]);
          }
        }


        //if(sum>400&&pmtEvent->trigType!=7) printf(" SSSSSSSSS event %lli trig %i sum %f   \n", jentry, pmtEvent->trigType, sum );
        hCounts[ipmt]->Fill(sum);
        // peak finding
        //printf("calling find peaks event %i pmt %i noise %f 68v %f base %f \n",event_number,ipmt,noise,sdigi[0.68*sdigi.size()], baselineMedian);
        std::vector<Int_t> peakTime = findPeaks(ddigi,THRESHOLDHIGH,THRESHOLDLOW);
        //std::vector<Int_t> peakTime = findMaxPeak(fdigi,8.0*noise,3.0*noise);
        Int_t nhits = findHits(ipmt,sum,peakTime,ddigi,ddigiUn,pmtEvent->trigType);
        // now fill in the time since last RF pulse
        hOcc->Fill(ipmt+1,nhits);
        hNHits[ipmt]->Fill(nhits);
        //printf(" event %i nhits %i \n", pmtEvent->event, pmtEvent->nhits );
        for (UInt_t ip = 0; ip < peakTime.size(); ip++) {
          Int_t bin = peakTime[ip];
          //printf(" ipmt %i ip %i bin %i v %f \n",ipmt,ip,bin,fdigi[bin]);
          hPeaks[ipmt]->SetBinContent(bin+1, hPeaks[ipmt]->GetBinContent(bin+1)+ddigi[bin]);
        }
        hQMax[ipmt]->Fill(qmax);
        //ntPmt->Fill(double(pmtEvent->trigType),
        //    double(ipmt),tmax,qmax,sum,tmaxUn,qmaxUn,sumUn,noise,baselineMedian-baselineNominal[ipmt],nhits,qrf);
        pmtEvent->qmax.push_back(qmax);
        pmtEvent->qsum.push_back(sum);

      } // channel loop 
      getTimeToRF(ib);
    } // board loop 
    /*** after filling hits, get prompt time ****/
    pmtEvent->tPrompt=getPromptTime();
    Double_t tof=pmtEvent->tPrompt*4.0-GAMMAPEAK;//in ns 
    if(pmtEvent->tPrompt!=-9999 && pmtEvent->trigType==TPmtEvent::TRIG111 && tof>0){
      pmtSummary->tprompt.push_back(pmtEvent->tPrompt*4.0);//in ns with respect to rf 
      cout<<jentry<<" tprompt = "<<pmtEvent->tPrompt*4.0<<endl;
      pmtSummary->tof.push_back(tof);//in ns
      pmtSummary->ke.push_back(nmass*(sqrt(1/((tof*clight/L)*(tof*clight/L)-1)+1)-1));//in MeV
    }
    else
    {
      pmtSummary->tprompt.push_back(-9999);//in ns
      pmtSummary->tof.push_back(-9999);//in ns
      pmtSummary->ke.push_back(-9999);//in MeV
    }//ysun
    //pmtSummary->vprompt.push_back(pmtEvent->tPrompt);
    
    // do each board separatly 
    hTPrompt->Fill(pmtEvent->tPrompt);
    pmtEvent->nhits= pmtEvent->hit.size();
    pmtTree->Fill();
    //if(jentry%1000==0) printf(" \t\t jentry %lli nhits = %d \n",jentry,pmtEvent->nhits);
    //if(jentry%1000==0) pmtEvent->print();
  }   // end loop over entries
  printf(" finised looping  %u pmtTree size %llu \n",nloop,pmtTree->GetEntries());
  // normalize
  for(Int_t ipmt=0; ipmt<NALLCH; ++ipmt) {
    //UInt_t sampleNorm = hSamples[ipmt]->GetEntries();
    for(int ibin=1; ibin<= hSamples[ipmt]->GetNbinsX()+1; ++ibin ){   
      hSamples[ipmt]->SetBinContent(ibin, hSamples[ipmt]->GetBinContent(ibin)/Double_t(nloop));
    }
  }
  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    for(int ibin=1; ibin<= hSamples[ipmt]->GetNbinsX()+1; ++ibin ){   
      hSamplesPDS[ipmt]->SetBinContent(ibin, hSamplesPDS[ipmt]->GetBinContent(ibin)/Double_t(nloop));
    }
  }


  // sum
  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    for(int ibin=1; ibin<= hSamples[ipmt]->GetNbinsX()+1; ++ibin ){   
      hSamplesSum->SetBinContent(ibin, hSamplesSum->GetBinContent(ibin) + hSamples[ipmt]->GetBinContent(ibin) );
      hSamplesPDSSum->SetBinContent(ibin, hSamplesPDSSum->GetBinContent(ibin) + hSamplesPDS[ipmt]->GetBinContent(ibin) );
    }
  }



  for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
    hBase->SetBinContent(ipmt+1, hBase->GetBinContent(ipmt+1)/Double_t(nloop));
    hBase->SetBinError(ipmt+1, hBase->GetBinError(ipmt+1)/Double_t(nloop));

    UInt_t peakNorm = hPeaks[ipmt]->GetEntries();
    for(int ibin=1; ibin<= hPeaks[ipmt]->GetNbinsX()+1; ++ibin ) {  
      hPeaks[ipmt]->SetBinContent(ibin, hPeaks[ipmt]->GetBinContent(ibin)/Double_t(peakNorm));
    }

    for(int ibin=1; ibin<=  hCounts[ipmt]->GetNbinsX()+1; ++ibin ) 
      hCounts[ipmt]->SetBinContent(ibin,  hCounts[ipmt]->GetBinContent(ibin)/Double_t(nloop));

  }
  for(int ibin=1; ibin<=  hNoise->GetNbinsX()+1; ++ibin ) hNoise->SetBinContent(ibin,  hNoise->GetBinContent(ibin)/Double_t(nloop));
  for(int ibin=1; ibin<=  hOcc->GetNbinsX()+1; ++ibin ) hOcc->SetBinContent(ibin,  hOcc->GetBinContent(ibin)/Double_t(nloop));

  return nloop;
}

Int_t lowAna::readGainConstants(TString fileName)
{
  TString filename("pmtGoodGains_07-31-1555_0.root");
  TFile *fgain = new TFile(filename,"READONLY");
  if(fgain->IsZombie()) { printf(" no gain file %s found \n",filename.Data()); return 0;}
  TTree *gtree=NULL;
  fgain->GetObject("gainsTree",gtree);
  if(!gtree) { printf(" no gainsTree not found \n"); return 0;}
  goodGains = new TPmtGains();
  gtree->SetBranchAddress("pmtGains",&goodGains);
  gtree->GetEntry(0);
  printf(" \n \t using gains from file %s\n",filename.Data()); 
  goodGains->print();
  for(int ipmt=1; ipmt<NPMT; ++ipmt) gain[ipmt]=goodGains->gain[ipmt]/goodGains->gain[0];
  gain[0]=1.0;
  printf(" using normalized gains \n");
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f ; ",ipmt,gain[ipmt]); 
  printf(" \n");

  return 21;
}

std::vector<Int_t> lowAna::findRFTimes(int ipmt, double& step) 
{
  std::vector<Int_t> rftimes;
  int ib; int ic;
  fromPmtNumber(ipmt,ib,ic);

  // find baseline
  std::vector<UShort_t> udigi; 
  for (UInt_t is=0; is<MAXSAMPLES; ++is) udigi.push_back(digitizer_waveforms[ib][ic][is]);
  std::sort(udigi.begin(), udigi.end());
  UShort_t baseline = udigi[0.5*double(MAXSAMPLES)];

  // looking for negative values.  
  UShort_t digiMin=MAXADC;
  for (UInt_t is=0; is<MAXSAMPLES; ++is) {
    digitizer_waveforms[ib][ic][is]=TMath::Min( baseline , digitizer_waveforms[ib][ic][is]);
    if(digitizer_waveforms[ib][ic][is]<digiMin) digiMin=digitizer_waveforms[ib][ic][is];
  }

  step = double(digiMin) - double(baseline);
  // return if step down is too small
  if(step>-500) return rftimes;
  // pick off start of rising edge
  bool isRF=false;
  for (UInt_t is=0; is<MAXSAMPLES; ++is){
    double digi = double(digitizer_waveforms[ib][ic][is]) - double(baseline);
    //histoDraw[iB][iC]->Fill(iS+0.5, ((1.*waveforms[iB][iC][iS]-baseline)*offsetstepADC/(1.*ADCrange+1.)+offset) );
    //int ADCrange = 4095;
    //offsetstepADC = 50.
    //double digi7 =  digi*50./4096.;
    hSamples[ipmt]->SetBinContent(int(is+1),hSamples[ipmt]->GetBinContent(int(is+1))+digi);
    if(digi<0.75*step&&!isRF) {
      rftimes.push_back(is);
      isRF=true;
    } else if(digi>0.75*step) {
      isRF=false;
    }
  }
  return rftimes;
}

std::vector<Int_t> lowAna::findMaxPeak(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of max digi peak
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());


  Double_t vmax=0;
  Int_t    imax=0;
  for(Int_t  ibin=0; ibin<= vsize; ++ibin ) {
    if( v[ibin]>vmax){
      vmax=v[ibin];
      imax=ibin;
    }
  }

  if( vmax<threshold) return peakTime;

  // consider this a "seed" and find full hit
  klow=imax;
  for(Int_t k=imax-1; k>=max(0,imax-maxHalfLength); --k) {
    if(v[k]<=sthreshold) break;
    klow=k;
  }
  khigh=imax;
  for(Int_t k=imax+1; k<min(imax+maxHalfLength,vsize); ++k) {
    if(v[k]<=sthreshold) break;
    khigh=k;
  }
  kover = khigh-klow+1;
  // found good pulse
  if(kover>minLength) { 
    for(Int_t k=klow ; k<= khigh; ++k) peakTime.push_back(k);
    //printf(" peakTime %i, %i ?  %i size %i .... \n ",klow,khigh,kover,peakTime.size());
    //for(UInt_t ih=0; ih<peakTime.size(); ++ih) printf(" \t  %i t= %i \n",ih,peakTime[ih]);
  }
  // skip to end of sthreshold search 

  return peakTime;
}


std::vector<Int_t> lowAna::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());

  //printf(" findPeaks \n");
  for( Int_t ibin=0; ibin< vsize; ++ibin ) {
    if( v[ibin]>threshold) {// starting possible new hit
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=max(0,ibin-maxHalfLength); --k) {
        if(v[k]<=sthreshold) break;
        klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<min(ibin+maxHalfLength,vsize); ++k) {
        if(v[k]<=sthreshold) break;
        khigh=k;
      }
      kover = khigh-klow+1;
      // found good pulse
      if(kover>minLength) {
        double qsum=0;
        for(Int_t k=klow ; k<= khigh; ++k) {
          peakTime.push_back(k);
          qsum +=v[k];
        }
        //printf(" peakTime  qsum %f sthreshod %f ibin %i klow %i khigh %i ", qsum,sthreshold,ibin,klow,khigh);
        //for(Int_t k=klow ; k<= khigh; ++k) printf(" v[%i]=%f ",k,v[k]);
        //printf("\n");
      }
      // skip to end of sthreshold search 
      ibin=khigh;
    }
  }

  //for(UInt_t ih=0; ih<peakTime.size(); ++ih) printf("  %i t= %i ADC= %f\n",ih,peakTime[ih],v[peakTime[ih]]);
  return peakTime;
}

Int_t lowAna::findHits(Int_t ipmt, Double_t sum, std::vector<Int_t> peakTime, std::vector<Double_t> ddigi, std::vector<Double_t> ddigiUn, Int_t type) 
{
  //printf(" findHits called with  peakTime size %i  \n",peakTime.size());

  if(peakTime.size()<1) {
    return 0;
  }
  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(UInt_t it=nlast; it>0; --it) {
    //printf(" .... %i %i hitTime.size %u hitList.size %u \n",it,peakTime[it],hitTime.size(),hitList.size());
    bool makeHit=false;
    if(peakTime[it]-peakTime[it-1]!=1||(it==1&&hitTime.size()>=minLength)) makeHit=true;

    if(makeHit) {
      hitTime.push_back(peakTime[it]);
      hitList.push_back(hitTime);
      //printf(" saving list %i size %i \n",hitList.size(),hitTime.size());
      if(hitTime.size()<minLength) printf(" WARNING:: saving list %zu size %zu \n",hitList.size(),hitTime.size());
      //for(UInt_t ih=0; ih<hitTime.size(); ++ih) printf(" \t\t %i t= %i digi %f \n",ih,hitTime[ih],ddigi[hitTime[ih]]);
      hitTime.clear();
      continue;
    }
    hitTime.push_back(peakTime[it]);
    //printf(" building list %i size %i \n",hitList.size(),hitTime.size());
  }

  //printf(" list of hits  %lu \n",hitList.size());
  Int_t nhits=0;
  for(UInt_t il=0; il<hitList.size(); ++il) {
    TPmtHit phit;
    hitTime=hitList[il];
    //printf(" il= %i hitTime.Size %lu \n ",il,hitTime.size());
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    Double_t qUnhit=0;
    Double_t qUnpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      //printf("  ih = %i time  %i ADC %f XXX ",ih,hitTime[ih],ddigi[hitTime[ih]]);
      phit.tsample.push_back(hitTime[ih]);	
      phit.qsample.push_back(ddigi[hitTime[ih]]);	
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
        qUnpeak = ddigiUn[hitTime[ih]];
      }
      qhit+=ddigi[hitTime[ih]];
      qUnhit+=ddigiUn[hitTime[ih]];
    }

    // fwhm
    phit.fwhm=0; 
    //printf("\n qhit %f qpeak %f samples size %i \n",qhit,qpeak,phit.qsample.size());
    for(UInt_t ih=0; ih<phit.qsample.size(); ++ih) {
      //printf(" %i %f fwhm %f \n ",ih,phit.qsample[ih],phit.fwhm);
      if( phit.qsample[ih] > qpeak/2.0 ) ++phit.fwhm;
    }

    //printf(" \t hit %i start %i stop %i  length %i qhit %f  \n",il, hitTime[hitTime.size()-1],hitTime[0],hitTime.size(),qhit);
    // fill in hit
    phit.ipmt=ipmt;
    phit.timeToRF=0;
    phit.tstart=hitTime[hitTime.size()-1];
    phit.tstop=hitTime[0];
    phit.qhit=qhit;
    phit.qUnhit=qUnhit;
    phit.qpeak=qpeak;
    phit.qUnpeak=qUnpeak;
    phit.ratio=phit.qpeak/phit.qhit;
    phit.peakTime=Double_t(peakt);
    phit.offset=0;
    phit.nsamples=phit.qsample.size();
    //
    Double_t length = TMath::Abs(phit.tstop-phit.tstart)+1;
    // time past latest RF pulse
    if(pmtEvent->trigType==TPmtEvent::TRIG000)  hRawQ[ipmt]->Fill(qUnpeak); 
    // hRawQ[ipmt]->Fill(qUnhit);
    if(qUnpeak<1) phit.print();
    hHitQ[ipmt]->Fill(qhit);
    bool bad=false;
    for(UInt_t ih=0; ih<phit.qsample.size(); ++ih) if(phit.qsample[ih]<=0) bad=true;
    if(bad) {
      printf(" \n\t !!pmtHit::fidHits WARNING l.e. zero sample in hit  \n");
      phit.print();
    }
    //phit.print();
    /* error if qhit < 0! */
    if(phit.qhit<0) { 
      printf("\n\t !!lowAna::findHits WARNING negative charge qhit!! qhit %f \n",phit.qhit) ; phit.print(); }
    if(phit.qUnhit<0) { 
      printf("\n\t !!lowAna::findHits WARNING negative charge qUnhit!! qUnhit %f \n",phit.qUnhit) ; phit.print(); }

    pmtEvent->hit.push_back(phit);
    ++nhits;
  }
  //printf(" findHits found %i \n",nhits);
  return  nhits;
}


TH1D* lowAna::FFTFilter(Int_t ipmt)
{
  int ib,ic;
  fromPmtNumber(ipmt,ib,ic);
  printf(" called FFTFilter pmt %i board %i channel %i \n",ipmt,ib,ic);
  for(int is =0; is<nFFTSize; ++is) {
    fFFT->SetPoint(is, digitizer_waveforms[ib][ic][is]);
  }

  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
  }

  fFFT->Transform();
  TH1D* hfft = new TH1D(Form("FFTPmt%i",ipmt),Form("FFT PMT %i",ipmt),nFFTSize/2,0,nFFTSize/2);

  // fill samples FFT histogram && elec response in time domain
  // printf(" created %s %s \n",hfft->GetName(),hfft->GetTitle());
  // skip first bin which is pedestal
  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im);
    hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
  } 
  return hfft;      
}
Int_t lowAna::triggerInfo()
{
  Int_t type = TPmtEvent::TRIGUNKNOWN; // unknosn 
  // RF channels 
  double s1,s2,s3;
  rftime21 = findRFTimes(21,s1);
  rftime22 = findRFTimes(22,s2);
  rftime23 = findRFTimes(23,s3);

  // calculate average RF time
  Double_t rft = 0;
  // average over 3 boards 
  Int_t nrftimes=0;

  if(rftime21.size()>0) {rft += double(rftime21[0]); ++nrftimes; }
  if(rftime22.size()>0) {rft += double(rftime22[0]); ++nrftimes; }
  if(rftime23.size()>0) {rft += double(rftime23[0]); ++nrftimes; }
  if( nrftimes>0) rft /= double(nrftimes);
  pmtEvent->tRFave=rft;
  for(int itime=0; itime<NB; ++itime)  pmtEvent->dtime[itime] = digitizer_time[itime];


  //UInt_t totalTimes = rftime21.size()+rftime21.size()+rftime21.size();
  double t1 = 0; if(rftime21.size()>0) t1 = rftime21[0];
  double t2 = 0; if(rftime22.size()>0) t2 = rftime22[0];
  double t3 = 0; if(rftime23.size()>0) t3 = rftime23[0];

  pmtSummary->vrf1.push_back(t1);
  pmtSummary->vrf2.push_back(t2);
  pmtSummary->vrf3.push_back(t3);

  Int_t r1 = Int_t(rftime21.size());
  Int_t r2 = Int_t(rftime22.size());
  Int_t r3 = Int_t(rftime23.size());

  // determine trigger type and count types 

  if(r1==0&&r2==0&r3==0) { // zero 
    type = TPmtEvent::TRIG000;
    ++pmtSummary->ntrig000;
  } else if ( r1==1 && r2==1 && r3==1 ) {
    type = TPmtEvent::TRIG111; 
    ++pmtSummary->ntrig111;
  } else if ( r1==1 || r2==1 || r3==1 ) {
    type = TPmtEvent::TRIG1XX; 
    ++pmtSummary->ntrig1xx;
  } else if (r1==0||r2==0||r3==0) { 
    type = TPmtEvent::TRIG0XX;
    ++pmtSummary->ntrig5xx; 
  } else if (r1==5&&r2==5&&r3==5) { //five 
    type = TPmtEvent::TRIG555; 
    ++pmtSummary->ntrig555;
  } else if ( r1==5 || r2==5 || r3==5 ) {
    type = TPmtEvent::TRIG5XX; 
    ++pmtSummary->ntrig5xx;
  } else if (r1==4&&r2==4&&r3==4) {  //four
    type = TPmtEvent::TRIG444; 
    ++pmtSummary->ntrig444;
  } else if ( r1==4 || r2==4 || r3==4 ) {
    type = TPmtEvent::TRIG4XX; 
    ++pmtSummary->ntrig4xx;
  }
  //if(type!= TPmtEvent::TRIG000 && type!=TPmtEvent::TRIG111)
      printf(" %s %i %i %i %lld %i %i %i %u %u %u  \n", pmtEvent->tag.c_str(),pmtEvent->run,pmtEvent->event,pmtEvent->compSec,pmtEvent->compNano,
          int(t1),int(t2),int(t3),pmtEvent->dtime[0],pmtEvent->dtime[1],pmtEvent->dtime[2]);
  ntTrig->Fill(float(pmtEvent->run),float(pmtEvent->event),float(r1),float(r2),float(r3),float(t1),float(t2),float(t3) );
  return type;
}

void lowAna::ADCFilter(int iB, int iC) 
{
  for (int is = 0; is<MAXSAMPLES; ++is) {
    if (digitizer_waveforms[iB][iC][is] > MAXADC) {
      if (is > 0) { digitizer_waveforms[iB][iC][is] = digitizer_waveforms[iB][iC][is-1];}
      else {
        int is2 = 0;
        while (digitizer_waveforms[iB][iC][is2] > MAXADC) {
          digitizer_waveforms[iB][iC][0] = digitizer_waveforms[iB][iC][is2+1];
          ++is2;
        }
      }
    }
  }
}


/******************************** auto generated stuff below. ****************************/
Int_t lowAna::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t lowAna::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void lowAna::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
  fChain->SetBranchAddress("computer_secIntoEpoch", &computer_secIntoEpoch, &b_computer_secIntoEpoch);
  fChain->SetBranchAddress("computer_nsIntoSec", &computer_nsIntoSec, &b_computer_nsIntoSec);
  fChain->SetBranchAddress("gps_nsIntoSec", &gps_nsIntoSec, &b_gps_nsIntoSec);
  fChain->SetBranchAddress("gps_secIntoDay", &gps_secIntoDay, &b_gps_secIntoDay);
  fChain->SetBranchAddress("gps_daysIntoYear", &gps_daysIntoYear, &b_gps_daysIntoYear);
  fChain->SetBranchAddress("gps_Year", &gps_Year, &b_gps_Year);
  fChain->SetBranchAddress("gps_ctrlFlag", &gps_ctrlFlag, &b_gps_ctrlFlag);
  fChain->SetBranchAddress("digitizer_size", digitizer_size, &b_digitizer_size);
  fChain->SetBranchAddress("digitizer_chMask", digitizer_chMask, &b_digitizer_chMask);
  fChain->SetBranchAddress("digitizer_evNum", digitizer_evNum, &b_digitizer_evNum);
  fChain->SetBranchAddress("digitizer_time", digitizer_time, &b_digitizer_time);
  fChain->SetBranchAddress("digitizer_waveforms", digitizer_waveforms, &b_digitizer_waveforms);
  fChain->SetBranchAddress("nDigitizers", &nDigitizers, &b_nDigitizers);
  fChain->SetBranchAddress("nChannels", &nChannels, &b_nChannels);
  fChain->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
  fChain->SetBranchAddress("nData", &nData, &b_nData);
  Notify();
}

Bool_t lowAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void lowAna::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t lowAna::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

Double_t lowAna::getPromptTime()
{
  // fill histogram to find peak bin in event.
  hTPromptEvent->Reset();
    
  std::vector<Int_t> rft;
  for(unsigned ihit=0; ihit< pmtEvent->hit.size(); ++ihit) {//ysun
    if(pmtEvent->hit[ihit].ipmt<7) rft = pmtEvent->rft21;//ysun
    else if(pmtEvent->hit[ihit].ipmt>=7 && pmtEvent->hit[ihit].ipmt<14) rft = pmtEvent->rft22;//ysun
    else if(pmtEvent->hit[ihit].ipmt>=14 && pmtEvent->hit[ihit].ipmt<21) rft = pmtEvent->rft23;//ysun
    if(rft.size()>0){
      for (int i=0;i<pmtEvent->hit[ihit].nsamples;i++) {//ysun
        hTPromptEvent->Fill(pmtEvent->hit[ihit].tsample[i]-rft[0],pmtEvent->hit[ihit].qsample[i]);//ysun
      }//ysun
    }
  }//ysun
  //return Double_t(hTPromptEvent->GetMaximumBin())-pmtEvent->tRFave; //ysun
  if(hTPromptEvent->GetEntries()>0) return Double_t(hTPromptEvent->GetMaximumBin()-MAXSAMPLES); //ysun
  else return -9999;//ysun
}

// nearest RF time to hit peak time
void lowAna::getTimeToRF(UInt_t board) 
{
  std::vector<Int_t> rft;
  if(board==0) rft = pmtEvent->rft21;//ysun
  else if(board==1) rft = pmtEvent->rft22;//ysun
  else if(board==2) rft = pmtEvent->rft23;//ysun
  for(unsigned ihit=0; ihit< pmtEvent->hit.size(); ++ihit) {
    Int_t time = MAXSAMPLES;
    Int_t hitTime = pmtEvent->hit[ihit].peakTime;
    for(unsigned i=0; i<rft.size() ; ++i) {//ysun
      Int_t tdiff = hitTime - rft[i];//ysun
      if(tdiff<0) break;
      if( tdiff<time ) time=tdiff;
    }
    pmtEvent->hit[ihit].timeToRF = time;
  }
}


// neils filter
std::vector<Double_t> lowAna::MovingAverageFilter(std::vector<Double_t> signal,Int_t aveN)
{
  Int_t N = aveN;
  if(aveN%2==0) ++N; 
  std::vector<Double_t> filter;
  Int_t N2 = std::floor(N/2);
  for(int i = N2; i < Int_t(signal.size())-N2; i++){
    Double_t sum = 0;
    for(int j = i-N2; j <= i+N2; j++){
      sum += signal[j];
    }
    sum /= N;
    filter.push_back(sum);
  }

  for(int i = 0; i < N2 ; i++){
    std::vector<Double_t>::iterator it;
    it = filter.begin();
    filter.insert(it,0.);
    filter.push_back(0);
  }
  return filter;
}

// summarize run quality
void lowAna::qualitySummary()
{

  /*
  Double_t x[NPMT], y[NPMT], z[NPMT],y2[NPMT],z2[NPMT],ex[NPMT], ey[NPMT], ez[NPMT];
  Double_t norm[NPMT]; 
  Double_t yun[NPMT], zun[NPMT],yun2[NPMT],zun2[NPMT], eyun[NPMT], ezun[NPMT];
  Double_t normun[NPMT]; 
  Double_t qRF[NPMT];

  for(Int_t j=0; j<NPMT; ++j) {
    x[j]=Double_t(j); ex[j]=0;  
    y[j]=0; z[j]=0; y2[j]=0; z2[j]=0; ey[j]=0; ez[j]=0;
    norm[j]=0;
    qRF[j]=0;
  }


  for (Int_t k=0 ;k<pmtEntries;k++){
    ntPmt->GetEntry(k);
    int ipmt = int(fpmt);
    // cosmic cut
    bool cut = tmax>440&&tmax<480&&sum<5000;
    qRF[ipmt]=Double_t(qrf);
    if(cut) { 
      y[ipmt]+=qmax;
      y2[ipmt]+=pow(qmax,2.);
      z[ipmt]+=sum;
      z2[ipmt]+=pow(sum,2.);
      norm[ipmt]+=1.0;
    }
    if(cut) { 
      yun[ipmt]+=qmaxUn;
      yun2[ipmt]+=pow(qmaxUn,2.);
      zun[ipmt]+=sumUn;
      zun2[ipmt]+=pow(sumUn,2.);
      normun[ipmt]+=1.0;
    }
  }
  for(Int_t j=0; j<NPMT; ++j) {
    y[j]/= norm[j]; z[j]/=norm[j]; y2[j]/=norm[j]; z2[j]/=norm[j];
    yun[j]/= normun[j]; zun[j]/=normun[j]; yun2[j]/=normun[j]; zun2[j]/=normun[j];
  }

  for(Int_t j=0; j<NPMT; ++j) {
    ey[j]= sqrt( (y2[j]-pow(y[j],2.))/norm[j]);
    ez[j]= sqrt( (z2[j]-pow(z[j],2.))/norm[j]);
    eyun[j]= sqrt( (yun2[j]-pow(yun[j],2.))/normun[j]);
    ezun[j]= sqrt( (zun2[j]-pow(zun[j],2.))/normun[j]);
  }

  // store averages in summary
  for(Int_t j=0; j<NPMT; ++j) {
    pmtSummary->norm[j]=norm[j];
    pmtSummary->qmax[j]=y[j];
    pmtSummary->eqmax[j]=ey[j];
    pmtSummary->qsum[j]=z[j];
    pmtSummary->eqsum[j]=ez[j];
    pmtSummary->qrf[j]=qRF[j];
  }
  */
 
  // fixed number
}



