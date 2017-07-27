#include "pmtAna.hh"

pmtAna::pmtAna(TString tag, Int_t maxLoop)
{
  fChain=NULL;
  TString fileName = TString("pdsData/PDSout_") + TString(tag) + TString(".root");
  printf(" looking for file %s\n",fileName.Data());
  TFile *f = new TFile(fileName,"readonly");
  if(!f) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  TTree *tree;
  f->GetObject("pmt_tree",tree);
  tree->ls();
  Init(tree);
  if(!fChain) return;

  
  // initicalize fft 
  nFFTSize = int(MAXSAMPLES);
  fFFT = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");

  // open ouput file and make some histograms
  TString outputFileName = TString("pdsOutput/pmtAna_")+tag+ TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outputFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);


  //ntuples
  ntDigi = new TNtuple("ntDigi"," digi  ","ipmt:idigi:digi"); 
  ntPmt = new TNtuple("ntPmt"," pmts ","ipmt:tmax:qmax:sum:noise:base:nhit");
  ntHit = new TNtuple("ntHit", " hits ","ipmt:sum:time:length:qpeak:qhit:fwhm:ratio");

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

  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("Samples_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Samples board%u channel%u pmt%u",ib,ic,ipmt);
      hSamples[ipmt] = new TH1D(hname,htitle,NS,0,NS);
      hSamples[ipmt]->SetXTitle(" sample number ");

      hname.Form("Peaks_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Peaks board%u channel%u pmt%u",ib,ic,ipmt);
      hPeaks[ipmt] = new TH1D(hname,htitle,NS,0,NS);
      hPeaks[ipmt]->SetXTitle(" sample number ");
      hPeaks[ipmt]->SetLineColor(kRed);
      hPeaks[ipmt]->SetMarkerColor(kRed);
      hPeaks[ipmt]->SetFillColor(kRed);
      hPeaks[ipmt]->SetFillStyle(3002);
      
      hname.Form("Counts_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" counts board%u channel%u pmt %u",ib,ic,ipmt);
      hCounts[ipmt] = new TH1D(hname,htitle,500,0,5000);
      hCounts[ipmt]->SetXTitle(" baseline subtracted summed ADC counts ");

      hname.Form("Baseline_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" Baseline board%u channel%u pmt %u",ib,ic,ipmt);
      hBaseline[ipmt] = new TH1D(hname,htitle,100,-50,50);
      hBaseline[ipmt]->SetXTitle(" baseline fluctuation (ADC counts) ");

      hname.Form("hitQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" hit charge board%u channel%u pmt %u",ib,ic,ipmt);
      hHitQ[ipmt] = new TH1D(hname,htitle,250,0,2500);
      hHitQ[ipmt]->SetXTitle(" hit charge in pulse time  (ADC counts) ");

      
      hname.Form("QMax_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" max ADC c board%u channel%u pmt %u",ib,ic,ipmt);
      hQMax[ipmt] = new TH1D(hname,htitle,100,0,500);
      hQMax[ipmt]->SetXTitle(" q max (ADC counts) ");

      hname.Form("NHits_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" Hits per event c board%u channel%u pmt %u",ib,ic,ipmt);
      hNHits[ipmt] = new TH1D(hname,htitle,20,0,20);
      hNHits[ipmt]->SetXTitle(" number of hits per event ");
      


    }
  }
  //gDirectory->ls();
  
  // loop over entries zero = all 
  UInt_t nLoop = Loop(maxLoop);

  outfile->Write();
  // do some plotting


  if(0) {
    TString canname;
    enum {NCAN=7};
    TCanvas *can1[NCAN];
    TCanvas *can2[NCAN];
    TCanvas *can3[NCAN];
    TCanvas *can4[NCAN];
    TCanvas *can5[NCAN];

    int ican=-1;
    int ip=0;
    for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
      if(ipmt%3==0) {
        ip=0;
        ++ican;
        canname.Form("FFT-set%i-run-%s",ican,tag.Data());
        can1[ican] = new TCanvas(canname,canname);
        can1[ican]->Divide(1,3);
        canname.Form("counts-set%i-run-%s",ican,tag.Data());
        can2[ican] = new TCanvas(canname,canname);
        can2[ican]->Divide(1,3);
        canname.Form("samples-set%i-run-%s",ican,tag.Data());
        can3[ican] = new TCanvas(canname,canname);
        can3[ican]->Divide(1,3);
        canname.Form("hitCharge-set%i-run-%s",ican,tag.Data());
        can4[ican] = new TCanvas(canname,canname);
        can4[ican]->Divide(1,3);
        canname.Form("qMax-set%i-run-%s",ican,tag.Data());
        can5[ican] = new TCanvas(canname,canname);
        can5[ican]->Divide(1,3);
      }
      can1[ican]->cd(ip+1); hFFT[ipmt]->Draw();
      can4[ican]->cd(ip+1); gPad->SetLogy(); hHitQ[ipmt]->Draw();
      can5[ican]->cd(ip+1); gPad->SetLogy(); hQMax[ipmt]->Draw();
      can3[ican]->cd(ip+1); 
      hPeaks[ipmt]->Draw();
      hSamples[ipmt]->Draw("sames");
      can2[ican]->cd(ip+1);  gPad->SetLogy(); hCounts[ipmt]->Draw();
      ++ip;
    }

    for(int ican=0; ican<NCAN; ++ican) {
      can1[ican]->Print(".pdf");
      can2[ican]->Print(".pdf");
      can3[ican]->Print(".pdf");
      can4[ican]->Print(".pdf");
      can5[ican]->Print(".pdf");
    }
  }

}


pmtAna::~pmtAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

UInt_t pmtAna::Loop(UInt_t nToLoop)
{
   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntriesFast();

   UInt_t nbytes=0;
   std::vector<Double_t> sdigi;  // source
   std::vector<Double_t> ddigi;  // baseline subtracted
   UInt_t nloop=nentries;
   if(nToLoop!=0) nloop = nToLoop;
   printf(" entries %lld looping %d \n",nentries,nloop);
  // loop over entries
   for (Long64_t jentry=0; jentry<nloop;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(jentry%100==0) printf(" \t.... %lld \n",jentry);
      nbytes += fChain->GetEntry(jentry);
      // clear the event
      pmtEvent->clear();

      // save event info 
      //pmtEvent.run;
      pmtEvent->event=event_number;
      //pmtEvent.tpcTrig;
      //pmtEvent.pdsTrig;
      pmtEvent->gpsYear=gps_Year;
      pmtEvent->gpsDay=gps_daysIntoYear;;
      pmtEvent->gpsSec=gps_secIntoDay;
      pmtEvent->gpsNs=gps_nsIntoSec;;

       // check if an RF trigger
       /*
        if (digitizer_waveforms[1][15][2000] < 8000) {
            double fulltime = computer_secIntoEpoch+computer_nsIntoSec/1.0E9;
            printf(" this is a RF trigger %20.9f\n", fulltime);
        }
        */
          
      for(UInt_t ib=0; ib<NB; ++ib) {
        UInt_t time = digitizer_time[ib];
        //printf(" board %u time %u \n",ib,time);
        for(UInt_t ic=0; ic<NC; ++ic) {
          // get pmt number
          int ipmt = toPmtNumber(ib,ic);
          if(ipmt<0||ipmt>=NPMT) continue;
          // make a vector of samples for sorting.
          sdigi.clear();
          ddigi.clear();
          double sum=0;

          // Find the sample median and it's "sigma".
          for (UInt_t is=0; is<NS; ++is) sdigi.push_back(double(digitizer_waveforms[ib][ic][is]));

   
          std::sort(sdigi.begin(), sdigi.end());
          double baselineMedian = sdigi[0.5*double(NS)];
          double baselineSigma = sdigi[0.16*double(NS)];
          baselineSigma = std::abs(baselineSigma-baselineMedian);

          //baselineSigma = std::abs(baselineSigma-baselineMedian);
          //noise = sdigi[0.68*sdigi.size()];
          hBase->SetBinContent(ipmt+1,hBase->GetBinContent(ipmt+1)+baselineMedian);
          hBase->SetBinError(ipmt+1,hBase->GetBinError(ipmt+1)+baselineSigma);
          double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);
          hNoise->SetBinContent(ipmt+1,hNoise->GetBinContent(ipmt+1)+noise);
          if(ientry==0) baselineNominal[ipmt]= baselineMedian;
          else hBaseline[ipmt]->Fill(baselineMedian-baselineNominal[ipmt]);
         
          if(ientry==0) hFFT[ipmt]=FFTFilter(ipmt);
          
          UInt_t tmax=0;
          double qmax=0;
          for(UInt_t is=0 ; is<NS; ++is) {
            double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])-baselineMedian);
            if(digi>qmax) {
              qmax=digi;
              tmax=is+1;
            }
            ddigi.push_back(digi);
            if(jentry%100==0)ntDigi->Fill(double(ipmt),double(is),digi);
            hSamples[ipmt]->SetBinContent(int(is+1),hSamples[ipmt]->GetBinContent(int(is+1))+digi);
            //if(digi>3.0*noise) sum+=digi;
            if(is>300&&is<1000) sum+=digi;
          }
          hCounts[ipmt]->Fill(sum);
          //if(sum>500) hOcc->Fill(ipmt+1,1);

          // peak finding
          std::vector<Int_t> peakTime = findPeaks(ddigi,20.0*noise,3.0*noise);
          //std::vector<Int_t> peakTime = findMaxPeak(ddigi,8.0*noise,3.0*noise);
          Int_t nhits = findHits(ipmt,sum,peakTime,ddigi);
          hOcc->Fill(ipmt+1,nhits);
          hNHits[ipmt]->Fill(nhits);
          //printf(" event %i nhits %i \n", pmtEvent->event, pmtEvent->nhits );
          for (UInt_t ip = 0; ip < peakTime.size(); ip++) {
            Int_t bin = peakTime[ip];
            //printf(" ipmt %i ip %i bin %i v %f \n",ipmt,ip,bin,ddigi[bin]);
            hPeaks[ipmt]->SetBinContent(bin+1, hPeaks[ipmt]->GetBinContent(bin+1)+ddigi[bin]);
          }
          hQMax[ipmt]->Fill(qmax);
          ntPmt->Fill(double(ipmt),tmax,qmax,sum,noise,baselineMedian-baselineNominal[ipmt],nhits);
          pmtEvent->qmax.push_back(qmax);
          pmtEvent->qsum.push_back(sum);
        }
      }
      pmtEvent->nhits= pmtEvent->hit.size();
      if(jentry%100==0) printf(" \t\t nhits = %ld \n",pmtEvent->nhits);

      pmtTree->Fill();
   }

   // normalize 
   for(Int_t ipmt=0; ipmt<NPMT; ++ipmt) {
     hBase->SetBinContent(ipmt+1, hBase->GetBinContent(ipmt+1)/Double_t(nloop));
     hBase->SetBinError(ipmt+1, hBase->GetBinError(ipmt+1)/Double_t(nloop));

     //UInt_t sampleNorm = hSamples[ipmt]->GetEntries();
     for(int ibin=1; ibin<= hSamples[ipmt]->GetNbinsX()+1; ++ibin ){   
       hSamples[ipmt]->SetBinContent(ibin, hSamples[ipmt]->GetBinContent(ibin)/Double_t(nloop));
     }

     UInt_t peakNorm = hPeaks[ipmt]->GetEntries();
     for(int ibin=1; ibin<= hPeaks[ipmt]->GetNbinsX()+1; ++ibin ) {  
       hPeaks[ipmt]->SetBinContent(ibin, hPeaks[ipmt]->GetBinContent(ibin)/Double_t(peakNorm));
     }
     
     for(int ibin=1; ibin<=  hCounts[ipmt]->GetNbinsX()+1; ++ibin ) 
       hCounts[ipmt]->SetBinContent(ibin,  hCounts[ipmt]->GetBinContent(ibin)/Double_t(nloop));
     
   }
   for(int ibin=1; ibin<=  hNoise->GetNbinsX()+1; ++ibin ) hNoise->SetBinContent(ibin,  hNoise->GetBinContent(ibin)/Double_t(nloop));
   for(int ibin=1; ibin<=  hOcc->GetNbinsX()+1; ++ibin ) hOcc->SetBinContent(ibin,  hOcc->GetBinContent(ibin)/Double_t(nloop));
   
 
   printf(" finised looping  %u pmtTree size %llu \n",nloop,pmtTree->GetEntries());
   return nloop;
}

std::vector<Int_t> pmtAna::findMaxPeak(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of max digi peak
  std::vector<Int_t> peakTime;
  Int_t minLength=3;
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
  for(Int_t k=imax-1; k>=0; --k) {
    if(v[k]<sthreshold) break;
    klow=k;
  }
  khigh=imax;
  for(Int_t k=imax+1; k<vsize; ++k) {
    if(v[k]<sthreshold) break;
    khigh=k;
  }
  kover = khigh-klow+1;
  // found good pulse
  if(kover>minLength) { 
    for(UInt_t k=klow ; k<= khigh; ++k) peakTime.push_back(k);
    //printf(" peakTime %i, %i ?  %i size %i .... \n ",klow,khigh,kover,peakTime.size());
    //for(UInt_t ih=0; ih<peakTime.size(); ++ih) printf(" \t  %i t= %i \n",ih,peakTime[ih]);
  }
  // skip to end of sthreshold search 

  return peakTime;
}


std::vector<Int_t> pmtAna::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t minLength=5;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());
  Int_t maxlength=5;

  //printf(" findPeaks \n");
  for(Int_t  ibin=0; ibin<= vsize; ++ibin ) {
    if( v[ibin]>threshold) {// starting possible new hit
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=max(0,ibin-maxlength); --k) {
        if(v[k]<sthreshold) break;
        klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<min(ibin+maxlength,vsize); ++k) {
        if(v[k]<sthreshold) break;
        khigh=k;
      }
      kover = khigh-klow+1;
      // found good pulse
      if(kover>minLength) { 
        for(UInt_t k=klow ; k<= khigh; ++k) peakTime.push_back(k);
        //printf(" peakTime %i, %i ?  %i size %i .... \n ",klow,khigh,kover,peakTime.size());
        //for(UInt_t ih=0; ih<peakTime.size(); ++ih) printf(" \t  %i t= %i \n",ih,peakTime[ih]);
      }
      // skip to end of sthreshold search 
      ibin=khigh;
    }
  }
   
  return peakTime;
}

Int_t pmtAna::findHits(Int_t ipmt, Double_t sum, std::vector<Int_t> peakTime, std::vector<Double_t> ddigi) 
{
  //printf(" findHits called with  peakTime size %i  \n",peakTime.size());

  if(peakTime.size()<1) {
    return 0;
  }
  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(UInt_t it=nlast; it>0; --it) {
    //printf(" .... %i %i %u %u \n",it,peakTime[it],hitTime.size(),hitList.size());
    bool makeHit=false;
    if(peakTime[it]-peakTime[it-1]!=1||it==1) makeHit=true;

    if(makeHit) {
      hitList.push_back(hitTime);
      //printf(" saving list %i size %i \n",hitList.size(),hitTime.size());
      //for(UInt_t ih=0; ih<hitTime.size(); ++ih) printf(" \t\t %i t= %i \n",ih,hitTime[ih]);
      hitTime.clear();
      continue;
    }
    hitTime.push_back(peakTime[it]);
    //printf(" building list %i size %i \n",hitList.size(),hitTime.size());
  }

  //printf(" list of hits  %u \n",hitList.size());
  Int_t nhits=0;
  for(UInt_t il=0; il<hitList.size(); ++il) {
    TPmtHit phit;
    hitTime=hitList[il];
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      phit.qsample.push_back(ddigi[hitTime[ih]]);	
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
      }
      qhit+=ddigi[hitTime[ih]];
    }
    // fwhm
    phit.fwhm=0; 
    for(UInt_t ih=0; ih<phit.qsample.size(); ++ih) if( phit.qsample[ih] > qpeak/2.0 ) ++phit.fwhm;

    // fill in hit
    phit.ipmt=ipmt;
    phit.time=0;
    phit.tstart=hitTime[0];
    phit.tstop=hitTime[hitTime.size()-1];
    phit.qhit=qhit;
    phit.qpeak=qpeak;
    phit.ratio=phit.qpeak/phit.qhit;
    phit.peakTime=peakt;
    phit.offset=0;
    phit.nsamples=phit.qsample.size();
    //
    Double_t length = TMath::Abs(phit.tstop-phit.tstart)+1;
    ntHit->Fill(ipmt,sum,peakt,length,qpeak,qhit,phit.fwhm,phit.ratio);
    hHitQ[ipmt]->Fill(qhit);
    //printf(" \t hit %i start %i end %i  length %i qhit %f  \n",il, hitTime[hitTime.size()-1],hitTime[0],hitTime.size(),qhit);
    //for(UInt_t ih=0; ih<times.size(); ++ih) printf(" \t\t %i t= %i \n",ih,times[ih]);
    //gg/phit.print();
    pmtEvent->hit.push_back(phit);
    ++nhits;
  }
  return  nhits;
}


TH1D* pmtAna::FFTFilter(Int_t ipmt)
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
  printf(" created %s %s \n",hfft->GetName(),hfft->GetTitle());
  // skip first bin which is pedestal
  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im);
    hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
  } 
  return hfft;      
}


/******************************** auto generated stuff below. ****************************/
Int_t pmtAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pmtAna::LoadTree(Long64_t entry)
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

void pmtAna::Init(TTree *tree)
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

Bool_t pmtAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pmtAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pmtAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
