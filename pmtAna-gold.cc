#include "pmtAna.hh"
#include <TF1.h>
#include <TPaveStats.h>
#include <TText.h>

pmtAna::pmtAna(TString tag, Int_t maxLoop, Int_t firstEntry)
{
  if(readGainConstants()==0) {
    printf(" cannot read gain constants file so abort \n");
    return;
  }
 

  fChain=NULL;
  TString fileName = TString("pdsData/PDSout_") + TString(tag) + TString(".root");
  printf(" looking for file %s\n",fileName.Data());
  TFile *f = new TFile(fileName,"readonly");
  if(f->IsZombie()) {
    printf(" couldnt open file %s so abort.\n",fileName.Data());
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

  TString summaryFileName = TString("pdsOutput/pmtSummary_")+tag+ TString(".root");
  summaryFile = new TFile(summaryFileName,"recreate");
  summaryFile->cd();
  printf(" opening summary file %s \n",summaryFileName.Data());
  TTree *summaryTree = new TTree("summaryTree","summaryTree");
  pmtSummary  = new TPmtSummary();
  summaryTree->Branch("pmtSummary",&pmtSummary);
  
  // open ouput file and make some histograms
  TString outputFileName = TString("pdsOutput/pmtAna_")+tag+ TString(".root");
  outFile = new TFile(outputFileName,"recreate");
  outFile->cd();
  printf(" opening output file %s \n",outputFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);



  //pmtTree->ls();
  //
  //ntuples
  ntDigi = new TNtuple("ntDigi"," digi  ","ipmt:idigi:digi"); 
  ntPmt = new TNtuple("ntPmt"," pmts ","trig:ipmt:tmax:qmax:sum:tmaxUn:qmaxUn:sumUn:noise:base:nhit");
  ntHit = new TNtuple("ntHit", " hits ","ipmt:sum:time:rftime:length:qpeak:qUnpeak:qhit:qUnhit:fwhm:ratio");

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
      if(ipmt<0) continue;
      hname.Form("Samples_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Samples board%u channel%u pmt%u",ib,ic,ipmt);
      hSamples[ipmt] = new TH1D(hname,htitle,NS,0,NS);
      hSamples[ipmt]->SetXTitle(" sample number ");
    }
  }
  hSamplesSum = new TH1D("SampleSum"," samples summed over PMTs",NS,0,NS);
  hSamplesSum->SetXTitle(" sample number ");

  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("SamplesPDS_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Samples board%u channel%u pmt%u",ib,ic,ipmt);
      hSamplesPDS[ipmt] = new TH1D(hname,htitle,NS,0,NS);
      hSamplesPDS[ipmt]->SetXTitle(" sample number ");
    }
  }
  hSamplesPDSSum = new TH1D("SampleSumPDS"," samples summed over PMTs",NS,0,NS);
  hSamplesPDSSum->SetXTitle(" sample number ");
  
  
  for(UInt_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      hname.Form("Peaks_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Peaks board%u channel%u pmt%u",ib,ic,ipmt);
      hPeaks[ipmt] = new TH1D(hname,htitle,NS,0,NS);
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

      hname.Form("hitQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" hit charge board%u channel%u pmt %u",ib,ic,ipmt);
      hHitQ[ipmt] = new TH1D(hname,htitle,25,0,25);
      hHitQ[ipmt]->SetXTitle(" hit charge in pulse time  (ADC counts) ");

      
      hname.Form("QMax_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" max ADC c board%u channel%u pmt %u",ib,ic,ipmt);
      hQMax[ipmt] = new TH1D(hname,htitle,50,0,100);
      hQMax[ipmt]->SetXTitle(" q max (ADC counts) ");

      hname.Form("NHits_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form(" Hits per event c board%u channel%u pmt %u",ib,ic,ipmt);
      hNHits[ipmt] = new TH1D(hname,htitle,20,0,20);
      hNHits[ipmt]->SetXTitle(" number of hits per event ");
      
      hname.Form("RawQ_b%u_ch%u_pmt%u",ib,ic,ipmt);
      //htitle.Form(" Raw Charge board%u channel%u pmt%u",ib,ic,ipmt);
      htitle.Form("");
      hQUnPeak[ipmt] = new TH1D(hname,htitle,60,0,60);
      hQUnPeak[ipmt]->SetXTitle(Form(" Raw Charge [b%u c%u pmt%u] (ADC)",ib,ic,ipmt));
      hQUnPeak[ipmt]->SetYTitle(" Entries / [ADC]");
    }
  }
  //gDirectory->ls();
  
  //pmtSummary->SetName(Form("summary-%s",TString(tag(0,10)).Data()));
  //pmtSummary->tag = tag;
  
  // loop over entries zero = all 
  UInt_t nLoop = Loop(maxLoop,firstEntry);
  qualitySummary(tag);
  
  outFile->Write();
  printf(" wrote output file %s \n",outFile->GetName());

  summaryTree->Fill();
  summaryFile->Write();
  printf(" wrote summary file %s \n",summaryFile->GetName());


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

UInt_t pmtAna::Loop(UInt_t nToLoop,UInt_t firstEntry)
{
   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes=0;
   std::vector<Double_t> sdigi;  // source
   std::vector<Double_t> ddigi;  // baseline subtracted

   // no gain applied
   std::vector<Double_t> sdigiUn;  // source
   std::vector<Double_t> ddigiUn;  // baseline subtracted

   UInt_t nloop=nentries;
   if(nToLoop!=0) nloop = nToLoop;
   printf(" entries %lld looping %d first %d \n",nentries,nloop,firstEntry);
  // loop over entries
   for (Long64_t jentry=firstEntry; jentry<nloop+firstEntry; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { printf(" load tree returns %lld\n",ientry); break;}
      nbytes += fChain->GetEntry(jentry);
      if(jentry%100==0) printf(" \t.... %lld nbytes %lld pmtTree entries %lld \n",jentry,nbytes,pmtTree->GetEntries());
      // clear the event
      pmtEvent->clear();
      // trigger type
      pmtEvent->trigType = triggerInfo();

      // save event info 
      //pmtEvent.run;
      pmtEvent->event=event_number;
      //pmtEvent.tpcTrig;
      //pmtEvent.pdsTrig;
      pmtEvent->gpsYear=gps_Year;
      pmtEvent->gpsDay=gps_daysIntoYear;;
      pmtEvent->gpsSec=gps_secIntoDay;
      pmtEvent->gpsNs=gps_nsIntoSec;;

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
          sdigiUn.clear();
          ddigiUn.clear();
          
          double sum=0;
          double sumUn=0;

          // Find the sample median and it's "sigma".
          for (UInt_t is=0; is<NS; ++is) {
            sdigi.push_back(double(digitizer_waveforms[ib][ic][is])/gain[ipmt]);
            sdigiUn.push_back(double(digitizer_waveforms[ib][ic][is]));
          }

          std::sort(sdigi.begin(), sdigi.end());
          double baselineMedian = sdigi[0.5*double(NS)];
          double baselineSigma = sdigi[0.16*double(NS)];
          baselineSigma = std::abs(baselineSigma-baselineMedian);

          std::sort(sdigiUn.begin(), sdigiUn.end());
          double baselineMedianUn = sdigiUn[0.5*double(NS)];
          double baselineSigmaUn = sdigiUn[0.16*double(NS)];
          baselineSigmaUn = std::abs(baselineSigmaUn-baselineMedianUn);
          

          //baselineSigma = std::abs(baselineSigma-baselineMedian);
          //noise = sdigi[0.68*sdigi.size()];/
          hBase->SetBinContent(ipmt+1,hBase->GetBinContent(ipmt+1)+baselineMedian);
          hBase->SetBinError(ipmt+1,hBase->GetBinError(ipmt+1)+baselineSigma);
          double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);
          hNoise->SetBinContent(ipmt+1,hNoise->GetBinContent(ipmt+1)+noise);
          if(ientry==0) baselineNominal[ipmt]= baselineMedian;
          else hBaseline[ipmt]->Fill(baselineMedian-baselineNominal[ipmt]);
         
          if(ientry==0) hFFT[ipmt]=FFTFilter(ipmt);
          
          UInt_t tmax=0;
          double qmax=0;
          UInt_t tmaxUn=0;
          double qmaxUn=0;
          for(UInt_t is=0 ; is<NS; ++is) {
            double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])/gain[ipmt]-baselineMedian);
            if(digi>qmax) {
              qmax=digi;
              tmax=is+1;
            }
            // witout gain
            double digiUn = -1.0*(double(digitizer_waveforms[ib][ic][is])-baselineMedianUn);
            if(digiUn>qmaxUn) {
              qmaxUn=digiUn;
              tmaxUn=is+1;
            }
            
            ddigi.push_back(digi);
	    ddigiUn.push_back(digiUn);
            if(jentry%100==0)ntDigi->Fill(double(ipmt),double(is),digi);
            if(pmtEvent->trigType == TPmtEvent::TRIG000) hSamplesPDS[ipmt]->SetBinContent(int(is+1),hSamplesPDS[ipmt]->GetBinContent(int(is+1))+digi);
            else hSamples[ipmt]->SetBinContent(int(is+1),hSamples[ipmt]->GetBinContent(int(is+1))+digi);
            // here I am not worrying about the difference between noise and gain-corrected noise.  just using gain-corrected noise
            if(digi>3.0*noise) { 
              //if(is>450&&is<470) 
              sum+=digi;
              sumUn+=digiUn;
            }
          } // loop over digitizations

          //if(sum>400&&pmtEvent->trigType!=7) printf(" SSSSSSSSS event %lli trig %i sum %f   \n", jentry, pmtEvent->trigType, sum );
          hCounts[ipmt]->Fill(sum);
          // peak finding
          std::vector<Int_t> peakTime = findPeaks(ddigi,4.0*noise,1.0*noise);
          //std::vector<Int_t> peakTime = findMaxPeak(ddigi,8.0*noise,3.0*noise);
          Int_t nhits = findHits(ipmt,sum,peakTime,ddigi,ddigiUn,pmtEvent->trigType);
          hOcc->Fill(ipmt+1,nhits);
          hNHits[ipmt]->Fill(nhits);
          //printf(" event %i nhits %i \n", pmtEvent->event, pmtEvent->nhits );
          for (UInt_t ip = 0; ip < peakTime.size(); ip++) {
            Int_t bin = peakTime[ip];
            //printf(" ipmt %i ip %i bin %i v %f \n",ipmt,ip,bin,ddigi[bin]);
            hPeaks[ipmt]->SetBinContent(bin+1, hPeaks[ipmt]->GetBinContent(bin+1)+ddigi[bin]);
          }
          hQMax[ipmt]->Fill(qmax);
          ntPmt->Fill(double(pmtEvent->trigType),double(ipmt),tmax,qmax,sum,tmaxUn,qmaxUn,sumUn,noise,baselineMedian-baselineNominal[ipmt],nhits);
          pmtEvent->qmax.push_back(qmax);
          pmtEvent->qsum.push_back(sum);
        } // channel loop 
      } // board loop 
      pmtEvent->nhits= pmtEvent->hit.size();
      if(jentry%100==0) printf(" \t\t jentry %lli nhits = %d \n",jentry,pmtEvent->nhits);
      pmtTree->Fill();
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
 
Int_t pmtAna::readGainConstants(TString fileName)
{
  ifstream in;
  in.open(fileName);
  Int_t ngains=0;
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return ngains;
  }
  Double_t rgain[NPMT];
  printf(" readGainConstants from file %s \n",fileName.Data());
  string line,label,type,sgain;
  while (in.good()) {
    in >> label >> type >> sgain;
    if(in.eof()) break;
    // look for comment or blank line
    if( label.find("%") != std::string::npos || label.size()<4 ) {
      getline(in,line); // throw away line
      continue;
    }
    if( label.size()<2) continue;
    //cout << label << "  " << type << "  " << sgain << endl;
    int b = atoi(&label[1]);
    int c = atoi(&label[3]);
    int ipmt = toPmtNumber(b,c);
    rgain[ipmt] = atof(sgain.c_str());
    ++ngains;
  }
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f \n",ipmt,rgain[ipmt]); 
  // normalize 
  for(int ipmt=0; ipmt<NPMT; ++ipmt) { gain[ipmt] = rgain[ipmt]/rgain[0];}
  printf(" normalized \n");
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f \n",ipmt,gain[ipmt]); 

  printf(" pmt gain\n");
  in.close();
  return ngains;
}

std::vector<Int_t> pmtAna::findRFTimes(int ipmt, double& step) 
{
  std::vector<Int_t> rftimes;
  int ib; int ic;
  fromPmtNumber(ipmt,ib,ic);

  // find baseline
  std::vector<UShort_t> udigi; 
  for (UInt_t is=0; is<NS; ++is) udigi.push_back(digitizer_waveforms[ib][ic][is]);
  std::sort(udigi.begin(), udigi.end());
  UShort_t baseline = udigi[0.5*double(NS)];
  
  // looking for negative values.  
  UShort_t digiMin=MAXADC;
  for (UInt_t is=0; is<NS; ++is) {
    digitizer_waveforms[ib][ic][is]=TMath::Min( baseline , digitizer_waveforms[ib][ic][is]);
    if(digitizer_waveforms[ib][ic][is]<digiMin) digiMin=digitizer_waveforms[ib][ic][is];
  }
  
  step = double(digiMin) - double(baseline);
  // return if step down is too small
  if(step>-500) return rftimes;
  // pick off start of rising edge
  bool isRF=false;
  for (UInt_t is=0; is<NS; ++is){
    double digi = double(digitizer_waveforms[ib][ic][is]) - double(baseline);
    //printf(" is %i digi %f,%f thresh %f base %f \n",is,double(digitizer_waveforms[ib][ic][is]),digi,0.75*step,double(baseline));
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

std::vector<Int_t> pmtAna::findMaxPeak(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
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
    if(v[k]<sthreshold) break;
    klow=k;
  }
  khigh=imax;
  for(Int_t k=imax+1; k<min(imax+maxHalfLength,vsize); ++k) {
    if(v[k]<sthreshold) break;
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


std::vector<Int_t> pmtAna::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());

  //printf(" findPeaks \n");
  for(Int_t  ibin=0; ibin<= vsize; ++ibin ) {
    if( v[ibin]>threshold) {// starting possible new hit
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=max(0,ibin-maxHalfLength); --k) {
        if(v[k]<sthreshold) break;
        klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<min(ibin+maxHalfLength,vsize); ++k) {
        if(v[k]<sthreshold) break;
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
      ibin=khigh;
    }
  }
   
  return peakTime;
}

Int_t pmtAna::findHits(Int_t ipmt, Double_t sum, std::vector<Int_t> peakTime, std::vector<Double_t> ddigi, std::vector<Double_t> ddigiUn, Int_t type) 
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
      //printf(" saving list %i size %i \n",hitList.size(),hitTime.size());
      hitList.push_back(hitTime);
      if(hitTime.size()<minLength) printf(" WARNING:: saving list %zu size %zu \n",hitList.size(),hitTime.size());
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
   // printf(" %i hitTime.Size %i \n ",il,hitTime.size());
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    Double_t qUnhit=0;
    Double_t qUnpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      //printf(" \t ih = %i time  %i sample %f  ",ih,hitTime[ih],ddigi[hitTime[ih]]);
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
   // printf("\n qhit %f qpeak %f samples size %i \n",qhit,qpeak,phit.qsample.size());
    for(UInt_t ih=0; ih<phit.qsample.size(); ++ih) {
      //printf(" %i %f fwhm %f \n ",ih,phit.qsample[ih],phit.fwhm);
      if( phit.qsample[ih] > qpeak/2.0 ) ++phit.fwhm;
    }

    //printf(" \t hit %i start %i stop %i  length %i qhit %f  \n",il, hitTime[hitTime.size()-1],hitTime[0],hitTime.size(),qhit);
    // fill in hit
    phit.ipmt=ipmt;
    phit.time=0;
    phit.tstart=hitTime[hitTime.size()-1];
    phit.tstop=hitTime[0];
    phit.qhit=qhit;
    phit.qhit=qUnhit;
    phit.qpeak=qpeak;
    phit.qUnpeak=qUnpeak;
    phit.ratio=phit.qpeak/phit.qhit;
    phit.peakTime=peakt;
    phit.offset=0;
    phit.nsamples=phit.qsample.size();
    //
    Double_t length = TMath::Abs(phit.tstop-phit.tstart)+1;
    // time past latest RF pulse
    Int_t rft = -99;
    //if(rfTime22.size()>0) rft=peakt%rfTime22[0];
    if(type==TPmtEvent::TRIG000) hQUnPeak[ipmt]->Fill(qUnpeak);
    ntHit->Fill(ipmt,sum,peakt,rft,length,qpeak,qUnpeak,qhit,qUnhit,phit.fwhm,phit.ratio);
    hHitQ[ipmt]->Fill(qhit);
    //for(UInt_t ih=0; ih<times.size(); ++ih) printf(" \t\t %i t= %i \n",ih,times[ih]);
    //phit.print();
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

// trigger information
Int_t pmtAna::triggerInfo()
{
  Int_t type = TPmtEvent::TRIGUNKNOWN; // unknosn 
  // RF channels 
  double s1,s2,s3;
  rftime21 = findRFTimes(21,s1);
  rftime22 = findRFTimes(22,s2);
  rftime23 = findRFTimes(23,s3);
  //UInt_t totalTimes = rftime21.size()+rftime21.size()+rftime21.size();
  double t1 = 0; if(rftime21.size()>0) t1 = rftime21[0];
  double t2 = 0; if(rftime22.size()>0) t2 = rftime22[0];
  double t3 = 0; if(rftime23.size()>0) t3 = rftime23[0];

  Int_t r1 = Int_t(rftime21.size());
  Int_t r2 = Int_t(rftime22.size());
  Int_t r3 = Int_t(rftime23.size());

  // determine trigger type and count types 
  
  if(r1==0&&r2==0&r3==0) { // zero 
    type = TPmtEvent::TRIG000;
    ++pmtSummary->ntrig000;
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
  return type;
 }

// summarize run quality
void pmtAna::qualitySummary(TString tag)
{
    Int_t pmtEntries = (Int_t) ntPmt->GetEntries();
    if(ntPmt->GetEntries()<1) return;
    //cout<<"Number of quality entries: "<<entries<< endl;
    Float_t ftrig,fpmt,tmax,qmax,sum,tmaxUn,qmaxUn,sumUn,noise,base,nhit;

    ntPmt->SetBranchAddress("trig",&ftrig);
    ntPmt->SetBranchAddress("ipmt",&fpmt);
    ntPmt->SetBranchAddress("tmax",&tmax);
    ntPmt->SetBranchAddress("qmax",&qmax);
    ntPmt->SetBranchAddress("sum",&sum);
    ntPmt->SetBranchAddress("tmaxUn",&tmaxUn);
    ntPmt->SetBranchAddress("qmaxUn",&qmaxUn);
    ntPmt->SetBranchAddress("sumUn",&sumUn);
    ntPmt->SetBranchAddress("noise",&noise);
    ntPmt->SetBranchAddress("base",&base);
    ntPmt->SetBranchAddress("nhit",&nhit);

    Double_t x[NPMT], y[NPMT], z[NPMT],y2[NPMT],z2[NPMT],ex[NPMT], ey[NPMT], ez[NPMT];
    Double_t norm[NPMT]; 
    Double_t yun[NPMT], zun[NPMT],yun2[NPMT],zun2[NPMT], eyun[NPMT], ezun[NPMT];
    Double_t normun[NPMT]; 
    
    for(Int_t j=0; j<NPMT; ++j) {
       x[j]=Double_t(j); ex[j]=0;  
       y[j]=0; z[j]=0; y2[j]=0; z2[j]=0; ey[j]=0; ez[j]=0;
       norm[j]=0;
       yun[j]=0; zun[j]=0; yun2[j]=0; zun2[j]=0; eyun[j]=0; ezun[j]=0;
       normun[j]=0;
    }
  
 
    for (Int_t k=0 ;k<pmtEntries;k++){
      ntPmt->GetEntry(k);
      int ipmt = int(fpmt);
      // cosmic cut
      bool cut = tmax>440&&tmax<480&&sum<5000;
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
    }
    
    gROOT->SetStyle("C43");
    gStyle->SetOptLogy(1);
    gStyle->SetOptStat(0000);
    
    int binmax;
    double histmax;

    TCanvas *myc1 = new TCanvas(Form("pmtRawQ-%s_1",tag.Data()),Form("pmtRawQ-%s_1",tag.Data()),0,0,1860,900);
    myc1->UseCurrentStyle();
    myc1->Divide(4,3);

    TCanvas *myc2 = new TCanvas(Form("pmtRawQ-%s_2",tag.Data()),Form("pmtRawQ-%s_2",tag.Data()),0,0,1860,900);
    myc2->UseCurrentStyle();
    myc2->Divide(3,3);

    //double gain[21],width[21];
    for (int i=0;i<21;i++)
    {
      //cout<<"111"<<endl;
      binmax = hQUnPeak[i]->GetMaximumBin();
      histmax = hQUnPeak[i]->GetBinContent(binmax);
      //cout<<"222"<<endl;
      TF1* f2 = new TF1("f2","[2]*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))+[5]*exp(-(x-[3])*(x-[3])/(2*[4]*[4]))",0,35);
      double par[6]={binmax,binmax/4.,histmax,15,15/2.,histmax/500.};
      f2->SetParLimits(1,0.50,1);
      f2->SetParLimits(3,10,25);
      f2->SetParLimits(4,5,10);
      f2->SetParameters(par);
      f2->SetLineColor(2);
      //cout<<"333"<<endl;
      if(i<12)myc1->cd(i+1);
      else    myc2->cd(i-11);
      hQUnPeak[i]->Draw();
      hQUnPeak[i]->Fit("f2","R");
      //cout<<"444"<<endl;
      TF1* f1 = new TF1("f1","[2]*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))+[5]*exp(-(x-[3])*(x-[3])/(2*[4]*[4]))",0,f2->GetParameter(3)+2*f2->GetParameter(4));
      double par1[6]={f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5)};
      f1->SetParLimits(4,5,10);
      f1->SetParameters(par1);
      f1->SetLineColor(2);
      hQUnPeak[i]->Fit("f1","R");
      hQUnPeak[i]->SetAxisRange(0,60);
      //gain[i]=f1->GetParameter(3);
      //width[i]=f1->GetParameter(4);
      pmtSummary->gain[i]=f1->GetParameter(3);
      pmtSummary->gain_e[i]=f1->GetParError(3);
      //cout<<"555"<<endl;
      TPaveStats *ptstats = new TPaveStats(0.65,0.5, 0.95,0.95,"brNDC");
      ptstats->SetName("stats");
      ptstats->SetBorderSize(1);
      ptstats->SetFillColor(0);
      ptstats->SetTextAlign(12);
	TText *text;// = ptstats->AddText(Form(""));
	text = ptstats->AddText(Form("Noise  = %0.3f",f1->GetParameter(0)));
	text = ptstats->AddText(Form("Noise_{#sigma}  = %0.3f",f1->GetParameter(1)));
	text = ptstats->AddText(Form("A_{noise}  = %0.0f",f1->GetParameter(2)));
	text = ptstats->AddText(Form("Gain  = %0.3f",f1->GetParameter(3)));
	text = ptstats->AddText(Form("#sigma  = %0.3f",f1->GetParameter(4)));
	text = ptstats->AddText(Form("A_{SPE}  = %0.0f",f1->GetParameter(5)));
      //cout<<"666"<<endl;
      ptstats->SetOptStat(0);
      ptstats->SetOptFit(9);
      ptstats->Draw();
      //cout<<"777"<<endl;
      delete f2;
      delete f1;
      //delete ptstats;
      //delete text;
      //cout<<"888"<<endl;
    }
    
    pmtSummary->print();
    myc1->Print(".pdf");
    myc2->Print(".pdf");
    printf(" \n uncorrected PMT averages \n");
    for(Int_t j=0; j<NPMT; ++j) 
      printf(" ipmt %i norm %i qmax %.2f +/- %.2f sum %.2f +/- %.2f gain %.2f +/- %.2f \n",j,int(normun[j]),yun[j],eyun[j],zun[j],ezun[j],pmtSummary->gain[j],pmtSummary->gain_e[j]);

     
    
    TGraphErrors* gr1 = new TGraphErrors(NPMT,x,y,ex,ey);
    TGraphErrors* grUn1 = new TGraphErrors(NPMT,x,yun,ex,eyun);
    TCanvas *c1 = new TCanvas(Form("pmtQMaxAverages-%s",tag.Data()),Form("pmt qmax average %s",tag.Data()));
    gr1->SetMarkerStyle(21);
    grUn1->SetMarkerStyle(22);
    grUn1->SetMarkerColor(kBlue);  grUn1->SetLineColor(kBlue); 
    gr1->SetMarkerColor(kGreen);  gr1->SetLineColor(kGreen); 
    grUn1->SetName(Form("uncorrected pmt qmax average %s",tag.Data()));
    grUn1->SetTitle(Form("uncorrected pmt qmax average %s",tag.Data()));
    gr1->SetName(Form("pmt qmax average %s",tag.Data()));
    gr1->SetTitle(Form("pmt qmax average %s",tag.Data()));
    gr1->GetXaxis()->SetTitle(" pmt number ");
    gr1->GetYaxis()->SetRangeUser(20,50);
    gr1->Draw("ap");
    grUn1->Draw("psame");
    c1->Print(".pdf");
    outFile->Append(gr1);
    outFile->Append(grUn1);
    TGraphErrors* gr2 = new TGraphErrors(NPMT,x,z,ex,ez);
    TGraphErrors* grUn2 = new TGraphErrors(NPMT,x,zun,ex,ezun);
    TCanvas *c2 = new TCanvas(Form("pmtPeakAverages-%s",tag.Data()),Form("pmt peak sum averages %s",tag.Data()));
    gr2->SetMarkerStyle(21);
    grUn2->SetMarkerStyle(22);
    grUn2->SetMarkerColor(kBlue);  grUn2->SetLineColor(kBlue); 
    gr2->SetMarkerColor(kGreen);  gr2->SetLineColor(kGreen); 
    gr2->SetName(Form("pmt peak sum average %s",tag.Data()));
    gr2->SetTitle(Form("pmt peak sum average %s",tag.Data()));
    grUn2->SetName(Form("uncorrected pmt peak sum average %s",tag.Data()));
    grUn2->SetTitle(Form("uncorrected pmt peak sum average %s",tag.Data()));
    gr2->GetXaxis()->SetTitle(" pmt number ");
    gr2->GetYaxis()->SetRangeUser(500,1500);
    gr2->Draw("ap");
    grUn2->Draw("psame");
    c2->Print(".pdf");
    outFile->Append(gr2);
    outFile->Append(grUn2);
}


void pmtAna::ADCFilter(int iB, int iC) 
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
