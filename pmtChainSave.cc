#include "pmtChain.hh"
#include <TF1.h>
#include <TPaveStats.h>
#include <TText.h>
/*
** remove bad events. adapted from pmtAna so has lots of unused, irrelevant code
*/
pmtChain::pmtChain(Int_t maxLoop, Int_t firstEntry)
{
  if(readGainConstants()==0) {
    printf(" cannot read gain constants file so abort \n");
    return;
  }
  TString tag("low_intensity");
  makeChain();  
  TTree *tree=NULL;
  if(!fChain) return;
  Init();
  maxLoop=10000;

  // initicalize fft 
  nFFTSize = int(MAXSAMPLES);
  fFFT = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");

  
  // open ouput file and make some histograms
  TString outputFileName = TString("pdsOutput/pmtChainLow") + TString(".root");
  outFile = new TFile(outputFileName,"recreate");
  promptDir = outFile->mkdir("promptDir");
  outFile->cd();
  printf(" opening output file %s \n",outputFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);


  /***  loop over entries zero = all ***/
  UInt_t nLoop = Loop(maxLoop,firstEntry);
  //qualitySummary(tag);

  outFile->Write();
  printf(" wrote output file %s \n",outFile->GetName());

}


pmtChain::~pmtChain()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

UInt_t pmtChain::Loop(UInt_t nToLoop,UInt_t firstEntry)
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
  // loop over entries
  for (Long64_t jentry=firstEntry; jentry<nloop+firstEntry; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) { printf(" load tree returns %lld\n",ientry); break;}
    nbytes += fChain->GetEntry(jentry);
    if(jentry%1000==0) printf(" \t.... %lld nbytes %lld pmtTree entries %lld \n",jentry,nbytes,pmtTree->GetEntries());
    // clear the event
    pmtEvent->clear();
    pmtEvent->tag=TString(stag.c_str());
    // trigger type
    // pmtEvent->trigType = triggerInfo();

    // save event info 
    //pmtEvent.run;
    pmtEvent->event=jentry;
    //pmtEvent.tpcTrig;
    //pmtEvent.pdsTrig;
    //pmtEvent->rft21=rftime21;
    //pmtEvent->rft22=rftime22;
    //pmtEvent->rft23=rftime23;
    pmtEvent->compSec=computer_secIntoEpoch;
    pmtEvent->compNano=computer_nsIntoSec;
   
    UInt_t rftime[3];
    rftime[0]=0; if(rftime21.size()>0) rftime[0]=UInt_t(rftime21[0]);
    rftime[1]=0; if(rftime22.size()>0) rftime[1]=UInt_t(rftime22[0]);
    rftime[2]=0; if(rftime23.size()>0) rftime[2]=UInt_t(rftime23[0]);

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
        for (UInt_t is=0; is<MAXSAMPLES; ++is) sdigi.push_back(double(digitizer_waveforms[ib][ic][is]));
        std::sort(sdigi.begin(), sdigi.end());
        double baselineMedian = sdigi[0.5*double(MAXSAMPLES)];
        double baselineSigma = sdigi[0.16*double(MAXSAMPLES)];
        baselineSigma = std::abs(baselineSigma-baselineMedian);
        double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);

        //if(ientry==0) hFFT[ipmt]=FFTFilter(ipmt);

        UInt_t tmax=0;
        double qmax=0;
        Double_t qrf=0;
        for(UInt_t is=0 ; is<MAXSAMPLES; ++is) {
          double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])/gain[ipmt]-baselineMedian);
          if(digi>qmax) {
            qmax=digi;
            tmax=is+1;
          } 
          ddigi.push_back(digi);
          if(digi>3.0*noise) { 
            sum+=digi;
          }
        } // loop over digitizations
        pmtEvent->qmax.push_back(qmax);
        pmtEvent->qsum.push_back(sum);
      } // channel loop 
    } // board loop 
    pmtTree->Fill();
  }   // end loop over entries
  printf(" finised looping  %u pmtTree size %llu \n",nloop,pmtTree->GetEntries());
  return nloop;
}

Int_t pmtChain::readGainConstants(TString fileName)
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

bool pmtChain::readAlignmentConstants(TString tag, TString fileName)
{
  bool found=0;
  cout << " reading alignment file " << fileName << " looking for " << tag  << endl;
  TFile*  fAlignIn = new TFile(fileName, "READ");
  if(fAlignIn->IsZombie()) {
    printf(" cannot read file %s\n",fileName.Data());
    return found;
  }
  TTree* atree = (TTree *)fAlignIn->Get("alignTree");
  if(!atree) {
    printf(" cannot find alignTree in file %s\n",fileName.Data());
    fAlignIn->Close();
    return found;
  }
  ULong64_t nEntry = atree->GetEntries();
  printf(" have %d alignments runs \n",int(nEntry));
  TPmtAlign* pmtAlign = new TPmtAlign();
  atree->SetBranchAddress("pmtAlign",&pmtAlign);
  // get alignments for this tag
  //align0.clear();
  //align1.clear();
  //align2.clear();
  
  for(ULong64_t entry=0; entry< nEntry; ++entry){
    atree->GetEntry(entry);
    if(pmtAlign->tag==tag) found=true;
    if(found) { // load constants from file

      Long64_t start0 = Long64_t(pmtAlign->start0);
      Long64_t start1 = Long64_t(pmtAlign->start1);
      Long64_t start2 = Long64_t(pmtAlign->start2); 
      addAlign1 = start0-start1;
      addAlign2 = start0-start2;
      printf(" readAlignmentConstants: %s  %lli  %lli  %lli  %lli  %lli \n", pmtAlign->tag.c_str(),start0,start1,start2,start0-start1,start0-start2);

      /*
      for(unsigned i=0; i < pmtAlign->align0.size(); ++i) {
        align0.push_back(pmtAlign->align0[i]);
        align1.push_back(pmtAlign->align1[i]);
        align2.push_back(pmtAlign->align2[i]);
      }
      */
    }
    if(found) break;
  }
  fAlignIn->Close();
  return found;
}

std::vector<Int_t> pmtChain::findRFTimes(int ipmt, double& step) 
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

std::vector<Int_t> pmtChain::findMaxPeak(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
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


std::vector<Int_t> pmtChain::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
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

Int_t pmtChain::findHits(Int_t ipmt, Double_t sum, std::vector<Int_t> peakTime, std::vector<Double_t> ddigi, std::vector<Double_t> ddigiUn, Int_t type) 
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
    ntHit->Fill(ipmt,sum,peakt,pmtEvent->tRFave,length,qpeak,qUnpeak,qhit,qUnhit,phit.fwhm,phit.ratio);
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
      printf("\n\t !!pmtChain::findHits WARNING negative charge qhit!! qhit %f \n",phit.qhit) ; phit.print(); }
    if(phit.qUnhit<0) { 
      printf("\n\t !!pmtChain::findHits WARNING negative charge qUnhit!! qUnhit %f \n",phit.qUnhit) ; phit.print(); }

    pmtEvent->hit.push_back(phit);
    ++nhits;
  }
  //printf(" findHits found %i \n",nhits);
  return  nhits;
}


TH1D* pmtChain::FFTFilter(Int_t ipmt)
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
Int_t pmtChain::triggerInfo()
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
  return type;
}

// summarize run quality
void pmtChain::qualitySummary(TString tag)
{
  Int_t pmtEntries = (Int_t) ntPmt->GetEntries();
  if(ntPmt->GetEntries()<1) return;
  //cout<<"Number of quality entries: "<<entries<< endl;
  Float_t ftrig,fpmt,tmax,qmax,sum,tmaxUn,qmaxUn,sumUn,noise,base,nhit,qrf;

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
  ntPmt->SetBranchAddress("qrf",&qrf);

  Double_t x[NPMT], y[NPMT], z[NPMT],y2[NPMT],z2[NPMT],ex[NPMT], ey[NPMT], ez[NPMT];
  Double_t norm[NPMT]; 
  Double_t yun[NPMT], zun[NPMT],yun2[NPMT],zun2[NPMT], eyun[NPMT], ezun[NPMT];
  Double_t normun[NPMT]; 
  Double_t qRF[NPMT];

  for(Int_t j=0; j<NPMT; ++j) {
    x[j]=Double_t(j); ex[j]=0;  
    y[j]=0; z[j]=0; y2[j]=0; z2[j]=0; ey[j]=0; ez[j]=0;
    norm[j]=0;
    yun[j]=0; zun[j]=0; yun2[j]=0; zun2[j]=0; eyun[j]=0; ezun[j]=0;
    normun[j]=0;
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
  for (int i=0;i<21;i++) {
    //cout<<"111"<<endl;
    binmax = hRawQ[i]->GetMaximumBin();
    double fbinmax = double(binmax);
    histmax = hRawQ[i]->GetBinContent(fbinmax);
    //cout<<"222"<<endl;
    TF1* f2 = new TF1("f2","[2]*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))+[5]*exp(-(x-[3])*(x-[3])/(2*[4]*[4]))",0,35);
    double par[6]={fbinmax,fbinmax/4.,histmax,15,15/2.,histmax/500.};
    f2->SetParLimits(1,0.50,1);
    f2->SetParLimits(3,10,25);
    f2->SetParLimits(4,5,10);
    f2->SetParameters(par);
    f2->SetLineColor(2);
    //cout<<"333"<<endl;
    if(i<12)myc1->cd(i+1);
    else    myc2->cd(i-11);
    hRawQ[i]->Draw();
    hRawQ[i]->Fit("f2","R");
    //cout<<"444"<<endl;
    TF1* f1 = new TF1("f1","[2]*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))+[5]*exp(-(x-[3])*(x-[3])/(2*[4]*[4]))",0,f2->GetParameter(3)+2*f2->GetParameter(4));
    double par1[6]={f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(2),f2->GetParameter(3),f2->GetParameter(4),f2->GetParameter(5)};
    f1->SetParLimits(4,5,10);
    f1->SetParameters(par1);
    f1->SetLineColor(2);
    hRawQ[i]->Fit("f1","R");
    hRawQ[i]->SetAxisRange(0,60);
    //gain[i]=f1->GetParameter(3);
    //width[i]=f1->GetParameter(4);
    pmtSummary->gain[i]=f1->GetParameter(3);
    pmtSummary->gain_e[i]=f1->GetParError(3);
    pmtGains->gain[i]=f1->GetParameter(3);
    pmtGains->egain[i]=f1->GetParError(3);

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

  myc1->Print(".pdf");
  myc2->Print(".pdf");


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

  // fit for tzero from this run 
  // find first rise !!this is only a first attempt should be improved!!
  Double_t zero=0;
  TF1 *gpfit[NB];
  TCanvas *cPromptFit[NB];
  for(int ib=0; ib<NB; ++ib) {
    for (int i=1;i<hTPrompt[ib]->GetNbinsX();i++){
      zero = hTPrompt[ib]->GetBinLowEdge(i);
      Double_t step = hTPrompt[ib]->GetBinContent(i+1)-hTPrompt[ib]->GetBinContent(i);
      if(step>5) break;
    }

    TString tnamePromptFit,tcanNamePromptFit, gpfitName;
    tnamePromptFit.Form("promptFit-%i-%s",ib,tag.Data());
    tcanNamePromptFit.Form("promptFit-%i",ib);
    gpfitName.Form("gPromptFit-%i",ib);
    cPromptFit[ib] = new TCanvas( tcanNamePromptFit.Data(), tcanNamePromptFit.Data());
    gpfit[ib]= new TF1(gpfitName.Data(),"gaus",-160,-156);
    gpfit[ib]->SetLineColor(2);
    printf("\t\t qualitySummary fitting to %s %s %f \n",gpfitName.Data(),hTPrompt[ib]->GetName(),hTPrompt[ib]->GetEntries());
    hTPrompt[ib]->Fit(gpfitName.Data(),"R");
    printf(" TPrompt fit board %i parameter = %f +/- %f low edge is %f \n",ib,gpfit[ib]->GetParameter(1),gpfit[ib]->GetParError(1),zero); 
    pmtSummary->tZero[ib] = gpfit[ib]->GetParameter(1);
    summaryFile->Append(cPromptFit[ib]);
    summaryFile->Append(gpfit[ib]);
  }
  // fill neutron spect
  printf(" qualitySummary calling fillNeutrons with  %zu \n",pmtSummary->vprompt1.size());
  pmtSummary->fillNeutrons();

  pmtGains->print();
  pmtSummary->print();
}


void pmtChain::ADCFilter(int iB, int iC) 
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
Int_t pmtChain::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t pmtChain::LoadTree(Long64_t entry)
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

void pmtChain::Init()
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!fChain) return;
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

Bool_t pmtChain::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  string fname = string(fChain->GetFile()->GetName());
  stag = fname.substr( fname.find_last_of("/")+8, fname.find(".") -1  - fname.find_last_of("/")-7); 
  cout << fCurrent << "  " << fChain->GetFile()->GetName() << " stag " << stag <<  endl;
  return kTRUE;
}

void pmtChain::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t pmtChain::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void pmtChain::getPromptTime()
{
  // fill histogram to find peak bin in event.
  for(int ib=0; ib<NB; ++ib) hTPromptEvent[ib]->Reset();

  for(unsigned ihit=0; ihit< pmtEvent->hit.size(); ++ihit) {
    Int_t hitTime = pmtEvent->hit[ihit].peakTime;
    Double_t qpeak = pmtEvent->hit[ihit].qpeak;
    Int_t ib,ic;
    fromPmtNumber(pmtEvent->hit[ihit].ipmt, ib,ic);
    hTPromptEvent[ib]->Fill(Int_t(hitTime),qpeak);
  }

  for(UInt_t ib=0; ib<NB; ++ib) pmtEvent->tPrompt[ib] = Double_t(hTPromptEvent[ib]->GetMaximumBin()) - pmtEvent->tRFave; 
}

// nearest RF time to hit peak time
void pmtChain::getTimeToRF() 
{
  for(unsigned ihit=0; ihit< pmtEvent->hit.size(); ++ihit) {
    Int_t time = MAXSAMPLES;
    Int_t hitTime = pmtEvent->hit[ihit].peakTime;
    for(unsigned i=0; i<pmtEvent->rft21.size() ; ++i) {
      Int_t tdiff = hitTime - pmtEvent->rft21[i];
      if(tdiff<0) break;
      if( tdiff<time ) time=tdiff;
    }
    for(unsigned i=0; i<pmtEvent->rft22.size() ; ++i) {
      Int_t tdiff = hitTime - pmtEvent->rft22[i];
      if(tdiff<0) break;
      if( tdiff<time ) time=tdiff;
    }
    for(unsigned i=0; i<pmtEvent->rft23.size() ; ++i) {
      Int_t tdiff = hitTime - pmtEvent->rft23[i];
      if(tdiff<0) break;
      if( tdiff<time ) time=tdiff;
    }
    pmtEvent->hit[ihit].timeToRF = time;
  }
}


// neils filter
std::vector<Double_t> pmtChain::MovingAverageFilter(std::vector<Double_t> signal,Int_t aveN)
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

void pmtChain::makeChain()
{
  fChain= new TChain("pmt_tree");
  // get list of files
  printf("Low intensity runs range from 0731_1518 to 0731_2130. They are all in one day. \n");
  //TString sumTag("low-intensity");
  TString dirname("/data1/gold/2017/PDS_beamtime_files/");
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles(); 
  
  if (!files) {
    cout << " makeFileList for "<<  dirname <<" returning with NULL directory pointer ! \n"; 
    return;
  }
  cout << dirname << " has " << files->GetSize() << " files " << endl;
   TSystemFile *file; 
   TIter next(files); 
   while ((file=(TSystemFile*)next())) { 
     string fname = string(file->GetName()); 
     //cout << fname << endl;
     if ( strstr(fname.c_str(), "PDSout_" )==NULL ) continue;
     if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
     //cout << fname << endl;
     getTag(fname); 
     Int_t month = getMonth();
     Int_t day =  getDay();
     Int_t min =  getMin();
     Int_t segment =  getSegment();
     Int_t time = min +60*24*day+60*24*30*(month-7)- 60*24*20;

     if( month==7 && day == 31 && min>1517 && min < 2131 ) {
       fileTime.push_back(time);
       timeMap.insert ( std::pair<int,std::string>(time,fname) );
     }
   }
   
  printf(" total of files in %s is %lu  \n",dirname.Data(),timeMap.size()-1);

  Int_t count =0;
  std::cout << "timeMap contains:\n";
  std::map<int,std::string>::iterator iter;
  for (iter=timeMap.begin(); iter!=timeMap.end(); ++iter) {
    TString addName = dirname + TString(iter->second.c_str());
    cout << "adding to tree " << ++count << " "  << iter->first << " => " << iter->second << "  file " << addName << endl;
    //cout << count++ << " "  << iter->first << " => " << iter->second << endl;
    fChain->Add(addName);
  }
  cout << " made chain pmt_tree with " <<  fChain->GetEntries() << " entries" << endl;
  // fChain->ls();
}


