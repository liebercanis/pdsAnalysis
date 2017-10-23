/*
** summary of PDS data 
*/
#include "TPdsSummary.hxx"
ClassImp(TPdsSummary)


TPdsSummary::TPdsSummary(TString theDirName): TNamed("TPdsSummary","TPdsSummary")
{
  dirName = theDirName;
  badFiles=0;
  goodFiles=0;
  // get list of files
  fullDirName= TString("2017/")+dirName;
  void *dirp = gSystem->OpenDirectory(fullDirName);
  cout << " TPdsSummary full directory name is " << fullDirName << endl;
  if (!dirp) {
    cout << " returning with NULL directory pointer ! \n"; 
    return;
  }
  char *direntry;
  Long_t id, size,flags,modtime;
  //loop on all entries of this directory
  fileList.clear();
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) { 
  printf(" total of files in %s is %lu  \n",dirName.Data(),fileList.size());
    //cout << direntry << endl;
    string fname = string(direntry);
    if ( strstr(fname.c_str(), "PDSout" )==NULL ) continue;
    if ( strstr(fname.c_str(), ".root" )== NULL ) continue;
    getTag(fname);
    Int_t month = getMonth();
    Int_t day =  getDay();
    Int_t min =  getMin();
    Int_t segment =  getSegment();
    if( month==7 && day < 20) {
      printf(" skipping early run formatted fill, tag is %s month %i day %i min %i segment %i \n",tag.c_str(),month,day,min,segment);
      continue;
    } 
    fileList.push_back(fname);
  }
  // why do I have to do this??
  isEmpty = fileList.empty();
  if(isEmpty) { 
    printf(" \t WARNINNG:: file list is empty \n");
    fileList.clear();
  }
  printf(" total of files in %s is %lu  \n",dirName.Data(),fileList.size());
}

TPdsSummary::~TPdsSummary(){}


void TPdsSummary::run(Int_t fFirst, Int_t maxFiles) 
{ 
  if(isEmpty) {
    printf(" run returning because fileList is empty\n");
    return;
  }
  // structure for holding pmt info 
  //open output file
  TString fileTag;
  fileTag.Form("files-%i-%i_",fFirst,maxFiles);
  TString summaryFileName = TString("pdsOutput/pdsSummary_") + fileTag +dirName + TString(".root");
  summaryFile = new TFile(summaryFileName,"recreate");
  summaryFile->cd();
  printf(" opening summary file %s \n",summaryFileName.Data());
  summaryTree = new TTree("summaryTree","summaryTree");
  pmtSummary  = new TPmtSummary();
  summaryTree->Branch("pmtSummary",&pmtSummary);
  summaryFile->ls();



  // loop over files
  unsigned ffirst =fFirst;
  unsigned fmax = ffirst+fileList.size();
  if(maxFiles>0) fmax=ffirst+UInt_t(maxFiles);
  printf(" now loop over files in %s is %lu reading from %u to %u \n",dirName.Data(),fileList.size(),ffirst,fmax);
  for( unsigned ifile =ffirst ; ifile < fmax ; ++ifile ) {
    printf(" %i %s \n",ifile,fileList[ifile].c_str());
    readFile(ifile);
    printf(" have written %i bytes \n",summaryTree->FlushBaskets());
  }
  summaryFile->Write();
  cout<< "  summary tree has   " << summaryTree->GetEntries() << " entries " << endl;
  printf(" number of good files %i \n number of bad files %i \n",goodFiles,badFiles);
  summaryFile->Close();
}



std::vector<Int_t> TPdsSummary::findRFTimes(int ipmt, double& step) 
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
    if(digi<0.75*step&&!isRF) {
      rftimes.push_back(is);
      isRF=true;
    } else if(digi>0.75*step) {
      isRF=false;
    }
  }
  return rftimes;
}


// trigger information
Int_t TPdsSummary::triggerInfo()
{
  Int_t type = TPmtEvent::TRIGUNKNOWN; // unknosn 
  // RF channels 
  double s1,s2,s3;

  std::vector<Int_t> rftime21 = findRFTimes(21,s1);
  std::vector<Int_t> rftime22 = findRFTimes(22,s2);
  std::vector<Int_t>  rftime23 = findRFTimes(23,s3);
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


void TPdsSummary::ADCFilter(int iB, int iC) 
{
  for (int is = 0; is<NS; ++is) {
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


void TPdsSummary::loop() 
{
 
  Long64_t entries = pmt_tree->GetEntries(); 
  cout << "Looping over tag " << tag << " pmt_tree has entries = " << entries << endl;
  pmtSummary->clear();
  // save the file tag
  pmtSummary->tag = tag;
  std::vector<Double_t> sdigi;  // source
  std::vector<Double_t> ddigi;  // baseline subtracted

  // collect info for summary 
  std::vector<Int_t> vpmt;
  std::vector<Int_t> vtmax;
  std::vector<Double_t> vqmax;
  std::vector<Double_t> vsum;
  std::vector<Double_t> vnoise;

  // loop over events 
  for(Long64_t ientry=0; ientry<entries; ++ientry) {
    pmt_tree->GetEntry(ientry);
    if(ientry%1000==0) printf(" ..... %lld \n",ientry);
    Int_t trigType = triggerInfo();
    // output based on trigger type
    if(trigType==TPmtEvent::TRIG000) noBeamTree->Fill();
    if(trigType==TPmtEvent::TRIG111) lowBeamTree->Fill();
    if(trigType==TPmtEvent::TRIG555 || trigType==TPmtEvent::TRIG444 ) highBeamTree->Fill();

    // store all the event time info
    pmtSummary->vtrig.push_back(trigType);
    pmtSummary->vevent.push_back(event);
    pmtSummary->ventry.push_back(ientry);    
    pmtSummary->vcompSec.push_back(compSec);
    pmtSummary->vcompNano.push_back(compNano);
    //pmtSummary->vgpsNs.push_back(gpsNs);
    //pmtSummary->vgpsSec.push_back(gpsSec);
    //pmtSummary->vgpsDay.push_back(gpsDay);
    //pmtSummary->vgpsYear.push_back(gpsYear);
    //pmtSummary->print();
    

    for(UInt_t ib=0; ib<NB; ++ib) {
      for(UInt_t ic=0; ic<NC; ++ic) {
          // filter waveforms for stuck bits
          ADCFilter(ib,ic);
          // get pmt number
          int ipmt = toPmtNumber(ib,ic);
          if(ipmt<0||ipmt>=NPMT) continue;
              
          // make a vector of samples for sorting.
          sdigi.clear();
          ddigi.clear();
          
          // Find the sample median and it's "sigma".
          for (UInt_t is=0; is<NS; ++is) {
            sdigi.push_back(double(digitizer_waveforms[ib][ic][is]));
          }

          std::sort(sdigi.begin(), sdigi.end());
          double baselineMedian = sdigi[0.5*double(NS)];
          double baselineSigma = sdigi[0.16*double(NS)];
          baselineSigma = std::abs(baselineSigma-baselineMedian);
          double noise = std::abs( sdigi[0.68*sdigi.size()] - baselineMedian);
          
          // loop over digitizations
          Int_t tmax=0;
          double qmax=0;
          double sum=0; 
           for(UInt_t is=0 ; is<NS; ++is) {
            double digi = -1.0*(double(digitizer_waveforms[ib][ic][is])-baselineMedian);
            if(digi>qmax) {
              qmax=digi;
              tmax=Int_t(is)+1;
            }
            ddigi.push_back(digi);
            if(digi>3.0*noise) { 
              //if(is>450&&is<470) 
              sum+=digi;
            }
          } // loop over digitizations
          //if(ientry%1000==0) printf(" \t ipmt %i trig %i tmax %i qmax %f sum %f noise %f \n",ipmt,trigType,tmax,qmax,sum,noise);
          vpmt.push_back(ipmt);
          vtmax.push_back(tmax);
          vqmax.push_back(qmax);
          vsum.push_back(sum);
          vnoise.push_back(noise);
      }
    }
  }
  // fill rest of summary
  // ugly utility arrays used in summary
  Double_t x[NPMT], y[NPMT], z[NPMT],y2[NPMT],z2[NPMT],ex[NPMT], ey[NPMT], ez[NPMT];
  Double_t norm[NPMT], aveNoise[NPMT]; 

  // zero utility arrays
  for(Int_t j=0; j<NPMT; ++j) {
    x[j]=Double_t(j); ex[j]=0;  
    y[j]=0; z[j]=0; y2[j]=0; z2[j]=0; ey[j]=0; ez[j]=0;
    norm[j]=0; aveNoise[j]=0;
  }
  
  printf(" \n \n \t calculating averages with %lu entries \n",vpmt.size());
 
  for (UInt_t k=0 ;k<vpmt.size(); k++) {
    int ipmt = vpmt[k];
    // cosmic cut
    bool cut = vtmax[k]>440&&vtmax[k]<480&&vsum[k]<5000;
    //if(k%100==0) printf(" \t ipmt %i norm %.0f trig %i tmax %i qmax %f sum %f noise %f \n",ipmt,norm[ipmt],vpmt[k],vtrig[k],vtmax[k],vqmax[k],vsum[k],vnoise[k]);
    if(cut) { 
      y[ipmt]+=vqmax[k];
      y2[ipmt]+=pow(vqmax[k],2.);
      z[ipmt]+=vsum[k];
      z2[ipmt]+=pow(vsum[k],2.);
      aveNoise[ipmt]+=vnoise[k];
      norm[ipmt]+=1.0;
    }
  }
  for(Int_t j=0; j<NPMT; ++j) {
    y[j]/= norm[j]; z[j]/=norm[j]; y2[j]/=norm[j]; z2[j]/=norm[j]; aveNoise[j]/=norm[j];
  }

  for(Int_t j=0; j<NPMT; ++j) {
    ey[j]= sqrt( (y2[j]-pow(y[j],2.))/norm[j]);
    ez[j]= sqrt( (z2[j]-pow(z[j],2.))/norm[j]);
  }

  //for(Int_t j=0; j<NPMT; ++j) 
  //  printf(" ipmt %i norm %i qmax %.2f +/- %.2f sum %.2f +/- %.2f \n",j,int(norm[j]),y[j],ey[j],z[j],ez[j]);

  // store averages in summary
  for(Int_t j=0; j<NPMT; ++j) {
    pmtSummary->norm[j]=norm[j];
    pmtSummary->qmax[j]=y[j];
    pmtSummary->eqmax[j]=ey[j];
    pmtSummary->qsum[j]=z[j];
    pmtSummary->eqsum[j]=ez[j];
    pmtSummary->noise[j]=aveNoise[j];
  }

  pmtSummary->print();
  summaryTree->Fill();
}

void TPdsSummary::readFile(UInt_t ifile) 
{
  TString fullName = fullDirName+TString("/")+TString(fileList[ifile].c_str());
  cout << " reading file " << fullName << endl;
  TFile*  fin = new TFile(fullName, "READ");
  if(fin->IsZombie()) {
    printf(" cannot read file %s\n",fullName.Data());
    ++badFiles;
    return;
  }
  pmt_tree = (TTree *)fin->Get("pmt_tree");
  if(!pmt_tree) {
    printf(" cannot find pmt_tree in file %s\n",fullName.Data());
    fin->Close();
    ++badFiles;
    return;
  }

  pmt_tree->SetBranchAddress("event_number", &event);
  pmt_tree->SetBranchAddress("computer_secIntoEpoch", &compSec);
  pmt_tree->SetBranchAddress("computer_nsIntoSec", &compNano);
  pmt_tree->SetBranchAddress("gps_nsIntoSec", &gpsNs);
  pmt_tree->SetBranchAddress("gps_secIntoDay", &gpsSec);
  pmt_tree->SetBranchAddress("gps_daysIntoYear", &gpsDay);
  pmt_tree->SetBranchAddress("gps_Year", &gpsYear);
  pmt_tree->SetBranchAddress("digitizer_waveforms", &digitizer_waveforms);
  
  getTag(fileList[ifile]); 
  ++goodFiles;

  // open strip output files
  TString noBeamFileName = TString("2017/") + dirName+ TString("/pdsNoBeam-")  + TString(tag.c_str()) + TString(".root");
  noBeamFile = new TFile(noBeamFileName,"recreate");
  noBeamFile->cd();
  pmt_tree->ls();
  noBeamTree = pmt_tree->CloneTree(0);
  printf(" opening pdsNoBeam file %s \n",noBeamFileName.Data());
  noBeamTree->ls();

  TString lowBeamFileName =  TString("2017/") + dirName+ TString("/pdsLowBeam-")  + TString(tag.c_str()) + TString(".root");
  lowBeamFile = new TFile(lowBeamFileName,"recreate");
  lowBeamFile->cd();
  lowBeamTree = pmt_tree->CloneTree(0);
  printf(" opening pdsNoBeam file %s \n",lowBeamFileName.Data());
  lowBeamTree->ls();

  TString highBeamFileName =  TString("2017/") + dirName+ TString("/pdsHighBeam-")  + TString(tag.c_str()) + TString(".root");
  highBeamFile = new TFile(highBeamFileName,"recreate");
  highBeamFile->cd();
  highBeamTree = pmt_tree->CloneTree(0);
  printf(" opening pdsNoBeam file %s \n",highBeamFileName.Data());
  highBeamTree->ls();

  loop();

  fin->Close();
  // write and close strip files
  printf(" stripped file %i,  %s noBeam events %lu lowBeam events %lu highBeam events %lu \n",  ifile, fileList[ifile].c_str(), 
      noBeamTree->GetEntries(), lowBeamTree->GetEntries(), highBeamTree->GetEntries());
  noBeamFile->Write(); noBeamFile->Close();
  lowBeamFile->Write(); lowBeamFile->Close();
  highBeamFile->Write();highBeamFile->Close();
  
  
}

