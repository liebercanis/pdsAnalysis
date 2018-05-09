#include "pmtChain.hh"
#include <TF1.h>
#include <TPaveStats.h>
#include <TText.h>
/*
** remove bad events. adapted from pmtAna so has lots of unused, irrelevant code
** this version does board fixes
*/
pmtChain::pmtChain(Int_t maxLoop, Long64_t firstEntry)
{
  //maxLoop = 5000;
  //firstEntry = 5000*12;
  TTree *tree=NULL;
  TString tag("low_intensity");
  makeChain();  
  if(!fChain) return;
  Init();
  //cout << " fChain " << endl; fChain->ls(); fChain->GetListOfBranches()->ls();
  
 
  // open ouput file and make some histograms
  TString spost;
  spost.Form("-fix4-%i-%lli",int(maxLoop),firstEntry);
  TString outputFileName = TString("pdsOutput/pmtChain") + spost + TString(".root");
  outFile = new TFile(outputFileName,"recreate");
  InitOutChain();
  cout << " outTree " << endl; outTree->GetListOfBranches()->ls();
  
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
  // zero event number buff 
  for(int ib=0; ib<NB; ++ib) eventNumberBuff[ib]=0;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes=0;
  UInt_t nloop=nentries;
 
  if(nToLoop!=0) nloop = nToLoop;
  float rawLast[4]={0,0,0,0};
  float rdiff[4]={0,0,0,0};
  float jump = 195.0E6;
 
  for(Int_t ir=0; ir< MAXRUN; ++ir) misAlignCount[ir]=0;
  printf(" entries %lld looping %d first %d \n",nentries,nloop,firstEntry);
  // loop over entries
  for (Long64_t jentry=firstEntry; jentry<nloop+firstEntry; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) { printf(" load tree returns %lld\n",ientry); break;}
    nbytes += fChain->GetEntry(jentry);
    UInt_t runEvent = UInt_t(jentry-5000*(jentry/5000));
    if(jentry%1000==0) printf(" \t.... %lld nbytes runEvent %u fCurrent %u \n",jentry,runEvent,fCurrent);
    // fill output chain
    copyChain();
    oevent_number=runEvent; // dont want UInt_t to overflow
    //cout << " runevent  " << runEvent << " fCurrent " << fCurrent << endl;
    ogps_nsIntoSec=minList[fCurrent];  
    ogps_secIntoDay=segList[fCurrent];
    ogps_ctrlFlag=UInt_t(fCurrent);

    ULong_t pdsCompTime = ocomputer_secIntoEpoch*1000000000 + ocomputer_nsIntoSec ;
    
   
    bool skipEvent = false;
    // fix the output stream 
    // NFIX is all 
    for(int ifix=0; ifix<NFIX; ++ifix) {
      if( skipMin[ifix] == Int_t(minList[fCurrent]) && abs(skipFirst[ifix]) == runEvent ) {
        if( skipFirst[ifix] > 0) {
          eventNumberBuff[skipBoard[ifix]] += skipLast[ifix]-skipFirst[ifix]+1;
          cout << " \t IFIX update eventNumber file " << fCurrent << " min "  << minList[fCurrent] << " runEvent " << runEvent 
           << " file event " << outTree->GetEntries() << " eventNumberBuff[" << skipBoard[ifix] << "]= "  << eventNumberBuff[skipBoard[ifix]] << endl;
        } 
        else {// all misses are single events 
          eventNumberBuff[skipBoard[ifix]] -= 1;
          cout << " \t IFIX update eventNumber file " << fCurrent << " min "  << minList[fCurrent] << " runEvent " << runEvent 
           << " file event " << outTree->GetEntries() << " eventNumberBuff[" << skipBoard[ifix] << "]= "  << eventNumberBuff[skipBoard[ifix]] << endl;
        }

      }
    }
    
    if(skipEvent) {
      cout << " \t IFIX skipEvent file " << fCurrent << " min "  << minList[fCurrent] << " runEvent " << runEvent << " file event " << outTree->GetEntries() << endl;
      continue;  // toss event
    }

    for(int ib=0; ib<NB; ++ib) if( eventNumberBuff[ib]!=0 ) {
      getBoard(jentry + eventNumberBuff[ib], ib);
      fillBoard(ib);
    }

    // alignment check
    for(int ib=0; ib<NB; ++ib) {
      rdiff[ib]= float(odigitizer_time[ib]) - rawLast[ib];
      rawLast[ib]=float(odigitizer_time[ib]);
    }
    rdiff[3]   = float(pdsCompTime) - rawLast[3];
    rawLast[3] =float(pdsCompTime);
    if(rdiff[0]>0&&rdiff[1]>0&&rdiff[2]>0) {
      if(rdiff[3] > jump && ( rdiff[0]*8<jump || rdiff[1]*8<jump || rdiff[2]*8<jump)){ 
        printf(" MISALIGN event %lld run %i total %i  %f (%f ,%f , %f) \n",jentry,fCurrent,++misAlignCount[fCurrent],rdiff[3],rdiff[0]*8,rdiff[1]*8,rdiff[2]*8); 
      }
    }
    //

    /*
     printf(" \t  ev %llu run ev %i run %i  (%u %u %u ) (%u %u %u )    \n",jentry,int(runEvent),int(fCurrent),
         digitizer_time[0],digitizer_time[1],digitizer_time[2],
         odigitizer_time[0],odigitizer_time[1],odigitizer_time[2]);
         */
    
    outTree->Fill();

  } // end loop over entries
  printf(" finished looping  %u outTree size %llu \n",nloop,outTree->GetEntries());
  for(Int_t ir=0; ir< MAXRUN; ++ir) printf(" run %i misAlignCount %i \n",ir,misAlignCount[ir]);
  return nloop;
}

void pmtChain::makeChain()
{
  minList.clear();
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
     Int_t time = 10*(min +60*24*day+60*24*30*(month-7)- 60*24*20) + segment;

     if( month==7 && day == 31  && min> 1517 && min < 2135 ) {
       fileTime.push_back(time);
       timeMap.insert ( std::pair<int,std::string>(time,fname) );
       //cout << tag << " " << timeMap.size() <<  " " << minList.size() << endl;
     }
   }
   
  printf(" total of files in %s is %lu (%lu)  \n",dirname.Data(),timeMap.size(), minList.size() );

  Int_t count =0;
  std::cout << "timeMap contains:\n";
  std::map<int,std::string>::iterator iter;
  for (iter=timeMap.begin(); iter!=timeMap.end(); ++iter) {
    TString addName = dirname + TString(iter->second.c_str());
    getTag(iter->second.c_str()); 
    Int_t month = getMonth();
    Int_t day =  getDay();
    Int_t min =  getMin();
    Int_t segment =  getSegment();
    //cout << "adding to tree " << ++count << " "  << iter->first << " => " << iter->second << "  file " << addName << endl;
    //cout << count++ << " "  << iter->first << " => " << iter->second << endl;
    minList.push_back(UInt_t(min));
    segList.push_back(UInt_t(segment));
    fChain->Add(addName);
  }
  cout << " made chain pmt_tree with " <<  fChain->GetEntries() << " entries" << endl;
  fChain->ls();
  cout << " list of files (minute) \n";
  for(unsigned it=0; it< minList.size() ;  ++it) cout << "    " << it << "  " << minList[it] << " , ";
  cout << endl;
}

void pmtChain::getBoard(Long64_t entry, Int_t ib)
{

  // initialize to zero
  sizeBuff[ib]=0;
  evNumBuff[ib]=0;
  timeBuff[ib]=0;
  for(int ic=0; ic<NC; ++ic) {
    chMaskBuff[ib][ic]=0;
    for(int is=0; is<MAXSAMPLES; ++is) waveBuff[ib][ic][is]=0;
  }

  // load local memory from fChain for event
  Long64_t ientry = LoadTree(entry);
  if (ientry < 0) { printf(" load tree returns %lld\n",ientry); return;}
  fChain->GetEntry(entry);

  // fill buffers
  sizeBuff[ib]=digitizer_size[ib];
  evNumBuff[ib]=digitizer_evNum[ib];
  timeBuff[ib]=digitizer_time[ib];
  for(int ic=0; ic<NC; ++ic) {
    chMaskBuff[ib][ic]=digitizer_chMask[ib][ic];
    for(int is=0; is<MAXSAMPLES; ++is) waveBuff[ib][ic][is]=digitizer_waveforms[ib][ic][is];
  }
  //printf(" calling getBoard %lli board %i time %u \n ",entry,ib,timeBuff[ib]);
} 

// fill output for this board from buff
void pmtChain::fillBoard(Int_t ib)
{
  odigitizer_size[ib]=sizeBuff[ib];
  odigitizer_evNum[ib]=evNumBuff[ib];
  odigitizer_time[ib]=timeBuff[ib];
  for(int ic=0; ic<NC; ++ic) {
    odigitizer_chMask[ib][ic]=chMaskBuff[ib][ic];
    for(int is=0; is<MAXSAMPLES; ++is) odigitizer_waveforms[ib][ic][is]=waveBuff[ib][ic][is];
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
void pmtChain::InitOutChain()
{
  outTree= new TTree("pdsTree"," aligned raw data ");
  outTree->Branch("event_number",&oevent_number,"event_number/i");
  outTree->Branch("computer_secIntoEpoch", &ocomputer_secIntoEpoch,"computer_secIntoEpoch/I");
  outTree->Branch("computer_nsIntoSec", &ocomputer_nsIntoSec,"computer_nsIntoSec/L");
  outTree->Branch("gps_nsIntoSec", &ogps_nsIntoSec,"gps_nsIntoSec/i");
  outTree->Branch("gps_secIntoDay", &ogps_secIntoDay,"gps_secIntoDay/i");
  outTree->Branch("gps_daysIntoYear", &ogps_daysIntoYear,"gps_daysIntoYear/s");
  outTree->Branch("gps_Year", &ogps_Year,"gps_Year/s");
  outTree->Branch("gps_ctrlFlag", &ogps_ctrlFlag,"gps_ctrlFlag/s");
  outTree->Branch("digitizer_size", odigitizer_size,"digitizer_size[3]/i");
  outTree->Branch("digitizer_chMask", odigitizer_chMask,"digitizer_chMask[3][8]/i");
  outTree->Branch("digitizer_evNum", odigitizer_evNum,"digitizer_evNum[3]/i");
  outTree->Branch("digitizer_time", odigitizer_time,"digitizer_time[3]/i");
  outTree->Branch("digitizer_waveforms", odigitizer_waveforms,"digitizer_waveforms[3][8][2100]/s");
  outTree->Branch("nDigitizers", &onDigitizers,"nDigitizers/i");
  outTree->Branch("nChannels", &onChannels,"nChannels/i");
  outTree->Branch("nSamples", &onSamples,"nSamples/i");
  outTree->Branch("nData", &onData,"nData/i");
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

