#include <vector>
#include "TPmtEvent.hxx"
  
enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
  //TH1D* hQHitTime[NPMT];
  TPmtEvent *ev;
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


void readFile(UInt_t ifile) 
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

  loop();

  fin->Close();
  // write and close strip files
  printf(" stripped file %i,  %s noBeam events %lu lowBeam events %lu highBeam events %lu \n",  ifile, fileList[ifile].c_str(), 
      noBeamTree->GetEntries(), lowBeamTree->GetEntries(), highBeamTree->GetEntries());
  noBeamFile->Write(); noBeamFile->Close();
  lowBeamFile->Write(); lowBeamFile->Close();
  highBeamFile->Write();highBeamFile->Close();
  
  
}



void tof(TString tag="PDS_beamtime_files")
{
  // get list of files
  fullDirName= TString("../2017")+tag;
  void *dirp = gSystem->OpenDirectory(fullDirName);
  cout << "  tof full directory name is " << fullDirName << endl;
  if (!dirp) {
    cout << " returning with NULL directory pointer ! \n"; 
    return;
  }
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
    if( month>=1581 && month <= 1609) fileList.push_back(fname);
  }


  // open ouput file and make some histograms
  TString outputFileName = TString("tof-")+tag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());


   // loop over files
  unsigned ffirst =fFirst;
  unsigned fmax = ffirst+fileList.size();
  if(maxFiles>0) fmax=ffirst+UInt_t(maxFiles);
  printf(" now loop over files in %s is %lu reading from %u to %u \n",dirName.Data(),fileList.size(),ffirst,fmax);
  for( unsigned ifile =ffirst ; ifile < fmax ; ++ifile ) {
    printf(" %i %s \n",ifile,fileList[ifile].c_str());
    readFile(ifile);
    //printf(" have written %i bytes \n",summaryTree->FlushBaskets());
  }
  //summaryFile->Write();
  //cout<< "  summary tree has   " << summaryTree->GetEntries() << " entries " << endl;
  //printf(" number of good files %i \n number of bad files %i \n",goodFiles,badFiles);
  //summaryFile->Close();
  // end of tof 
  outfile->Write();
}
