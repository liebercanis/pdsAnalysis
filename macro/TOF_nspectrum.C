#include <iostream>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <string>
#include <sstream>

//Loop through low intensity runs to get TOF and Neutron spectrum
//require pmtAna_07-31-????_0.root
int low_intensity_start=1518;//1518 low intensity starts
int low_intensity_end=1609;//2130 low intensity ends
double L=23.2;//m
double c=0.299792458;//m/ns
double nmass=939.565;//MeV
double En;
int sample_size=2100;
void TOF_nspectrum()
{
  int threshold;
  int trigType,run,event,tpcTrig,pdsTrig,qmax_i,qsum,nhits;
  unsigned short gpsYear,gpsDay;
  int gpsSec,gpsNs;
  
  float ipmt_ntHit,ipmt_ntPmt,sum,time,rftime,length,qpeak,qUnpeak,qhit,qUnhit,fwhm,ratio;
  float trig,tmax,qmax,tmaxUn,qmaxUn,sumUn,noise,base,nhit;
  vector <double> run_number, trign[9];
  

    struct dirent *entry;
    DIR *dp;
    string path="/home/sunyujing/pdsf/pdsAnalysis_20171012/pdsOutput/";//directory of pmtAna_07-31-xxxx_0.root files
    dp = opendir(path.c_str());
     if (dp == NULL) {
	perror("opendir: Path does not exist or could not be read.");
	return -1;
     }
     TH1F *promptt_rft_h = new TH1F("promptt_rft_h","Prompt Time - RFtime;Time (ns);Number of entries",2*sample_size,-4*sample_size,4*sample_size); 
     TH1F *tof_h = new TH1F("tof_h","TOF (Assuming L = 23.2 m);Time (ns);Number of entries",2*sample_size,-4*sample_size,4*sample_size); 
     TH1F *nspectrum_h = new TH1F("neutron_spectrum_h",";Neutron E_{Kin} (MeV);Frac of triggers",400,0,4000);
     //TH1F *tof_h = new TH1F("tof","tof",4200,-2100,2100);
    while ((entry = readdir(dp))) {
     string filename=entry->d_name;
     std::size_t found=std::string::npos;
     found=filename.find("pmtAna_07-31-");
     if (found!=std::string::npos) {
       string runn = filename.substr (13,4);
       istringstream buffer(runn);
       int value;
       buffer >> value;
       if(value>=low_intensity_start && value<=low_intensity_end)
       {
	  cout<<"analysing run pmtAna_07-31-"<<value<<".root"<<endl;
	  run_number.push_back(double (value));
	  TString outputFileName = path+filename;
	  TFile *fi = new TFile(outputFileName);
    
	  TTree *pmtTree = (TTree*) fi->Get("pmtTree");
	  Long64_t aSize=0;
	  if(pmtTree) aSize=pmtTree->GetEntriesFast();

	  Long64_t aSize=0;
	  pmtTree = (TTree*) fi->Get("pmtTree");
	  if(pmtTree) aSize=pmtTree->GetEntriesFast();

	  TNtuple* ntHit = (TNtuple*)(fi->Get("ntHit"));
  
	  ntHit->SetBranchAddress("ipmt", &ipmt_ntHit);
	  ntHit->SetBranchAddress("sum", &sum);
	  ntHit->SetBranchAddress("time", &time);
	  ntHit->SetBranchAddress("rftime", &rftime);
	  ntHit->SetBranchAddress("length", &length);
	  ntHit->SetBranchAddress("qpeak", &qpeak);
	  ntHit->SetBranchAddress("qUnpeak", &qUnpeak);
	  ntHit->SetBranchAddress("qhit", &qhit);
	  ntHit->SetBranchAddress("qUnhit", &qUnhit);
	  ntHit->SetBranchAddress("fwhm", &fwhm);
	  ntHit->SetBranchAddress("ratio", &ratio);  

	  TNtuple* ntPmt = (TNtuple*)(fi->Get("ntPmt"));
	  
	  ntPmt->SetBranchAddress("trig", &trig);
	  ntPmt->SetBranchAddress("ipmt", &ipmt_ntPmt);
	  ntPmt->SetBranchAddress("tmax", &tmax);
	  ntPmt->SetBranchAddress("qmax", &qmax);
	  ntPmt->SetBranchAddress("tmaxUn", &tmaxUn);
	  ntPmt->SetBranchAddress("qmaxUn", &qmaxUn);
	  ntPmt->SetBranchAddress("sumUn",   &sumUn);
	  ntPmt->SetBranchAddress("noise",&noise);
	  ntPmt->SetBranchAddress("base",   &base);
	  ntPmt->SetBranchAddress("nhit", &nhit);
  
	  //read ntPmt
	  int trigger;
	  int nentries_ntPmt = ntPmt->GetEntries();
	  double totalqUn;

	  if(aSize==0) return;

	  TPmtEvent *ev = new TPmtEvent();
	  pmtTree->SetBranchAddress("pmtEvent",&ev);
	  TFile *outputfile = new TFile(Form("TOF_nspectrum_%d_%d.root",low_intensity_start,low_intensity_end),"recreate");
	  std::vector<TPmtHit> hit;
		  
	  int counter=0;
	  for(unsigned ientry =0; ientry < aSize; ++ientry ) {
	    pmtTree->GetEntry(ientry);
	    totalqUn=0;

	    hit.clear();
	    hit = ev->hit;
	    TH1F *sum_h = new TH1F("sum","sum",2100,0,2100);// temperarily sum all the peaks and take the maximum as prompt time; eventually will use summed channel

	    ntPmt->GetEntry(ientry*21);
	    if(trig==5) {
		for(int ihit =0; ihit < hit.size(); ++ihit) {
		      TPmtHit* phit = &(ev->hit[ihit]);
		      sum_h->Fill(phit->peakTime,phit->qpeak);
		      ntHit->GetEntry(counter);
		      totalqUn+=qUnhit;
		      counter++;
		    }
		if(totalqUn>threshold)promptt_rft_h->Fill((sum_h->GetMaximumBin()-rftime)*4);
		if(totalqUn>threshold && (sum_h->GetMaximumBin()-rftime)*4<-800)cout<<ientry<<endl;
	      }
	    else {
		counter+=hit.size();
	      }
	    delete sum_h;
	  }
        }
      }
    }//end of while loop
    
    
  TF1 *g1 = new TF1("g1","gaus",-640,-620);
  g1->SetLineColor(2);
  promptt_rft_h->Fit("g1","R");
  //promptt_rft_h->SetAxisRange(-700,-200);
  
  cout<<g1->GetParameter(1)<<endl;
  for (int i=1;i<=4200;i++)
  {
    if(i>1945)tof_h->SetBinContent(i-g1->GetParameter(1)/4.+L/c/4.,promptt_rft_h->GetBinContent(i));
  }

  for (int i=0;i<1000000;i++)
  {
     double t=tof_h->GetRandom(); 
     En=nmass*(sqrt(1/((t*c/L)**2-1)+1)-1);
     nspectrum_h->Fill(En);
  }  
  nspectrum_h->Scale(1./1000000);
   
  outputfile->cd();
  promptt_rft_h->Write();
  tof_h->Write();
  nspectrum_h->Write();
  outputfile->Close();
  delete outputfile;


}