#include <iostream>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <dirent.h>
#include <TFile.h>
#include <TTree.h>

//Loop through low intensity runs to get TOF and Neutron spectrum
//require pmtAna_07-31-????_0.root
int low_intensity_start=1600;//1518 low intensity starts
int low_intensity_end=1600;//2130 low intensity ends
double L=23.2;//m
double c=0.299792458;//m/ns
double nmass=939.565;//MeV
double En;
int sample_size=2100;
double gain[21],gaine[21],rgain[21];
int NPMT=21;
double bins[481];
double promptt_rft;

Int_t readGainConstants()
{
  TString filename("pmtGoodGains_07-31-1555_0.root");
  TFile *fgain = new TFile(filename,"READONLY");
  if(fgain->IsZombie()) { printf(" no gain file %s found \n",filename.Data()); return 0;}
  TTree *gtree=NULL;
  fgain->GetObject("gainsTree",gtree);
  if(!gtree) { printf(" no gainsTree not found \n"); return 0;}
  TPmtGains *goodGains = new TPmtGains();
  gtree->SetBranchAddress("pmtGains",&goodGains);
  gtree->GetEntry(0);
  printf(" \n \t using gains from file %s\n",filename.Data()); 
  goodGains->print();
  for(int ipmt=0; ipmt<NPMT; ++ipmt) {gain[ipmt]=goodGains->gain[ipmt];}
  for(int ipmt=1; ipmt<NPMT; ++ipmt) {rgain[ipmt]=goodGains->gain[ipmt]/goodGains->gain[0];}
  rgain[0]=1.0;
  printf(" using normalized gains \n");
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f ; ",ipmt,gain[ipmt]); 
  printf(" \n");
  //for(int ipmt=1; ipmt<NPMT; ++ipmt) cout<<gain[ipmt]<<endl;
  return 21;
}

void TOF_nspectrum()
{
  
  for (int tick=500;tick>=20;tick--)
  {
  bins[500-tick]=nmass*(sqrt(1/((tick*4*c/L)*(tick*4*c/L)-1)+1)-1);
  //cout<<bins[500-tick]<<endl;
  }
  
  if(readGainConstants()==0) {
    printf(" cannot read gain constants file so abort \n");
    return;
  }
  int rft;
  int threshold;
  int trigType,run,event,tpcTrig,pdsTrig,qmax_i,qsum,nhits;
  unsigned short gpsYear,gpsDay;
  int gpsSec,gpsNs;
  
  float ipmt_ntHit,ipmt_ntPmt,sum,time,rftime,length,qpeak,qUnpeak,qhit,qUnhit,fwhm,ratio;
  float trig,tmax,qmax,tmaxUn,qmaxUn,sumUn,noise,base,nhit;
  vector <double> run_number, trign[9];
  
  TFile *outputfile = new TFile(Form("TOF_nspectrum_pulse_%d_%d.root",low_intensity_start,low_intensity_end),"recreate");
  
    struct dirent *entry;
    DIR *dp;
    string path="/home/sunyujing/pdsf/pdsAnalysis/pdsOutput/";//directory of pmtAna_07-31-xxxx_0.root files
    dp = opendir(path.c_str());
     if (dp == NULL) {
	perror("opendir: Path does not exist or could not be read.");
	return -1;
     }
     TH1F *promptt_rft_h = new TH1F("promptt_rft_h","Prompt Time - RFtime;Time (ns);Number of entries",2*sample_size,-4*sample_size,4*sample_size); 
     TH1F *tof_h = new TH1F("tof_h","TOF (Assuming L = 23.2 m);Time (ns);Number of entries",2*sample_size,-4*sample_size,4*sample_size); 
     TH1F *nspectrum_h = new TH1F("neutron_spectrum_h",";Neutron E_{Kin} (MeV);Frac of triggers",480,bins);
     
     TH2F *height_integral_ADC_h = new TH2F("height_integral_ADC_h",";PMT Pulse Height (ADC);PMT Pulse Integral (ADC)",500,0,10000,1500,0,30000);
     TH2F *height_integral_PE_h = new TH2F("height_integral_PE_h",";PMT Pulse Height (PE);PMT Pulse Integral (ADC)",500,0,1000,1500,0,30000);
     
     TH2F *lightyield_nenergy_ADC_h = new TH2F("lightyield_nenergy_ADC_h","Light Yield VS n-Energy;n Energy (MeV);Light Yield (ADC)",480,bins,500,0,10000);
     TH2F *lightyield_nenergy_PE_h = new TH2F("lightyield_nenergy_PE_h","Light Yield VS n-Energy;n Energy (MeV);Light Yield (PE)",480,bins,500,0,1000);
     
     TH1F *prompt_height_total_ratio_ADC_h = new TH1F("prompt_height_total_ratio_ADC_h","Prompt/Total;Ratio;Number of entries",100,0,1); 
     TH1F *prompt_height_total_ratio_PE_h = new TH1F("prompt_height_total_ratio_PE_h","Prompt/Total;Ratio;Number of entries",100,0,1);
     TH1F *prompt_integral_total_ratio_h = new TH1F("prompt_integral_total_ratio_h","Prompt/Total;Ratio;Number of entries",100,0,1);

     TH1F *light_timing_h = new TH1F("light_timing_h",";Time (ns);Number of entries",10000,-10000,30000);
     //TH1F *tof_h = new TH1F("tof","tof",4200,-2100,2100);
     
    while ((entry = readdir(dp))) {
     string filename=entry->d_name;
     std::size_t found=std::string::npos;
     found=filename.find("pmtAna_07-31-");
     if (found!=std::string::npos) {
       string runn = filename.substr (13,4);
       //cout<<runn<<endl;
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

	  pmtTree = (TTree*) fi->Get("pmtTree");
	  if(pmtTree) aSize=pmtTree->GetEntriesFast();

	  TNtuple* ntHit = (TNtuple*)(fi->Get("ntHit"));
  
	  ntHit->SetBranchAddress("ipmt", &ipmt_ntHit);
	  ntHit->SetBranchAddress("sum", &sum);
	  ntHit->SetBranchAddress("time", &time);
	  //ntHit->SetBranchAddress("rftime", &rftime);
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
	  
	  //ev->print();
	  //return;
	  std::vector<TPmtHit> hit;
		  
	  int counter=0;
	  //aSize=100;
	  for(unsigned ientry =0; ientry < aSize; ++ientry ) {
	    
	    pmtTree->GetEntry(ientry);
	    totalqUn=0;

	    hit.clear();
	    hit = ev->hit;
	    //cout<<"####################  RF pulses #####################"<<endl;
	    //cout<<"ev->rft21: "<<ev->run<<endl;
	    //cout<<"ev->rft22: "<<ev->event<<endl;
	    ntPmt->GetEntry(ientry*21);

	    //cout<<endl;
	    TH1F *sum_h = new TH1F("sum","sum",4200,-2100,2100);// temperarily sum all the peaks and take the maximum as prompt time; eventually will use summed channel
	    //TCanvas *myc = new TCanvas("myc","histograms",1);

   
	    //cout<<"event "<<ientry<<endl;
	    double ly=0,lyPE=0;
	    if(trig==5) {
	      	    rft=min(ev->rft21[0], ev->rft22[0]);
		    rft=min(rft, ev->rft23[0]); 
		    //cout<<hit.size()<<endl;
		for(int ihit =0; ihit < hit.size(); ++ihit) {
		      TPmtHit* phit = &(ev->hit[ihit]);

		      int shiftt;
		      if(phit->ipmt<7)shiftt=ev->rft21[0]-rft;
		      else if(phit->ipmt<14)shiftt=ev->rft22[0]-rft;
		      else if(phit->ipmt<21)shiftt=ev->rft23[0]-rft;
		      
		      double pulse_height=0;
		      double pulse_integral=0;
		      for (int i=0;i<phit->nsamples;i++) {
		      sum_h->Fill(phit->tsample[i]-shiftt,phit->qsample[i]*rgain[phit->ipmt]);// sum all the hits and take the maximum as prompt time; Time are shifted according to RF signal
		      //cout<<phit->tsample[i]<<"	"<<shiftt<<"	"<<phit->qsample[i]*rgain[phit->ipmt]<<endl;
		      pulse_height=max(pulse_height,phit->qsample[i]*rgain[phit->ipmt]);
		      pulse_integral+=phit->qsample[i]*rgain[phit->ipmt];
		      }
		      if(pulse_height-phit->qUnpeak>0.1)
		      {
			cout<<"wrong !!!"<<endl;
			cout<<pulse_height<<" 	"<<phit->qUnpeak<<endl;
		      }
		      /*if(pulse_integral/pulse_height>100) 
		      {
			cout<<"check event "<<ientry <<" pmt "<<phit->ipmt <<" time from "<<phit->tsample[phit->nsamples-1]<<" to "<<phit->tsample[0]<<endl;
			cout<<pulse_height<<" 	"<<pulse_integral<<endl;
		      }*/
		      ly+=pulse_height;
		      lyPE+=pulse_height/float (gain[phit->ipmt]);
		      
		      height_integral_ADC_h->Fill(pulse_height,pulse_integral);
		      height_integral_PE_h->Fill(pulse_height/float (gain[phit->ipmt]),pulse_integral);
		      ntHit->GetEntry(counter);
		      totalqUn+=qUnhit;
		      counter++;
		      
		    }


  
		    //sum_h->Draw();
		    //myc->Update();sleep(2);
		    //cout<<"haha"<<endl;
		promptt_rft=sum_h->GetMaximumBin()-rft-2100;
		//cout<<ientry<<" tprompt = "<<promptt_rft*4<<endl;
		if(totalqUn>threshold){
		lightyield_nenergy_ADC_h->Fill(nmass*(sqrt(1/(((promptt_rft*4+625)*c/L)*((promptt_rft*4+625)*c/L)-1)+1)-1),ly);  
		lightyield_nenergy_PE_h->Fill(nmass*(sqrt(1/(((promptt_rft*4+625)*c/L)*((promptt_rft*4+625)*c/L)-1)+1)-1),lyPE);
		promptt_rft_h->Fill(promptt_rft*4);
		cout<<ientry<<" tprompt ... = "<<promptt_rft*4<<endl;
		}
		double promptpulse_height=0;
		double promptpulse_height_PE=0;
		double promptpulse_integral=0;
		double delayedpulse_height=0;
		double delayedpulse_integral=0;
		for(int ihit =0; ihit < hit.size(); ++ihit) {
		      TPmtHit* phit = &(ev->hit[ihit]);

		      int shiftt;
		      if(phit->ipmt<7)shiftt=ev->rft21[0]-rft;
		      else if(phit->ipmt<14)shiftt=ev->rft22[0]-rft;
		      else if(phit->ipmt<21)shiftt=ev->rft23[0]-rft;
		      
		      double temppulse_height=0;
		      double temppulse_integral=0;
		      
		      for (int i=0;i<phit->nsamples;i++) {
			//sum_h->Fill(phit->tsample[i]-shiftt,phit->qsample[i]*rgain[phit->ipmt]);
			temppulse_height=max(temppulse_height,phit->qsample[i]*rgain[phit->ipmt]);
			temppulse_integral+=phit->qsample[i]*rgain[phit->ipmt];
		      }
		      //cout<<"*** "<<phit->tsample[0]-shiftt-rft<<"--"<<phit->tsample[phit->nsamples-1]-shiftt-rft<<" < "<<promptt_rft+2<<endl;
		      if( (phit->peakTime-shiftt-rft) <= promptt_rft+5 )//Unit tick(=4ns)
		      {
			promptpulse_height+=temppulse_height;
			promptpulse_height_PE+=temppulse_height/float (gain[phit->ipmt]);
			promptpulse_integral+=temppulse_integral;
		      }
		      else
		      {
			delayedpulse_height+=temppulse_height;
			delayedpulse_integral+=temppulse_integral;
		      }
	      light_timing_h->Fill((phit->peakTime-shiftt-sum_h->GetMaximumBin()+2100)*4);
		  }
		  //cout<<ev->event<<"	"<<promptt_rft<<"	"<<rft<<"	"<<promptpulse_height<<"	"<<delayedpulse_height<<"	"<<ly<<endl;
		if(delayedpulse_height>0)
		{
		  prompt_height_total_ratio_ADC_h->Fill(promptpulse_height/ly);
		  prompt_height_total_ratio_PE_h->Fill(promptpulse_height_PE/lyPE);
		}
		if(delayedpulse_integral>0)prompt_integral_total_ratio_h->Fill(promptpulse_integral/(promptpulse_integral+delayedpulse_integral));
		//if(totalqUn>threshold && (sum_h->GetMaximumBin()-rft-2100)*4<-800)cout<<ientry<<endl;
	      }      
	    else {
		counter+=hit.size();
	      }
	    delete sum_h;
	  }
        }
      }
    }//end of while loop
    
   
  TF1 *g1 = new TF1("g1","gaus",-640,-616);
  g1->SetLineColor(2);
  promptt_rft_h->Fit("g1","R");
  //promptt_rft_h->SetAxisRange(-700,-200);
  
  cout<<g1->GetParameter(1)<<endl;
  for (int i=1;i<=4200;i++)
  {
    if(i>1945)tof_h->SetBinContent(i-g1->GetParameter(1)/4.+L/c/4.,promptt_rft_h->GetBinContent(i));
  }

  for (int tick=2599;tick>=2119;tick--)
  {
     double height=tof_h->GetBinContent(tick); 
     //cout<<height<<endl;
     //En=nmass*(sqrt(1/((t*c/L)*(t*c/L)-1)+1)-1);
     nspectrum_h->SetBinContent(2600-tick,height);
  }
  //nspectrum_h->Scale(1./1000000);
   
  outputfile->cd();
  promptt_rft_h->Write();
  tof_h->Write();
  nspectrum_h->Write();
  height_integral_ADC_h->Write();
  height_integral_PE_h->Write();
  prompt_height_total_ratio_ADC_h->Write();
  prompt_height_total_ratio_PE_h->Write();
  prompt_integral_total_ratio_h->Write();
  lightyield_nenergy_ADC_h->Write();
  lightyield_nenergy_PE_h->Write();
  light_timing_h->Write();
  
  outputfile->Close();
  
  delete outputfile;
}


