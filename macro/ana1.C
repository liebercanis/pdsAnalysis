#include <vector>
const Double_t GAMMAPEAK=-628.089;//ysun 
const double L=23.2;//m
const double clight=0.299792458;//m/ns
const double nmass=939.565;//MeV

enum {MAXSAMPLES=2100};
enum {NB=3,NCPMT=7,NC=NCPMT+1,NS=MAXSAMPLES};
enum {NPMT=NB*NCPMT};
  TH1D *hTimeToRF;
  TH1D *hOcc;
  TH1D *hNoise;
  TH1D *hBase;
  TH1D* hSamples[NPMT];
  TH1D* hPeaks[NPMT];
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hRawQ[NPMT];
  TH1D* hNHits[NPMT];
  TH1D* hQMax[NPMT];
  TH1D* hCounts[NPMT];
  TH1D* hBaseline[NPMT];

  TH1D* hQinHit; 
  TH1D* hQHitBeam[NPMT];
  TH1D* hQHitNoBeam[NPMT];

  TH1D *hTPrompt;
  TH1D *hTPromptEvent;
  TH1D *tof_h;
  TH1D *nspectrum_h;

  //TH1D* hQHitTime[NPMT];
  TPmtEvent *pmtEvent;
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





void ana1(Long64_t max=0)
{
  TString fileTag("pmtAnaLow_07-31-1728_0");
  TString inputFileName = TString("../pdsOutput/")+fileTag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  TTree *pmtTree=NULL;

  // tree has to be in file
  Long64_t aSize=0;
  pmtTree = (TTree*) infile->Get("pmtTree");
  if(pmtTree) aSize=pmtTree->GetEntriesFast();
  printf(" pmtTree with %i entries \n",int(aSize));

  Long64_t maxLoop = aSize;
  if(max>0) maxLoop = max;

  if(aSize==0) return;

  // open ouput file and make some histograms
  TString outputFileName = TString("ana-")+fileTag+TString(".root");
  TFile *outfile = new TFile(outputFileName,"recreate");
  printf(" opening output file %s \n",outputFileName.Data());
  hTimeToRF = new TH1D("TimeToRF"," hit time to RF time",MAXSAMPLES/2,0,MAXSAMPLES+1);
  TNtuple* ntAnaHit = new TNtuple("ntAnaHit", " ana hits ","ipmt:rft:peakt:length:ratio:q:qp");
  TNtuple* ntAnaEv = new TNtuple("ntAnaEv"," event info ","ev:dt1:dt2:dt3:prompt:rft:tof:ke:qsum0:qsum1:qsum2");

  TString hname;
  TString htitle;

  for(Int_t ib=0; ib<NB; ++ib) {
    for(UInt_t ic=0; ic<NC; ++ic) {
      int ipmt = toPmtNumber(ib,ic);
      if(ipmt<0||ipmt>=NPMT) continue;
      
      hname.Form("QhitNoBeam_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit between RF board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitNoBeam[ipmt] = new TH1D(hname,htitle,200,0,200);
      hQHitNoBeam[ipmt]->SetXTitle(" ADC counts ");
      hQHitNoBeam[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

      hname.Form("QhitBeam_b%u_ch%u_pmt%u",ib,ic,ipmt);
      htitle.Form("Qhit beam trig board%u channel%u pmt%u",ib,ic,ipmt);
      hQHitBeam[ipmt] = new TH1D(hname,htitle,200,0,200);
      hQHitBeam[ipmt]->SetXTitle(" ADC counts ");
      hQHitBeam[ipmt]->SetYTitle(Form(" # hits %i ",ipmt));

     }
  }
  


  hQinHit = new TH1D("QinHit"," ADC counts of bins in hit ",2000,0,20);
  outfile->ls();
    
  TH1F *sum_h = new TH1F("sum_h","sum",2100,0,2100);// temperarily sum all the peaks and take the maximum as prompt time; eventually will use summed channel
  TH1F *qun_h = new TH1F("qun_h","sum",3000,0,30000);
  
  // for neutron spectrum
  tof_h = new TH1D("tof_h","TOF (Assuming L = 23.2 m);Time (ns);Number of entries",2*MAXSAMPLES,-4*MAXSAMPLES,4*MAXSAMPLES); 
  nspectrum_h = new TH1D("neutron_spectrum_h",";Neutron E_{Kin} (MeV);Frac of triggers",400,0,4000);
  
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent",&pmtEvent);
  std::vector<TPmtHit> hit;

  Double_t qboard[NB];

  for(unsigned entry =0; entry < maxLoop; ++entry ) {
    pmtTree->GetEntry(entry);
    Double_t tpromptToRF = pmtEvent->tPromptToRF;

    for(Int_t ib=0; ib<NB; ++ib) qboard[ib]=0;

    Double_t promptt=pmtEvent->tPromptToRF*4.0;//in ns
    Double_t tof=pmtEvent->tPromptToRF*4.0-GAMMAPEAK+L/clight;//in ns 
    Double_t ke=-9999;
    Double_t tpromptNs = -9999;
    if(pmtEvent->tPrompt!=-9999 && pmtEvent->trigType==TPmtEvent::TRIG111 && tof>0){
      tpromptNs = pmtEvent->tPrompt*4.0;
      double beta = L/tof/clight;
      double gamma2 = 1./(1.-beta*beta);
      if(gamma2>0) ke = nmass*(sqrt(gamma2)-1.0);
    }//ysun


       
    hit.clear();
    hit = pmtEvent->hit;
    /*
       if( pmtEvent->compSec >= 1501543798 && pmtEvent->compSec <= 1501543799) 
       printf("...entry  %i sec %i nano %ll %i hits %lu trig type %i (%lu,%lu,%lu) \n",
       entry,pmtEvent->compSec, pmtEvent->compNano, hit.size(), pmtEvent->trigType,
       pmtEvent->rft21.size(),pmtEvent->rft23.size(),pmtEvent->rft23.size())
    */

    for(int ihit =0; ihit < hit.size(); ++ihit) {
      TPmtHit* phit = &(pmtEvent->hit[ihit]);
      for(int is=0; is<phit->nsamples; ++is) hQinHit->Fill(phit->qsample[is]);
      Int_t timeToRF = phit->timeToRF; 
      if(pmtEvent->trigType==TPmtEvent::TRIG111||pmtEvent->trigType==TPmtEvent::TRIG444||pmtEvent->trigType==TPmtEvent::TRIG555) 
        hTimeToRF->Fill(timeToRF);
      //if(timeToRF<MAXSAMPLES) printf(" trig %i hit %i pmt %i time %i \n",pmtEvent->trigType, ihit, phit->ipmt,timeToRF);
      int length = TMath::Abs(phit->tstop-phit->tstart)+1;
      //printf(" \t %i %i ipmt %i length %i qhit %f \n",entry,ihit,phit->ipmt,length,phit->qhit);
      // cut on time since RF pulse
   
      if(pmtEvent->trigType==TPmtEvent::TRIG000)  hQHitNoBeam[phit->ipmt]->Fill(phit->qhit); 
      else  hQHitBeam[phit->ipmt]->Fill(phit->qhit); 
 
      ntAnaHit->Fill(phit->ipmt,phit->timeToRF, phit->peakTime, length, phit->ratio, phit->qhit, phit->qpeak);
      
      Int_t iboard, ichan;
      fromPmtNumber(phit->ipmt,iboard,ichan);
      qboard[iboard] += phit->qhit;

    }
    ntAnaEv->Fill(double(entry),pmtEvent->dtime[0],pmtEvent->dtime[1],pmtEvent->dtime[2],pmtEvent->tPrompt,pmtEvent->tRFave,tof,ke,qboard[0],qboard[1],qboard[2]);
  }

  // end of ana 
  //newCanPlots(tag);
  outfile->Write();

}
