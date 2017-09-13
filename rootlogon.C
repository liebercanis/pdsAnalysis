{
  printf("\n this is rootlogon for pdsAnalysis \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  int iload = gSystem->Load("./obj/libPmtRoot.so");
  printf(" loaded libPmtRoot = %i zero is success! \n",iload);
  gSystem->AddIncludePath(" -I./obj ");
  //printf(" dynamic path %s \n\n",gSystem->GetDynamicPath());
  //printf(" include path %s \n\n",gSystem->GetIncludePath());
  //printf(" libaries %s \n\n",gSystem->GetLibraries());
  
  TStyle *C43Style= new TStyle("C43","C43 approved plots style");

  // use plain black on white colors
  C43Style->SetFrameBorderMode(0);
  C43Style->SetCanvasBorderMode(0);
  C43Style->SetPadBorderMode(0);
  C43Style->SetCanvasColor(0);
  C43Style->SetPadColor(0);
  C43Style->SetTitleFillColor(0);
  C43Style->SetStatColor(0);
  C43Style->SetFillColor(0);
  // set the paper & margin sizes
  C43Style->SetPadTopMargin(0.05);
  C43Style->SetPadRightMargin(0.05);
  C43Style->SetPadBottomMargin(0.15);
  C43Style->SetPadLeftMargin(0.1);
  // use large Times-Roman fonts
  C43Style->SetTitleX(0.25); //title X location
  C43Style->SetTitleY(0.95); //title Y location
  C43Style->SetTitleW(0.5); //title width
  C43Style->SetTitleH(0.1); //title height
  C43Style->SetStatX(0.9); //title X location
  C43Style->SetStatY(0.95); //title Y location
  C43Style->SetStatW(0.3); //title width
  C43Style->SetStatH(0.5); //title height
  C43Style->SetLabelSize(0.1,"x");
  C43Style->SetTitleSize(0.1,"x");
  C43Style->SetTitleOffset(0.8,"x");
  C43Style->SetLabelSize(0.1,"y");
  C43Style->SetTitleSize(0.1,"y");
  C43Style->SetTitleOffset(0.8,"y");
  C43Style->SetLegendBorderSize(0);
  C43Style->SetCanvasBorderMode(0);
  C43Style->SetCanvasBorderSize(0);
  C43Style->SetFrameBorderMode(0);
  C43Style->SetFrameBorderSize(0);
  C43Style->SetMarkerStyle(20);
  C43Style->SetHistLineWidth(2);
  
}

void browse() 
{
 TBrowser *browser = new TBrowser();
 gSystem->Load("libTreeViewer");
}

void dir() 
{
  gDirectory->ls("-m");
}


