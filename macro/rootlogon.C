{
  printf("\n this is rootlogon for pdsAnalysis \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  int iload = gSystem->Load("../obj/libPmtRoot.so");
  printf(" loaded libPmtRoot = %i zero is success! \n",iload);
  gSystem->AddIncludePath(" -I../obj ");
  //printf(" dynamic path %s \n\n",gSystem->GetDynamicPath());
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  //printf(" libaries %s \n\n",gSystem->GetLibraries());
}

void browse() 
{
 TBrowser* browser = new TBrowser();
 gSystem->Load("libTreeViewer");
}

void dir() 
{
  gDirectory->ls("-m");
}



