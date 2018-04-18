{
  printf("\n this is rootlogon for pdsAnalysis/macro \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  int iload = gSystem->Load("../obj/libPmtRoot.so");
  printf(" loaded libPmtRoot = %i zero is success! \n",iload);
  gSystem->AddIncludePath(" -I../obj ");
  //printf(" dynamic path %s \n\n",gSystem->GetDynamicPath());
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  //printf(" libaries %s \n\n",gSystem->GetLibraries());
  gROOT->LoadMacro("util.C");
  
}


