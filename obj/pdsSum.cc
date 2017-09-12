#include "TPdsSummary.hxx"
using namespace std;
int main(int argc, char **argv) {
  cout << " executing " << argv[0] << endl;
  //for(int i=1; i<argc; ++i) printf(" %i %s ",i,argv[i]);
  TString dirName;
  Int_t fmax=0;
  if(argc>1)  fmax = atoi(argv[1]);
  if(argc>2) dirName = TString(argv[2]);
  else dirName = TString("PDS_beamtime_files");
  cout << " \n now creating TPdsSummary with directory " << dirName << " will read " << fmax << " files " << endl;
  TPdsSummary *psum = new TPdsSummary(dirName);
  psum->printFiles();
  cout << " now run over files " <<  endl;
  psum->run(fmax);
}

