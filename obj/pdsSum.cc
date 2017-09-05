#include "TPdsSummary.hxx"
using namespace std;
int main(int argc, char **argv) {
  cout << " executing " << argv[0] << endl;
  for(int i=1; i<argc; ++i) printf(" %i %s ",i,argv[i]);
  TString dirName;
  if(argc>1) dirName = TString(argv[1]);
  else dirName = TString("PDS_beamtime_files");
  cout << " now creating TPdsSummary with directory " << dirName << endl;
  TPdsSummary *psum = new TPdsSummary(dirName);
  psum->printFiles();
  cout << " now run over files " <<  endl;
  psum->run();
}

