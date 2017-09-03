#include "TPdsSummary.hxx"
using namespace std;
int main(int argc, char **argv) {
  cout << " executing " << argv[0] << endl;
  for(int i=1; i<argc; ++i) printf(" %i %s ",i,argv[i]);
  if(argc>1) new TPdsSummary(argv[0]);
  else new TPdsSummary();
}

