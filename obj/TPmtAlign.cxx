#include "TPmtAlign.hxx"
ClassImp(TPmtAlign)

TPmtAlign::TPmtAlign(): TNamed("TPmtAlign","TPmtAlign")
{
  clear();
}


TPmtAlign::~TPmtAlign(){}

void TPmtAlign::clear()
{
  tag.clear();
  trig.clear();
  align0.clear();
  align1.clear();
  align2.clear();
  rf0.clear();
  rf1.clear();
  rf2.clear();
  start0=0;
  start1=0;
  start2=0;
  eventList.clear();
  eventList.resize(3);
}

void TPmtAlign::print()
{
  printf(" \n\t AAAAAAAAAAA PMT Align AAAAAAAAAA run %s start times: board 0 %u board 1 %u board 2 %u triggers %u \n",tag.c_str(),start0,start1,start2,trig.size());
}


void TPmtAlign::makeEventLists(std::map<int,std::string> timeMap)
{
  // fill lists
  std::map<int,std::string>::iterator iter;
  Int_t count=0;
  cout << " timeMap size " << timeMap.size() << endl;
  for (iter=timeMap.begin(); iter!=timeMap.end(); ++iter) {
    string fname = string(iter->second);
    Int_t segment = atoi(  fname.substr(fname.find(".")-6,4).c_str());
    for(int i=0; i<NFIX; ++i) {
      if(segment == skipSegment[i]) { 
        for(unsigned ib = 0; ib<3 ; ++ib) { 
          if( skipBoard[i] == ib ) eventList[ib].push_back( eventList[ib].size() + skipLast[i] - skipFirst[i]) ;
          else  eventList[ib].push_back(eventList[ib].size());
        }
      } else {
        for(unsigned ib = 0; ib<3 ; ++ib) eventList[ib].push_back(eventList[ib].size()); 
      }
      cout << i << " count " << ++count << " "  << iter->first << " => " << iter->second << "  " << segment << endl;
    }
  }
};

