/**
** MG, Jan 12 2018
**/
#ifndef TPMTALIGN_DEFINED
#define TPMTALIGN_DEFINED
#include <iostream>
#include <string>
#include <map>
#include <TNamed.h>

using namespace std;

// class to store info for the data file labeld by tag 

class TPmtAlign: public TNamed {
	public:
    /* yujing 
       1614   2032-2035            B2 has 4 extra events.
       1646   3561-3562            B0 has 2 extra events.
       1714   2058-2059            B2 has 2 extra events.
       1742   2096                     B2 has 1 extra event.    
       1756   347-348                B0 misses 1 event between 347 & 348
       1853   4492-4497            B1 has 6 extra events than B0.
              4493-4496            B1 has 4 extra events than B2.
       1858   1307                     B1 has 1 extra event than B0.
              1309-1310            B1 has 2 extra events than B2.
       1910   4110                      B2 has 1 extra events.
       2020   1066-1067            B1 misses 1 event between 1066 & 1067
              2635-2636            B1 misses 1 event between 2635 & 2636
       2025   467                   B1 has 1 extra event than B0.
              471-472                B2 has 2 extra events than B0.
       2054   210-212               B1 misses 1 event between 210-212.
       2059   852                    B2 has 1 extra event
       */
    enum {NPMT=21};
    enum {NFIX=16};
    enum {MAXEVENT=5000};
    Int_t skipSegment[NFIX]={1614,1646,1714,1742,1756,1853,1853,1858,1858,1910,2020,2020,2025,2025,2054,2059};
    Int_t skipFirst[NFIX]={2032,3561,2058,2096,-347,4492,4493,1307,1309,4110,-1066,-2635,467,471,-210,852};
    Int_t skipLast[NFIX]= {2035,3562,2056,2096,-348,4497,4496,1307,1310,4110,-1067,-2636,467,472,-212,852};
    Int_t skipBoard[NFIX]={2,0,2,2,0,1,1,1,1,2,1,1,1,2,1,2};

    TPmtAlign();
		~TPmtAlign();
		void clear();
		void print();
    UInt_t event[3];
    
		ClassDef(TPmtAlign,5)
};
#endif

