# captain pdsAnalysis
Author: M. Gold
July 24, 2017


code is on GitHub, 

https://github.com/liebercanis/pdsAnalysis

and on pds
/project/projectdirs/captain/users/mgold/pdsAnalysis/

instructions

add to your ~/.bash_profile.ext
export ROOTSYS=/project/projectdirs/captain/releases/LCGCMT/2.0.4/LCG_Settings/../EXTERNALS/ROOT/5.34.34/x86_64-linux-gcc44-opt
source $ROOTSYS/bin/thisroot.sh

in directory pdsAnalysis add soft link
ln -s /global/homes/m/mgold/data/2017/PDS_beamtime_alternate_runs pdsData

build the library:
	cd obj; make; cd ../

pdsf6 $ runAna
 <tag> (e.g. 07-12-1900_0) <max events> 

for example, test of 100 events

runAna 07-22-1429_0 100

to submit batch job, edit bsub and execute it.

for now, the directory /global/homes/m/mgold/mgold/pdsAnalysis/pdsData is explicitly written in pmtAna.cc code.

here is list of commands I did to get started.

   >git clone https://github.com/liebercanis/pdsAnalysis.git
   >cd pdsAnalysis/obj
   >make
   >cd ..
   // soft links or mkdir to create directories pdsData for input data and pdsOutput for root outputfiles 
   >ln -s  ~/captain/pds/pdsAnalysis/pdsData/
   >ln -s  ~/captain/pds/pdsAnalysis/pdsObj
   >root pmtAna.cc+

