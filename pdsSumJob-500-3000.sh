#!/bin/bash -l
## for batch,  qsub -l h_vmem=2G pdsSumJob-500-3000.sh 
export ROOTSYS=/project/projectdirs/captain/releases/LCGCMT/2.0.4/LCG_Settings/../EXTERNALS/ROOT/5.34.34/x86_64-linux-gcc44-opt
source /project/projectdirs/captain/releases/LCGCMT/2.0.4/EXTERNALS/ROOT/5.34.34/x86_64-linux-gcc44-opt/bin/thisroot.sh
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/global/homes/m/mgold/mgold/pdsAnalysis/obj/
/global/homes/m/mgold/mgold/pdsAnalysis/obj/pdsSum 500 3000 PDS_beamtime_files
