#!/bin/bash

cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/

START_TIME=`/bin/date`
echo "started at $START_TIME"
 echo "started at $START_TIME on ${HOSTNAME}"

# setup software environment at UMD
#
source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh
echo "ran setup"
source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh
echo "ran thisdd4hep"
#
# run 
#


root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_10GeV_e-.root","./output/out_DualTestBeam_10GeV_pi-.root","./output/out_DRFSCEPonly_10GeV_e-.root",10,0,1,1,1,3,"hists_DualTestBeam_3_10GeV.root","DRSNoSegment","DRSNoSegment")' >& 10GeV.log


root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_20GeV_e-.root","./output/out_DualTestBeam_20GeV_pi-.root","./output/out_DRFSCEPonly_20GeV_e-.root",20,0,1,1,1,3,"hists_DualTestBeam_3_20GeV.root","DRSNoSegment","DRSNoSegment")' >& 20GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_30GeV_e-.root","./output/out_DualTestBeam_30GeV_pi-.root","./output/out_DRFSCEPonly_30GeV_e-.root",30,0,1,1,1,3,"hists_DualTestBeam_3_30GeV.root","DRSNoSegment","DRSNoSegment")' >& 30GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_40GeV_e-.root","./output/out_DualTestBeam_40GeV_pi-.root","./output/out_DRFSCEPonly_40GeV_e-.root",40,0,1,1,1,3,"hists_DualTestBeam_3_40GeV.root","DRSNoSegment","DRSNoSegment")' >& 40GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_50GeV_e-.root","./output/out_DualTestBeam_50GeV_pi-.root","./output/out_DRFSCEPonly_50GeV_e-.root",50,0,1,1,1,3,"hists_DualTestBeam_3_50GeV.root","DRSNoSegment","DRSNoSegment")' >& 50GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_DualTestBeam_100GeV_e-.root","./output/out_DualTestBeam_100GeV_pi-.root","./output/out_DRFSCEPonly_100GeV_e-.root",100,0,1,1,1,3,"hists_DualTestBeam_3_100GeV.root","DRSNoSegment","DRSNoSegment")' >& 100GeV.log





exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
