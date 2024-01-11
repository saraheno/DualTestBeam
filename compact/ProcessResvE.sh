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


root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_10GeV_e-.root","./output/out_FSCEPonly_10GeV_pi-.root","./output/out_FSCEPonly_10GeV_pi-.root",10,0,1,0,1,4,"hists_10GeV.root","DRFNoSegment","DRFNoSegment")' > 10GeV.log


root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_20GeV_e-.root","./output/out_FSCEPonly_20GeV_pi-.root","./output/out_FSCEPonly_20GeV_pi-.root",20,0,1,0,1,4,"hists_20GeV.root","DRFNoSegment","DRFNoSegment")' >20GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_30GeV_e-.root","./output/out_FSCEPonly_30GeV_pi-.root","./output/out_FSCEPonly_30GeV_pi-.root",30,0,1,0,1,4,"hists_30GeV.root","DRFNoSegment","DRFNoSegment")' >30GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_40GeV_e-.root","./output/out_FSCEPonly_40GeV_pi-.root","./output/out_FSCEPonly_40GeV_pi-.root",40,0,1,0,1,4,"hists_40GeV.root","DRFNoSegment","DRFNoSegment")' >40GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_50GeV_e-.root","./output/out_FSCEPonly_50GeV_pi-.root","./output/out_FSCEPonly_50GeV_pi-.root",50,0,1,0,1,4,"hists_50GeV.root","DRFNoSegment","DRFNoSegment")' >50GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_100GeV_e-.root","./output/out_FSCEPonly_100GeV_pi-.root","./output/out_FSCEPonly_100GeV_pi-.root",100,0,1,0,4,1,"hists_100GeV.root","DRFNoSegment","DRFNoSegment")' >100GeV.log





exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
