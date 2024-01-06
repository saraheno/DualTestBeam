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


root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_15GeV_e-.root","./output/out_FSCEPonly_15GeV_pi-.root","./output/out_FSCEPonly_15GeV_pi-.root",15,0,1,0,1,1,"hists_10GeV.root","DRFNoSegment","DRFNoSegment")' > 15GeV.log


root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_25GeV_e-.root","./output/out_FSCEPonly_25GeV_pi-.root","./output/out_FSCEPonly_25GeV_pi-.root",25,0,1,0,1,1,"hists_25GeV.root","DRFNoSegment","DRFNoSegment")' >25GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_35GeV_e-.root","./output/out_FSCEPonly_35GeV_pi-.root","./output/out_FSCEPonly_35GeV_pi-.root",35,0,1,0,1,1,"hists_35GeV.root","DRFNoSegment","DRFNoSegment")' >35GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_45GeV_e-.root","./output/out_FSCEPonly_45GeV_pi-.root","./output/out_FSCEPonly_45GeV_pi-.root",45,0,1,0,1,1,"hists_45GeV.root","DRFNoSegment","DRFNoSegment")' >45GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_FSCEPonly_55GeV_e-.root","./output/out_FSCEPonly_55GeV_pi-.root","./output/out_FSCEPonly_55GeV_pi-.root",55,0,1,0,1,1,"hists_55GeV.root","DRFNoSegment","DRFNoSegment")' >55GeV.log






exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
