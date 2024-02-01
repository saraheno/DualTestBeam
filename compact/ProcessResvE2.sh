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


root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_15GeV_e-.root","./output/out_SampOnly_15GeV_pi-.root","./output/out_SampOnly_15GeV_pi-.root",15,0,1,1,1,3,"hists_SampOnly_3_15GeV.root","DRSNoSegment","DRSNoSegment")' > 15GeV.log


root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_25GeV_e-.root","./output/out_SampOnly_25GeV_pi-.root","./output/out_SampOnly_25GeV_pi-.root",25,0,1,1,1,3,"hists_SampOnly_3_25GeV.root","DRSNoSegment","DRSNoSegment")' >25GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_35GeV_e-.root","./output/out_SampOnly_35GeV_pi-.root","./output/out_SampOnly_35GeV_pi-.root",35,0,1,1,1,3,"hists_SampOnly_3_35GeV.root","DRSNoSegment","DRSNoSegment")' >35GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_45GeV_e-.root","./output/out_SampOnly_45GeV_pi-.root","./output/out_SampOnly_45GeV_pi-.root",45,0,1,1,1,3,"hists_SampOnly_3_45GeV.root","DRSNoSegment","DRSNoSegment")' >45GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_55GeV_e-.root","./output/out_SampOnly_55GeV_pi-.root","./output/out_SampOnly_55GeV_pi-.root",55,0,1,1,1,3,"hists_SampOnly_3_55GeV.root","DRSNoSegment","DRSNoSegment")' >55GeV.log






exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
