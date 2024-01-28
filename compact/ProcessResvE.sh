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


root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_10GeV_e-.root","./output/out_SampOnly_10GeV_pi-.root","./output/out_SampOnly_10GeV_pi-.root",10,0,1,1,1,3,"hists_10GeV.root","DRSNoSegment","DRSNoSegment")' > 10GeV.log


root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_20GeV_e-.root","./output/out_SampOnly_20GeV_pi-.root","./output/out_SampOnly_20GeV_pi-.root",20,0,1,1,1,3,"hists_20GeV.root","DRSNoSegment","DRSNoSegment")' >20GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_30GeV_e-.root","./output/out_SampOnly_30GeV_pi-.root","./output/out_SampOnly_30GeV_pi-.root",30,0,1,1,1,3,"hists_30GeV.root","DRSNoSegment","DRSNoSegment")' >30GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_40GeV_e-.root","./output/out_SampOnly_40GeV_pi-.root","./output/out_SampOnly_40GeV_pi-.root",40,0,1,1,1,3,"hists_40GeV.root","DRSNoSegment","DRSNoSegment")' >40GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_50GeV_e-.root","./output/out_SampOnly_50GeV_pi-.root","./output/out_SampOnly_50GeV_pi-.root",50,0,1,1,1,3,"hists_50GeV.root","DRSNoSegment","DRSNoSegment")' >50GeV.log

root -b -l -q 'Resolution.C(500,"./output/out_SampOnly_100GeV_e-.root","./output/out_SampOnly_100GeV_pi-.root","./output/out_SampOnly_100GeV_pi-.root",100,1,1,0,3,1,"hists_100GeV.root","DRSNoSegment","DRSNoSegment")' >100GeV.log





exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
