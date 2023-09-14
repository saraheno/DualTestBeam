#!/bin/bash

cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/

START_TIME=`/bin/date`
echo "started at $START_TIME"


# setup software environment at UMD
#
source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh
echo "ran setup"
source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh
echo "ran thisdd4hep"
#
# run 
#






 ddsim --compactFile=/home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/DRFonly.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/SCEPCALsteering.py --outputFile=./output/out_fonly_20GeV_e-_100.root --part.userParticleHandler='' -G --gun.position="0.05,0.,-210*cm" --gun.direction "0 0 1" --gun.energy "20*GeV" --gun.particle="e-" -N 1000 >& ./output/sce_e_fonly_20.log





exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
