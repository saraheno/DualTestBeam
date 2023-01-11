#!/bin/bash

cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/

START_TIME=`/bin/date`
echo "started at $START_TIME"


# setup software environment at UMD
#
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh
echo "ran setup"
source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh
echo "ran thisdd4hep"
#
# run 
#


ddsim --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/SCEPCALsteering.py --compact /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/DRDualTestBeam.xml --runType batch --part.userParticleHandler='' -G --gun.position="0.,0.,-150." --gun.direction "0 0 1" --gun.energy "30*GeV" --gun.particle="pi-" --gun.distribution=uniform -N 50 -O out.root >& haha.log

exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
