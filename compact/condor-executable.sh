#!/bin/bash

cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/SingleDualCrystal/compact/

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


ddsim --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/SingleDualCrystal/compact/SCEPCALsteering.py --compact /home/eno/CalVision/dd4hep/DD4hep/examples/SingleDualCrystal/compact/DRSingleCrystal.xml --runType batch --part.userParticleHandler='' -G --gun.position="0.,10.,0." --gun.direction "0 -1 0" --gun.energy "1*GeV" --gun.particle="mu-" --gun.distribution=uniform -N 50 -O out.root >& haha.log

exitcode=$?


#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
