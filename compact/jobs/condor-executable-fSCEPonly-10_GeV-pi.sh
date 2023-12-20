#!/bin/bash
cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/
START_TIME=`/bin/date`
echo "started at $START_TIME"
echo "started at $START_TIME on ${HOSTNAME}"
source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh
echo "ran setup"
source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh
echo "ran thisdd4hep"
ddsim --compactFile=/home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/DRFSCEPonly.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/SCEPCALsteering.py --outputFile=./output/out_fSCEPonly_20GeV_pi-_100.root --part.userParticleHandler= -G --gun.position="0.,0.,-210*cm" --gun.direction "0 0.05 0.99875" --gun.energy "10*GeV" --gun.particle="pi-" -N 100 >& ./output/sce_pi_fSCEPonly_10.log
exitcode=$?
echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
