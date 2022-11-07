This is a simulation of a dual readout crystal calorimeter (currently the code is work in progress).  See https://iopscience.iop.org/article/10.1088/1748-0221/15/11/P11005 for the concept.

Some instructions on how to use:


Each subsequent time
$ cd dd4hep/DD4hep
$ source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

$ source bin/thisdd4hep.sh

To make changes
When make change, go to DD4hep/build and do “make make install”
To visualize the geometry
cd examples/SingleDualCrystal/compact
geoDisplay DRSingleCrystal.xml

to run
cd examples/SingleDualCrystal/compact
ddsim --steeringFile SCEPCALsteering.py --compact ./DRSingleCrystal.xml --runType batch --part.userParticleHandler='' -G --gun.position="0.,10.,0." --gun.direction "0 -1 0" --gun.energy "1*GeV" --gun.particle="mu-" --gun.distribution=uniform -N 1 -O out.root

to run interactively
DOES NOT WORK
ddsim --compactFile=./DRSingleCrystal.xml --runType=vis -G --steeringFile SCEPCALsteering.py --outputFile=testSCEPCAL.root --part.userParticleHandler='' -G --gun.position="0.,10.,0." --gun.direction "0 -1 0" --gun.energy "1*GeV" --gun.particle="mu-" --gun.distribution=uniform

/control/execute vis.mac

/run/beamOn 1

On the window that pops up, choose “Miscellany” and “Exit to G4Vis >”

Then do typical GEANT4 visualization commands like:

/vis/viewer/refresh

/vis/viewer/zoomTo 10

/vis/viewer/pan -100 200 cm

exit


  




