This is a simulation of a dual readout crystal calorimeter (currently the code is work in progress).  See https://iopscience.iop.org/article/10.1088/1748-0221/15/11/P11005 for the concept.

## Build instructions

```
git clone ssh://git@gitlab.cern.ch:7999/calvisionsimulation/DualTestBeam.git 
cd DualTestBeam
```

On alma 8, run
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt/setup.sh
```
On alma 9, run
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh
```

(At Baylor University, also run `module load gl_fix`)


```
mkdir build install
cd build
cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DROOT_DIR=$ROOTSYS -D CMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -D DD4HEP_USE_EDM4HEP=ON ..
make -j4
make install
```

## Running instructions

Go to the DualTestBeam directory (the directory where this README lives).
Run the same `source /cvmfs...` line as before, then
```
source install/bin/thisDualTestBeam.sh
cd compact
ddsim --compactFile=DRConly.xml --runType=batch -G --steeringFile SCEPCALsteering.py --outputFile=junk.root --part.userParticleHandler= -G --gun.position="0. 0. -1*cm" --gun.direction "0. 0. 1." --gun.energy "20*GeV" --gun.particle="pi-" --outputFile=junk.root -v VERBOSE -N 1 >& crash.txt
```

## Development instructions

When you have made a change, go to `DualTestBeam/build` and run
```
make -j4
make install
```

## Other instructions

These instructions may not be updated.

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


  




