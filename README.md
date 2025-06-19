This is a simulation of a dual readout crystal calorimeter (currently the code is work in progress).
See https://iopscience.iop.org/article/10.1088/1748-0221/15/11/P11005 for the concept.

## If you are not on alma9-like OS, but can use singularity
```
singularity run -B /cvmfs:/cvmfs -B /data:/data docker://gitlab-registry.cern.ch/sft/docker/alma9-core:latest
# at Baylor
# singularity run -B /cvmfs:/cvmfs -B /cms/data:/cms/data docker://gitlab-registry.cern.ch/sft/docker/alma9-core:latest
```

## All times:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh
```

## First time only:
```
# setup directory
mkdir stuff4stuff
cd stuff4stuff

# git clone, compile, install
git clone ssh://git@gitlab.cern.ch:7999/calvisionsimulation/DualTestBeam.git
# or try sarah's version
# git clone git@github.com:saraheno/DualTestBeam.git
cd DualTestBeam
mkdir build
mkdir install
cd build
cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DROOT_DIR=$ROOTSYS -D CMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -D DD4HEP_USE_EDM4HEP=ON ..
make -j4
make install
```

## All times
```
cd to stuff4stuff/DualTestBeam
source ./install/bin/thisDualTestBeam.sh
cd compact
```

## running in batch mode

Look in [massjobs.py](https://gitlab.cern.ch/calvisionsimulation/DualTestBeam/-/blob/master/compact/massjobs.py) to see how to run it
to analyze the output, look at [massjobs_s2.py](https://gitlab.cern.ch/calvisionsimulation/DualTestBeam/-/blob/master/compact/massjobs_s2.py).
For Baylor users, see [massjobs_pbsarray.py](https://gitlab.cern.ch/calvisionsimulation/DualTestBeam/-/blob/master/compact/massjobs_pbsarray.py)

or see examples in CI (continuous integration) yaml file for running `ddsim` and `Resolution.C` in
[.gitlab-ci.yml](https://gitlab.cern.ch/calvisionsimulation/DualTestBeam/-/blob/master/.gitlab-ci.yml)

## running interactively
Change `--runType=batc` above to `--runType=vis`.
Then
```
/control/execute vis.mac
/run/beamOn 1
```
On the window that pops up, choose “Miscellany” and “Exit to G4Vis >”
Then do typical GEANT4 visualization commands such as:
```
/vis/viewer/refresh
/vis/viewer/zoomTo 10
/vis/viewer/pan -100 200 cm
exit
```
