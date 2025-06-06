##  all times:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh
```

## First time only:
```
# setup directory
mkdir stuff4stuff
cd stuff4stuff
# git clone, compile, install
git clone git@github.com:saraheno/DualTestBeam.git
cd DualTestBeam
mkdir build
mkdir install
cd build
cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DROOT_DIR=$ROOTSYS -D CMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -D DD4HEP _USE_EDM4HEP=ON 
make -j4
make install
```

## all times
```
cd to stuff4stuff/DualTestBeam	
source ./install/bin/thisDualTestBeam.sh
cd compact
```

## running
```
look in massjobs.py to see how to run it
to analyze the output, look at massjobs_s2.py
```
