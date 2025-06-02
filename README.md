all times:

   source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt/setup.sh


First time only:

   mkdir stuff4stuff

   cd stuff4stuff
 
   git clone git@githubNOSPAM.com:saraheno/DualTestBeam.git

   cd DualTestBeam

   mkdir build

   mkdir install

   cd build

   cmake -DDD4HEP_USE_GEANT4=ON -DBoost_NO_BOOST_CMAKE=ON -DDD4HEP_USE_LCIO=ON -DROOT_DIR=$ROOTSYS -D CMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install -D DD4HEP _USE_EDM4HEP=ON ..

   make -j4

  make install


all times

   cd to stuff4stuff/DualTestBeam	

   source ./install/bin/thisDualTestBeam.sh

   cd compact

ddsim --compactFile=DRConly.xml --runType=batch -G --steeringFile SCEPCALsteering.py --outputFile=junk.root --part.userParticleHandler= -G --gun.position="0. 0. -1*cm" --gun.direction "0. 0. 1." --gun.energy "20*GeV" --gun.particle="pi-" --outputFile=junk.root -v VERBOSE -N 1 >& crash.txt