DualRead out calorimeter

Setup the code using instruction in: https://foswiki.web.cern.ch/Calvision/DualCrystalDD4hep#Setup_the_working_area 

To compare rifractive index values with the latest values in DRDualTestBeam.xml file (from Mekhala Single Crystal Cosmic Ray analysis):

`python RI_refrac.py -f1 DRSingleCrystalCosmicRay.xml -f2 DRDualTestBeam.xml`

To run condor jobs:

`python massjobs.py -g FSCEPonly -p e-`

change particle type to by replacing `-e` with `pi-`.

To run Resolution do interactively:

`root -l -b -q 'Resolution.C(nevt,"electron_rootfile.root", "pion_rootfile.root","hcalonly_rootfile.root",energy,doecal,dohcal,doedge,gendet,hcaltype,"output.root","ECALleaf","HCALleaf")'`

the parameters to run Resolution.C file are described on top of Resolution.C file.

Take the output file from last step and run ResvE.C file:

`root -l -b -q 'ResvE.C()’>resve.txt`
