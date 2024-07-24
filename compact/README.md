DualReadOut Calorimeter using crystal ECAL and fiber HCAL

To setup your environment use instruction on `https://foswiki.web.cern.ch/Calvision/DualCrystalDD4hep#Simulation_with_dual_readout_calorimetry`

To check the rifractive indices in the xml file (e.g. DRDualTestBeam.xml) with the latest file: DRSingleCrystalCosmicRay.xml (from Mekhala):

`python RI_refrac.py -f1 DRSingleCrystalCosmicRay.xml -f2 DRDualTestBeam.xml`

To run ddsim with condor:

`python massjobs.py -g <geometry> -p <particle>`

`bash massjobs.sh`

select geometry from: DRConly.xml (crystal ecal), DRFSCEPonly.xml (fiber hcal), DRDualTestBeam.xml (crystal ecal + fiber hcal)
run the above code for electron (`e-`) and pion (pi-) gun

Run Resolution.C code using:

`root -l -b -q 'Resolution.C(nevt,"electron_rootfile.root", "pion_rootfile.root","hcalonly_rootfile.root",energy,doecal,dohcal,doedge,gendet,hcaltype,"output.root","ECALleaf","HCALleaf")'`

Options described on the top of Resolution.C code

To run ResvE.C:

`root -l -b -q 'ResvE.C()'`

