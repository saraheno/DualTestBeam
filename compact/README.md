To display FSCEPonly geometry: 

```bash
ddsim --compactFile=EventDisplay_DRFSCEPonly.xml --runType=vis -G --steeringFile SCEPCALsteering.py --outputFile=out.root --part.userParticleHandler= -G --gun.position="0.,0.,0.0*cm" --gun.direction "0.,0.,1.0" --gun.energy "20 * GeV" --gun.particle="e-" -N 1
```

DualReadOut Calorimeter using crystal ECAL and fiber HCAL

To setup your environment use instruction on `https://foswiki.web.cern.ch/Calvision/DualCrystalDD4hep#Simulation_with_dual_readout_calorimetry`

xml files used in ddsim commands:

Crystal ECAL: DRConly.xml
Fiber HCAL: DRFSCEPonly.xml
Dual (crystal + fiber): DRDualTestBeam.xml

To check the rifractive indices in the xml file (e.g. DRDualTestBeam.xml) with the latest file: DRSingleCrystalCosmicRay.xml (from Mekhala):

```bash
python RI_refrac.py -f1 DRSingleCrystalCosmicRay.xml -f2 DRDualTestBeam.xml
```

To run ddsim with condor:

```bash
python massjobs.py -g <geometry> -p <particle>
bash massjobs.sh
```

select geometry from: DRConly.xml (crystal ecal), DRFSCEPonly.xml (fiber hcal), DRDualTestBeam.xml (crystal ecal + fiber hcal)
run the above code for electron (`e-`) and pion (pi-) gun

Run Resolution.C code using:

```bash
root -l -b -q 'Resolution.C(nevt,"electron_rootfile.root", "pion_rootfile.root","hcalonly_rootfile.root",energy,doecal,dohcal,doedge,gendet,hcaltype,"output.root","ECALleaf","HCALleaf")'
```
To run Resolution.C in condor:

```bash
python massjobs_res.py -g DualTestBeam -ho=1 -hc=1 -ec=1 -ed=1 -gd=3 -e=1
```

Options described at the top of Resolution.C code

To run ResvE.C:

```bash
root -l -b -q 'ResvE.C()'
```
