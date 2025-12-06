To access the root file produced by the instructions:
```python
import ROOT
```
Load the dictionaries:
```python
ROOT.gSystem.Load("libDDG4.so")
ROOT.gSystem.Load("libDDG4Plugins.so")
ROOT.gSystem.Load("./lib/libDualTestBeam.so")
```
Open the file:
```python
file = ROOT.TFile("junk.root")
```
Explore the content:
```python
>>> file.ls()
TFile**		junk.root	dd4hep Simulation data
 TFile*		junk.root	dd4hep Simulation data
  KEY: TTree	EVENT;1	Geant4 EVENT information
```
Or make a RDataFrame:
```python
df = ROOT.RDataFrame("EVENT", "junk.root")
```

```python
>>> df.Describe()
Dataframe from TChain EVENT in file junk.root

Property                Value
--------                -----
Columns in total            3
Columns from defines        0
Event loops run             0
Processing slots            1

Column                  Type                                                    Origin
------                  ----                                                    ------
DRCNoSegment            ROOT::VecOps::RVec<CalVision::DualCrysCalorimeterHit*>  Dataset
EdgeDetNoSegment        ROOT::VecOps::RVec<CalVision::DualCrysCalorimeterHit*>  Dataset
MCParticles             ROOT::VecOps::RVec<dd4hep::sim::Geant4Particle*>        Dataset
```
Declare a C++ function:
```python
converter_cpp = """
... using namespace ROOT::VecOps;
... ROOT::RVecD getEnergyDeposit(const RVec<CalVision::DualCrysCalorimeterHit> &hit_collection)
... {
...   auto energy = [](const CalVision::DualCrysCalorimeterHit& hit){return hit.energyDeposit;};
...   return Map(hit_collection, energy);
... """
>>> ROOT.gInterpreter.Declare(converter_cpp)
```
