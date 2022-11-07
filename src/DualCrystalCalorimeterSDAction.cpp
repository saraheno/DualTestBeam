//==========================================================================
//  AIDA Detector description implementation 
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================

// Framework include files
#include "DualCrystalCalorimeterHit.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"









namespace CalVision {

  G4double fromEvToNm(G4double energy)
  {
    // there is some bug somewhere.  shouldn't need this factor
    return 1239.84187 / energy*1000.;
  }


  int SCECOUNT=0;

  class DualCrystalCalorimeterSD {
  public:
    typedef DualCrystalCalorimeterHit Hit;
    // If we need special data to personalize the action, be put it here
    //int mumDeposits = 0;
    //double integratedDeposit = 0;
  };
}

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {
  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim   {

    using namespace CalVision;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //               Geant4SensitiveAction<MyTrackerSD>
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /** \addtogroup Geant4SDActionPlugin
     *
     * @{
     * \package DualCrystalCalorimeterSDAction
     *
     * @}
     */

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DualCrystalCalorimeterSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<CalVision::DualCrystalCalorimeterSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DualCrystalCalorimeterSD>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {


      bool SCEPRINT=(SCECOUNT<10000);
      //if(SCEPRINT) std::cout<<"scecount is "<<SCECOUNT<<" print is "<<SCEPRINT<<std::endl;



      G4StepPoint *thePrePoint = step->GetPreStepPoint();
      G4StepPoint *thePostPoint = step->GetPostStepPoint();
      //      const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
      //const G4ThreeVector &thePostPosition = thePostPoint->GetPosition();
      G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
      G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
      G4String thePrePVName = "";
      if (thePrePV)
	thePrePVName = thePrePV->GetName();
      G4String thePostPVName = "";
      if (thePostPV)
	thePostPVName = thePostPV->GetName();
      //G4Track *theTrack = step->GetTrack();
      //G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

      //      if(thePrePVName.contains("slice")==0) {
      //std::cout<<"entering DualCrystalAction"<<std::endl;
      //  std::cout<<" prevolume is "<<thePrePVName<<std::endl;
      //  std::cout<<" postvolume is "<<thePostPVName<<std::endl;
      //  std::cout<<" pid is "<<TrPDGid<<std::endl;
	  //}


      Geant4StepHandler h(step);
      HitContribution contrib = DualCrystalCalorimeterHit::extractContribution(step);

      Geant4HitCollection*  coll    = collection(m_collectionID);
      VolumeID cell = 0;

      try {
        cell = cellID(step);
      } catch(std::runtime_error &e) {
	std::stringstream out;
        out << std::setprecision(20) << std::scientific;
        out << "ERROR: " << e.what()  << std::endl;
        out << "Position: "
            << "Pre (" << std::setw(24) << step->GetPreStepPoint()->GetPosition() << ") "
            << "Post (" << std::setw(24) << step->GetPostStepPoint()->GetPosition() << ") "
            << std::endl;
        out << "Momentum: "
            << " Pre (" <<std::setw(24) << step->GetPreStepPoint() ->GetMomentum()  << ") "
            << " Post (" <<std::setw(24) << step->GetPostStepPoint()->GetMomentum() << ") "
            << std::endl;

	std::cout << out.str();

        return true;
      }


      DualCrystalCalorimeterHit* hit = coll->findByKey<DualCrystalCalorimeterHit>(cell);
      if ( !hit ) {
        Geant4TouchableHandler handler(step);
	DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new DualCrystalCalorimeterHit(global);
        hit->cellID = cell;
        coll->add(cell, hit);
        printM2("CREATE hit with deposit:%e MeV  Pos:%8.2f %8.2f %8.2f  %s",
                contrib.deposit,pos.X,pos.Y,pos.Z,handler.path().c_str());
	std::cout<<"DRcalo deposit "<<contrib.deposit<<" position ("<<pos.X<<","<<pos.Y<<","<<pos.Z<<") string "<<handler.path().c_str()<<" volume id "<<cell<<std::endl;

        if ( 0 == hit->cellID )  { // for debugging only!
          hit->cellID = cellID(step);
          except("+++ Invalid CELL ID for hit!");
        }
      } else {
	//	std::cout<<"updating old hit"<<std::endl;
      }


      G4Track * track =  step->GetTrack();

      if(SCEPRINT) std::cout<< (track->GetDefinition())->GetParticleName()<<std::endl;

      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	if(SCEPRINT) std::cout<<"     in volume ID "<<cell<<std::endl;

	SCECOUNT+=1;
	//	if(SCEPRINT) std::cout<<"optical photon"<<std::endl;

	bool OptN = (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")||(track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys");

	//if(track->GetParentID()!=1) SCEPRINT=1;
	if( (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "CerenkovPhys")&&(track->GetCreatorProcess()->G4VProcess::GetProcessName() != "ScintillationPhys")  ) SCEPRINT=1;

	if(SCEPRINT) {
	  std::cout<<"     SCECOUNT="<<SCECOUNT<<std::endl;
	
	  std::cout<<"     will robinson have photon "<<track->GetCreatorProcess()->G4VProcess::GetProcessName() <<std::endl;
	  std::cout<<"     photon mother is "<<track->GetParentID()<<std::endl;
	  std::cout<<"     photon material is "<<(track->GetMaterial())->GetName()<<std::endl;
	  std::cout<<"     photon creator process is "<<(track->GetCreatorProcess())->GetProcessName()<<std::endl;
	  std::cout<<"     photon  process  type is "<<(track->GetCreatorProcess())->GetProcessType()<<std::endl;
	  std::cout<<"     photon sub process is "<<(track->GetCreatorProcess())->GetProcessSubType()<<std::endl;
	  std::cout<<"     photon current step number is "<<track->GetCurrentStepNumber()<<std::endl;
	//(track->GetCreatorProcess())->DumpInfo();
	  std::cout<<"     photon energy is "<<track->GetTotalEnergy()/eV<<std::endl;
	  std::cout<<"     photon wavelength is "<<fromEvToNm(track->GetTotalEnergy()/eV)<<std::endl;
	  std::cout<<"     number of cherenkov is  is "<<hit->ncerenkov<<std::endl;
	  std::cout<<"     number of scintillation is  is "<<hit->nscintillator<<std::endl;
	}

	float wavelength=fromEvToNm(track->GetTotalEnergy()/eV);
	int ibin=-1;
	float binsize=(hit->wavelenmax-hit->wavelenmin)/hit->nbin;
	ibin = (wavelength-hit->wavelenmin)/binsize;
	

	if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")  {
	  if(SCEPRINT) std::cout<<" found cerenkov photon"<<std::endl;
	  if(((track->GetMaterial())->GetName())=="killMedia") 
	    { 
	      hit->ncerenkov+=1;
	      if(ibin>-1&&ibin<hit->nbin) ((hit->ncerwave).at(ibin))+=1;
	      track->SetTrackStatus(fStopAndKill);}
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1)  ) hit->ncerenkov+=1;
	    if( (track->GetCurrentStepNumber()==1)  ) hit->ncerenkov+=1;
	  }
          //return false;
        } 
	else if (  track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys"  ) {
          if(SCEPRINT) std::cout<<"     scintillation photon"<<std::endl;
          if(((track->GetMaterial())->GetName())=="killMedia") 
	    {
	      hit->nscintillator+=1;
	      if((ibin>-1)&&(ibin<hit->nbin)) ((hit->nscintwave).at(ibin))+=1;
	      track->SetTrackStatus(fStopAndKill);}
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1) ) hit->nscintillator+=1; 
	    if( (track->GetCurrentStepNumber()==1) ) hit->nscintillator+=1; 
	  }

          //return false;
        }
	else {
          if(SCEPRINT) std::cout<<"      other photon"<<std::endl;
          //track->SetTrackStatus(fStopAndKill);
          //return false;
	}

      }

    else {   // particles other than optical photons
	
      //if(SCEPRINT) std::cout<<"NOT optical photon"<<std::endl;

        hit->energyDeposit += contrib.deposit;
        hit->truth.emplace_back(contrib);

        mark(h.track);
        //return true;
      }

	
      return true;

    }
	
  }
} // end namespace calvision







namespace dd4hep { namespace sim {

    using namespace CalVision;

    struct WavelengthMinimumCut : public dd4hep::sim::Geant4Filter  {
  /// Energy cut value
      double m_wavelengthCut;
    public:
  /// Constructor.
      WavelengthMinimumCut(dd4hep::sim::Geant4Context* c, const std::string& n);
  /// Standard destructor
      virtual ~WavelengthMinimumCut();
  /// Filter action. Return true if hits should be processed
      virtual bool operator()(const G4Step* step) const  override  final  {
	bool test=true;
	G4Track *theTrack = step->GetTrack();
	if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
	  float energy=theTrack->GetTotalEnergy()/eV;
	  float wave=fromEvToNm(energy);
	  if(wave < m_wavelengthCut) {
	    theTrack->SetTrackStatus(fStopAndKill);
	    test=false;}
	}
	return test;
      }
      virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
	return true;
      }
    };

  /// Constructor.
    WavelengthMinimumCut::WavelengthMinimumCut(Geant4Context* c, const std::string& n)
      : Geant4Filter(c,n) {
      InstanceCount::increment(this);
      declareProperty("Cut",m_wavelengthCut=0.0);
    }

  /// Standard destructor
    WavelengthMinimumCut::~WavelengthMinimumCut() {
      InstanceCount::decrement(this);
    }



    struct WavelengthnmwindCut : public dd4hep::sim::Geant4Filter  {
  /// Energy cut value
      double m_wavelengthstart;
    public:
  /// Constructor.
      WavelengthnmwindCut(dd4hep::sim::Geant4Context* c, const std::string& n);
  /// Standard destructor
      virtual ~WavelengthnmwindCut();
  /// Filter action. Return true if hits should be processed
      virtual bool operator()(const G4Step* step) const  override  final  {
	bool test=true;
	G4Track *theTrack = step->GetTrack();
	if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
	  float energy=theTrack->GetTotalEnergy()/eV;
	  float wave=fromEvToNm(energy);
	  if((wave<m_wavelengthstart) || (wave > m_wavelengthstart+0.5) ) {
	    theTrack->SetTrackStatus(fStopAndKill);
	    test=false;}
	}
	return test;
      }
      virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
	return true;
      }
    };

  /// Constructor.
    WavelengthnmwindCut::WavelengthnmwindCut(Geant4Context* c, const std::string& n)
      : Geant4Filter(c,n) {
      InstanceCount::increment(this);
      declareProperty("Cut",m_wavelengthstart=0.0);
    }

  /// Standard destructor
    WavelengthnmwindCut::~WavelengthnmwindCut() {
      InstanceCount::decrement(this);
    }


  }}  // end using namespace





//--- Factory declaration
namespace dd4hep { namespace sim {
    typedef Geant4SensitiveAction<DualCrystalCalorimeterSD> DualCrystalCalorimeterSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DualCrystalCalorimeterSDAction)
DECLARE_GEANT4ACTION(WavelengthMinimumCut)
DECLARE_GEANT4ACTION(WavelengthnmwindCut)
