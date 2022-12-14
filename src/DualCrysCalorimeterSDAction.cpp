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
#include <CLHEP/Units/PhysicalConstants.h>
#include "DualCrysCalorimeterHit.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "DD4hep/InstanceCount.h"
#include "DDG4/Geant4Random.h"


// way too slow if track all photons for now
// so randomly delete photons after creation according to this fraction
double dialCher= 0.00001;
double dialScint=1.;
int printlimitSCE=100;

// this doesn't work.  it is in mks
//double conversioneVnm=2.*CLHEP::pi*CLHEP::hbarc;



namespace CalVision {

  G4double fromEvToNm(G4double energy)
  {
    // there is some bug somewhere.  shouldn't need this facto
    //return  conversioneVnm/ energy*1000.;
    return 1239.84187 / energy*1000.;

  }


  int SCECOUNT=0;
  int SCECOUNT2=0;

  class DualCrysCalorimeterSD {
  public:
    typedef DualCrysCalorimeterHit Hit;
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
     * \package DualCrysCalorimeterSDAction
     *
     * @}
     */

    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DualCrysCalorimeterSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<CalVision::DualCrysCalorimeterSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DualCrysCalorimeterSD>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {


      if(SCECOUNT==1) {
	std::cout<<"DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!"<<std::endl;
	std::cout<<"dialCher is "<<dialCher<<std::endl;
	std::cout<<"dialScint is "<<dialScint<<std::endl;
	std::cout<<" you need to use this to interpret your results"<<std::endl;
      }

      bool SCEPRINT=(SCECOUNT<printlimitSCE);
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
      //std::cout<<"entering DualCrysAction"<<std::endl;
      //  std::cout<<" prevolume is "<<thePrePVName<<std::endl;
      //  std::cout<<" postvolume is "<<thePostPVName<<std::endl;
      //  std::cout<<" pid is "<<TrPDGid<<std::endl;
	  //}


      Geant4StepHandler h(step);
      HitContribution contrib = DualCrysCalorimeterHit::extractContribution(step);

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


      DualCrysCalorimeterHit* hit = coll->findByKey<DualCrysCalorimeterHit>(cell);
      if ( !hit ) {
        Geant4TouchableHandler handler(step);
	DDSegmentation::Vector3D pos = m_segmentation.position(cell);
        Position global = h.localToGlobal(pos);
        hit = new DualCrysCalorimeterHit(global);
        hit->cellID = cell;
        coll->add(cell, hit);
        printM2("CREATE hit with deposit:%e MeV  Pos:%8.2f %8.2f %8.2f  %s",
                contrib.deposit,pos.X,pos.Y,pos.Z,handler.path().c_str());
	if(SCECOUNT2<printlimitSCE) {
	  std::cout<<"DRcalo deposit "<<contrib.deposit<<" position ("<<pos.X<<","<<pos.Y<<","<<pos.Z<<") string "<<handler.path().c_str()<<" volume id "<<cell<<std::endl;
	  SCECOUNT2+=1;
	}

        if ( 0 == hit->cellID )  { // for debugging only!
          hit->cellID = cellID(step);
          except("+++ Invalid CELL ID for hit!");
        }
      } else {
	//	std::cout<<"updating old hit"<<std::endl;
      }



      G4Track * track =  step->GetTrack();
      std::string amedia = ((track->GetMaterial())->GetName());
      if(SCEPRINT) std::cout<< (track->GetDefinition())->GetParticleName()<<std::endl;

      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	if(SCEPRINT) std::cout<<"     in volume ID "<<cell<<std::endl;

	SCECOUNT+=1;
	//	if(SCEPRINT) std::cout<<"optical photon"<<std::endl;

	bool OptN = (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")||(track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys");

	//if(track->GetParentID()!=1) SCEPRINT=1;
	if( (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "CerenkovPhys")&&(track->GetCreatorProcess()->G4VProcess::GetProcessName() != "ScintillationPhys")  ) SCEPRINT=1;  // print if less than scecount (set earlier) or if photon has weird source


	float wavelength=fromEvToNm(track->GetTotalEnergy()/eV);
	int ibin=-1;
	float binsize=(hit->wavelenmax-hit->wavelenmin)/hit->nwlbin;
	ibin = (wavelength-hit->wavelenmin)/binsize;
	int phstep = track->GetCurrentStepNumber();
	
	
	if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")  {
	  if(SCEPRINT) std::cout<<" found Cerenkov photon"<<std::endl;

	  if(amedia.find("kill")!=std::string::npos) 
	    { 
	      if(SCEPRINT) std::cout<<"killing photon"<<std::endl;
	      //	      SCEPRINT=1;
	      if(phstep>1) {  // don't count photons created in kill media
		hit->ncerenkov+=1;
		if(ibin>-1&&ibin<hit->nwlbin) ((hit->ncerwave).at(ibin))+=1;
	      }
	      track->SetTrackStatus(fStopAndKill);
	    }
	  else if(amedia.find("BlackHole")!=std::string::npos) {
	      if(phstep>1) {  // don't count photons created in kill media
		hit->ncerenkov+=1;
		if(ibin>-1&&ibin<hit->nwlbin) ((hit->ncerwave).at(ibin))+=1;
	      }
	      track->SetTrackStatus(fStopAndKill);
	  }
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1)  ) hit->ncerenkov+=1;
	    if( (phstep==1)  ) {
	      hit->ncerenkov+=1;
	      Geant4Event&  evt = context()->event();
	      dd4hep::sim::Geant4Random& rnd = evt.random();
	      if(rnd.rndm()>dialCher) track->SetTrackStatus(fStopAndKill);
	    }
	  }
	}
	else if (  track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys"  ) {
          if(SCEPRINT) std::cout<<"     scintillation photon"<<std::endl;
	  std::string amedia = ((track->GetMaterial())->GetName());
	  if(amedia.find("kill")!=std::string::npos) 
	    //          if(((track->GetMaterial())->GetName())=="killMedia") 
	    {
	      if(SCEPRINT) std::cout<<"killing photon"<<std::endl;
	      if(phstep>1) {
		hit->nscintillator+=1;
		if((ibin>-1)&&(ibin<hit->nwlbin)) ((hit->nscintwave).at(ibin))+=1;
	      }
	      track->SetTrackStatus(fStopAndKill);}
	  else if(amedia.find("BlackHole")!=std::string::npos) {
	      if(phstep>1) {  // don't count photons created in kill media
		hit->nscintillator+=1;
		if(ibin>-1&&ibin<hit->nwlbin) ((hit->nscintwave).at(ibin))+=1;
	      }
	      track->SetTrackStatus(fStopAndKill);
	  }
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1) ) hit->nscintillator+=1; 
	    if( (phstep==1) ) {
	      hit->nscintillator+=1; 
	      Geant4Event&  evt = context()->event();
	      dd4hep::sim::Geant4Random& rnd = evt.random();
	      if(rnd.rndm()>dialScint) track->SetTrackStatus(fStopAndKill);	    
	    }
	  }

          //return false;
        }
	else {
          if(SCEPRINT) std::cout<<"      other photon"<<std::endl;
          //track->SetTrackStatus(fStopAndKill);
          //return false;
	}

	if(SCEPRINT) {
	  std::cout<<"     SCECOUNT="<<SCECOUNT<<std::endl;
	
	  std::cout<<"     will robinson have photon "<<track->GetCreatorProcess()->G4VProcess::GetProcessName() <<std::endl;
	  std::cout<<"     photon mother is "<<track->GetParentID()<<std::endl;
	  std::cout<<"     photon material is "<<(track->GetMaterial())->GetName()<<std::endl;
	  std::cout<<"     photon creator process is "<<(track->GetCreatorProcess())->GetProcessName()<<std::endl;
	  std::cout<<"     photon  process  type is "<<(track->GetCreatorProcess())->GetProcessType()<<std::endl;
	  std::cout<<"     photon sub process is "<<(track->GetCreatorProcess())->GetProcessSubType()<<std::endl;
	  std::cout<<"     photon current step number is "<<track->GetCurrentStepNumber()<<std::endl;
	  std::cout<<"     the pre volume name is "<<thePrePVName<<std::endl;
	  std::cout<<"     the post volume name is "<<thePostPVName<<std::endl;
	//(track->GetCreatorProcess())->DumpInfo();
	  std::cout<<"     photon energy is "<<track->GetTotalEnergy()/eV<<std::endl;
	  std::cout<<"     photon wavelength is "<<fromEvToNm(track->GetTotalEnergy()/eV)<<std::endl;
	  std::cout<<"     number of cherenkov is  is "<<hit->ncerenkov<<std::endl;
	  std::cout<<"     number of scintillation is  is "<<hit->nscintillator<<std::endl;
	}



      }
      else {   // particles other than optical photons
	
      //if(SCEPRINT) std::cout<<"NOT optical photon"<<std::endl;

	if(amedia.find("BlackHole")!=std::string::npos) {
	  hit->energyDeposit += track->GetKineticEnergy();
	  track->SetTrackStatus(fStopAndKill);
	} else {
      //add information about each contribution to the hit
      hit->truth.emplace_back(contrib);

	  hit->energyDeposit += contrib.deposit;
	  //	  hit->contribBeta.emplace_back(track->GetVelocity()/CLHEP::c_light*10000.);
	  hit->contribBeta.emplace_back(track->GetVelocity()/CLHEP::c_light);
	  hit->contribCharge.emplace_back((track->GetParticleDefinition())->GetPDGCharge());
	  //int tsize = (hit->truth).size();
	  //if((tsize<hit->ntruthbin)&&(tsize>0)) {
	  //hit->contribBeta[tsize-1]=track->GetVelocity();
	  //hit->contribCharge[tsize-1]=(track->GetParticleDefinition())->GetPDGCharge();
	    //	    std::cout<<"tsize is "<<tsize<<" ntruthbin is "<<hit->ntruthbin<<std::endl;
	  //} else {
	    //std::cout<<"tsize is "<<tsize<<" ntruthbin is "<<hit->ntruthbin<<std::endl;

	  //}

	}


	// comment for this routine says Mark the track to be kept for MC truth propagation during hit processing


        //return true;
      }

      mark(h.track);
	
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
    typedef Geant4SensitiveAction<DualCrysCalorimeterSD> DualCrysCalorimeterSDAction;
  }}
DECLARE_GEANT4SENSITIVE(DualCrysCalorimeterSDAction)
DECLARE_GEANT4ACTION(WavelengthMinimumCut)
DECLARE_GEANT4ACTION(WavelengthnmwindCut)
