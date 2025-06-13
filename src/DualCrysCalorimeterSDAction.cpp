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
#include <algorithm>
#include <CLHEP/Units/PhysicalConstants.h>
#include "DualCrysCalorimeterHit.h"
#include "DDG4/EventParameters.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "DDG4/Factories.h"
#include "DD4hep/InstanceCount.h"
#include "DDG4/Geant4Random.h"
#include <G4Event.hh>
using namespace std;

// way too slow if track all photons for now
// so randomly delete photons after creation according to this fraction
//   dialScint=1.0, dialCer=1.0 to keep all photons 
double m_dialCherC= 10./800000.;
double m_dialScintC=100./20000000.;
double m_dialCherO= 100./800000.;
double m_dialScintO=1./20000000.;
float m_betarel=1/1.544;
int m_printlimitSCE=10;
int m_MAXEVENTSCE=10;

namespace CalVision {

  G4double fromEvToNm(G4double energy)
  {
    // there is some bug somewhere.  shouldn't need this facto
    //return  conversioneVnm/ energy*1000.;
    return 1239.84187 / energy*1000.;

  }


  int SCECOUNT=0;
  int SCECOUNT2=0;
  int OLDEVENTNUMBER=-1;


  
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

template <> void Geant4SensitiveAction<DualCrysCalorimeterSD>::initialize() {
      declareProperty("dialCherC", m_dialCherC);
      declareProperty("dialScintC", m_dialScintC);
      declareProperty("dialCherO", m_dialCherO);
      declareProperty("dialScintO", m_dialScintO);
      declareProperty("betarel", m_betarel);
      declareProperty("printlimitSCE", m_printlimitSCE);
      declareProperty("MAXEVENTSCE", m_MAXEVENTSCE);
    }

    
    /// Define collections created by this sensitivie action object
    template <> void Geant4SensitiveAction<DualCrysCalorimeterSD>::defineCollections()    {
      m_collectionID = declareReadoutFilteredCollection<CalVision::DualCrysCalorimeterSD::Hit>();
    }

    /// Method for generating hit(s) using the information of G4Step object.
    template <> bool Geant4SensitiveAction<DualCrysCalorimeterSD>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {



      Geant4Event&  evt = context()->event();



      int eventNumber = static_cast<G4Event const&>(context()->event()).GetEventID();


      if(eventNumber != OLDEVENTNUMBER) {
	if(eventNumber<m_MAXEVENTSCE) {
        std::cout<<" SDAction event number is "<<eventNumber<<std::endl;
        OLDEVENTNUMBER=eventNumber;
        SCECOUNT=0;
        SCECOUNT2=0;
	}
      }



      //std::cout<<" in SD action event number is "<<eventNumber<<std::endl;


      if((eventNumber==1)&&(SCECOUNT==1)) {
	std::cout<<"event number is "<<eventNumber<<std::endl;
	std::cout<<"DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!"<<std::endl;
	std::cout<<"dialCher for PbWO4 and BGO is "<<m_dialCherC<<std::endl;
	std::cout<<"dialScint for PbWO4 and BGO is "<<m_dialScintC<<std::endl;
	std::cout<<"dialCher for others is "<<m_dialCherO<<std::endl;
	std::cout<<"dialScint for others is "<<m_dialScintO<<std::endl;
	std::cout<<" you need to use this to interpret your results"<<std::endl;
	std::cout<<" also counting as relativiistic particles with beta >"<<m_betarel<<" you should adjust according to the index of your media"<<std::endl;
      }

      bool SCEPRINT=(SCECOUNT<m_printlimitSCE);
      //if(SCEPRINT) std::cout<<"scecount is "<<SCECOUNT<<" print is "<<SCEPRINT<<std::endl;



      G4StepPoint *thePrePoint = step->GetPreStepPoint();
      G4StepPoint *thePostPoint = step->GetPostStepPoint();
      //      const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
      //const G4ThreeVector &thePostPosition = thePostPoint->GetPosition();
      G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
      G4double pretime = thePrePoint->GetGlobalTime();
      G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
      G4double posttime = thePostPoint->GetGlobalTime();
      G4String thePrePVName = "";
      if (thePrePV)
	thePrePVName = thePrePV->GetName();
      G4String thePostPVName = "";
      if (thePostPV)
	thePostPVName = thePostPV->GetName();
      G4Track *theTrack = step->GetTrack();
      G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

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
	if(SCECOUNT2<m_printlimitSCE) {
	  std::cout<<"DRcalo deposit "<<contrib.deposit<<" position ("<<pos.X<<","<<pos.Y<<","<<pos.Z<<") string "<<handler.path().c_str()<<" volume id "<<cell<<" event "<<eventNumber<<std::endl;
	  if(SCECOUNT2<m_printlimitSCE+1) SCECOUNT2+=1;
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
      //      if(SCEPRINT) std::cout<<"track name is "<< (track->GetDefinition())->GetParticleName()<<std::endl;

      float avearrival=(pretime+posttime)/2.;
      int jbin=-1;
      float tbinsize=(hit->timemax-hit->timemin)/hit->nfinebin;
      jbin = (avearrival-hit->timemin)/tbinsize;
      jbin = std::min(jbin,(hit->nfinebin)-1);
      int jbinz=-1;
      float tbinsizez=(hit->timemaxz-hit->timemin)/hit->nfinebin;
      jbinz = (avearrival-hit->timemin)/tbinsizez;
      jbinz = std::min(jbinz,(hit->nfinebin)-1);


      if(thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("Inelast")){
	hit->n_inelastic+=1;
      }

      //photons
      if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
	if(SCEPRINT) std::cout<<"     in volume ID "<<cell<<std::endl;

	if(SCECOUNT<m_printlimitSCE+1) SCECOUNT+=1;
	//	if(SCEPRINT) std::cout<<"optical photon"<<std::endl;

	bool OptN = (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")||(track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys");

	//if(track->GetParentID()!=1) SCEPRINT=1;
	if( (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "CerenkovPhys")&&(track->GetCreatorProcess()->G4VProcess::GetProcessName() != "ScintillationPhys")  ) SCEPRINT=1;  // print if less than scecount (set earlier) or if photon has weird source


        float wavelength=fromEvToNm(track->GetTotalEnergy()/eV);
        int ibin=-1;
        float binsize=(hit->wavelenmax-hit->wavelenmin)/hit->nfinebin;
        ibin = (wavelength-hit->wavelenmin)/binsize;
	ibin = std::min(ibin,(hit->nfinebin)-1);




        int xbin=-1;
        float xbinsize=(hit->xmax-hit->xmin)/hit->ncoarsebin;
        xbin = (thePrePoint->GetPosition().x()-hit->xmin)/xbinsize;


        int ybin=-1;
        float ybinsize=(hit->ymax-hit->ymin)/hit->ncoarsebin;
        ybin = (thePrePoint->GetPosition().y()-hit->ymin)/ybinsize;


	int phstep = track->GetCurrentStepNumber();
	
	
	if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")  {
	  if(SCEPRINT) std::cout<<" found Cerenkov photon"<<std::endl;

	  if(amedia.find("kill")!=std::string::npos) 
	    { 
	      if(SCEPRINT) std::cout<<"killing photon"<<std::endl;
	      //std::cout<<"killing photon"<<std::endl;
	      //	      SCEPRINT=1;
	      if(phstep>1) {  // don't count photons created in kill media
		hit->ncerenkov+=1;
                //if(ibin>-1&&ibin<hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
                if(jbin>-1&&jbin<hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
		if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->ncertimez).at(jbinz))+=1;
		//                if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->cerhitpos).at(xbin).at(ybin))+=1;
                if(SCEPRINT) std::cout<<" cer photon kill pre time "<<pretime<<std::endl;
                if(SCEPRINT) std::cout<<" cer photon kill post time "<<posttime<<std::endl;
	      }
	      track->SetTrackStatus(fStopAndKill);
	    }
	  else if(amedia.find("BlackHole")!=std::string::npos) {
	      if(phstep>1) {  // don't count photons created in kill media
		hit->ncerenkov+=1;
		//if(ibin>-1&&ibin<hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
              if(jbin>-1&&jbin<hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
	      if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->ncertimez).at(jbinz))+=1;
	      //              if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->cerhitpos).at(xbin).at(ybin))+=1;
	      }
	      track->SetTrackStatus(fStopAndKill);
	  }
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1)  ) hit->ncerenkov+=1;
	    if( (phstep==1)  ) {
	      hit->ncerenkov+=1;
	      //if(ibin>-1&&ibin<hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
              if(jbin>-1&&jbin<hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
	      if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->ncertimez).at(jbinz))+=1;
	      //              if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->cerhitpos).at(xbin).at(ybin))+=1;

	      //Geant4Event&  evt = context()->event();
	      dd4hep::sim::Geant4Random& rnd = evt.random();
	      if((amedia.find("E_PbWO4")!=std::string::npos)||(amedia.find("BGO")!=std::string::npos)) {
		if(rnd.rndm()>m_dialCherC) track->SetTrackStatus(fStopAndKill);
	      } else{
		if(rnd.rndm()>m_dialCherO) track->SetTrackStatus(fStopAndKill);
	      }
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
	      //std::cout<<"killing photon"<<std::endl;
	      if(phstep>1) {
		hit->nscintillator+=1;
		//if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
                if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
		if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->nscinttimez).at(jbinz))+=1;
		//                if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->scinthitpos).at(xbin).at(ybin))+=1;

                if(SCEPRINT) std::cout<<" scint photon kill pre time "<<pretime<<std::endl;
                if(SCEPRINT) std::cout<<" scint photon kill post time "<<posttime<<std::endl;
	      }
	      track->SetTrackStatus(fStopAndKill);}
	  else if(amedia.find("BlackHole")!=std::string::npos) {
	      if(phstep>1) {  // don't count photons created in kill media
		hit->nscintillator+=1;
		//if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
                if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
		if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->nscinttimez).at(jbinz))+=1;
		//                if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->scinthitpos).at(xbin).at(ybin))+=1;

	      }
	      track->SetTrackStatus(fStopAndKill);
	  }
	  else {
	    //	    if( (track->GetParentID()==1)&&(track->GetCurrentStepNumber()==1) ) hit->nscintillator+=1; 
	    if( (phstep==1) ) {
	      hit->nscintillator+=1; 
	      //if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
                if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
		if(jbinz>-1&&jbinz<hit->nfinebin) ((hit->nscinttimez).at(jbinz))+=1;
		//                if(( xbin>-1&&xbin<hit->ncoarsebin)&&(ybin>-1&&ybin<hit->ncoarsebin)) ((hit->scinthitpos).at(xbin).at(ybin))+=1;


		//Geant4Event&  evt = context()->event();
	      dd4hep::sim::Geant4Random& rnd = evt.random();
	      if((amedia.find("E_PbWO4")!=std::string::npos)||(amedia.find("BGO")!=std::string::npos)) {
		if(rnd.rndm()>m_dialScintC) track->SetTrackStatus(fStopAndKill);
	      } else {
		if(rnd.rndm()>m_dialScintO) track->SetTrackStatus(fStopAndKill);
	      }
	    }
	  }

          //return false;
        }
	else {
          if(SCEPRINT) std::cout<<"      other photon"<<std::endl;
          track->SetTrackStatus(fStopAndKill);
          //return false;
	}

	if(SCEPRINT) {
	  std::cout<<"     SCECOUNT="<<SCECOUNT<<std::endl;
	  std::cout<<"     event number is "<<eventNumber<<std::endl;
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

	if(amedia.find("BlackHole")!=std::string::npos) {  //edge detector
	  // note this is the edge detector, not the kill media of photodetector
	  hit->energyDeposit += track->GetKineticEnergy();
	  if((abs(TrPDGid)==11)||(abs(TrPDGid)==22)) {
	    hit->edepepgam+=contrib.deposit;	  
	  }
	  track->SetTrackStatus(fStopAndKill);
	} else {
      //add information about each contribution to the hit
	  hit->truth.emplace_back(contrib);

	  hit->energyDeposit += contrib.deposit;
	  //	  hit->contribBeta.emplace_back(track->GetVelocity()/CLHEP::c_light*10000.);
	  float aabeta=track->GetVelocity()/CLHEP::c_light;
	  hit->contribBeta.emplace_back(aabeta);
	  hit->contribCharge.emplace_back((track->GetParticleDefinition())->GetPDGCharge());
	  if(jbin>-1&&jbin<hit->nfinebin) ((hit->edeptime).at(jbin))+=contrib.deposit;
	  if(aabeta>m_betarel) {
	    hit->edeprelativistic+=contrib.deposit;	  
	    if(jbin>-1&&jbin<hit->nfinebin) ((hit->ereldeptime).at(jbin))+=contrib.deposit;
	  }
	  if((abs(TrPDGid)==11)||(abs(TrPDGid)==22)) {
	    hit->edepepgam+=contrib.deposit;	  
	  }

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

namespace dd4hep { 
	namespace sim {
		using namespace CalVision;
		struct WavelengthMinimumCut : public dd4hep::sim::Geant4Filter  {// Energy cut value
			double m_wavelengthCut;
			public: // Constructor.
			WavelengthMinimumCut(dd4hep::sim::Geant4Context* c, const std::string& n);
			virtual ~WavelengthMinimumCut(); // Standard destructor
			
			// Filter action. Return true if hits should be processed
			virtual bool operator()(const G4Step* step) const  override  final  {
				bool test=true;
				G4Track *theTrack = step->GetTrack();
				if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
					float energy = theTrack->GetTotalEnergy()/eV;
					float wave   = fromEvToNm(energy);
					if(wave < m_wavelengthCut) {
						theTrack->SetTrackStatus(fStopAndKill);
						test=false;
					}
				}
				return test;
			}
			virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
				return true;
			}
		};
		WavelengthMinimumCut::WavelengthMinimumCut(Geant4Context* c, const std::string& n) //Constructor
			: Geant4Filter(c,n) {
				InstanceCount::increment(this);
				declareProperty("Cut",m_wavelengthCut=0.0);
			}
		WavelengthMinimumCut::~WavelengthMinimumCut() { // Standard destructor
			InstanceCount::decrement(this);
		}
		
		struct WavelengthnmwindCut : public dd4hep::sim::Geant4Filter  { // Energy cut value
			double m_wavelengthstart;
			public: // Constructor.
			WavelengthnmwindCut(dd4hep::sim::Geant4Context* c, const std::string& n);
			virtual ~WavelengthnmwindCut(); // Standard destructor
			// Filter action. Return true if hits should be processed
			virtual bool operator()(const G4Step* step) const  override  final  {
				bool test=true;
				G4Track *theTrack = step->GetTrack();
				if(theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
					float energy = theTrack->GetTotalEnergy()/eV;
					float wave   = fromEvToNm(energy);
					if((wave < m_wavelengthstart) || (wave > m_wavelengthstart+0.5) ) {
						theTrack->SetTrackStatus(fStopAndKill);
						test=false;
					}
				}
				return test;
			}
			virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
				return true;
			}
		};
		WavelengthnmwindCut::WavelengthnmwindCut(Geant4Context* c, const std::string& n) // Constructor.
			: Geant4Filter(c,n) {
				InstanceCount::increment(this);
				declareProperty("Cut",m_wavelengthstart=0.0);
			}
		WavelengthnmwindCut::~WavelengthnmwindCut() { // Standard destructor
			InstanceCount::decrement(this);
		}
	}
}  // end using namespace

//--- Factory declaration
namespace dd4hep { 
	namespace sim {
		typedef Geant4SensitiveAction<DualCrysCalorimeterSD> DualCrysCalorimeterSDAction;
	}
}

DECLARE_GEANT4SENSITIVE(DualCrysCalorimeterSDAction)
DECLARE_GEANT4ACTION(WavelengthMinimumCut)
DECLARE_GEANT4ACTION(WavelengthnmwindCut)
