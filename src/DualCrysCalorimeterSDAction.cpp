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
//   so randomly delete photons after creation according to this fraction
//   dialScint=1.0, dialCer=1.0 to keep all photons 
double dialCherC  = 10./800000.;
double dialScintC = 100./20000000.;
double dialCherO  = 100./800000.;
double dialScintO = 1./20000000.;
float betarel=1/1.544;	//depends on the media refractive index
int printlimitSCE=10;
int MAXEVENT=10;

namespace CalVision {
	G4double fromEvToNm(G4double energy){
		// there is some bug somewhere.  shouldn't need this factor
		//return  conversioneVnm/ energy*1000.;
		return 1239.84187 / energy*1000.;
	}	
	int SCECOUNT=0;
	int SCECOUNT2=0;
	int OLDEVENTNUMBER=-1;
	class DualCrysCalorimeterSD {
		public:
			typedef DualCrysCalorimeterHit Hit;
			// If we need special data to personalize the action, put it here
			//int mumDeposits = 0;
			//double integratedDeposit = 0;
			};
}

// Namespace for the AIDA detector description toolkit
namespace dd4hep {  // Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
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
		
		// Method for generating hit(s) using the information of G4Step object.
		template <> bool Geant4SensitiveAction<DualCrysCalorimeterSD>::process(const G4Step* step,G4TouchableHistory* /*hist*/ ) {
			Geant4Event&  evt = context()->event();
			int eventNumber = static_cast<G4Event const&>(context()->event()).GetEventID();
			if(eventNumber != OLDEVENTNUMBER) {
				if(eventNumber<MAXEVENT) {
					cout<< " SDAction event number is " << eventNumber << endl;
					OLDEVENTNUMBER=eventNumber;
					SCECOUNT=0;
					SCECOUNT2=0;
				}
			}	
			if(( eventNumber == 1 ) && ( SCECOUNT == 1 ) ) {
				cout << "Info to interpret results:"		 << endl;
				cout << "	event number is " << eventNumber << endl;
				cout << "	DANGER DANGER "   		 << endl;
				cout << "	dialCher  for PbWO4 and BGO = "	 << dialCherC  << ", for others" << dialCherO  << endl;
				cout << "	dialScint for PbWO4 and BGO = "  << dialScintC << ", for others" << dialScintO << endl;
				cout << "	relativistic particles count with beta = " << betarel << " adjust with media RI" << endl;
			}
			
			bool SCEPRINT = ( SCECOUNT < printlimitSCE);
			G4StepPoint *thePrePoint  = step->GetPreStepPoint();
			G4StepPoint *thePostPoint = step->GetPostStepPoint();
			
			G4VPhysicalVolume *thePrePV 	= thePrePoint->GetPhysicalVolume();
			G4double pretime 		= thePrePoint->GetGlobalTime();
			G4VPhysicalVolume *thePostPV 	= thePostPoint->GetPhysicalVolume();
			G4double posttime 		= thePostPoint->GetGlobalTime();
			
			G4String thePrePVName  = "";
			if (thePrePV)	thePrePVName = thePrePV->GetName();
			G4String thePostPVName = "";
			if (thePostPV)	thePostPVName = thePostPV->GetName();
			G4Track *theTrack = step->GetTrack();
			G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

			Geant4StepHandler h(step);
			HitContribution contrib = DualCrysCalorimeterHit::extractContribution(step);
			Geant4HitCollection*  coll    = collection(m_collectionID);
			VolumeID cell = 0;	
			try {
				cell = cellID(step);
			} catch(std::runtime_error &e) {
				stringstream out;
				out << setprecision(20) << scientific;
				out << "ERROR: " << e.what()  << endl;
				out << "Position: "
					<< "Pre (" <<  setw(24) << step->GetPreStepPoint()->GetPosition() << ") "
					<< "Post (" << setw(24) << step->GetPostStepPoint()->GetPosition() << ") "
					<< endl;
				out << "Momentum: "
					<< " Pre (" <<setw(24) << step->GetPreStepPoint() ->GetMomentum()  << ") "
					<< " Post (" <<setw(24) << step->GetPostStepPoint()->GetMomentum() << ") "
					<< endl;
				cout << out.str();
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
				  cout << "DRcalo deposit " << contrib.deposit << " position (" <<
					  pos.X << "," << pos.Y << "," << pos.Z << ") string "  <<
					  handler.path().c_str() << " volume id " << cell << " event " << eventNumber << endl;
				  if(SCECOUNT2<printlimitSCE+1) SCECOUNT2+=1;
				}
				if ( 0 == hit->cellID )  { // for debugging only!
					hit->cellID = cellID(step);
					except("+++ Invalid CELL ID for hit!");
				}
			}
			G4Track *track =  step->GetTrack();
			string amedia = ((track->GetMaterial())->GetName());
			float avearrival=(pretime+posttime)/2.;
			int jbin=-1;
			float tbinsize=(hit->timemax-hit->timemin)/hit->nfinebin;
			jbin = (avearrival-hit->timemin)/tbinsize;
			jbin = std::min(jbin,(hit->nfinebin)-1);
			if(thePostPoint->GetProcessDefinedStep()->GetProcessName().contains("Inelast")) hit->n_inelastic+=1;
			
			//photons
			if( track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() )  {
				if(SCEPRINT)	cout << "     in volume ID " << cell << endl;
				if(SCECOUNT<printlimitSCE+1)	SCECOUNT+=1;
				
				bool OptN = (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")
					||  (track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys");
				// print if less than scecount (set earlier) or if photon has weird source
				if( (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "CerenkovPhys")
				&&  (track->GetCreatorProcess()->G4VProcess::GetProcessName() != "ScintillationPhys")  ) SCEPRINT=1;

				float wavelength = fromEvToNm(track->GetTotalEnergy()/eV);
				int ibin=-1;
				float binsize    = (hit->wavelenmax - hit->wavelenmin) / hit->nfinebin;
				ibin = (wavelength - hit->wavelenmin) / binsize;
				ibin = std::min(ibin, (hit->nfinebin) - 1);
				
				int xbin=-1;
				float xbinsize   = (hit->xmax - hit->xmin) / hit->ncoarsebin;
				xbin = (thePrePoint->GetPosition().x() - hit->xmin) / xbinsize;

				int ybin=-1;
				float ybinsize=(hit->ymax - hit->ymin) / hit->ncoarsebin;
				ybin = (thePrePoint->GetPosition().y() - hit->ymin) / ybinsize;

				int phstep = track->GetCurrentStepNumber();
				if ( track->GetCreatorProcess()->G4VProcess::GetProcessName() == "CerenkovPhys")  {
					if(SCEPRINT) cout << " found Cerenkov photon" << endl;
					if(amedia.find("kill")!=string::npos){
						if(SCEPRINT) cout << "	killing photon" << endl;
						if(phstep>1) {  // don't count photons created in kill media
							hit->ncerenkov+=1;
							if(ibin > -1 && ibin < hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
							if(jbin > -1 && jbin < hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
							if(SCEPRINT) cout<<" cer photon kill pre time "  << pretime  << endl;
							if(SCEPRINT) cout<<" cer photon kill post time " << posttime << endl;
						}
						track->SetTrackStatus(fStopAndKill);
					}
					else if(amedia.find("BlackHole")!=string::npos) {
						if(phstep>1) {  // don't count photons created in kill media
							hit->ncerenkov += 1;
							if(ibin > -1 && ibin < hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
							if(jbin > -1 && jbin < hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
						}
						track->SetTrackStatus(fStopAndKill);
					}
					else {
						if( (phstep==1)  ) {
							hit->ncerenkov += 1;
							if(ibin > -1 && ibin < hit->nfinebin) ((hit->ncerwave).at(ibin))+=1;
							if(jbin > -1 && jbin < hit->nfinebin) ((hit->ncertime).at(jbin))+=1;
							dd4hep::sim::Geant4Random& rnd = evt.random();
							if((amedia.find("E_PbWO4")!=std::string::npos)||(amedia.find("BGO")!=std::string::npos)) {
								if(rnd.rndm()>dialCherC) track->SetTrackStatus(fStopAndKill);
							} else{
								if(rnd.rndm()>dialCherO) track->SetTrackStatus(fStopAndKill);
							}
						}
					}
				}
				else if (  track->GetCreatorProcess()->G4VProcess::GetProcessName() == "ScintillationPhys"  ) {
					if(SCEPRINT) cout << "found Scintillation photon" << endl;
					std::string amedia = ((track->GetMaterial())->GetName());
					if(amedia.find("kill")!=std::string::npos) {
						if(SCEPRINT) cout << "	killing photon " << endl;
						if(phstep>1) {
							hit->nscintillator += 1;
							if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
							if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
							if(SCEPRINT) cout<<"		scint photon kill pre time "  << pretime  << endl;
							if(SCEPRINT) cout<<"		scint photon kill post time " << posttime << endl;
						}
						track->SetTrackStatus(fStopAndKill);
					}
					else if(amedia.find("BlackHole")!=std::string::npos) {
						if(phstep>1) {  // don't count photons created in kill media
							hit->nscintillator += 1;
							if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
							if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
						}
						track->SetTrackStatus(fStopAndKill);
					}
					else {
						if( (phstep==1) ) {
							hit->nscintillator += 1; 
							if((ibin>-1)&&(ibin<hit->nfinebin)) ((hit->nscintwave).at(ibin))+=1;
							if(jbin>-1&&jbin<hit->nfinebin) ((hit->nscinttime).at(jbin))+=1;
							//Geant4Event&  evt = context()->event();
							dd4hep::sim::Geant4Random& rnd = evt.random();
							if((amedia.find("E_PbWO4")!=std::string::npos)||(amedia.find("BGO")!=std::string::npos)) {
								if(rnd.rndm()>dialScintC) track->SetTrackStatus(fStopAndKill);
							} else {
								if(rnd.rndm()>dialScintO) track->SetTrackStatus(fStopAndKill);
							}
						}
					}
				}
				else {
					if(SCEPRINT) cout << "other photon" << endl;
					track->SetTrackStatus(fStopAndKill);
				}
				
				if(SCEPRINT) {
					cout<< " SCECOUNT="<<SCECOUNT<<std::endl;
					cout<< " event number is " << eventNumber<<std::endl;
					cout<< " Found photon:   " << track->GetCreatorProcess()->G4VProcess::GetProcessName() 		<< endl;
					cout<< "   PreVolume name:  		" << thePrePVName					<< endl;
                                        cout<< "   PostVolume name: 		" << thePostPVName					<< endl;
					cout<< "   photon mother:   		" << track->GetParentID()				<< endl;
					cout<< "   photon material: 		" << (track->GetMaterial())->GetName()			<< endl;
					cout<< "   photon creator process: 	" << (track->GetCreatorProcess())->GetProcessName()	<< endl;
					cout<< "   photon  process type: 	" << (track->GetCreatorProcess())->GetProcessType()	<< endl;
					cout<< "   photon sub process: 		" << (track->GetCreatorProcess())->GetProcessSubType()	<< endl;
					cout<< "   photon current step number: 	" << track->GetCurrentStepNumber()			<< endl;
					cout<< "   photon energy: 		" << track->GetTotalEnergy() / eV			<< endl;
					cout<< "   photon wavelength: 		" << fromEvToNm(track->GetTotalEnergy() / eV)		<< endl;
					cout<< "   number of cherenkov:		" << hit->ncerenkov					<< endl;
					cout<< "   number of scintillation: 	" << hit->nscintillator					<< endl;
				}
			}
			else {   // particles other than optical photons
				if(SCEPRINT) cout << " NOT optical photon" << endl;
				if(amedia.find("BlackHole")!=std::string::npos) {  //edge detector, not the kill PD media
					hit->energyDeposit += track->GetKineticEnergy();
					track->SetTrackStatus(fStopAndKill);
				} 
				else { //add information about each contribution to the hit
					hit->truth.emplace_back(contrib);
					hit->energyDeposit += contrib.deposit;
					
					float aabeta=track->GetVelocity()/CLHEP::c_light;
					hit->contribBeta.emplace_back(aabeta);
					hit->contribCharge.emplace_back((track->GetParticleDefinition())->GetPDGCharge());
					if(jbin>-1&&jbin<hit->nfinebin) ((hit->edeptime).at(jbin))+=contrib.deposit;
					if(aabeta>betarel) {
						hit->edeprelativistic+=contrib.deposit;	  
						if(jbin>-1&&jbin<hit->nfinebin) ((hit->ereldeptime).at(jbin))+=contrib.deposit;
					}
					if((abs(TrPDGid)==11)||(abs(TrPDGid)==22)) {
						hit->edepepgam+=contrib.deposit;
					}
				}// comment for this routine says mark the track to be kept for MC truth propagation during hit processing
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
