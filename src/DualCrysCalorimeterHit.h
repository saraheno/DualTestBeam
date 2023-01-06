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
#ifndef EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H
#define EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H

/// Framework include files
#include "DDG4/Geant4Data.h"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"


typedef ROOT::Math::XYZVector Position;
typedef ROOT::Math::XYZVector Direction;



namespace CalVision {

    const int wavelengthnbin=1000;
    const int truthnbin=1000;

  /// This is the hit definition.
  /** I took here the same definition of the default Geant4Tracker class,
   *  (see DDG4/Geant4Data.h)  but it could be anything else as well.
   *
   *  Please note:
   *  ============
   *  The MC truth handling as implemented in the Geant4ParticleHandler
   *  will not work with this class if the object(s) are saved with 
   *  the standard Geant4Output2ROOT event action. If the hit is 
   *  specialized, the output writing also must be specialized if
   *  MC truth handling should be supported.
   *  Otherwise it is sufficient to provide a ROOT dictionary as long as the
   *  base class dd4hep::sim::Geant4HitData is kept.
   *
   *  \author  M.Frank
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */
  class DualCrysCalorimeterHit : public dd4hep::sim::Geant4Calorimeter::Hit   {

  public:
    int ncerenkov,nscintillator;
    float wavelenmin=300;
    float wavelenmax=1000;
    int nwlbin=wavelengthnbin;
    std::array<int,wavelengthnbin>  ncerwave;
    std::array<int,wavelengthnbin> nscintwave;

    int ntruthbin=truthnbin;
    std::array<float,truthnbin> contribBeta;
    std::array<float,truthnbin> contribCharge;


  public:
    /// Default constructor
    DualCrysCalorimeterHit() = default;
    /// Initializing constructor
  DualCrysCalorimeterHit(const Position& cell_pos):dd4hep::sim::Geant4Calorimeter::Hit(cell_pos),ncerenkov(0),nscintillator(0) {

      for( int i=0;i<nwlbin;i++){
	ncerwave[i]=0;
	nscintwave[i]=0;
      }
      for( int i=0;i<ntruthbin;i++) {
	contribBeta[i]=0.;
	contribCharge[i]=0.;
      }

}

    /// Default destructor
    virtual ~DualCrysCalorimeterHit() = default;
    /// Assignment operator
    //DualCrysCalorimeterHit& operator=(const DualCrysCalorimeterHit& c);
  };

  /// Helper to dump data file
  /**
   *  Usage:  
   *  $> root.exe
   *  ....
   *  root [0] gSystem->Load("libDDG4Plugins.so");
   *  root [1] gSystem->Load("libDDG4_MySensDet.so");
   *  root [2] CalVision::Dump::dumpData(<num-ebents>,<file-name>);
   *
   */
  class Dump   {
  public:
    /// Standalone function to dump data from a root file
    static int DualCrysCalorimeterdumpData(int num_evts, const char* file_name);
  };
}

// CINT configuration
#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

/// Define namespaces
#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::sim;
#pragma link C++ namespace CalVision;
#pragma link C++ class     CalVision::DualCrysCalorimeterHit+;
#pragma link C++ class     CalVision::DualCrysCalorimeterDump;
#endif

#endif // EXAMPLES_DDDualCrys_SRC_DualCrysCalorimeterHIT_H
