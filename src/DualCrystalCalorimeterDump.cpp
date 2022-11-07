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
#include "DualCrystalCalorimeterHit.h"
#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>

#include <iostream>
#include <vector>
#include <cerrno>

using namespace std;

namespace CalVision   {

  void DualCrystalCalorimeterdump_branch(int num_evts, TBranch* b)    {
    typedef vector<dd4hep::sim::Geant4Calorimeter::Hit*> CalHits;
    int num_evt;
    num_evt = b->GetEntries();
    if ( num_evt > 0 )   {
      CalHits*  hits = new CalHits();
      b->SetAddress(&hits);
      for(int ievt=0; ievt<num_evts; ++ievt )   {
	int nbyte = b->GetEntry(ievt);
	if ( nbyte > 0 )    {
	  cout << "Cal hits: " << nbyte << " bytes " << hits->size()  << endl;
	  for(size_t i=0; i<min(hits->size(),10UL); ++i)    {
	    DualCrystalCalorimeterHit* hit = (DualCrystalCalorimeterHit*)hits->at(i);
	    cout << b->GetName() << " Event " << ievt
		 << " Hit " << (int)i
		 << " type: "    << typeid(*hit).name()
		 << " deposit:"  << hit->energyDeposit
		 << " ncerenkov:"  << hit->ncerenkov
		 << endl;
	  }
	}
	else   {
	  cout << b->GetName() << " Event " << ievt << " NO DATA." << endl;
	}
      }
    }
  }

  /// Standalone function to dump data from a root file
  int Dump::DualCrystalCalorimeterdumpData(int num_evts, const char* file_name)    {
    if ( !file_name )    {
      cout << "Illegal file name: Input file cannot be opened!" << endl;
      return ENOENT;
    }
    TFile* f = TFile::Open(file_name);
    if ( !f )    {
      cout << "File " << file_name << " cannot be opened!" << endl;
      return ENOENT;
    }
    TTree* t = (TTree*)f->Get("EVENT");
    if ( !t )    {
      f->Close();
      cout << "File " << file_name << " cannot be read!" << endl;
      cout << "No branch with name EVENT!" << endl;
      return ENOENT;
    }
    t->Print();
    DualCrystalCalorimeterdump_branch(num_evts,t->GetBranch("CrystalEcalBarrelHits"));

    f->Close();
    delete f;
    return 0;
  }
}

	
