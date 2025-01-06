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
//
// Specialized generic detector constructor
// 
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "DD4hep/OpticalSurfaces.h"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  std::cout<<"Creating DRUpS "<<std::endl;

  static double tol = 0.001;

  // material to underly it all
  Material      air       = description.air();

  // get stuff from xml file
  xml_det_t     x_det     = e;
  Layering      layering (e);


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is .. "<<det_name<<std::endl;


  DetElement    sdet      (det_name, det_id);
  std::cout<<"detector name is "<<det_name<<std::endl;




  std::cout<<"exiting DRUpS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DRUpS,create_detector)


