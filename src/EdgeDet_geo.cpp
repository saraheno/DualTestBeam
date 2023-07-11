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


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {

  std::cout<<"Creating EdgeDet"<<std::endl;

  static double tol = 0.0;


  xml_det_t     x_det     = e;


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is .. "<<det_name<<std::endl;




  // material to underly it all
  Material      air       = description.air();





  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_dim     = x_det.dimensions();
  double        hzlength   = x_dim.z_length()/2.;
  double        hxlength   = x_dim.height()/2.;
  double        hylength   = x_dim.width()/2.;
  double        thickness   = x_dim.thickness();

  std::cout<<" half zlength "<<hzlength<<
    " half height "<<hxlength<<
    " half width "<<hylength<<
" thickness "<<thickness<<std::endl;




  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_side( fX_struct.child( _Unicode(side) ) );




  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);



  //PolyhedraRegular hedra  (nphi,inner_r,outer_r+tol*2e0,zmaxt);

  Position a_pos(0.,0.,0.);


  dd4hep::Box abox1   (hxlength,hylength,hzlength);
  dd4hep::Box abox2   (hxlength-2*thickness,hylength-2*thickness,hzlength-2*thickness);
  dd4hep::Solid tmps = dd4hep::SubtractionSolid(abox1,abox2,a_pos);
  Volume        envelope  (det_name,tmps,air);
  std::cout<<"setting envelope visibility to "<<x_det.visStr()<<std::endl;
  envelope.setVisAttributes(description,x_det.visStr());

  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);
  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
  sens.setType("calorimeter");
  envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());






  dd4hep::Box absbox1(hxlength-tol,hylength-tol,hzlength-tol);
  dd4hep::Box absbox2(hxlength-tol-thickness,hylength-tol-thickness,hzlength-tol-thickness);
  dd4hep::Solid tmpsolid = dd4hep::SubtractionSolid(absbox1,absbox2,a_pos);
  Volume        Vside  ("side",tmpsolid,description.material(fX_side.materialStr()));
  Vside.setAttributes(description,fX_side.regionStr(),fX_side.limitsStr(),fX_side.visStr());
  if ( fX_side.isSensitive() ) {
    std::cout<<"setting EdgeDet side sensitive "<<std::endl;
    Vside.setSensitiveDetector(sens);
  }
  DetElement edge_det("side",det_id);


  PlacedVolume  edge_phv   = envelope.placeVolume(Vside,a_pos);
  edge_phv.addPhysVolID("side",1);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);










  std::cout<<"exiting DRFiber creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_EdgeDet,create_detector)

DECLARE_DEPRECATED_DETELEMENT(EdgeDet,create_detector)
