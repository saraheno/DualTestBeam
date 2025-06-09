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

  // an air box
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





  // make 6 sides of the edge detector
  DetElement edge_det("side",det_id);
  PlacedVolume  edge_phv;
  Transform3D tr;


  dd4hep::Box absbox1(hxlength-tol,hylength-tol,thickness);
  Volume        Vside1  ("side",absbox1,description.material(fX_side.materialStr()));
  Vside1.setAttributes(description,fX_side.regionStr(),fX_side.limitsStr(),fX_side.visStr());
  if ( fX_side.isSensitive() ) {
    std::cout<<"setting EdgeDet side sensitive "<<std::endl;
    Vside1.setSensitiveDetector(sens);
  }

  tr={RotationZYX(0.,0.,0.),Position(0.,0.,-hzlength+tol)};
  edge_phv   = envelope.placeVolume(Vside1,tr);
  edge_phv.addPhysVolID("side",1);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);

  tr={RotationZYX(0.,0.,0.),Position(0.,0.,hzlength-tol)};
  edge_phv   = envelope.placeVolume(Vside1,tr);
  edge_phv.addPhysVolID("side",2);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);



  dd4hep::Box absbox2(hxlength-tol,hzlength-tol,thickness);
  Volume        Vside2  ("side",absbox2,description.material(fX_side.materialStr()));
  Vside2.setAttributes(description,fX_side.regionStr(),fX_side.limitsStr(),fX_side.visStr());
  if ( fX_side.isSensitive() ) {
    std::cout<<"setting EdgeDet side sensitive "<<std::endl;
    Vside2.setSensitiveDetector(sens);
  }

  tr={RotationZYX(0.,0.,M_PI/2.),Position(0.,-hylength+tol,0.)};
  edge_phv   = envelope.placeVolume(Vside2,tr);
  edge_phv.addPhysVolID("side",3);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);

  tr={RotationZYX(0.,0.,M_PI/2.),Position(0.,hylength-tol,0.)};
  edge_phv   = envelope.placeVolume(Vside2,tr);
  edge_phv.addPhysVolID("side",4);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);




  dd4hep::Box absbox3(thickness,hylength-tol,hzlength-tol);
  Volume        Vside3  ("side",absbox3,description.material(fX_side.materialStr()));
  Vside3.setAttributes(description,fX_side.regionStr(),fX_side.limitsStr(),fX_side.visStr());
  if ( fX_side.isSensitive() ) {
    std::cout<<"setting EdgeDet side sensitive "<<std::endl;
    Vside3.setSensitiveDetector(sens);
  }

  tr={RotationZYX(0.,0.,0.),Position(-hxlength+tol,0.,0.)};
  edge_phv   = envelope.placeVolume(Vside3,tr);
  edge_phv.addPhysVolID("side",5);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);

  tr={RotationZYX(0.,0.,0.),Position(hxlength-tol,0.,0.)};
  edge_phv   = envelope.placeVolume(Vside3,tr);
  edge_phv.addPhysVolID("side",6);
  edge_det.setPlacement(edge_phv);
  edge_phv.addPhysVolID("system",det_id);











  std::cout<<"exiting EdgeDetector creator"<<std::endl;
  std::cout<<"**********************************************************************************"<<endl;

  return sdet;
}

DECLARE_DETELEMENT(EdgeDet,create_detector)


