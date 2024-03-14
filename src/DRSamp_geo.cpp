//
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
// Specialized generic detector constructor for fiber tube calorimeter.
// See https://github.com/AIDASoft/DD4hep/issues/1173 for details
//
// Detector geometry structure
//                <system>  <layer>  <tube>  <hole> <type>
// /world_volume/FiberCalo/rowtube_1/brass_1/hole_1/quartz_1
//                                  /brass_2/hole_1/scintillator_1    brass_1 == brass_2 == brass....n
//                                  /brass_3/hole_1/quartz_1         
//                                  /brass_4/hole_1/scintillator_1
//                                  /brass_5/hole_1/quartz_1
// brass_1/quartz_1   Volume(brass) / hole
//                    Volume(hole)  / quartz
//              alt:  Volume(hole)  / scintillator
// 
// /world_volume/FiberCalo/rowtube_1/brass/quartz
//
// Dump using:
// $> geoPluginRun -input examples/ClientTests/compact/FiberTubeCalorimeter.xml 
//              -volmgr -print 3 -plugin DD4hep_DetectorCheck 
//              -name FiberTubeCalorimeter -geometry -structure -sensitive -volmgr -ignore
//
// Display using:
// $> geoDisplay examples/ClientTests/compact/FiberTubeCalorimeter.xml
//
//==========================================================================
#include <DD4hep/DetFactoryHelper.h>
#include "DD4hep/OpticalSurfaces.h"
#include <DD4hep/Printout.h>


#include <iomanip>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
  std::cout<<"create sampling calorimeter "<<std::endl;
  //
  constexpr double tol = 0.0;
  xml_det_t     x_det = e;

  // for volume tags in detector
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  Material      air       = description.air();

  DetElement    sdet      (det_name, det_id);
  std::cout<<"detector name is "<<det_name<<std::endl;

  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_dim     = x_det.dimensions();

  double hwidth   = x_dim.width()/2.;
  int Ncount = x_dim.dim_x();
  int Nlayers = x_dim.dim_z();
  double zoffset = x_dim.z1();
  double agap = x_dim.gap();
  double shave=x_dim.z2();

  std::cout<<"half width tower is "<<hwidth<<std::endl;
  std::cout<<" array size is "<<Ncount<<std::endl;
  std::cout<<" number of layers is "<<Nlayers<<std::endl;
  std::cout<<" z offset is "<<zoffset<<std::endl;
  std::cout<<" gap between array elements is "<<agap<<std::endl;

  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_absorb( fX_struct.child( _Unicode(absorb) ) );
  xml_comp_t fX_kill1( fX_struct.child( _Unicode(kill1) ) );
  xml_comp_t fX_scint( fX_struct.child( _Unicode(scint) ) );
  xml_comp_t fX_kill2( fX_struct.child( _Unicode(kill2) ) );
  xml_comp_t fX_separate1( fX_struct.child( _Unicode(separate1) ) );
  xml_comp_t fX_kill3( fX_struct.child( _Unicode(kill3) ) );
  xml_comp_t fX_ceren( fX_struct.child( _Unicode(ceren) ) );
  xml_comp_t fX_kill4( fX_struct.child( _Unicode(kill4) ) );
  xml_comp_t fX_separate2( fX_struct.child( _Unicode(separate2) ) );
  xml_comp_t fX_unitbox( fX_struct.child( _Unicode(unitbox) ) );
  xml_comp_t fX_tower( fX_struct.child( _Unicode(tower) ) );


  std::cout<<" got fXs "<<std::endl;

  std::cout<<"making optical surface"<<std::endl;
  OpticalSurfaceManager surfMgr = description.surfaceManager();
  OpticalSurface cryS  = surfMgr.opticalSurface("/world/DRSamp#mirrorSurface");

  std::cout<<" surface made"<<std::endl;

  Material mat;
  Transform3D trafo;
  PlacedVolume pv;
  Solid sol;

  sens.setType("calorimeter");

  // absorber
  std::cout<<" absorber made of "<<fX_absorb.materialStr()<<std::endl;
  sol = Box(hwidth,hwidth,fX_absorb.thickness()/2.);
  mat = description.material(fX_absorb.materialStr());
  Volume absorb_vol(fX_absorb.nameStr(), sol, mat);
  absorb_vol.setAttributes(description,fX_absorb.regionStr(),fX_absorb.limitsStr(),fX_absorb.visStr());
  if ( fX_absorb.isSensitive() ) {
    absorb_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << absorb_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_absorb.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_absorb.isSensitive()) << endl;


  // kill1
  sol = Box(hwidth,hwidth,fX_kill1.thickness()/2.);
  mat = description.material(fX_kill1.materialStr());
  Volume kill1_vol(fX_kill1.nameStr(), sol, mat);
  kill1_vol.setAttributes(description,fX_kill1.regionStr(),fX_kill1.limitsStr(),fX_kill1.visStr());
  if ( fX_kill1.isSensitive() ) {
    kill1_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << kill1_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_kill1.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_kill1.isSensitive()) << endl;


  // scintillator
  sol = Box(hwidth,hwidth,fX_scint.thickness()/2.-shave/2.);
  mat = description.material(fX_scint.materialStr());
  Volume scint_vol(fX_scint.nameStr(), sol, mat);
  scint_vol.setAttributes(description,fX_scint.regionStr(),fX_scint.limitsStr(),fX_scint.visStr());
  if ( fX_scint.isSensitive() ) {
    scint_vol.setSensitiveDetector(sens);
  }
  // this does work, but need to think where to put the kill media before I implement
  //SkinSurface mirror = SkinSurface(description,sdet,"HallCrys", cryS,scint_vol);
  //mirror.isValid();
  cout << setw(28) << left << scint_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_scint.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_scint.isSensitive()) << endl;


  // kill2
  sol = Box(hwidth,hwidth,fX_kill2.thickness()/2.);
  mat = description.material(fX_kill2.materialStr());
  Volume kill2_vol(fX_kill2.nameStr(), sol, mat);
  kill2_vol.setAttributes(description,fX_kill2.regionStr(),fX_kill2.limitsStr(),fX_kill2.visStr());
  if ( fX_kill2.isSensitive() ) {
    kill2_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << kill2_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_kill2.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_kill2.isSensitive()) << endl;


  //  separate
  sol = Box(hwidth,hwidth,fX_separate1.thickness()/2.);
  mat = description.material(fX_separate1.materialStr());
  Volume separate1_vol(fX_separate1.nameStr(), sol, mat);
  separate1_vol.setAttributes(description,fX_separate1.regionStr(),fX_separate1.limitsStr(),fX_separate1.visStr());
  if ( fX_separate1.isSensitive() ) {
    separate1_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << separate1_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_separate1.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_separate1.isSensitive()) << endl;


  // kill3
  sol = Box(hwidth,hwidth,fX_kill3.thickness()/2.);
  mat = description.material(fX_kill3.materialStr());
  Volume kill3_vol(fX_kill3.nameStr(), sol, mat);
  kill3_vol.setAttributes(description,fX_kill3.regionStr(),fX_kill3.limitsStr(),fX_kill3.visStr());
  if ( fX_kill3.isSensitive() ) {
    kill3_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << kill3_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_kill3.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_kill3.isSensitive()) << endl;


  // cerenkov
  sol = Box(hwidth,hwidth,fX_ceren.thickness()/2.-shave/2.);
  mat = description.material(fX_ceren.materialStr());
  Volume ceren_vol(fX_ceren.nameStr(), sol, mat);
  ceren_vol.setAttributes(description,fX_ceren.regionStr(),fX_ceren.limitsStr(),fX_ceren.visStr());
  if ( fX_ceren.isSensitive() ) {
    ceren_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << ceren_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_ceren.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_ceren.isSensitive()) << endl;


  // kill4
  sol = Box(hwidth,hwidth,fX_kill4.thickness()/2.);
  mat = description.material(fX_kill4.materialStr());
  Volume kill4_vol(fX_kill4.nameStr(), sol, mat);
  kill4_vol.setAttributes(description,fX_kill4.regionStr(),fX_kill4.limitsStr(),fX_kill4.visStr());
  if ( fX_kill4.isSensitive() ) {
    kill4_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << kill4_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_kill4.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_kill4.isSensitive()) << endl;

  //  separate
  sol = Box(hwidth,hwidth,fX_separate2.thickness()/2.);
  mat = description.material(fX_separate2.materialStr());
  Volume separate2_vol(fX_separate2.nameStr(), sol, mat);
  separate2_vol.setAttributes(description,fX_separate2.regionStr(),fX_separate2.limitsStr(),fX_separate2.visStr());
  if ( fX_separate2.isSensitive() ) {
    separate2_vol.setSensitiveDetector(sens);
  }
  cout << setw(28) << left << separate2_vol.name()
       << " mat: "   << setw(15) << left << mat.name()
       << " vis: "   << setw(15) << left<< fX_separate2.visStr()
       << " solid: " << setw(20) << left << sol.type()
       << " sensitive: " << yes_no(fX_separate2.isSensitive()) << endl;


  //*********************************************************************************
  // unit cell volume


  float cellthick=fX_absorb.thickness()+fX_kill1.thickness()+fX_scint.thickness()+fX_kill2.thickness()+fX_separate1.thickness()+fX_kill3.thickness()+fX_ceren.thickness()+fX_kill4.thickness()+fX_separate2.thickness()+agap;
  float towerthick=Nlayers*cellthick;
  sol = Box(hwidth+tol,hwidth+tol,cellthick/2.+tol);
  std::cout<<"cell thick tower thick are "<<cellthick<<" "<<towerthick<<std::endl;
  mat = description.material(fX_unitbox.materialStr());
  Volume unit_box_vol( "unit_box", sol,mat);
  unit_box_vol.setSensitiveDetector(sens);
  unit_box_vol.setAttributes(description, fX_unitbox.regionStr(), fX_unitbox.limitsStr(), fX_unitbox.visStr());
  if ( fX_unitbox.isSensitive() ) {
    unit_box_vol.setSensitiveDetector(sens);
  }


  sol = Box(hwidth+2*tol,hwidth+2*tol,cellthick/2+2.*tol);
  mat = description.material(fX_unitbox.materialStr());
  Volume unit_box_vol2( "unit_box2", sol,mat);
  unit_box_vol2.setSensitiveDetector(sens);
  unit_box_vol2.setAttributes(description, fX_unitbox.regionStr(), fX_unitbox.limitsStr(), fX_unitbox.visStr());
  if ( fX_unitbox.isSensitive() ) {
    unit_box_vol2.setSensitiveDetector(sens);
  }


  //absorber
  float ahere = -cellthick/2.+agap/2.+fX_absorb.thickness()/2.;
  trafo = Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv    = unit_box_vol.placeVolume(absorb_vol, trafo);
  pv.addPhysVolID("type",1); 
  //kill1
  ahere+=fX_absorb.thickness()/2.+fX_kill1.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(kill1_vol, trafo);
  pv.addPhysVolID("type",2);
  //scintillator
  ahere+=fX_kill1.thickness()/2.+fX_scint.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  PlacedVolume pv_scint = unit_box_vol.placeVolume(scint_vol, trafo);
  pv_scint.addPhysVolID("type",3);


  //kill2
  ahere+=fX_scint.thickness()/2.+fX_kill2.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(kill2_vol, trafo);
  pv.addPhysVolID("type",4);

  //separate
  ahere+=fX_kill2.thickness()/2.+fX_separate1.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(separate1_vol, trafo);
  pv.addPhysVolID("type",8);


  //kill 3
  ahere+=fX_separate1.thickness()/2.+fX_kill3.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(kill3_vol, trafo);
  //quartz
  pv.addPhysVolID("type",5);
  ahere+=fX_kill3.thickness()/2.+fX_ceren.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(ceren_vol, trafo);
  pv.addPhysVolID("type",6);
  //kill4
  ahere+=fX_ceren.thickness()/2.+fX_kill4.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(kill4_vol, trafo);
  pv.addPhysVolID("type",7);


  //separate2
  ahere+=fX_kill4.thickness()/2.+fX_separate2.thickness()/2.;
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,ahere));
  pv = unit_box_vol.placeVolume(separate2_vol, trafo);
  pv.addPhysVolID("type",9);


  //kluge to use border surface
  trafo =  Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,0.));
  pv = unit_box_vol2.placeVolume(unit_box_vol, trafo);
  pv.addPhysVolID("box2",1);
  // this also works but need to think where to put photodetector
  //BorderSurface mirror = BorderSurface(description,sdet,"HallCrys", cryS,pv_scint,pv);
  //mirror.isValid();




  // tower
  Box    tower_box(hwidth+3.*tol+agap/2., hwidth+3.*tol+agap/2., towerthick/2+3.*tol);
  mat = description.material(fX_unitbox.materialStr());
  Volume tower_vol("tower", tower_box, mat);
  tower_vol.setVisAttributes(description, fX_tower.visStr());
  tower_vol.setSensitiveDetector(sens);
  for (int i=0;i<Nlayers;i++) {
    double zppp=-towerthick/2.+cellthick/2.+i*cellthick;
    std::cout<<" placing cell in tower at z of "<<zppp<<std::endl;
    pv=tower_vol.placeVolume(unit_box_vol2,Position(0.,0.,zppp));
    pv.addPhysVolID("layer",i+1);
  }

  // setup the volumes with the shapes and properties in one horixontal layer
  double dx = 2*(Ncount + Ncount+1)/2e0 * (hwidth+agap/2.) + tol;
  double dy = hwidth + agap/2. + tol;
  double dz = towerthick/2. + 2*tol;
  Box    tube_row_box(dx, dy, dz);
  Volume tube_row_vol("row", tube_row_box, air);
  tube_row_vol.setVisAttributes(description, fX_tower.visStr());
  tube_row_vol.setSensitiveDetector(sens);
  cout << tube_row_vol.name()
       << " dx: " << tube_row_box.x()
       << " dy: " << tube_row_box.y()
       << " dz: " << tube_row_box.z() << endl;
  tube_row_vol.setVisAttributes(description, "layerVis");


  // set towers in a row
  for (int ijk=0; ijk<2*Ncount+1; ijk++) {
    double mod_x_off = -dx+ (hwidth + agap/2.) + ijk * 2.*(hwidth+agap/2.);
    int towernum = ijk + 1;
    pv = tube_row_vol.placeVolume(tower_vol, Position(mod_x_off,0.,0.));
    pv.addPhysVolID("ix", towernum);
    //Box bounding_box = pv.volume().solid().GetBoundingBox();
    cout << "Placing row  "    << setw(5) << right << ijk
         << " x-offset: "      << setw(7) << right << mod_x_off
         << " volume of type " << pv.volume().name()
         << endl;
  }

  dx = 2*(Ncount + Ncount+1)/2e0 * (hwidth+agap/2.) + tol;
  dy = 2*(Ncount + Ncount+1)/2e0 * (hwidth+agap/2.) + tol;



  Box           env_box   (dx+tol, dy+tol, dz+tol);
  Volume        envelopeVol  (det_name, env_box, air);
  envelopeVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  // Now stack multiple horizontal layers to form the final box
  for (int ijk=0; ijk<2*Ncount+1; ijk++) {
    double mod_y_off = -dy+ (hwidth + agap/2.) + ijk * 2.*(hwidth+agap/2.);
    Transform3D tr(RotationZYX(0.,0.,0.), Position(0.,mod_y_off,0.));
    pv = envelopeVol.placeVolume(tube_row_vol, tr);
    pv.addPhysVolID("iy", ijk+1);

    DetElement de_layer(_toString(Ncount+ijk, "layer_%d"), det_id);
    de_layer.setPlacement(pv);
    sdet.add(de_layer);
    cout << "Placing " << setw(28) << left << de_layer.name()
         << " y-offset: "  << setw(7) << right << mod_y_off
         << " volume of type " << pv.volume().name()
         << endl;
  }

  // detector element for entire detector.  
  Volume        motherVol = description.pickMotherVolume(sdet);
  pv = motherVol.placeVolume(envelopeVol, Position(0.,0.,(zoffset+towerthick+tol)));
  pv.addPhysVolID("system", det_id);
  sdet.setPlacement(pv);  // associate the placed volume to the detector element


  return sdet;
}

DECLARE_DETELEMENT(DRSamp,create_detector)
