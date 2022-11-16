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

  std::cout<<"Creating DRFiber"<<std::endl;

  static double tol = 0.1;


  xml_det_t     x_det     = e;



  // material to underly it all
  Material      air       = description.air();


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is "<<det_name<<std::endl;


  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_dim     = x_det.dimensions();
  double        thick   = x_dim.thickness();
  double        zlength   = x_dim.z_length();
  double        zextra   = x_dim.dz();
  int Ncount  = x_dim.numsides();

  std::cout<<" thick "<<thick<<" zlength "<<zlength<<" Ncount "<<Ncount<<std::endl;
  std::cout<<" extra fiber length is "<<zextra<<std::endl;



  int agap = x_dim.gap();
  std::cout<<" gap is "<<agap<<std::endl;


  double azmin = x_dim.zmin();
  std::cout<<" placing at zmin "<<azmin<<std::endl;

  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_absorb( x_det.child( _Unicode(absorb) ) );
  xml_comp_t fX_core( fX_struct.child( _Unicode(core) ) );
  xml_comp_t fX_hole( fX_struct.child( _Unicode(hole) ) );



  // three structures, volumes, placedvolumes, and detelements
  // volumes need a setVisAttribute
  // DetElements. you can have volumes inherit attrivutesby setting them here
  //              instead of in the volumes
  // placed volumes need a physvolid, a setplacement in a detelement,
  //                and are created with a mother.placevolume
  
  // detector element for entire detector.  
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);


  //PolyhedraRegular hedra  (nphi,inner_r,outer_r+tol*2e0,zmaxt);
  dd4hep::Box abox   ((2*Ncount+1)*(thick+agap+zextra+tol),(2*Ncount+1)*(thick+agap+zextra+tol), zlength+tol);
  Volume        envelope  (det_name,abox,air);
  Position a_pos(0.,0.,azmin+zlength);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);

  env_phv.addPhysVolID("system",det_id);

  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element

  sens.setType("calorimeter");



    // tower envelope
  dd4hep::Box towertrap(thick,thick,zlength);
  dd4hep::Volume towerVol( "tower", towertrap, air);
  towerVol.setVisAttributes(description, x_det.visStr());


  int itower=0;
  string t_name = _toString(itower,"towerp%d");
  DetElement tower_det(t_name,det_id);  



    //passive
  dd4hep::Box absstrap(thick,thick,zlength);
  dd4hep::Volume absVol( "towerAbsorber", absstrap, description.material(fX_absorb.materialStr()) );
  std::cout<<"    material is "<<fX_absorb.materialStr()<<std::endl;
  absVol.setVisAttributes(description, fX_absorb.visStr());
  string a_name = _toString(itower,"absorberm%d");
  DetElement abs_det(tower_det,a_name,det_id);  // detector element for absorber



    // hole for fiber
  dd4hep::Tube fiberhole = dd4hep::Tube(0.,fX_hole.rmax(),zlength);
  dd4hep::Volume holeVol("hole", fiberhole, description.material(fX_hole.materialStr()));
  holeVol.setVisAttributes(description, fX_hole.visStr());
  string h_name = _toString(itower,"holem%d");
  DetElement hole_det(abs_det,h_name,det_id);  // detector element for absorber


    // fiber
  dd4hep::Tube fiber = dd4hep::Tube(0.,fX_core.rmax(),(zlength+zextra));
  std::cout<<" making fiber from "<<fX_core.materialStr()<<std::endl;
  dd4hep::Volume fiberVol("fiber", fiber, description.material(fX_core.materialStr()));
  std::cout<<"fX_core.isSensitive is "<<fX_core.isSensitive()<<std::endl;
  if ( fX_core.isSensitive() ) {
    std::cout<<"setting DRFiber fiber sensitive "<<std::endl;
    fiberVol.setSensitiveDetector(sens);
  }
  string f_name = _toString(itower,"fiberm%d");
  DetElement fiber_det(abs_det,f_name,det_id);  // detector element for absorber
  fiberVol.setAttributes(description,fX_core.regionStr(),fX_core.limitsStr(),fX_core.visStr());


  Position tra(0.,0.,zlength/2.);
  Position tra2(0.,0.,0.);
  PlacedVolume fiber_phv = holeVol.placeVolume( fiberVol, tra2);
  PlacedVolume hole_phv = absVol.placeVolume( holeVol, tra2);
  PlacedVolume abs_phv = towerVol.placeVolume( absVol, tra);


  fiber_phv.addPhysVolID("fiber",1);
    //hole_phv.addPhysVolID("hole",1);
    //abs_phv.addPhysVolID("abs",1);

  fiber_det.setPlacement(fiber_phv);
  hole_det.setPlacement(hole_phv);
  abs_det.setPlacement(abs_phv);


    
  for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
      double mod_x_off = (ijk1)*2*(thick+agap);
      double mod_y_off = (ijk2)*2*(thick+agap);
      std::cout<<"placing fiber cell at ("<<mod_x_off<<","<<mod_y_off<<")"<<std::endl;


      Transform3D tr(RotationZYX(0.,0.,0.),Position(mod_x_off,mod_y_off,0.));
      PlacedVolume pv = envelope.placeVolume(towerVol,tr);

      pv.addPhysVolID("system",det_id);
      pv.addPhysVolID("ix",ijk1);
      pv.addPhysVolID("iy",ijk2);


      int towernum = (ijk1+2)*(2*Ncount+1)+(ijk2+2);
      std::cout<<"placing tower "<<towernum<<std::endl;
      string t_name = _toString(towernum,"0%d");
      DetElement sd = tower_det.clone(t_name,det_id);

      sd.setPlacement(pv);
      //      sdet.add(sd);

    }
  }


    //          }



  // Set envelope volume attributes.
  envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());








  std::cout<<"exiting DRFiber creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRFiber,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRFiber,create_detector)
