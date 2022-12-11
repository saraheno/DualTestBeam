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

  static double tol = 0.01;


  xml_det_t     x_det     = e;


  // for volume tags in detector
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is .. "<<det_name<<std::endl;




  // material to underly it all
  Material      air       = description.air();




  OpticalSurfaceManager surfMgr = description.surfaceManager();
  OpticalSurface cryS  = surfMgr.opticalSurface("/world/"+det_name+"#mirrorSurface");





  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_dim     = x_det.dimensions();
  double        thick   = x_dim.thickness();
  double        zlength   = x_dim.z_length();
  double        zextra   = x_dim.dz();
  double        zph = x_dim.z1();
  int Ncount  = x_dim.numsides();

  std::cout<<" thick "<<thick<<" zlength "<<zlength<<" Ncount "<<Ncount<<std::endl;
  std::cout<<" extra fiber length is "<<zextra<<std::endl;



  int agap = x_dim.gap();
  std::cout<<" gap is "<<agap<<std::endl;


  double azmin = x_dim.zmin();
  std::cout<<" placing at zmin "<<azmin<<std::endl;

  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_absorb( fX_struct.child( _Unicode(absorb) ) );
  xml_comp_t fX_core( fX_struct.child( _Unicode(core) ) );
  xml_comp_t fX_hole( fX_struct.child( _Unicode(hole) ) );
  xml_comp_t fX_phdet( fX_struct.child( _Unicode(phdet) ) );



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
  dd4hep::Box abox   ((2*Ncount+1)*(thick+agap+tol),(2*Ncount+1)*(thick+agap+tol), (zlength+zextra+zph+tol));
  Volume        envelope  (det_name,abox,air);
  std::cout<<"setting envelope visibility to "<<x_det.visStr()<<std::endl;
  envelope.setVisAttributes(description,x_det.visStr());
  Position a_pos(0.,0.,azmin+zlength+zextra+zph+tol);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);

  env_phv.addPhysVolID("system",det_id);

  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element

  sens.setType("calorimeter");



    // tower envelope
  dd4hep::Box towertrap(thick,thick,zlength+zextra+zph);
  dd4hep::Volume towerVol( "tower", towertrap, air);
  towerVol.setVisAttributes(description, x_det.visStr());
  towerVol.setSensitiveDetector(sens);

  int itower=0;
  string t_name = _toString(itower,"towerp%d");
  DetElement tower_det(t_name,det_id);  
 
  SkinSurface haha = SkinSurface(description,sdet, "HallCrys", cryS, towerVol);
  haha.isValid();




    //absorber
  dd4hep::Box absbox(thick,thick,zlength);
  dd4hep::Tube fiberhole = dd4hep::Tube(0.,fX_hole.rmax(),zlength);
  dd4hep::Solid tmpsolid = dd4hep::SubtractionSolid(absbox,fiberhole);
  dd4hep::Volume absVol( "towerAbsorber", tmpsolid, description.material(fX_absorb.materialStr()) );
  std::cout<<"    material is "<<fX_absorb.materialStr()<<std::endl;
  absVol.setAttributes(description,fX_absorb.regionStr(),fX_absorb.limitsStr(),fX_absorb.visStr());
    // hole for fiber

  dd4hep::Volume holeVol("hole", fiberhole, description.material(fX_hole.materialStr()));




  string a_name = _toString(itower,"absorberm%d");
  DetElement abs_det(tower_det,a_name,det_id);  // detector element for absorber
  if ( fX_absorb.isSensitive() ) {
    std::cout<<"setting DRFiber absorber sensitive "<<std::endl;
    absVol.setSensitiveDetector(sens);
  }





    // fiber
  dd4hep::Tube fiber = dd4hep::Tube(0.,fX_core.rmax(),(zlength+zextra));
  std::cout<<" making fiber from "<<fX_core.materialStr()<<std::endl;
  dd4hep::Volume fiberVol("fiber", fiber, description.material(fX_core.materialStr()));
  fiberVol.setAttributes(description,fX_core.regionStr(),fX_core.limitsStr(),fX_core.visStr());
    string f_name = _toString(itower,"fiberm%d");
  DetElement fiber_det(tower_det,f_name,det_id);  // detector element for absorber

  if ( fX_core.isSensitive() ) {
    std::cout<<"setting DRFiber fiber sensitive "<<std::endl;
    fiberVol.setSensitiveDetector(sens);
  }



    // photodeector
  dd4hep::Tube photod = dd4hep::Tube(0.,fX_phdet.rmax(),(zph));
  dd4hep::Volume photodVol("phdet", photod, description.material(fX_phdet.materialStr()));
  photodVol.setAttributes(description,fX_phdet.regionStr(),fX_phdet.limitsStr(),fX_phdet.visStr());
    string ph_name = _toString(itower,"phdetm%d");
  DetElement photod_det(tower_det,ph_name,det_id);  // detector element for absorber

  if ( fX_phdet.isSensitive() ) {
    std::cout<<"setting DRFiber photodetector sensitive "<<std::endl;
    photodVol.setSensitiveDetector(sens);
  }



  std::cout<<"placing DRFiber tower subvolumes"<<std::endl;
  Position tra(0.,0.,-(zextra+zph));
  Position tra3(0.,0.,-(zph));
  Position tra4(0.,0.,(zextra+zlength));

  //
  PlacedVolume fiber_phv = towerVol.placeVolume( fiberVol, tra3);
  fiber_phv.addPhysVolID("fiber",1);
  fiber_phv.addPhysVolID("abs",0);
  fiber_phv.addPhysVolID("phdet",0);
  fiber_det.setPlacement(fiber_phv);



  PlacedVolume photod_phv = towerVol.placeVolume( photodVol, tra4);
  photod_phv.addPhysVolID("fiber",0);
  photod_phv.addPhysVolID("abs",0);
  photod_phv.addPhysVolID("phdet",1);
  photod_det.setPlacement(fiber_phv);



  PlacedVolume abs_phv = towerVol.placeVolume( absVol, tra);
  abs_phv.addPhysVolID("fiber",0);
  abs_phv.addPhysVolID("abs",1);
  abs_phv.addPhysVolID("phdet",0);
  abs_det.setPlacement(abs_phv);


    
    for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
  //int ijk1=0;
  //int ijk2=0;
      double mod_x_off = (ijk1)*2*(thick+agap);
      double mod_y_off = (ijk2)*2*(thick+agap);
      std::cout<<"placing fiber cell "<<ijk1<<","<<ijk2<<" at "<<" ("<<mod_x_off<<","<<mod_y_off<<")"<<std::endl;


      Transform3D tr(RotationZYX(0.,0.,0.),Position(mod_x_off,mod_y_off,0.));
      PlacedVolume pv = envelope.placeVolume(towerVol,tr);

      pv.addPhysVolID("system",det_id);
      pv.addPhysVolID("ix",ijk1);
      pv.addPhysVolID("iy",ijk2);


      int towernum = (ijk1+2)*(2*Ncount+1)+(ijk2+2);
      std::cout<<"placing tower "<<towernum<<std::endl;
      string t_name2 = _toString(towernum,"0%d");
      DetElement sd = tower_det.clone(t_name2,det_id);

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
