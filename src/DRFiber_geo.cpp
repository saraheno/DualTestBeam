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
  std::cout<<"DANGER DANGER WILL ROBINSON:  if the beam is aimed directly at a fiber, you will get lots of punch through.  need to add some tilt to this geometry to prevent this"<<std::endl;

  static double tol = 0.0;


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
  double        hthick   = x_dim.thickness()/2.;
  double        hzlength   = x_dim.z_length()/2.;
  double        hzph = x_dim.z1();
  int Ncount  = x_dim.numsides();

  std::cout<<"half  thick "<<hthick<<" half zlength "<<hzlength<<" Ncount "<<Ncount<<std::endl;
  std::cout<<" half kill media length is "<<hzph<<std::endl;



  double agap = x_dim.gap();
  std::cout<<" gap is "<<agap<<std::endl;


  double azmin = x_dim.zmin();
  std::cout<<" placing at zmin "<<azmin<<std::endl;

  // these refer to different fields in the xml file for this detector
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_absorb( fX_struct.child( _Unicode(absorb) ) );
  xml_comp_t fX_core1( fX_struct.child( _Unicode(core1) ) );
  xml_comp_t fX_core2( fX_struct.child( _Unicode(core2) ) );
  xml_comp_t fX_hole( fX_struct.child( _Unicode(hole) ) );
  xml_comp_t fX_phdet1( fX_struct.child( _Unicode(phdet1) ) );
  xml_comp_t fX_phdet2( fX_struct.child( _Unicode(phdet2) ) );



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
  dd4hep::Box abox   ((2*Ncount+1)*(hthick+agap+tol),(2*Ncount+1)*(hthick+agap+tol), (hzlength+hzph+tol));
  Volume        envelope  (det_name,abox,air);
  std::cout<<"setting envelope visibility to "<<x_det.visStr()<<std::endl;
  envelope.setVisAttributes(description,x_det.visStr());
  Position a_pos(0.,0.,azmin+hzlength+hzph+tol);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);

  env_phv.addPhysVolID("system",det_id);

  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element

  sens.setType("calorimeter");



  int itower=0;

    // tower envelope for scint fiber
  dd4hep::Box tower1trap(hthick,hthick,hzlength+hzph);
  dd4hep::Volume tower1Vol( "tower1", tower1trap, air);
  tower1Vol.setVisAttributes(description, x_det.visStr());
  tower1Vol.setSensitiveDetector(sens);
  string t1_name = _toString(itower,"tower1p%d");
  DetElement tower1_det(t1_name,det_id);  
  //  SkinSurface haha1 = SkinSurface(description,sdet, "HallSCINT", cryS, tower1Vol);
  //haha1.isValid();

    // tower envelope for quartz fiber
  dd4hep::Box tower2trap(hthick,hthick,hzlength+hzph);
  dd4hep::Volume tower2Vol( "tower2", tower2trap, air);
  tower2Vol.setVisAttributes(description, x_det.visStr());
  tower2Vol.setSensitiveDetector(sens);
  string t2_name = _toString(itower,"tower2p%d");
  DetElement tower2_det(t2_name,det_id);  
  //SkinSurface haha2 = SkinSurface(description,sdet, "HallCeren", cryS, tower2Vol);
  //haha2.isValid();




    //absorber
  dd4hep::Box absbox(hthick,hthick,hzlength);
  dd4hep::Volume abs1Vol( "towerAbsorber", absbox, description.material(fX_absorb.materialStr()) );
  dd4hep::Volume abs2Vol( "towerAbsorber", absbox, description.material(fX_absorb.materialStr()) );
  std::cout<<"    material is "<<fX_absorb.materialStr()<<std::endl;
  abs1Vol.setAttributes(description,fX_absorb.regionStr(),fX_absorb.limitsStr(),fX_absorb.visStr());
  abs2Vol.setAttributes(description,fX_absorb.regionStr(),fX_absorb.limitsStr(),fX_absorb.visStr());
  if ( fX_absorb.isSensitive() ) {
    std::cout<<"setting DRFiber absorber sensitive "<<std::endl;
    abs1Vol.setSensitiveDetector(sens);
    abs2Vol.setSensitiveDetector(sens);
  }
  string a_name1 = _toString(itower,"absorber1m%d");
  DetElement abs1_det(tower1_det,a_name1,det_id);  // detector element for absorber
  string a_name2 = _toString(itower,"absorber2m%d");
  DetElement abs2_det(tower2_det,a_name2,det_id);  // detector element for absorber

    //hole is absorber
  dd4hep::Tube fiberhole = dd4hep::Tube(0.,fX_hole.rmax(),hzlength);
  dd4hep::Volume absh1Vol( "towerAbsorberHole", fiberhole, air);
  dd4hep::Volume absh2Vol( "towerAbsorberHole", fiberhole, air);
  std::cout<<"    material is hardcoded to air"<<std::endl;
  absh1Vol.setAttributes(description,fX_hole.regionStr(),fX_hole.limitsStr(),fX_hole.visStr());
  absh2Vol.setAttributes(description,fX_hole.regionStr(),fX_hole.limitsStr(),fX_hole.visStr());
  if ( fX_hole.isSensitive() ) {
    std::cout<<"setting DRFiber absorber hole sensitive "<<std::endl;
    absh1Vol.setSensitiveDetector(sens);
    absh2Vol.setSensitiveDetector(sens);
  }
  string a_name3= _toString(itower,"absorberhole1m%d");
  DetElement absh1_det(tower1_det,a_name3,det_id);  // detector element for absorber
  string a_name4 = _toString(itower,"absorberhole2m%d");
  DetElement absh2_det(tower2_det,a_name4,det_id);  // detector element for absorber


    // scint fiber
  dd4hep::Tube fiber1 = dd4hep::Tube(0.,fX_core1.rmax(),(hzlength));
  std::cout<<" making fiber1 from "<<fX_core1.materialStr()<<std::endl;
  dd4hep::Volume fiber1Vol("fiber1", fiber1, description.material(fX_core1.materialStr()));
  std::cout<<"    core1 material is "<<fX_core1.materialStr()<<std::endl;
  fiber1Vol.setAttributes(description,fX_core1.regionStr(),fX_core1.limitsStr(),fX_core1.visStr());
  if ( fX_core1.isSensitive() ) {
    std::cout<<"setting DRFiber fiber1 sensitive "<<std::endl;
    fiber1Vol.setSensitiveDetector(sens);
  }
  string f1_name = _toString(itower,"fiber1m%d");
  DetElement fiber1_det(tower1_det,f1_name,det_id);  



    // quartz fiber
  dd4hep::Tube fiber2 = dd4hep::Tube(0.,fX_core2.rmax(),(hzlength));
  std::cout<<" making fiber2 from "<<fX_core2.materialStr()<<std::endl;
  dd4hep::Volume fiber2Vol("fiber2", fiber2, description.material(fX_core2.materialStr()));
  std::cout<<"    core 2material is "<<fX_core2.materialStr()<<std::endl;
  fiber2Vol.setAttributes(description,fX_core2.regionStr(),fX_core2.limitsStr(),fX_core2.visStr());
  std::cout<<"fiber2 vis is "<<fX_core2.visStr()<<std::endl;
  if ( fX_core2.isSensitive() ) {
    std::cout<<"setting DRFiber fiber2 sensitive "<<std::endl;
    fiber2Vol.setSensitiveDetector(sens);
  }
  string f2_name = _toString(itower,"fiber2m%d");
  DetElement fiber2_det(tower2_det,f2_name,det_id);  



    // photodeectors
  dd4hep::Tube photod1 = dd4hep::Tube(0.,fX_phdet1.rmax(),(hzph));
  dd4hep::Volume photod1Vol("phdet1", photod1, description.material(fX_phdet1.materialStr()));
  photod1Vol.setAttributes(description,fX_phdet1.regionStr(),fX_phdet1.limitsStr(),fX_phdet1.visStr());
  std::cout<<"   ph 1 material is "<<fX_phdet1.materialStr()<<std::endl;
  if ( fX_phdet1.isSensitive() ) {
    std::cout<<"setting DRFiber photodetector1 sensitive "<<std::endl;
    photod1Vol.setSensitiveDetector(sens);
  }

  dd4hep::Tube photod2 = dd4hep::Tube(0.,fX_phdet2.rmax(),(hzph));
  dd4hep::Volume photod2Vol("phdet2", photod2, description.material(fX_phdet2.materialStr()));
  photod2Vol.setAttributes(description,fX_phdet2.regionStr(),fX_phdet2.limitsStr(),fX_phdet2.visStr());
  std::cout<<"    ph2 material is "<<fX_phdet2.materialStr()<<std::endl;
  if ( fX_phdet2.isSensitive() ) {
    std::cout<<"setting DRFiber photodetector2 sensitive "<<std::endl;
    photod2Vol.setSensitiveDetector(sens);
  }


  string ph1_name = _toString(itower,"phdet1m%d");
  DetElement photod1_det(tower1_det,ph1_name,det_id);  
  string ph2_name = _toString(itower,"phdet2m%d");
  DetElement photod2_det(tower2_det,ph2_name,det_id);  




  std::cout<<"placing DRFiber tower subvolumes"<<std::endl;
  Transform3D tra(RotationZYX(0.,0.,0.),Position(0.,0.,-(hzph)));
  Transform3D tra2(RotationZYX(0.,0.,0.),Position(0.,0.,0.));
  Transform3D tra3(RotationZYX(0.,0.,0.),Position(0.,0.,-(hzph)));
  Transform3D tra4(RotationZYX(0.,0.,0.),Position(0.,0.,(hzlength)));


  //  scint
  PlacedVolume abs1_phv = tower1Vol.placeVolume( abs1Vol, tra);
  abs1_phv.addPhysVolID("fiber",0);
  abs1_phv.addPhysVolID("abs",1);
  abs1_phv.addPhysVolID("phdet",0);
  abs1_det.setPlacement(abs1_phv);
  PlacedVolume absh1_phv = abs1Vol.placeVolume( absh1Vol, tra2);
  absh1_phv.addPhysVolID("fiber",0);
  absh1_phv.addPhysVolID("abs",2);
  absh1_phv.addPhysVolID("phdet",0);
  absh1_det.setPlacement(absh1_phv);
  PlacedVolume fiber1_phv = absh1Vol.placeVolume( fiber1Vol, tra2);
  fiber1_phv.addPhysVolID("fiber",1);
  fiber1_phv.addPhysVolID("abs",0);
  fiber1_phv.addPhysVolID("phdet",0);
  fiber1_det.setPlacement(fiber1_phv);
  PlacedVolume photod1_phv = tower1Vol.placeVolume( photod1Vol, tra4);
  photod1_phv.addPhysVolID("fiber",0);
  photod1_phv.addPhysVolID("abs",0);
  photod1_phv.addPhysVolID("phdet",1);
  photod1_det.setPlacement(photod1_phv);



  //  quartz
  PlacedVolume abs2_phv = tower2Vol.placeVolume( abs2Vol, tra);
  abs2_phv.addPhysVolID("fiber",0);
  abs2_phv.addPhysVolID("abs",1);
  abs2_phv.addPhysVolID("phdet",0);
  abs2_det.setPlacement(abs2_phv);
  PlacedVolume absh2_phv = abs2Vol.placeVolume( absh2Vol, tra2);
  absh2_phv.addPhysVolID("fiber",0);
  absh2_phv.addPhysVolID("abs",2);
  absh2_phv.addPhysVolID("phdet",0);
  absh2_det.setPlacement(absh2_phv);
  PlacedVolume fiber2_phv = absh2Vol.placeVolume( fiber2Vol, tra2);
  fiber2_phv.addPhysVolID("fiber",2);
  fiber2_phv.addPhysVolID("abs",0);
  fiber2_phv.addPhysVolID("phdet",0);
  fiber2_det.setPlacement(fiber2_phv);
  PlacedVolume photod2_phv = tower1Vol.placeVolume( photod2Vol, tra4);
  photod2_phv.addPhysVolID("fiber",0);
  photod2_phv.addPhysVolID("abs",0);
  photod2_phv.addPhysVolID("phdet",2);
  photod2_det.setPlacement(photod2_phv);


  for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
      // int ijk1=0;
      // int ijk2=0;
      double mod_x_off = (ijk1)*2*(hthick+agap);
      double mod_y_off = (ijk2)*2*(hthick+agap);
      std::cout<<"  placing fiber cell "<<ijk1<<","<<ijk2<<" at "<<" ("<<mod_x_off<<","<<mod_y_off<<")"<<std::endl;
      Transform3D tr(RotationZYX(0.,0.,0.),Position(mod_x_off,mod_y_off,0.));


      DetElement sd;
      PlacedVolume pv;
      int towernum = (ijk1+2)*(2*Ncount+1)+(ijk2+2);
      std::cout<<"giving it tower number "<<towernum<<std::endl;
      if(towernum%2==0) {
	pv = envelope.placeVolume(tower1Vol,tr);
	pv.addPhysVolID("system",det_id);
	pv.addPhysVolID("ix",ijk1);
	pv.addPhysVolID("iy",ijk2);
	string t_name3 = _toString(towernum,"0%d");
	std::cout<<"    scint fiber"<<std::endl;
	sd = tower1_det.clone(t_name3,det_id);
      } else {
	pv = envelope.placeVolume(tower2Vol,tr);
	pv.addPhysVolID("system",det_id);
	pv.addPhysVolID("ix",ijk1);
	pv.addPhysVolID("iy",ijk2);
	string t_name3 = _toString(towernum,"0%d");
	std::cout<<"    quartz fiber"<<std::endl;
	sd = tower2_det.clone(t_name3,det_id);
      }

      string tts_name = _toString(towernum,"HallTower%d");
      BorderSurface haha = BorderSurface(description,sdet, tts_name, cryS, pv,env_phv);
      haha.isValid();



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
