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
  std::cout<<"Creating DRCrys"<<std::endl;

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


  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_towers  = x_det.staves();
  xml_comp_t    x_dim     = x_det.dimensions();



  double hwidth   = x_dim.width();
  double hzmax = x_dim.z_length();
  int Ncount = x_dim.repeat();


  std::cout<<"half width zmax are "<<hwidth<<" "<<hzmax<<std::endl;
  std::cout<<" repeat is "<<Ncount<<std::endl;



  int agap = x_dim.gap();
  std::cout<<" gap is "<<agap<<std::endl;



  OpticalSurfaceManager surfMgr = description.surfaceManager();
  OpticalSurface cryS  = surfMgr.opticalSurface("/world/"+det_name+"#mirrorSurface");



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
  //set containment area for whole calorimeter
  Box abox  ((2*Ncount+1)*(hwidth+agap+tol),(2*Ncount+1)*(hwidth+agap+tol), hzmax+tol);
  Volume        envelope  (det_name,abox,air);
  Position a_pos(0.,0.,hzmax);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);

  env_phv.addPhysVolID("system",det_id);

  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
  sens.setType("calorimeter");



  // create towers to put into the calorimeter

  // a tower may have several different patterns that repeat.
  //   for example, there may be 10 layers with one thickness of Pb and scint and 20 with anothe r set of thicknesses.
  //   each of these repeating things is a "layer". (so in this example, two "layers")
  //within a layer is a slice os the Pb and scint are slices
  //the assembled tower is a Stave


    // tower envelope
  dd4hep::Box towertrap(hwidth,hwidth,hzmax);
  dd4hep::Volume towerVol( "tower", towertrap, air);
  std::cout<<"   tower visstr is "<<x_towers.visStr()<<std::endl;
  towerVol.setVisAttributes(description, x_towers.visStr());
  int itower=0; 
  string t_name = _toString(itower,"towerp%d") ;
  DetElement tower_det(t_name,det_id);  // detector element for a tower
  towerVol.setSensitiveDetector(sens);




  //SkinSurface haha = SkinSurface(description,sdet, "HallCrys", cryS, towerVol);
  //haha.isValid();



    // Loop over the sets of layer elements in the detector.
  double z_bottoml  = -hzmax;
  int l_num = 1;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"DRCrys layer "<<li<<std::endl;
    xml_comp_t x_layer = li;
    int repeat = x_layer.repeat();
      // Loop over number of repeats for this layer.
    for (int j=0; j<repeat; j++)    {
      std::cout<<"DRCrys layer "<<li<<" repeat "<<j<<std::endl;
      string l_name = _toString(l_num,"layer%d");
      double l_hzthick = layering.layer(l_num-1)->thickness();  // Layer's thickness.
      std::cout<<"half  thickness is "<<l_hzthick<<std::endl;

	// find top and bottom lengths at this position and center
        // relative to tower bottom
      double z_topl=z_bottoml + 2.*l_hzthick;
      double z_midl=z_bottoml + l_hzthick;
      Position   l_pos(0.,0.,z_midl);      // Position of the layer.
      std::cout<<" placed at z of "<<z_midl<<std::endl;

      dd4hep::Box l_box(hwidth-tol,hwidth-tol,l_hzthick-tol);
      dd4hep::Volume     l_vol(l_name,l_box,air);
      DetElement layer(tower_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
      int s_num = 1;

      double z_bottoms2=-l_hzthick;  
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	std::cout<<" with slice "<<si<<std::endl;
	xml_comp_t x_slice = si;
	string     s_name  = _toString(s_num,"slice%d");
	double     s_hzthick = x_slice.thickness();
	std::cout<<" with half  thickness "<<s_hzthick<<std::endl;

	      // this is relative to tower bottom, not layer bottom

	double z_mids2 = z_bottoms2+s_hzthick;
	      

	Position   s_pos(0.,0.,z_mids2);      // Position of the layer.
	std::cout<<" placed at "<<z_mids2<<std::endl;
	dd4hep::Box s_box(hwidth-tol,hwidth-tol,s_hzthick);


	dd4hep::Volume     s_vol(s_name,s_box,description.material(x_slice.materialStr()));
	DetElement slice(layer,s_name,det_id);

	if ( x_slice.isSensitive() ) {
		s_vol.setSensitiveDetector(sens);
          }
	std::cout<<"          slice visstr is "<<x_slice.visStr()<<std::endl;
	slice.setAttributes(description,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());

          // Slice placement.
	PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);
	slice_phv.addPhysVolID("slice", s_num);
	slice.setPlacement(slice_phv);
          // Increment Z position of slice.

	z_bottoms2 += 2.*s_hzthick;

          // Increment slice number.
	      ++s_num;
      }


      // place the layer
        // Set region, limitset, and vis of layer.
      std::cout<<" layer visstr is "<<x_layer.visStr()<<std::endl;
      layer.setAttributes(description,l_vol,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

      PlacedVolume layer_phv = towerVol.placeVolume(l_vol,l_pos);
      layer_phv.addPhysVolID("layer", l_num);
      layer.setPlacement(layer_phv);
        // Increment to next layer Z position.


      z_bottoml=z_bottoml+2.*l_hzthick;

      ++l_num;
    }
  }
      
  // now that you put the layers and slices into the tower, place it  

    // now do the placement



  for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
      double mod_x_off = (ijk1)*2*(hwidth+agap);             
      double mod_y_off = (ijk2)*2*(hwidth+agap);
      std::cout<<"placing crystal at ("<<mod_x_off<<","<<mod_y_off<<")"<<std::endl;


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


      BorderSurface haha = BorderSurface(description,sdet, "HallCrys", cryS, pv,env_phv);
      haha.isValid();



  }
  }

  // Set envelope volume attributes.
      std::cout<<" envelope visstr is "<<x_det.visStr()<<std::endl;
    envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());




  std::cout<<"exiting DRCrys creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_DRCrys,create_detector)

DECLARE_DEPRECATED_DETELEMENT(DRCrys,create_detector)
