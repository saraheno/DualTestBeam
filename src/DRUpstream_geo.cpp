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


  // pointer to finding dimensins text in xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_comp_t    x_towers  = x_det.staves();
  xml_comp_t    x_dim     = x_det.dimensions();
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_bbb(  fX_struct.child( _Unicode(bbb) ) );

  double hwidth   = x_dim.width()/2.;
  double zoffset = x_dim.z1();
  std::cout<<"half width offset are "<<hwidth<<" "<<zoffset<<std::endl;

  // bbb thickness
  double hthick = fX_bbb.thickness();
  std::cout<<"bbb thickness is "<<hthick<<std::endl;

  // detector element for entire detector.  
  DetElement    sdet      (det_name,det_id);
  Volume        motherVol = description.pickMotherVolume(sdet);


  //PolyhedraRegular hedra  (nphi,inner_r,outer_r+tol*2e0,zmaxt);
  //set containment area for whole calorimeter
  Box abox  (hwidth+tol,hwidth+tol,fabs(zoffset));
  Volume        envelope  (det_name,abox,air);
  Position a_pos(0.,0.,zoffset);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope,a_pos);

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
  sens.setType("calorimeter");


  string t_name1 = _toString(0,"UpSp%d") ;
  DetElement UPS_det(t_name1,det_id);  // detector element for a tower

  
  // place bbb into the volume

  Position b_pos(0.,0.,0.);

  // an air box
  dd4hep::Box abox1   (hwidth,hwidth,hthick);
  Volume bbb_vol  (det_name,abox1,description.material(fX_bbb.materialStr()));
  bbb_vol.setVisAttributes(description,fX_bbb.visStr());
  DetElement shell(UPS_det, "bbb", det_id);
  PlacedVolume bbb_phv = env_vol.placeVolume(bbb_vol,b_pos);
  bbb_phv.addPhysVolID("layer", 0);
  shell.setPlacement(bbb_phv);

    // Loop over the sets of layer elements in the detector.
  double z_bottoml  = -zoffset;
  int l_num = 1;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"DRUpS layer (layers contain slices of material)"<<l_num<<std::endl;
    xml_comp_t x_layer = li;
    int repeat = x_layer.repeat();
      // Loop over number of repeats for this layer.
    for (int j=0; j<repeat; j++)    {
      std::cout<<"DRUpS layer "<<li<<" repeat "<<j<<std::endl;
      string l_name = _toString(l_num,"layer%d");
      double l_hzthick = layering.layer(l_num-1)->thickness()/2.;  // Layer's thickness.
      std::cout<<"half  thickness is "<<l_hzthick<<std::endl;

	// find top and bottom lengths at this position and center
        // relative to tower bottom
      double z_topl=z_bottoml + 2.*l_hzthick;
      double z_midl=z_bottoml + l_hzthick;
      Position   l_pos(0.,0.,z_midl);      // Position of the layer.
      std::cout<<" placed at z of "<<z_midl<<std::endl;

      dd4hep::Box l_box(hwidth,hwidth,l_hzthick);
      dd4hep::Volume     l_vol(l_name,l_box,air);
      DetElement layer(UPS_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
      int s_num = 1;

      double z_bottoms2=-l_hzthick;  
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	std::cout<<" with slice "<<s_num<<std::endl;
	xml_comp_t x_slice = si;
	string     s_name  = _toString(s_num,"slice%d");
	double     s_hzthick = x_slice.thickness()/2.;
	std::cout<<" with half  thickness "<<s_hzthick<<" and material "<<x_slice.materialStr()<<std::endl;

	      // this is relative to tower bottom, not layer bottom

	double z_mids2 = z_bottoms2+s_hzthick;
	      

	Position   s_pos(0.,0.,z_mids2);      // Position of the layer.
	std::cout<<" placed at "<<z_mids2<<std::endl;
	dd4hep::Box s_box(hwidth,hwidth,s_hzthick);


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

      PlacedVolume layer_phv = env_vol.placeVolume(l_vol,l_pos);
      layer_phv.addPhysVolID("layer", l_num);
      layer.setPlacement(layer_phv);
        // Increment to next layer Z position.


      z_bottoml=z_bottoml+2.*l_hzthick;

      ++l_num;
    }
  }
  // Set envelope volume attributes.
      std::cout<<" envelope visstr is "<<x_det.visStr()<<std::endl;
    envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());




  std::cout<<"exiting DRUpS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DRUpS,create_detector)


