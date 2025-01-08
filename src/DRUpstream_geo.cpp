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
#include "DD4hep/Printout.h"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {

  std::cout<<"Creating DRUpS "<<std::endl;
  static double tol = 0.001;
  Material air = description.air();
  Material mat;
  Transform3D trafo;
  PlacedVolume pv;
  Solid sol;

  // get stuff from xml file
  // look in DDCore/include/Parsers/detail/Dimension.h for arguments
  xml_det_t     x_det     = e;
  int           det_id    = x_det.id();
  std::cout<<" det_id is "<<det_id<<std::endl;
  string        det_name  = x_det.nameStr();
  std::cout<<" det_name is .. "<<det_name<<std::endl;


  // overall detector offset
  xml_comp_t    x_dimg     = x_det.dimensions();
  double zoffset = x_dimg.z1();
  std::cout<<" z offset is "<<zoffset<<std::endl;

  
  // calculate size of entire detector and create envelop
  std::cout<<std::endl;
  std::cout<<"calculating detector thickness "<<std::endl;
  Layering      layering (e);
  double detectorhthickness=0.;
  double detectorhwidth=0.;
  int l_num=0;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"  DRCrys layer "<<l_num<<std::endl;
    xml_comp_t x_layer = li;
    xml_comp_t x_dim = x_layer.child(_U(dimensions));
    double hwidth   = x_dim.width()/2.;
    int repeat = x_layer.repeat();
    double hthickness=repeat*layering.layer(l_num)->thickness()/2.;
    detectorhthickness+=hthickness;
    if(hwidth>detectorhwidth) detectorhwidth=hwidth;
    std::cout<<"    hwidth hthickness repeat "<<hwidth<<" "<<hthickness<<" "<<repeat<<std::endl;
    l_num++;
  }
  std::cout<<"half width and thickness are "<<detectorhwidth<<" "<<detectorhthickness<<std::endl;
  DetElement    sdet      (det_name, det_id);
  Box           env_box   (detectorhwidth+tol, detectorhwidth+tol, detectorhthickness+tol);
  Volume        envelopeVol  (det_name, env_box, air);
  envelopeVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  std::cout<<"done"<<std::endl;
  std::cout<<std::endl;

  // put the whole thing into the mother volume
  std::cout<<"putting calorimeter in mother volume"<<std::endl;

  Volume        motherVol = description.pickMotherVolume(sdet);
  PlacedVolume env_phv  = motherVol.placeVolume(envelopeVol, Position(0.,0.,zoffset+detectorhthickness));
  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
  sens.setType("calorimeter");
  std::cout<<" done"<<std::endl;
  std::cout<<std::endl;
  

  // Loop over the sets of layer elements in the detector.

  std::cout<<" starting to build layers "<<std::endl;
  l_num = 0;
  double z_bottoml=0.;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"DRCrys layer "<<l_num<<std::endl;
    xml_comp_t x_layer = li;
    xml_comp_t x_dim = x_layer.child(_U(dimensions));

    double hwidth   = x_dim.width()/2.;
    int repeat = x_layer.repeat();  // how many times slice pattern repeats in layer
    double hthickness=repeat*layering.layer(l_num)->thickness()/2.;
    double hlayerthick=layering.layer(l_num)->thickness()/2.;
    std::cout<<"  hwidth repeat hthickness "<<hwidth<<" "<<repeat<<" "<<hthickness<<std::endl;
    if(l_num<1) z_bottoml= -hthickness;
    
    // make a layer box volume 
    dd4hep::Box LayerBox(detectorhwidth,detectorhwidth,hthickness+tol);
    string lbox_name = _toString(l_num,"layerbox%d");
    dd4hep::Volume LayerBoxVol(lbox_name, LayerBox, air);
    LayerBoxVol.setAttributes(description,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());
    if ( x_layer.isSensitive() ) {
      LayerBoxVol.setSensitiveDetector(sens);
    }
    cout <<"  "<< setw(28) << left << LayerBoxVol.name()
	 << " mat: "   << setw(15) << left << air.name()
	 << " vis: "   << setw(15) << left<< x_layer.visStr()
	 << " solid: " << setw(20) << left << LayerBox.type()
	 << " sensitive: " << yes_no(x_layer.isSensitive()) << endl;
    double zzz=-hthickness+hlayerthick;
    // Loop over number of repeats for this layer.
    for (int j=0; j<repeat; j++)    {
      std::cout<<"  DRCrys layer "<<l_num<<" repeat "<<j<<std::endl;
      string l_name = _toString(j,"layer%d");
      double l_hzthick = layering.layer(l_num)->thickness()/2.;  // Layer's thickness.
      std::cout<<"    half  thickness is "<<l_hzthick<<std::endl;
      dd4hep::Box l_box(hwidth,hwidth,l_hzthick);
      dd4hep::Volume     l_vol(l_name,l_box,air);
      std::cout<<"    layer visstr is "<<x_layer.visStr()<<std::endl;
      l_vol.setAttributes(description,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

      
        // Loop overf the sublayers or slices for this layer.
      std::cout<<"    starting slices"<<std::endl;
      int s_num = 1;
      double z_bottoms2=-l_hzthick;  
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	std::cout<<"    with slice "<<s_num<<std::endl;
	xml_comp_t x_slice = si;
	string     s_name  = _toString(s_num,"slice%d");
	double     s_hzthick = x_slice.thickness()/2.;
	std::cout<<"    with half  thickness "<<s_hzthick<<" and material "<<x_slice.materialStr()<<std::endl;
	dd4hep::Box s_box(hwidth,hwidth,s_hzthick);
	dd4hep::Volume     s_vol(s_name,s_box,description.material(x_slice.materialStr()));
	if ( x_slice.isSensitive() ) {
	  s_vol.setSensitiveDetector(sens);
	}
	std::cout<<"    slice visstr is "<<x_slice.visStr()<<std::endl;
	s_vol.setAttributes(description,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	// Slice placement.
	double z_mids2 = z_bottoms2+s_hzthick;
	Position   s_pos(0.,0.,z_mids2);      // Position of the layer.
	std::cout<<"    placed at "<<z_mids2<<std::endl;
	PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);
	slice_phv.addPhysVolID("slice", s_num);
	// Increment Z position of slice.
	z_bottoms2 += 2.*s_hzthick;
	// Increment slice number.
	++s_num;
      }  // end of loop over slices

      // place the layer into the layerbox
      // Set region, limitset, and vis of layer.

      Position   l_pos(0.,0.,zzz);      // Position of the layer.
      std::cout<<" placed at z of "<<zzz<<std::endl;
      PlacedVolume layer_phv = LayerBoxVol.placeVolume(l_vol,l_pos);
      layer_phv.addPhysVolID("repeat", j+1);
      zzz+=2*hlayerthick;
      
    }  //end of repeat for this layer
    
      
    // place layerbox in envelope
    trafo=Transform3D(RotationZYX(0.,0.,0.),Position(0.,0.,z_bottoml+hthickness));
    pv = envelopeVol.placeVolume(LayerBoxVol,trafo);
    pv.addPhysVolID("layer",l_num);
    DetElement de_layer(_toString(l_num, "layer_%d"), det_id);
    de_layer.setPlacement(pv);
    sdet.add(de_layer);
    cout << "Placing " << setw(28) << left << de_layer.name()
         << " z-offset: "  << setw(7) << right << z_bottoml+hthickness
         << " volume of type " << pv.volume().name()
         << endl;


    
    // Increment to next layer Z position.
    z_bottoml=z_bottoml+2.*hthickness;
    ++l_num;
  }  //end of loop over layers






  
  std::cout<<"exiting DRUpS creator"<<std::endl;

  return sdet;
}

DECLARE_DETELEMENT(DRUpS,create_detector)


