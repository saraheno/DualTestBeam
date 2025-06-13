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
#include "DD4hep/Printout.h"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {

  std::cout<<"Creating DRCrys "<<std::endl;
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
  xml_comp_t    x_towers  = x_det.staves();

  // overall detector offset
  xml_comp_t    x_dimg     = x_det.dimensions();
  double zoffset = x_dimg.z1();
  std::cout<<" z offset is "<<zoffset<<std::endl;

  // honeycomb thickness
  xml_comp_t fX_struct( x_det.child( _Unicode(structure) ) );
  xml_comp_t fX_honey(  fX_struct.child( _Unicode(honey) ) );
  double hthick = fX_honey.thickness();
  std::cout<<"honeycomb thickness is "<<hthick<<std::endl;

  // calculate size of entire detector and create envelop
  Layering      layering (e);
  double detectorhthickness=0.;
  double detectorhwidth=0.;
  int l_num=0;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"DRCrys layer "<<l_num<<std::endl;
    xml_comp_t x_layer = li;
    xml_comp_t x_dim = x_layer.child(_U(dimensions));
    int Ncount = x_dim.repeat();
    double hwidth   = x_dim.width()/2.;
    double agap=x_dim.gap()/2.;
    int repeat = x_layer.repeat();
    double hthickness=repeat*layering.layer(l_num)->thickness()/2.;
    double hnwidth = (2*Ncount+1)*(hwidth+agap);
    detectorhthickness+=hthickness;
    if(hnwidth>detectorhwidth) detectorhwidth=hnwidth;
    std::cout<<" ncount hwidth repeat hthickness "<<Ncount<<" "<<hnwidth<<" "<<repeat<<" "<<hthickness<<std::endl;
    l_num++;
  }
  DetElement    sdet      (det_name, det_id);
  Box           env_box   (detectorhwidth+tol, detectorhwidth+tol, detectorhthickness+tol);
  Volume        envelopeVol  (det_name, env_box, air);
  envelopeVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  // put the whole thing into the mother volume
  std::cout<<"putting calorimeter in mother volume"<<std::endl;

  Volume        motherVol = description.pickMotherVolume(sdet);
  PlacedVolume env_phv  = motherVol.placeVolume(envelopeVol, Position(0.,0.,zoffset+detectorhthickness));
  env_phv.addPhysVolID("system", det_id);
  sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
  sens.setType("calorimeter");
  std::cout<<" done"<<std::endl;


  // crystal wrappings
  std::cout<<"setting up optical surfaces"<<std::endl;
  OpticalSurfaceManager surfMgr = description.surfaceManager();
  OpticalSurface cryS  = surfMgr.opticalSurface("/world/"+det_name+"#mirrorSurface");
  std::cout<<" done"<<std::endl;

  // Loop over the sets of layer elements in the detector.

  std::cout<<" starting to build layers "<<std::endl;
  int opt_num=0;
  l_num = 0;
  double z_bottoml=0.;
  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    std::cout<<"DRCrys layer "<<l_num<<std::endl;
    xml_comp_t x_layer = li;
    xml_comp_t x_dim = x_layer.child(_U(dimensions));

    int Ncount = x_dim.repeat();
    double hwidth   = x_dim.width()/2.;
    double agap=x_dim.gap()/2.;
    int repeat = x_layer.repeat();  // how many times slice pattern repeats in layer
    double hthickness=repeat*layering.layer(l_num)->thickness()/2.;
    std::cout<<" ncount hwidth repeat hthickness "<<Ncount<<" "<<hwidth<<" "<<repeat<<" "<<hthickness<<std::endl;
    if(l_num<1) z_bottoml= -hthickness;

    // make a layer box volume and a tower volume
    dd4hep::Box LayerBox(detectorhwidth,detectorhwidth,hthickness+tol);
    string lbox_name = _toString(l_num,"layerbox%d");
    dd4hep::Volume LayerBoxVol(lbox_name, LayerBox, air);
    LayerBoxVol.setAttributes(description,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());
    if ( x_layer.isSensitive() ) {
      LayerBoxVol.setSensitiveDetector(sens);
    }
    cout << setw(28) << left << LayerBoxVol.name()
	 << " mat: "   << setw(15) << left << air.name()
	 << " vis: "   << setw(15) << left<< x_layer.visStr()
	 << " solid: " << setw(20) << left << LayerBox.type()
	 << " sensitive: " << yes_no(x_layer.isSensitive()) << endl;


    // tower which has layers, each of which has slices.
    // a grid of these will be placed in the layer box for this set of layers
    dd4hep::Box towertrap(hwidth+agap,hwidth+agap,hthickness+tol);
    dd4hep::Volume towerVol( "tower", towertrap, air);
    towerVol.setAttributes(description,x_towers.regionStr(),x_towers.limitsStr(),x_towers.visStr());
    if ( x_towers.isSensitive() ) {
      towerVol.setSensitiveDetector(sens);
    }
    cout << setw(28) << left << towerVol.name()
	 << " mat: "   << setw(15) << left << air.name()
	 << " vis: "   << setw(15) << left<< x_towers.visStr()
	 << " solid: " << setw(20) << left << towertrap.type()
	 << " sensitive: " << yes_no(x_towers.isSensitive()) << endl;

    // place the honeycomb into the tower
    Position b_pos(0.,0.,0.);
    dd4hep::Box abox1   (hwidth+0.5*(agap-hthick)+hthick,hwidth+0.5*(agap-hthick)+hthick,hthickness);
    dd4hep::Box abox2   (hwidth+0.5*(agap-hthick),hwidth+0.5*(agap-hthick),hthickness);
    dd4hep::Solid tmps = dd4hep::SubtractionSolid(abox1,abox2,b_pos);
    Volume  honeycomb  (det_name,tmps,description.material(fX_honey.materialStr()));
    honeycomb.setAttributes(description, fX_honey.regionStr(), fX_honey.limitsStr(), fX_honey.visStr());
    PlacedVolume honeycomb_phv = towerVol.placeVolume(honeycomb,b_pos);
    honeycomb_phv.addPhysVolID("wc", 0);

    // Loop over number of repeats for this layer.

    for (int j=0; j<repeat; j++)    {
      std::cout<<"  DRCrys layer "<<l_num<<" repeat "<<j<<std::endl;
      string l_name = _toString(j,"layer%d");
      double l_hzthick = layering.layer(l_num)->thickness()/2.;  // Layer's thickness.
      std::cout<<"  half  thickness is "<<l_hzthick<<std::endl;
      dd4hep::Box l_box(hwidth,hwidth,l_hzthick);
      dd4hep::Volume     l_vol(l_name,l_box,air);
      std::cout<<" layer visstr is "<<x_layer.visStr()<<std::endl;
      l_vol.setAttributes(description,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());


      // Loop over the sublayers or slices for this layer.
      int s_num = 1;
      double z_bottoms2=-l_hzthick;
      for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
	std::cout<<" with slice "<<s_num<<std::endl;
	xml_comp_t x_slice = si;
	string     s_name  = _toString(s_num,"slice%d");
	double     s_hzthick = x_slice.thickness()/2.;
	std::cout<<" with half  thickness "<<s_hzthick<<" and material "<<x_slice.materialStr()<<std::endl;
	dd4hep::Box s_box(hwidth,hwidth,s_hzthick);
	dd4hep::Volume     s_vol(s_name,s_box,description.material(x_slice.materialStr()));
	if ( x_slice.isSensitive() ) {
	  s_vol.setSensitiveDetector(sens);
	}
	std::cout<<"          slice visstr is "<<x_slice.visStr()<<std::endl;
	s_vol.setAttributes(description,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	// Slice placement.
	double z_mids2 = z_bottoms2+s_hzthick;
	Position   s_pos(0.,0.,z_mids2);      // Position of the layer.
	std::cout<<" placed at "<<z_mids2<<std::endl;
	PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);
	slice_phv.addPhysVolID("slice", s_num);
	// Increment Z position of slice.
	z_bottoms2 += 2.*s_hzthick;
	// Increment slice number.
	++s_num;
      }  // end of loop over slices

      // place the layer into the tower
      // Set region, limitset, and vis of layer.
      double z_midl=0.;
      Position   l_pos(0.,0.,z_midl);      // Position of the layer.
      std::cout<<" placed at z of "<<z_midl<<std::endl;
      PlacedVolume layer_phv = towerVol.placeVolume(l_vol,l_pos);
      layer_phv.addPhysVolID("wc", j+1);
      string tt_name = _toString(opt_num,"HallCrys%d");
      BorderSurface haha = BorderSurface(description,sdet, tt_name, cryS, layer_phv,env_phv);
      haha.isValid();
      opt_num++;

    }  //end of repeat for this layer

    //place towers into a row
    double dx = 2*(Ncount + Ncount+1)/2e0 * (hwidth+agap) + tol;
    double dy = hwidth + tol;
    double dz = hthickness + tol;
    Box    tube_row_box(dx, dy, dz);
    Volume tube_row_vol("layer", tube_row_box, air);
    tube_row_vol.setVisAttributes(description, x_det.visStr());
    tube_row_vol.setSensitiveDetector(sens);
    cout << tube_row_vol.name()
	 << " dx: " << tube_row_box.x()
	 << " dy: " << tube_row_box.y()
	 << " dz: " << tube_row_box.z() << endl;
    tube_row_vol.setVisAttributes(description, "layerVis");


    for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
	double mod_x_off = (ijk1)*2*(hwidth+tol+agap);
	std::cout<<"placing crystal at ("<<mod_x_off<<")"<<std::endl;
	trafo= Transform3D(RotationZYX(0.,0.,0.),Position(mod_x_off,0.,0.));
	pv = tube_row_vol.placeVolume(towerVol,trafo);
	int towernum = Ncount + ijk1 + 1;
	pv.addPhysVolID("ix",towernum);
	std::cout<<"placing tower "<<towernum<<std::endl;
	//BorderSurface haha = BorderSurface(description,sdet, tt_name, cryS, pv,env_phv);
	//haha.isValid();
    }  //end of placing towers in layer envelope


    //place the rows into the layer box
    for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
	double mod_y_off = (ijk2)*2*(hwidth+tol+agap);
	std::cout<<"placing crystal at ("<<mod_y_off<<")"<<std::endl;
	trafo= Transform3D(RotationZYX(0.,0.,0.),Position(0.,mod_y_off,0.));
	pv = LayerBoxVol.placeVolume(tube_row_vol,trafo);
	int towernum = Ncount + ijk2 + 1;
	pv.addPhysVolID("iy",towernum);
	std::cout<<"placing tower "<<towernum<<std::endl;
	//BorderSurface haha = BorderSurface(description,sdet, tt_name, cryS, pv,env_phv);
	//haha.isValid();
    }  //end of placing towers in layer envelope

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

  std::cout<<"exiting DRCrys creator"<<std::endl;
  return sdet;
}

DECLARE_DETELEMENT(DRCrys,create_detector)
