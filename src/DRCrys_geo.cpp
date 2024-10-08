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


using namespace dd4hep;
using namespace dd4hep::detail;
using namespace std;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens)  {
	cout << "Creating DRCrys " << endl;
	static double tol = 0.001;
	dd4hep::Material      air       = description.air(); // material to underly it all
	
	// shortcuts to get XML structures; from: DDCore/include/XML/XML.h 
	// 	xml_comp_t 	components from xml file
	// 	xml_det_t	detector parts in xml file
	// 	xml_coll_t	collection from xml file
	xml_det_t     x_det     = e;
	Layering      layering (e);
	
	// for volume tags in detector
	int           det_id    = x_det.id();
	string        det_name  = x_det.nameStr();

	cout << "   det_id is " << det_id << " det_name is .. "	<< det_name << endl;
	
	// pointer to finding dimensins text in xml file
	// look in DDCore/include/Parsers/detail/Dimension.h for arguments
	xml_comp_t    x_towers  = x_det.staves();
	xml_comp_t    x_dim     = x_det.dimensions();
	
	double 	hwidth   = x_dim.width()/2.;
	double 	hzmax    = x_dim.z_length()/2.;
	int 	Ncount   = x_dim.repeat();
	double 	zoffset  = x_dim.z1();
	double 	agap 	 = x_dim.gap();
	
	cout << "   half width = " << hwidth  <<" , zmax= " << hzmax;
	cout << " , array size = " << Ncount;
	cout << " , z offset  = "  << zoffset;
	cout << " , gap between array elements = " << agap << endl;
	
	OpticalSurfaceManager 	surfMgr = description.surfaceManager();
	OpticalSurface 		cryS  	= surfMgr.opticalSurface("/world/"+det_name+"#mirrorSurface");
	
	// detector element for entire detector.  
	dd4hep::DetElement    sdet      (det_name,det_id);
	
	//set containment area for whole calorimeter
	dd4hep::Box		abox  (( 2 * Ncount + 1 ) * ( hwidth + agap + tol) , (2 * Ncount + 1) * (hwidth + agap + tol), hzmax + tol);
	dd4hep::Box 		towertrap(hwidth+tol,hwidth+tol,hzmax+tol);

	dd4hep::Volume        	motherVol = description.pickMotherVolume(sdet);
	dd4hep::Volume          envelope  (det_name,abox,air);
        dd4hep::Volume          towerVol( "tower", towertrap, air);
        towerVol.setVisAttributes(description, x_towers.visStr());
        towerVol.setSensitiveDetector(sens);
        cout << "   tower visStr: " << x_towers.visStr() << endl;

	dd4hep::Position	a_pos(0.,0.,(hzmax+zoffset));

	dd4hep::PlacedVolume	env_phv   = motherVol.placeVolume(envelope,a_pos);
	env_phv.addPhysVolID("system",det_id);
	sdet.setPlacement(env_phv);  // associate the placed volume to the detector element
	sens.setType("calorimeter");
	
	// create towers of calorimeter
	// 	tower may have different patterns that repeat,
	//	e.g. there may be 10 layers with one thickness of Pb and scint and 20 with another set of thicknesses.
	//	each of these repeating patterns is a "layer". (so in this example, two "layers")
	//	within each layer there is a slice, e.g. Pb and scint are slices
	//	the assembled tower is a Stave

	// tower envelope
	int itower=0;
	string t_name1 = _toString(itower,"towerp%d") ;
	dd4hep::DetElement tower_det(t_name1,det_id);  // detector element for a tower
	
	// Loop over the sets of layer elements in the detector.
	double z_bottoml  = -hzmax;
	int l_num = 1;
	for(xml_coll_t li(x_det,_U(layer)); li; ++li)  { // _U User-defined suffix to define custom calo parts, e.g. layer, slice, ...
		xml_comp_t x_layer 	= li;
		int repeat 		= x_layer.repeat();
		cout << setw(15) << " DRCrys layer (with slices of material) num = " << l_num << " , meterial name = " << x_layer.visStr() <<endl;
		
		// Loop over number of repeats for this layer.
		for (int j=0; j<repeat; j++)    {
			cout << setw(25) << "   DRCrys layer " << li << " , repeat " << j;
			string l_name = _toString(l_num,"layer%d");
			double l_hzthick = layering.layer(l_num-1)->thickness()/2.;  // Layer's thickness.
			cout << "   half thickness = " << l_hzthick;
			
			// find top and bottom lengths at this position and center
			// relative to tower bottom
			double z_topl = z_bottoml + 2.*l_hzthick;
			double z_midl = z_bottoml + l_hzthick;
			
			dd4hep::Position   l_pos(0.,0.,z_midl);      // Position of the layer.
			cout<< ", z = " << z_midl << endl;
			
			dd4hep::Box		l_box(hwidth,hwidth,l_hzthick);
			dd4hep::Volume		l_vol(l_name,l_box,air);
			dd4hep::DetElement	layer(tower_det, l_name, det_id);
			
			// Loop over the sublayers or slices for this layer.
			int s_num = 1;
			double z_bottoms2=-l_hzthick; 
			
			for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
				cout << setw(30) << " slice = " << s_num;
				xml_comp_t x_slice = si;
				string     s_name  = _toString(s_num,"slice%d");
				double     s_hzthick = x_slice.thickness()/2.;
				cout << " , half  thickness = " << s_hzthick << setw(15) << " , material: " << x_slice.materialStr() << setw(15);
				
				// this is relative to tower bottom, not layer bottom
				double z_mids2 = z_bottoms2+s_hzthick;
				cout << " , placed at " << z_mids2 << setw(15);
				dd4hep::Position   	s_pos(0.,0.,z_mids2);      // Position of the layer.
				dd4hep::Box		s_box(hwidth,hwidth,s_hzthick);
				dd4hep::Volume		s_vol(s_name,s_box,description.material(x_slice.materialStr()));
				dd4hep::DetElement 	slice(layer,s_name,det_id);
				
				if ( x_slice.isSensitive() ) s_vol.setSensitiveDetector(sens);
				cout << " visStr: " << x_slice.visStr() << endl;
				slice.setAttributes(description,s_vol,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
				
				// Slice placement.
				dd4hep::PlacedVolume slice_phv = l_vol.placeVolume(s_vol,s_pos);
				slice_phv.addPhysVolID("slice", s_num);
				slice.setPlacement(slice_phv);

				z_bottoms2 += 2.*s_hzthick;// Increment Z position of slice.
				++s_num; // Increment slice number.
			}
			// place the layer; Set region, limitset, and vis of layer.
			cout << setw(25) << " layer visStr is " << x_layer.visStr() << endl;
			layer.setAttributes(description,l_vol,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());

			dd4hep::PlacedVolume layer_phv = towerVol.placeVolume(l_vol,l_pos);
			layer_phv.addPhysVolID("layer", l_num);
			layer.setPlacement(layer_phv);
			z_bottoml=z_bottoml+2.*l_hzthick;
			++l_num; // Increment to next layer Z position.
		}
	}
	// now that you put the layers and slices into the tower do the placement	
	int towernum=-1;
	for (int ijk1=-Ncount; ijk1<Ncount+1; ijk1++) {
		for (int ijk2=-Ncount; ijk2<Ncount+1; ijk2++) {
			double mod_x_off = (ijk1)*2*(hwidth+tol+agap);
			double mod_y_off = (ijk2)*2*(hwidth+tol+agap);
			cout << setw(35) << "placing crystal at (" << mod_x_off << "," << mod_y_off << ")" << setw(20);
			dd4hep::Transform3D tr(dd4hep::RotationZYX(0.,0.,0.),dd4hep::Position(mod_x_off,mod_y_off,0.));
			dd4hep::PlacedVolume pv = envelope.placeVolume(towerVol,tr);
			pv.addPhysVolID("system",det_id);
			pv.addPhysVolID("ix",ijk1);
			pv.addPhysVolID("iy",ijk2);

			towernum +=1;
			cout<< ", placing tower = " << towernum << endl;
			string t_name2 = _toString(towernum,"0%d");
			dd4hep::DetElement sd = tower_det.clone(t_name2,det_id);
			sd.setPlacement(pv);

			string tt_name = _toString(towernum,"HallCrys%d");
			dd4hep::BorderSurface bordersurf = BorderSurface(description,sdet, tt_name, cryS, pv,env_phv);
			bordersurf.isValid();    //      sdet.add(sd);
		}
	}
	// Set envelope volume attributes.
	cout << setw(10) << "   envelope visstr is " << x_det.visStr() << endl;
	envelope.setAttributes(description,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
	cout << "   exiting DRCrys creator" << endl;
	cout << "**********************************************************************************" << endl;
	return sdet;
}
DECLARE_DETELEMENT(DRCrys,create_detector)
