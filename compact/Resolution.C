// force update
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../src/DualCrysCalorimeterHit.h"

#include <vector>
#include <functional>
#include <map>
#include <algorithm>

#include <iostream>
#include <sstream> // for ostringstream
#include <string>
using namespace std;


// This file is modified to accound for two-segmented ECAL crystals
// to run in batch copy:
//     root -l -b -q 'Resolution.C(nevt,"electron_rootfile.root", "pion_rootfile.root","hcalonly_rootfile.root",energy,doecal,dohcal,hcaltype,doedge,gendet,"output.root","ECALleaf","HCALleaf")'
//
//     nevt				number of events to process, e.g. 50;
//
//     electron_rootfile.root		ddsim input with electron gun, e.g. out_DualTestBeam-dial_e-10_.root;
//     pion_rootfile.root		ddsim input with electron gun, e.g. out_DualTestBeam-dial_e-10_.root;
//     hcalonly_rootfile.root		e- filename for hcal calibration if both ecal + hcal, e.g. " ";
//
//     energy				particle gun energy, e.g. 20;
//
//     doecal & dohcal			include ECAL & HCAL;
//
//     hcaltype
//     						fiber=0; -->
//     						sampling=1;
//
//     doedge				include Edge detector;
//
//     gendet				photon in:
//						active media (e.g. ecal crystal) = 1;
//						photodetector                    = 2;
//						energy deposit                   = 3; -->
//						debug                            = 4;
//     output.root			output root file with all histograms --> to be used in ResvE.C
//     ECALleaf   			ECAL ttree in electron_rootfile.root and pion_rootfile.root, e.g. DRCNoSegment
//     HCALleaf   			HCAL ttree in electron_rootfile.root and pion_rootfile.root, e.g. DRFNoSegment



// IGNORE THIS!!!!!!!!!!!!!!
// LOOK AT THE DECLARATION JUST AT THE crystalana def!!!!!!!!!!!!!!!!!!!!
const int nsliceecal = 6;
std::string nameecalslice[nsliceecal] = {"air","PD1","crystal1","gap","crystal2","PD2"};
int SCECOUNT=5;
bool dodualcorr=1;
bool doplots=1;
float timecut=1000;
const int finenbin=40;
const float timemin=0.;
const float timemax=400.;
const float timemaxz=40.;
const float timebinsize=(timemax-timemin)/float(finenbin);
const float timebinsizez=(timemaxz-timemin)/float(finenbin);

// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
// const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry
typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


void SCEDraw1   (TCanvas* canv,const char* name,TH1F* h1, const char* outfile, bool logy);
void SCEDraw1tp (TCanvas* canv, const char* name, TProfile* h1,			const char* outfile);
void SCEDraw1_2D(TCanvas* canv, const char* name, TH2F* h1,			const char* outfile,	float eohS,float eohC);
void SCEDraw2   (TCanvas* canv, const char* name, TH1F* h1,TH1F* h2,		const char* outfile,	bool logy);
void SCEDraw3   (TCanvas* canv, const char* name, TH1F* h1,TH1F* h2,TH1F* h3,	const char* outfile,	bool logy);
//void setCanvas  (TCanvas* canv, const char* name, TH1F* h1);
void setCanvas  (TCanvas* canv,const char* name, bool hist, bool tprofile, TH1F* h1, TProfile* t1);


void getStuff(map<string,int> mapecalslice, map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal,
		int hcaltype,  bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge, CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
		float  &eesum,float &eesumcal,float &eesumem,float &eesumair,float &eesumdead, float &eesumcrystal,float &eesumPDe,
		float &eesumfiber1, float &eesumfiber2,float &eesumabs,float &eesumPDh,
		float &eesumedge,float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal, float &timecut, float &eecaltimecut,
		float &ehcaltimecut, TH1F* eecaltime,  TH1F* ehcaltime, int &nine, int &ninh,
		TH1F *ecalpd1scint,  TH1F *ecalpd1cer, TH1F *ecalpd2scint,TH1F *ecalpd2cer,
		TH1F *hcalpd1scint,TH1F *hcalpd1cer,TH1F *hcalpd2scint,TH1F *hcalpd2cer,
	      	TH1F *ecalpd1scintz, TH1F *ecalpd1cerz,TH1F *ecalpd2scintz,TH1F *ecalpd2cerz,
	        TH1F *hcalpd1scintz,TH1F *hcalpd1cerz,TH1F *hcalpd2scintz,TH1F *hcalpd2cerz,
	        TH2F *eecal2d, TH2F *ehcal2d);


void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaecal, float kappahcal,
		float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype,
		TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, float &EEcal, float &EHcal,
		float &timecut, float &eecaltimecut, float &ehcaltimecut); //

void getMeanPhot(map<string, int> mapecalslice, map<string, int> mapsampcalslice,  //input
		int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input
		TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, //input
		float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,float &timecut, float &eecaltimecut, float &ehcaltimecut); // output


map<string,int> mapecalslice;
map<string,int>::iterator eii0;
map<string,int>::iterator eii1;
map<string,int>::iterator eii2;
map<string,int>::iterator eii3;
map<string,int>::iterator eii4;
map<string,int>::iterator eii5;
map<string,int>::iterator eii6;
map<string,int>::iterator eii7;

map<string,int> mapsampcalslice;
map<string,int>::iterator sii1;
map<string,int>::iterator sii2;
map<string,int>::iterator sii3;
map<string,int>::iterator sii4;
map<string,int>::iterator sii5;
map<string,int>::iterator sii6;
map<string,int>::iterator sii7;
map<string,int>::iterator sii8;
map<string,int>::iterator sii9;

void Resolution(int num_evtsmax, const char* einputfilename, const char* piinputfilename, const char* hcalonlyefilename,
                const float beamEE, bool doecal, bool dohcal, int hcaltype, bool doedge, int gendet, const char* outputfilename,
                const char* ECALleaf, const char* HCALleaf){
        mapecalslice["air"]=0;
	mapecalslice["PD1"]=1;
	mapecalslice["crystal1"]=2;
	mapecalslice["gap1"]=3;
	mapecalslice["middlemat"]=4;
	mapecalslice["gap2"]=5;
	mapecalslice["crystal2"]=6;
	mapecalslice["PD2"]=7;

	eii0 = mapecalslice.find("air");
	eii1 = mapecalslice.find("PD1");
	eii2 = mapecalslice.find("crystal1");
	eii3 = mapecalslice.find("gap1");
	eii4 = mapecalslice.find("middlemat");
	eii5 = mapecalslice.find("gap2");
	eii6 = mapecalslice.find("crystal2");
	eii7 = mapecalslice.find("PD2");

	mapsampcalslice["air"]		=0;
	mapsampcalslice["Iron"]		=1;
	mapsampcalslice["PD1"]		=2;
	mapsampcalslice["PS"]		=3;
	mapsampcalslice["PD2"]		=4;
	mapsampcalslice["PD3"]		=5;
	mapsampcalslice["Quartz"]	=6;
	mapsampcalslice["PD4"]		=7;
	mapsampcalslice["Sep1"]		=8;
	mapsampcalslice["Sep2"]		=9;

	sii1 = mapsampcalslice.find("Iron");
	sii2 = mapsampcalslice.find("PD1");
	sii3 = mapsampcalslice.find("PS");
	sii4 = mapsampcalslice.find("PD2");
	sii5 = mapsampcalslice.find("PD3");
	sii6 = mapsampcalslice.find("Quartz");
	sii7 = mapsampcalslice.find("PD4");
	sii8 = mapsampcalslice.find("Sep1");
	sii9 = mapsampcalslice.find("Sep2");

	float beamE=beamEE*1000.;  // convert to MeV
	// read in libraries that define the classes
	Long_t result;
	char text[1024];
	const char* dd4hep = gSystem->Getenv("DD4hepINSTALL");
	snprintf(text,sizeof(text)," -I%s/include -D__DD4HEP_DDEVE_EXCLUSIVE__ -Wno-shadow -g -O0",dd4hep);
	gSystem->AddIncludePath(text);
	TString fname = "libDDG4IO";
	const char* io_lib = gSystem->FindDynamicLibrary(fname,kTRUE);
	result = gSystem->Load("libDDG4IO");
	result = gSystem->Load("libDDEvePlugins");
	result = gSystem->Load("libDDEvePlugins");
	result = gSystem->Load("libDualTestBeam");
	result = gSystem->Load("libDDG4Plugins");
	int num_evt_mc,num_evt;
	TBranch* b_mc;
	TBranch* b_ecal;
	TBranch* b_hcal;
	TBranch* b_edge;

	const int num_histograms = 6;
	vector<TH1F*> eenergy;
	vector<TH1F*> pienergy;
	vector<TH1F*> ephoton;
	vector<TH1F*> piphoton;
	//vector<const char*> energy_lowbin = {}
	vector<double> energy_hbin = {0.5,0.5,0.5,0.5,1.5,1.5};
	vector<const char*> energy_title = {"sum(crysE)/beamE","sum(fiberE)/beamE","sum(scintfiberE)/beamE","sum(cerfiberE)/beamE","sum(edgeE)/beamE","sum(beamE-edgeE)/beamE"};
	vector<const char*> nphoton_title = {"ecal: ntotcer/emeancer", "hcal: ntotcer/emeancer", "ecal: ntotscint/emeanscint", "hcal: ntotscint/emeanscint"};
	for (size_t i = 0; i < num_histograms; ++i) {
	  eenergy.push_back( new TH1F(Form("elenergy_%zu", i), Form("e-: %s", energy_title[i]), 100, 0.0, energy_hbin[i]));
	  pienergy.push_back(new TH1F(Form("pienergy_%zu", i), Form("pi-: %s",energy_title[i]), 100, 0.0, energy_hbin[i]));
	  if (i>3) { continue;}
	  ephoton.push_back(new  TH1F(Form("elnphoton_%zu",i), Form("e-: %s", nphoton_title[i]), 100, 0.0, 1.5));
	  piphoton.push_back(new TH1F(Form("pinphoton_%zu",i), Form("pi-: %s",nphoton_title[i]), 100, 0.0, 1.5));
	}

	TH1F *phcEcalcorr	= new TH1F("phcEcalcorr","ecal: dualCorrelation", 600,0.,1.5);
	TH1F *phcHcalcorr	= new TH1F("phcHcalcorr","hcal: dualCorrelation", 100,0.,1.5);

	TH2F *ehcEcalNsNc	= new TH2F("ehcEcalNsNc","e-  ecal ncer vs nscint",500,0.,1.5,500,0.,1.5);
	TH2F *phcEcalNsNc	= new TH2F("phcEcalNsNc","pi- ecal ncer vs nscint",500,0.,1.5,500,0.,1.5);
	TH2F *ehcHcalNsNc	= new TH2F("ehcHcalNsNc","e-  hcal ncer vs nscint",500,0.,1.5,500,0.,1.5);
	TH2F *phcHcalNsNc	= new TH2F("phcHcalNsNc","pi- hcal ncer vs nscint",500,0.,1.5,500,0.,1.5);

	TH2F *mes1Ecal          = new TH2F("mes1Ecal","Ecal missing energy versus em frac",500,0.,1.5,500,0.,1.5);
        TH2F *mes2Ecal          = new TH2F("mes2Ecal","Ecal missing energy versus num nucl int",500,0.,1.5,500,0.,1000);
        TH2F *mes3Ecal          = new TH2F("mes3Ecal","Ecal missing energy versus ncer",500,0.,1.5,500,0.,1.5);
	TH2F *mes4Ecal          = new TH2F("mes4Ecal","Ecal missing energy versus nscint",500,0.,1.5,500,0.,1.5);
	TH2F *mes1Hcal          = new TH2F("mes1Hcal","Hcal missing energy versus em frac",500,0.,1.5,500,0.,1.5);
	TH2F *mes2Hcal          = new TH2F("mes2Hcal","Hcal missing energy versus num nucl int",500,0.,1.5,500,0.,1000);
	TH2F *mes3Hcal          = new TH2F("mes3Hcal","Hcal missing energy versus ncer",500,0.,1.5,500,0.,1.5);
	TH2F *mes4Hcal          = new TH2F("mes4Hcal","Hcal missing energy versus nscint",500,0.,1.5,500,0.,1.5);

	TH2F *ehcEcalNsNctc     = new TH2F("ehcEcalNsNctc","ecal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);
	TH2F *phcEcalNsNctc     = new TH2F("phcEcalNsNctc","ecal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);

	TH2F *ehcEcalMarco	= new TH2F("ehcEcalMarco","e-  ecal c v c/s",500,0.,1.5,500,0.,1.5);
	TH2F *phcEcalMarco	= new TH2F("phcEcalMarco","pi- ecal c v c/s",500,0.,1.5,500,0.,1.5);
	TH2F *ehcHcalMarco	= new TH2F("ehcHcalMarco","e-  hcal c v c/s",500,0.,1.5,500,0.,1.5);
	TH2F *phcHcalMarco	= new TH2F("phcHcalMarco","pi- hcal c v c/s",500,0.,1.5,500,0.,1.5);

	TH2F *ehcHcalf1f2	= new TH2F("ehcHcalf1f2","e-  hcal-fiber scint vs cer",500,0.,2.5,500,0.,2.5);
	TH2F *phcHcalf1f2	= new TH2F("phcHcalf1f2","pi- hcal-fiber scint vs cer",500,0.,2.5,500,0.,2.5);

	TH1F *ehaphcal		= new TH1F("ehaphcal","e-  hcal-fiber (scint + cer)/ (tot_energy)" ,100,0.,0.2);
	TH1F *phaphcal		= new TH1F("phaphcal","pi- hcal-fiber (scint + cer)/ (tot_energy)" ,100,0.,0.2);

	TH1F *eheest		= new TH1F("eheest","e-  : (crys+fiber) / beamE",500,0.,1.25);
	TH1F *pheest		= new TH1F("pheest","pi- : (crys+fiber) / beamE",500,0.,1.25);

	TH1F *ehetrue		= new TH1F("ehetrue","e-  : alldepositedE / beamE",500,0.,1.25);
	TH1F *phetrue		= new TH1F("phetrue","pi- : alldepositedE / beamE",500,0.,1.25);

	TH1F *hedepcal		= new TH1F("hedepcal","e-  (alldepositE-edgeE) / beamE",100,0.6,1.1);
	TH1F *hpdepcal		= new TH1F("hpdepcal","pi- (alldepositE-edgeE) / beamE",100,-0.5,1.1);

	TH1F *eecaltime		= new TH1F("eecaltime","e- ecal: time",100,0.,40.);
	TH1F *ehcaltime		= new TH1F("ehcaltime","e- hcal: time",100,0.,40.);
	TH1F *piecaltime	= new TH1F("piecaltime","pi- ecal: time",100,0.,40.);
	TH1F *pihcaltime	= new TH1F("pihcaltime","pi- hcal: time",100,0.,40.);

	TH2F *enscvni		= new TH2F("enscvni", "e-  (alldepositE-edgeE) vs num_inelastic",500,0.,1.2,100,0.,500);
	TH2F *pinalvni          = new TH2F("pinalvni","pion all E versus number inelastic",100,0.,1.2,100,0.,1000);
	TH2F *pinscvni		= new TH2F("pinscvni","pi- (alldepositE-edgeE) vs num_inelastic",500,0.,1.2,100,0.,500);
	TH2F *pincevni          = new TH2F("pincevni","pion cherenkov  versus number inelastic",100,0.,1.2,100,0.,1000);

	TH1F *heesumcal 	= new TH1F("heesumcal",		"electron energy in calorimeter",		400,0.,1.5);
	TH1F *heesumemcal 	= new TH1F("heesumemcal",	"electron relativistic energy in calorimeter",	400,0.,1.5);
	TH1F *hefff 		= new TH1F("hefff",		"electron relativistic fraction in calorimeter",400,0.,1.5);
	TH1F *hpesumcal 	= new TH1F("hpesumcal",		"pion energy in calorimeter",			400,0.,1.5);
	TH1F *hpesumemcal 	= new TH1F("hpesumemcal",	"pion relativistic energy in calorimeter",	400,0.,1.5);
	TH1F *hpfff 		= new TH1F("hpfff",		"pion relativistic fraction in calorimeter",	400,0.,1.5);

	TH1F *eecalpd1scint     = new TH1F("eecalpd1scint","electron scint photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
	TH1F *eecalpd1cer       = new TH1F("eecalpd1cer","electron cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
	TH1F *pecalpd1scint     = new TH1F("pecalpd1scint","pion scint photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
	TH1F *pecalpd1cer       = new TH1F("pecalpd1cer","pion cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
	TH1F *eecalpd2scint     = new TH1F("eecalpd2scint","electron scint photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
	TH1F *eecalpd2cer       = new TH1F("eecalpd2cer","electron cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
	TH1F *pecalpd2scint     = new TH1F("pecalpd2scint","pion scint photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
	TH1F *pecalpd2cer       = new TH1F("pecalpd2cer","pion cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
	TH1F *ehcalpd1scint     = new TH1F("ehcalpd1scint","elec scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
	TH1F *ehcalpd1cer       = new TH1F("ehcalpd1cer","elec cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
	TH1F *phcalpd1scint     = new TH1F("phcalpd1scint","pion scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
	TH1F *phcalpd1cer       = new TH1F("phcalpd1cer","pion cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
	TH1F *ehcalpd2scint     = new TH1F("ehcalpd2scint","elec scint photon arrival time ns HCAL quartz fiber",finenbin,timemin,timemax);
	TH1F *ehcalpd2cer       = new TH1F("ehcalpd2cer","elec cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemax);
	TH1F *phcalpd2scint     = new TH1F("phcalpd2scint","pion scint photon arrival time ns quartz fiber",finenbin,timemin,timemax);
	TH1F *phcalpd2cer       = new TH1F("phcalpd2cer","pion cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemax);

	TH1F *eecalpd1scintz    = new TH1F("eecalpd1scintz","electron scint photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
	TH1F *eecalpd1cerz      = new TH1F("eecalpd1cerz","electron cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
	TH1F *pecalpd1scintz    = new TH1F("pecalpd1scintz","pion scint photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
	TH1F *pecalpd1cerz      = new TH1F("pecalpd1cerz","pion cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
	TH1F *eecalpd2scintz    = new TH1F("eecalpd2scintz","electron scint photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
	TH1F *eecalpd2cerz      = new TH1F("eecalpd2cerz","electron cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
	TH1F *pecalpd2scintz    = new TH1F("pecalpd2scintz","pion scint photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
	TH1F *pecalpd2cerz      = new TH1F("pecalpd2cerz","pion cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
	TH1F *ehcalpd1scintz    = new TH1F("ehcalpd1scintz","elec scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
	TH1F *ehcalpd1cerz      = new TH1F("ehcalpd1cerz","elec cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
	TH1F *phcalpd1scintz    = new TH1F("phcalpd1scintz","pion scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
	TH1F *phcalpd1cerz      = new TH1F("phcalpd1cerz","pion cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
	TH1F *ehcalpd2scintz    = new TH1F("ehcalpd2scintz","elec scint photon arrival time ns HCAL quartz fiber",finenbin,timemin,timemaxz);
	TH1F *ehcalpd2cerz      = new TH1F("ehcalpd2cerz","elec cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);
	TH1F *phcalpd2scintz    = new TH1F("phcalpd2scintz","pion scint photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);
	TH1F *phcalpd2cerz      = new TH1F("phcalpd2cerz","pion cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);

	TH2F *eecal2d = new TH2F("eecal2d","lego of ecal elec", 41,-20.,20.,41,-20.,20.);
	TH2F *ehcal2d = new TH2F("ehcal2d","lego of hcal elec", 41,-20.,20.,41,-20.,20.);

	TH1F *eleEcalncer = new TH1F("eleEcalncer","# ecal cerenkov / mean (ele)",  600,0.,1.5);
	TH1F *eleEcalnscint = new TH1F("eleEcalnscint","# ecal scintillation / mean (ele)", 600,0.,1.5);

	//****************************************************************************************************************************
	// process electrons

	TFile* ef = TFile::Open(einputfilename);
	TTree* et = (TTree*)ef->Get("EVENT;1");

	b_mc= et->GetBranch("MCParticles");
	if(doecal) b_ecal = et->GetBranch(ECALleaf);
	if(dohcal) b_hcal = et->GetBranch(HCALleaf);
	if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");

	num_evt_mc = b_mc->GetEntries();
	num_evt= min(num_evt_mc,num_evtsmax);
	cout<<"num_evt="<<num_evt<<", MC events="<<num_evt_mc<<endl;

	float meanscinEcal(0),meanscinHcal(0),meancerEcal(0),meancerHcal(0);
	float meaneecaltimecut(0),meanehcaltimecut(0);
	if(num_evt>0) { // loop over events
		CalHits* ecalhits = new CalHits();
		if(doecal) b_ecal->SetAddress(&ecalhits);
		CalHits* hcalhits = new CalHits();
		if(dohcal) b_hcal->SetAddress(&hcalhits);
		CalHits* edgehits = new CalHits();
		if(doedge) b_edge->SetAddress(&edgehits);
		for(int ievt=0;ievt<num_evt; ++ievt) {// first pass through file get hcal & ecal: meanCer and meanScint + timing info
			getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, doecal, dohcal, hcaltype, b_ecal,b_hcal, ecalhits, hcalhits, //input
					meanscinEcal, meanscinHcal, meancerEcal, meancerHcal,timecut, meaneecaltimecut, meanehcaltimecut); // output:
		}
		cout<<"done with getMeanPhot"<<endl;
		meanscinEcal=meanscinEcal/num_evt;
		meanscinHcal=meanscinHcal/num_evt;
		meancerEcal =meancerEcal/num_evt;
		meancerHcal =meancerHcal/num_evt;

		cout<<"ECAL: meanscint="<<meanscinEcal<<", meancer="<<meancerEcal<<endl;
		cout<<"HCAL: meanscint="<<meanscinHcal<<", meancer="<<meancerHcal<<endl;

		// second pass through file
		for(int ievt=0;ievt<num_evt; ++ievt) {
			float eesum(0), eesumcal(0.),eesumem(0.), eesumair(0), eesumdead(0), eesumcrystal(0), eesumPDe(0);
			float eesumfiber1(0), eesumfiber2(0), eesumabs(0), eesumPDh(0), eesumedge(0);
			float necertotecal(0), nescinttotecal(0), necertothcal(0), nescinttothcal(0), eecaltimecut(0), ehcaltimecut(0);
			int nine(0), ninh(0);
			getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,
					edgehits,eesum,eesumcal,eesumem,eesumair,eesumdead,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,
					eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal,timecut, eecaltimecut,
					ehcaltimecut,  eecaltime,ehcaltime,nine,ninh,
					eecalpd1scint, eecalpd1cer, eecalpd2scint, eecalpd2cer,
					ehcalpd1scint, ehcalpd1cer, ehcalpd2scint, ehcalpd2cer,
					eecalpd1scintz,eecalpd1cerz,eecalpd2scintz,eecalpd2cerz,
				        ehcalpd1scintz,ehcalpd1cerz,ehcalpd2scintz,ehcalpd2cerz,
				 	eecal2d, ehcal2d);

			vector<float> eesums   = {eesumcrystal, eesumfiber1 + eesumfiber2, eesumfiber1, eesumfiber2, eesumedge, beamE - eesumedge};
			vector<float> nphotons = {necertotecal/meancerEcal, necertothcal/meancerHcal, nescinttotecal/meanscinEcal, nescinttothcal/meanscinHcal};

			auto eit = eesums.begin();
			for (auto hist : eenergy) {
				hist->Fill(*eit/beamE);
				++eit;
			}
			auto nit = nphotons.begin();
			for (auto hist: ephoton) {
				hist->Fill(*nit);
				++nit;
			}

			heesumcal->Fill(eesumcal/beamE);
			heesumemcal->Fill(eesumem/beamE);
			if(eesumem>0) hefff->Fill(eesumem/eesumcal);
			else cout << "esumemcal is zero " << endl;

			ehcHcalf1f2->Fill(eesumfiber1/1000.,eesumfiber2/1000.);
			if((eesumfiber1+eesumfiber2)>0) ehaphcal->Fill((eesumfiber1+eesumfiber2)/(eesumabs+eesumfiber1+eesumfiber2));
			eheest->Fill((eesumcrystal+(eesumfiber1+eesumfiber2))/beamE);

			float ttt  = nescinttotecal / meanscinEcal;
			float ttt2 = necertotecal   / meancerEcal;
			float tty  = nescinttothcal / meanscinHcal;
			float tty2 = necertothcal   / meancerHcal;
			if( (ttt>0.1) && (ttt2>0.2) ) ehcEcalNsNc->Fill(ttt,ttt2);
			if( (tty>0.1) && (tty2>0.2) ) ehcHcalNsNc->Fill(tty,tty2);
			ehcEcalMarco->Fill(ttt2/ttt,ttt2);
			ehcHcalMarco->Fill(tty2/tty,tty2);

			float eachecks=eesumair+eesumdead+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge;
			ehetrue->Fill(eachecks/beamE);
			float eedepcal=eesumair+eesumdead+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh;
			hedepcal->Fill(eedepcal/beamE);
			enscvni->Fill(eedepcal/beamE,nine+ninh);

			eleEcalncer->Fill(ttt2);
			eleEcalnscint->Fill(ttt);

                        cout<<"******************************"<<endl;
			cout<<"GETSTUFF electrons"<<endl;
			cout<<"	hits Edeposit="<<eesum/1000.<<", beamE="<<beamE/1000.<<endl;
			cout<<" cal total energy deposit "<<eesumcal/1000.<<endl;
			cout<<" cal EM total energy deposit "<<eesumem/1000.<<endl;
			cout<<"	EDeposit: air="<<eesumair/1000.<<", ecalPD="<<eesumPDe/1000.<<", crys="<<eesumcrystal/1000.<<endl;
			cout<<"	          scintfiber="<<eesumfiber1/1000.<<", cerfiber="<<eesumfiber2/1000.<<", hcalPD="<<eesumPDh/1000.<<endl;
			cout<<"           absorber="<<eesumabs/1000.<<", edgeE="<<eesumedge/1000.<<endl;
			cout<<" sum EDeposit="<<eachecks/1000.<<", sum EDeposit/beamE="<<eachecks/beamE<<endl;
			cout<<"ecal, totncer="<<necertotecal<<", totnscint="<<nescinttotecal<<endl;
			cout<<"hcal, totncer="<<necertothcal<<", totnscint="<<nescinttothcal<<endl;
			cout<<"ehcaltimecut= "<<ehcaltimecut/1000.<<", num_inelastic="<<nine+ninh<<endl;
		}  //end loop 2nd pass through events
	}  // end process electron if no events
	ef->Close();
	cout<<"done with getstuff electrons"<<endl;
	//****************************************************************************************************************************
	if(doecal&&dohcal ) {//start calibration
		TFile* efa = TFile::Open(hcalonlyefilename);
		TTree* eta = (TTree*)efa->Get("EVENT;1");
		b_mc= eta->GetBranch("MCParticles");
		b_hcal = eta->GetBranch(HCALleaf);
		cout<<"Calibration **********************"<<endl;
		cout<<"mc branch="<<b_mc<<", hcal branch="<<b_hcal<<endl;
		num_evt_mc = b_mc->GetEntries();
		num_evt    = min(num_evt,num_evtsmax);
		cout<<"num_evt_ele hcalcalib file="<<num_evt<<endl;
		float meanscinEcal(0),meancerEcal(0),meanscinHcal(0),meancerHcal(0),meaneecaltimecut(0),meanehcaltimecut(0);
		if(num_evt>0) {
		      	CalHits* ecalhitsa = new CalHits();
		      	CalHits* hcalhitsa = new CalHits();
			b_hcal->SetAddress(&hcalhitsa);
			cout<<" branches set"<<endl;
			// first pass through file
			for(int ievt=0;ievt<num_evt; ++ievt) {
				getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, 0, dohcal, hcaltype, b_ecal,b_hcal, ecalhitsa, hcalhitsa,
						meanscinEcal, meanscinHcal, meancerEcal, meancerHcal,timecut,meaneecaltimecut,meanehcaltimecut);
			}
			meanscinHcal=meanscinHcal/num_evt;
			meancerHcal=meancerHcal/num_evt;
			cout<<"HCAL meanscint="<<meanscinHcal<<endl;
			cout<<"HCAL meancer="  <<meancerHcal<<endl;
		}
		efa->Close();
		cout<<"END calibration **********************"<<endl;
	} // end if(doecal&&dohcal)
	//
	//****************************************************************************************************************************
	// process pions
	TFile* pif = TFile::Open(piinputfilename);
	TTree* pit = (TTree*)pif->Get("EVENT;1");

	if(pif==0) cout<<" no file " <<endl;
	if(pit==0) cout<<" no event "<<endl;
	cout<<"PION ************************"<<endl;
	b_mc= pit->GetBranch("MCParticles");
	if(doecal) b_ecal = pit->GetBranch(ECALleaf);
	if(dohcal) b_hcal = pit->GetBranch(HCALleaf);
	if(doedge) b_edge = pit->GetBranch("EdgeDetNoSegment");

	float b1Ecal=0.;float m1Ecal=1.;
	float b1Hcal=0.;float m1Hcal=1.;

	num_evt_mc = b_mc->GetEntries();
	num_evt    = std::min(num_evt_mc,num_evtsmax);
	cout<<"num_evt_pi="<<num_evt<<endl;
	// loop over events
	if(num_evt>0) {
		CalHits* ecalhits = new CalHits();
		if(doecal) b_ecal->SetAddress(&ecalhits);
		CalHits* hcalhits = new CalHits();
		if(dohcal) b_hcal->SetAddress(&hcalhits);
		CalHits* edgehits = new CalHits();
		if(doedge) b_edge->SetAddress(&edgehits);
		for(int ievt=0;ievt<num_evt; ++ievt) {
			float pesum(0), pesumcal(0.),pesumem(0.), pesumair(0), pesumdead(0), pesumcrystal(0), pesumPDe(0);
			float pesumfiber1(0), pesumfiber2(0), pesumabs(0), pesumPDh(0), pesumedge(0);
			float npcertotecal(0), npscinttotecal(0), npcertothcal(0), npscinttothcal(0), eecaltimecut(0), ehcaltimecut(0);
			int nine(0),ninh(0);

			getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,
					ecalhits,hcalhits,edgehits,pesum,pesumcal,pesumem,pesumair,pesumdead, pesumcrystal,pesumPDe,
					pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumedge,
					npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,
					timecut, eecaltimecut, ehcaltimecut,piecaltime,pihcaltime,nine,ninh,
					pecalpd1scint,pecalpd1cer,pecalpd2scint,pecalpd2cer,phcalpd1scint,phcalpd1cer,phcalpd2scint,phcalpd2cer,
				        pecalpd1scintz,pecalpd1cerz,pecalpd2scintz,pecalpd2cerz,phcalpd1scintz,phcalpd1cerz,phcalpd2scintz,phcalpd2cerz,
				        eecal2d, ehcal2d);

			vector<float> piesums   = {pesumcrystal, pesumfiber1 + pesumfiber2, pesumfiber1, pesumfiber2, pesumedge, beamE-pesumedge};
			vector<float> npiphotons = {npcertotecal/meancerEcal, npcertothcal/meancerHcal, npscinttotecal/meanscinEcal, npscinttothcal/meanscinHcal};

			mes1Ecal->Fill(pesumedge/beamE,pesumem/beamE);
			mes2Ecal->Fill(pesumedge/beamE,nine);
			mes3Ecal->Fill(pesumedge/beamE,npcertotecal/meancerEcal);
			mes4Ecal->Fill(pesumedge/beamE,npscinttotecal/meanscinEcal);

			mes1Hcal->Fill(pesumedge/beamE,pesumem/beamE);
			mes2Hcal->Fill(pesumedge/beamE,ninh);
			mes3Hcal->Fill(pesumedge/beamE,npcertothcal/meancerHcal);
			mes4Hcal->Fill(pesumedge/beamE,npscinttothcal/meanscinHcal);

			auto pit = piesums.begin();
			for (auto hist : pienergy) {
				hist->Fill(*pit/beamE);
				++pit;
			}
			auto nit = npiphotons.begin();
			for (auto hist : piphoton) {
				hist->Fill(*nit);
				++nit;
			}

			hpesumcal->Fill(pesumcal/beamE);
			hpesumemcal->Fill(pesumem/beamE);
			if(pesumem>0) hpfff->Fill(pesumem/pesumcal);
			else cout << "psumemcal is zero " << endl;

			phcHcalf1f2->Fill(pesumfiber1/1000.,pesumfiber2/1000.);
			if((pesumfiber1+pesumfiber2)>0) phaphcal->Fill((pesumfiber1+pesumfiber2)/(pesumabs+pesumfiber1+pesumfiber2));
			pheest->Fill((pesumcrystal+(pesumfiber1+pesumfiber2))/beamE);

			float rrr=npscinttotecal/meanscinEcal;
			float rrr2=npcertotecal/meancerEcal;
			float rrx=npscinttothcal/meanscinHcal;
			float rrx2=npcertothcal/meancerHcal;
			if( (rrr>0.1)&&(rrr2>0.2) ) phcEcalNsNc->Fill(rrr,rrr2);
			if( (rrx>0.1)&&(rrx2>0.2) ) phcHcalNsNc->Fill(rrx,rrx2);

			phcEcalMarco->Fill(rrr2/rrr,rrr2);
			phcHcalMarco->Fill(rrx2/rrx,rrx2);

			float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumedge;
			float pedepcal=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh;
			hpdepcal->Fill(pedepcal/beamE);
			phetrue->Fill(pachecks/beamE);
			pinalvni->Fill(pedepcal/beamE,nine+ninh);
			pinscvni->Fill(rrr+rrx,  nine+ninh);
			pincevni->Fill(rrr2+rrx2,nine+ninh);

			cout<<"******************************"<<endl;
			cout<<"GETSTUFF pions"<<endl;
			cout<<" hits Edeposit="<<pesum/1000.<<", beamE="<<beamE/1000.<<endl;
			cout<<" cal total energy deposit "<<pesumcal/1000.<<std::endl;
			cout<<" cal EM total energy deposit "<<pesumem/1000.<<std::endl;
                        cout<<" EDeposit: air="<<pesumair/1000.<<", ecalPD="<<pesumPDe/1000.<<", crys="<<pesumcrystal/1000.<<endl;
                        cout<<"           scintfiber="<<pesumfiber1/1000.<<", cerfiber="<<pesumfiber2/1000.<<", hcalPD="<<pesumPDh/1000.<<endl;
                        cout<<"           absorber="<<pesumabs/1000.<<", edgeE="<<pesumedge/1000.<<endl;
                        cout<<" sum EDeposit="<<pachecks/1000.<<", sum EDeposit/beamE="<<pachecks/beamE<<endl;
                               cout<<"ecal, totnum_cer="<<npcertotecal<<", totnum_scint="<<npscinttotecal<<endl;
                        cout<<"hcal, totnum_cer="<<npcertothcal<<", totnum_scint="<<npscinttothcal<<endl;
                               cout<<"ehcaltimecut= "<<ehcaltimecut/1000.<<", num_inelastic="<<nine+ninh<<endl;
                               cout<<"******************************"<<endl;

		}  //end loop over events
		if(dodualcorr) {
			cout<<" starting fits"<<endl;
			//** fits
			TF1 *gEcale = new TF1("gEcale","[0]*(x-1.)+1",0.,1.);
			TF1 *gHcale = new TF1("gHcale","[0]*(x-1.)+1",0.,1.);
			TF1 *gEcalp = new TF1("gEcalp","[0]*(x-1.)+1",0.,1.);
			TF1 *gHcalp = new TF1("gHcalp","[0]*(x-1.)+1",0.,1.);
			// fit to get e/h
		       	if(doecal) {
				TProfile* phcEcalNsNc_pfx = phcEcalNsNc->ProfileX();
				phcEcalNsNc_pfx->Fit("gEcalp","W");
				m1Ecal=gEcalp->GetParameter(0);
				b1Ecal=1-m1Ecal;
				cout<<"ECAL b="<<b1Ecal<<", m="<<m1Ecal;
				//TCanvas* ceeee;
				//SCEDraw1(ceeee,"name",1,0,phcEdgeR,phcEcalNsNc_pfx,"outfile",0);
			}
			double kappaEcal = 1+(b1Ecal/m1Ecal);
			cout<<", k=1+b/m="<<kappaEcal<<endl;
			double hovereecalscint=piphoton[2]->GetMean();
			double hovereecalcer=piphoton[0]->GetMean();
			kappaEcal= (1-hovereecalscint)/(1.-hovereecalcer);
			cout<<"	h/e: scint="<<hovereecalscint<<", cer="<<hovereecalcer<<", k=(1-h/e(scinr))/(1-h/e(cer))="<<kappaEcal<<endl;
			if(dohcal) {
				TProfile* phcHcalNsNc_pfx = phcHcalNsNc->ProfileX();
				phcHcalNsNc_pfx->Fit("gHcalp","W");
				m1Hcal=gHcalp->GetParameter(0);
				b1Hcal=1-m1Hcal;
				cout<<"HCAL b="<<b1Hcal<<", m="<<m1Hcal<<", k=1+b/m="<<endl;
			}
			double kappaHcal = 1+(b1Hcal/m1Hcal);
			cout<<kappaHcal<<endl;
			double hoverehcalscint=piphoton[3]->GetMean();
			double hoverehcalcer=piphoton[1]->GetMean();
			kappaHcal= (1-hoverehcalscint)/(1.-hoverehcalcer);
			cout<<"    h/e: scint="<<hoverehcalscint<<", cer="<<hoverehcalcer<<", k=scint(1-h/e)/cer(1-h/e)="<<kappaHcal<<endl;

			for(int ievt=0;ievt<num_evt; ++ievt) {
				if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) cout<<"num_evt_pi="<<ievt<<endl;
				float EcorEcal(0),EcorHcal(0),eecaltimecutcor(0),ehcaltimecutcor(0);
				getStuffDualCorr(mapecalslice, mapsampcalslice, gendet, kappaEcal, kappaHcal, meanscinEcal, meancerEcal, meanscinHcal, meancerHcal,
						ievt,doecal,dohcal, hcaltype, b_ecal,b_hcal, ecalhits,hcalhits, EcorEcal, EcorHcal,
						timecut, eecaltimecutcor, ehcaltimecutcor);
				phcEcalcorr->Fill(EcorEcal);
				phcHcalcorr->Fill(EcorHcal);
			}  //end loop over events
		} // end if dualcorr
	}  // end if no events
	// close pion file
	pif->Close();
	cout<<"END PION ************************"<<endl;

	cout<<"********* Start plotting some plots ***********"<<endl;
	if(doplots) {
		int ncanv = 7;
		vector<TCanvas*> canv;
		for (int i = 0; i < ncanv; ++i) {canv.push_back(new TCanvas(Form("canvas%d", i), Form("Canvas %d", i), 800, 600));}
		SCEDraw2(canv[0],"allDepE: pi- vs e-",ehetrue,phetrue,"pivse_allE.png",0);
		SCEDraw2(canv[1],"EdgeE(escaping): pi- vs e-",eenergy[4],pienergy[4],"pivse_edgeE.png",0);
		SCEDraw2(canv[2],"BeamE-EdgeE: pi- vs e-",eenergy[5],pienergy[5],"pivse_allcaloE.png",0);
		SCEDraw1(canv[3],"e-: allDepE-EdgeE",hedepcal,"e_allcaloE.png",0);
		SCEDraw1(canv[4],"pi-: allDepE-EdgeE",hpdepcal,"pi_allcaloE.png",0);
		SCEDraw1_2D(canv[5],"e-: allDepE-EdgeE vs num_inelastic",enscvni,"e_allcaloE-vs-inel.png",0.,0.);
                SCEDraw1_2D(canv[6],"pi-: allDepE-EdgeE vs num_inelastic",pinscvni,"pi_allcaloE-vs-inel.png",0.,0.);
		if(doecal) {
			int necanv = 8;
			vector<TCanvas*> ecanv;
			for (int i = 0; i < necanv; ++i) {ecanv.push_back(new TCanvas(Form("ecanvas%d", i), Form("eCanvas %d", i), 800, 600));}
			SCEDraw2(ecanv[0],"pi- vs e-: Crystal Energy",eenergy[0],pienergy[0],"pivse_crysE.png",0);
			SCEDraw2(ecanv[1],"ECAL e-: num_scint vs num_cer",ephoton[0],ephoton[2],"eecal_scintvscer.png",0);
			SCEDraw3(ecanv[2],"ECAL pi-: num_cer, num_scint, dualcorr",piphoton[0],piphoton[2],phcEcalcorr,"piecal_ncerscintcorr.png",0);
			SCEDraw1_2D(ecanv[3],"ECAL e-: num_scint vs num_cer (TH2F)",ehcEcalNsNc,"eecal_cervsscint.png",0.,0.);
			SCEDraw1_2D(ecanv[4],"ECAL e- (Marco): num_scint vs num_cer",ehcEcalMarco,"eecal_scintvscerMarcoFit.png",0.,0.);
			SCEDraw1_2D(ecanv[5],"ECAL pi-: num_scint vs num_cer",phcEcalNsNc,"piecal_scintvscerFit.png",0.,0.);
			SCEDraw1_2D(ecanv[6],"ECAL pi- (Marco): num_scint vs num_cer",phcEcalMarco,"piecal_scintvscerMarcoFit.png",0.,0.);
			SCEDraw2(ecanv[7],"ECAL Time: pi- vs e-",eecaltime,piecaltime,"pivse_ecaltime.png",1);
		}
		if(dohcal) {
			int nhcanv = 12;
			vector<TCanvas*> hcanv;
			for (int i = 0; i < nhcanv; ++i) {hcanv.push_back(new TCanvas(Form("hcanvas%d", i), Form("hCanvas %d", i), 800, 600));}
			SCEDraw2(hcanv[0],"pi- vs e-: quatz+scint fiberE",eenergy[3],pienergy[3],"pivse_allfiber.png",0);
			SCEDraw2(hcanv[1],"pi- vs e-: scint fiberE",ephoton[1],piphoton[1],"pivse_scintfiber.png",0);
			SCEDraw2(hcanv[2],"pi- vs e-: quatz fiberE",eenergy[2],pienergy[2],"pivse_quartzfiber.png",0);
			SCEDraw2(hcanv[3],"HCAL e-: cer_num vs scint_num",ephoton[1],ephoton[3],"ehcal_cervsscint.png",0);
			SCEDraw3(hcanv[4],"HCAL pi-: num_cer, num_scint, dualcorr",piphoton[1],piphoton[3],phcHcalcorr,"pihcal_ncerscintcorr.png",0);
			SCEDraw1_2D(hcanv[5],"HCAL e-: scint_num vs cer_num",ehcHcalNsNc,"ehcal_scintvscerFit.png",0.,0.);
			SCEDraw1_2D(hcanv[6],"HCAL e- (Marco): num_scint vs num_cer",ehcHcalMarco,"ehcal_scintvscerMarcoFit.png",0.,0.);
			SCEDraw1_2D(hcanv[7],"HCAL pi-: num_scint vs num_cer",phcHcalNsNc,"pihcal_scintvscerFit.png",-b1Hcal/m1Hcal,0.);
			SCEDraw1_2D(hcanv[8],"HCAL pi-(Marco): num_scint vs num_cer",phcHcalMarco,"pihcal_scintvscerMarcoFit.png",-b1Hcal/m1Hcal,0.);
			SCEDraw1_2D(hcanv[9],"HCAL e-: scint vs quartz fiberE",ehcHcalf1f2,"ehcal_scintvsquartz.png",0.,0.);
			SCEDraw1_2D(hcanv[10],"HCAL pi-: scint vs quartz fiberE",phcHcalf1f2,"pihcal_scintvsquartz.png",0.,0.);
			SCEDraw2(hcanv[11],"HCAL Time: pi- vs e-",ehcaltime,pihcaltime,"pivse_hcaltime.png",1);
		}
	}

	//***********************************************************************************************************
	TFile * out = new TFile(outputfilename,"RECREATE");
	for (auto ee : eenergy)     {ee->Write();}
	for (auto pie : pienergy)   {pie->Write();}
	for (auto enph : ephoton)   {enph->Write();}
	for (auto pinph : piphoton) {pinph->Write();}
	phcEcalcorr->Write();
	phcHcalcorr->Write();
	ehaphcal->Write();
	phaphcal->Write();
	eheest->Write();
	pheest->Write();
	ehetrue->Write();
	phetrue->Write();
	hedepcal->Write();
	hpdepcal->Write();
	ehcEcalNsNc->Write();
	phcEcalNsNc->Write();
	ehcHcalNsNc->Write();
	phcHcalNsNc->Write();
	ehcEcalMarco->Write();
	phcEcalMarco->Write();
	ehcHcalMarco->Write();
	phcHcalMarco->Write();
	ehcHcalf1f2->Write();
	phcHcalf1f2->Write();
	eecaltime->Write();
	ehcaltime->Write();
	piecaltime->Write();
	pihcaltime->Write();
	enscvni->Write();
	pinscvni->Write();
	pincevni->Write();
	pinalvni->Write();
	heesumcal->Write();
	heesumemcal->Write();
	hefff->Write();
	hpesumcal->Write();
	hpesumemcal->Write();
	hpfff->Write();
	mes1Ecal->Write();
	mes2Ecal->Write();
	mes3Ecal->Write();
	mes4Ecal->Write();
	mes1Hcal->Write();
	mes2Hcal->Write();
	mes3Hcal->Write();
	mes4Hcal->Write();

	eecalpd1scint->Write();
	eecalpd1cer->Write();
  	pecalpd1scint->Write();
  	pecalpd1cer->Write();
  	eecalpd2scint->Write();
  	eecalpd2cer->Write();
  	pecalpd2scint->Write();
  	pecalpd2cer->Write();
  	ehcalpd1scint->Write();
  	ehcalpd1cer->Write();
  	phcalpd1scint->Write();
  	phcalpd1cer->Write();
  	ehcalpd2scint->Write();
  	ehcalpd2cer->Write();
  	phcalpd2scint->Write();
  	phcalpd2cer->Write();
      	eecalpd1scintz->Write();
       	eecalpd1cerz->Write();
        pecalpd1scintz->Write();
        pecalpd1cerz->Write();
        eecalpd2scintz->Write();
        eecalpd2cerz->Write();
        pecalpd2scintz->Write();
        pecalpd2cerz->Write();
        ehcalpd1scintz->Write();
        ehcalpd1cerz->Write();
        phcalpd1scintz->Write();
        phcalpd1cerz->Write();
        ehcalpd2scintz->Write();
        ehcalpd2cerz->Write();
        phcalpd2scintz->Write();
        phcalpd2cerz->Write();

	eecal2d->Write();
	ehcal2d->Write();

	eleEcalncer->Write();
	eleEcalnscint->Write();

	out->Close();
} // end of Resolution

void setCanvas(TCanvas* canv,const char* name) {
	canv= new TCanvas(name,name,200,10,700,500);
	canv->SetFillColor(0);
	canv->SetBorderMode(0);
	canv->SetFrameFillStyle(0);
	canv->SetFrameBorderMode(0);
	canv->SetTickx(0);
	canv->SetTicky(0);
	return;
}

void SCEDraw1(TCanvas* canv,const char* name,TH1F* h1, const char* outfile, bool logy){
	setCanvas(canv, name);
	canv->cd();
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(3);
	h1->SetStats(111111);
	h1->Draw("HIST");
	if(logy) canv->SetLogy();
	canv->Print(outfile);
	canv->Update();
	return;
}

void SCEDraw1tp(TCanvas* canv,const char* name,TProfile* t1, const char* outfile) {
	setCanvas(canv, name);
	canv->cd();
	t1->SetMarkerSize(20);
	t1->SetMarkerStyle(4);
	t1->SetMarkerColor(3);
	t1->SetLineColor(kGreen);
	t1->SetLineWidth(3);
	t1->SetStats(111111);
	t1->Draw("HIST");
	canv->Print(outfile,".png");
	canv->Update();
	return;
}

void SCEDrawp (TCanvas* canv,const char* name,TProfile* t1,const char* outfile) {
	setCanvas(canv, name);
	canv->cd();
	gStyle->SetOptFit();
	t1->SetLineColor(kGreen);
	t1->SetLineWidth(kGreen);
	t1->SetStats(111111);
	t1->SetMarkerSize(20);
	t1->SetMarkerStyle(4);
	t1->SetMarkerColor(3);
	gStyle->SetOptFit();
	t1->Draw("*");
	canv->Print(outfile,".png");
	canv->Update();
	return;
}

void SCEDraw1_2D(TCanvas* canv,const char* name,TH2F* h1,const char* outfile,float eohS,float eohC) {
	setCanvas(canv, name);
	canv->cd();
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(kGreen);
	h1->SetStats(111111);
	h1->Draw("");
	TLine line = TLine(eohS,eohC,1.,1.);
	line.SetLineColor(kYellow);
	line.SetLineWidth(2);
	line.Draw("same");
	canv->Print(outfile,".png");
	canv->Update();
	return;
}

void SCEDraw2 (TCanvas* canv,const char* name,TH1F* h1,TH1F* h2,const char* outfile, bool logy) {
	setCanvas(canv, name);
	canv->cd();
	if(logy) canv->SetLogy();
	float max = std::max(h1->GetMaximum(),h2->GetMaximum());
	h1->SetMaximum(max*1.3);
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(3);
	h1->SetStats(111111);
	h1->Draw("HIST");
	h2->SetLineColor(kRed);
	h2->SetLineWidth(3);
	h2->SetStats(111111);
	h2->Draw("HIST same");
	canv->Print(outfile,".png");
	canv->Update();
	return;
}

void SCEDraw3 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile,bool logy) {
	  setCanvas(canv, name);
	  canv->cd();
	  if(logy) canv->SetLogy();
	  float mx = max(h1->GetMaximum(),h3->GetMaximum());
	  float mn = min(h1->GetMinimum(),h3->GetMinimum());
	  h1->SetMaximum(mx*1.3);
	  h1->SetMinimum(mn*1.3);
	  h2->SetLineColor(kRed);
	  h2->SetLineWidth(3);
	  h2->SetStats(111111);
	  h2->Draw("HIST same");
	  h3->SetLineColor(kBlue);
	  h3->SetLineWidth(3);
	  h3->SetStats(111111);
	  h3->Draw("HIST same");
	  canv->Print(outfile,".png");
	  canv->Update();
	  return;
  }

void getMeanPhot(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, //input
		int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input
		TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, //input
		float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal, //output
		float &timecut, float &eecaltimecut, float &ehcaltimecut){//output
	int nbyteecal, nbytehcal, nbyteedge;
	if(doecal) {
		nbyteecal = b_ecal->GetEntry(ievt);
		// ecal hits
		if(ievt<SCECOUNT) cout<<" #ecal hits "<<ecalhits->size()<<endl;
		eecaltimecut=0.;
		for(size_t i=0;i<ecalhits->size(); ++i) { //loop over ecalhits
			CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
			long long int ihitchan=aecalhit->cellID;
			int idet = (ihitchan) & 0x07;
			int ix = (ihitchan >>3) & 0x3F ;  // is this right?
			if(ix>32) ix=ix-64;
			int iy =(ihitchan >>10) & 0x3F ; // is this right?
			if(iy>32) iy=iy-64;
			int  islice = (ihitchan >>17) & 0x07;
			int  ilayer = (ihitchan>> 20) & 0x07;
			Contributions zxzz=aecalhit->truth;
			if(gendet==1) {   // use photons as generated in optical material
				if((islice==(*eii2).second)||(islice==(*eii6).second) ) {  // crystal
					meancerEcal+=aecalhit->ncerenkov;
					meanscinEcal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==2) {
				if( (islice==(*eii1).second)||(islice==(*eii7).second) ) { // either photo detector
					meancerEcal+=aecalhit->ncerenkov;
					meanscinEcal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==3||gendet==4){
				if(idet==5) {
					if((islice==(*eii2).second)||(islice==(*eii6).second) ) {  // crystal
						meanscinEcal+=aecalhit->energyDeposit;
						if(gendet==3) meancerEcal+=aecalhit->edeprelativistic;
						if(gendet==4) meancerEcal+=aecalhit->energyDeposit;
						for(size_t j=0;j<zxzz.size(); j++) {
							if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
						}
					}
				}
			}
		} // end loop over ecalhits
	} // if doecal

	if(dohcal) { // hcal hits
		nbytehcal = b_hcal->GetEntry(ievt);
		if(ievt<SCECOUNT) cout<<" hcalhits->size()="<<hcalhits->size()<<endl;
		ehcaltimecut=0.;
		for(size_t i=0;i<hcalhits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
			long long int ihitchan=ahcalhit->cellID;
			Contributions zxzz=ahcalhit->truth;
			if(hcaltype==0) { // fiber
				int idet = (ihitchan) & 0xFF;
				int ilayer = (ihitchan >>8) & 0xFFF;
				int itube = (ihitchan >>20) & 0xFFF;
				int iair = (ihitchan >>32) & 0x7;
				int itype = (ihitchan >>35) & 0x7;
				int ifiber=0; int iabs=0; int iphdet=0;  int ihole=0;
				int ix=0; int iy=0;
				if((itype==0)&&(iair==0)&&(itube!=0)) iabs=1;
				if(itype==1) ifiber=1; // scint
				if(itype==2) ifiber=2; // quartz
				if(itype==3) iphdet=1; //scint pt
				if(itype==4) iphdet=2; // quartz pt
				if(((iair==1)||(iair==2))&&(itype==0)) ihole=1;
				if(itube==0) ihole=1;
				ix=itube;
				iy=ilayer;
				if(gendet==1) {  // take light as generated in fiber
					if(ifiber==1) {  // scintillating fibers
						meanscinHcal+=ahcalhit->nscintillator;
					}
					if(ifiber==2) {  // quartz fibers
						meancerHcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==2) {
					if(iphdet==1) {  // take light that hits photodetectors
						meanscinHcal+=ahcalhit->nscintillator;
					}
					if(iphdet==2) {  // take light that hits photodetectors
						meancerHcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if(ifiber==1) {
							meanscinHcal+=ahcalhit->energyDeposit;
						}
						if(ifiber==2) {
							if(gendet==3) meancerHcal+=ahcalhit->edeprelativistic;
							if(gendet==4) meancerHcal+=ahcalhit->energyDeposit;
						}
						for(size_t j=0;j<zxzz.size(); j++) {
							if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
						}
					}
				}
			} // end if hcaltype == fiber
			else {  // sampling
				int idet = (ihitchan) & 0x07;
				int iy = (ihitchan >>3) & 0xFFF;
				int ix = (ihitchan >>15) & 0xFFF;
				int ilayer = (ihitchan >>27) & 0xFFF;
				int ibox2 = (ihitchan >> 39) & 0x03;
				int islice = (ihitchan >>41) & 0xF;
				if(gendet==1) {  // take light as generated in media
					if(islice==(*sii3).second) {
						meanscinHcal+=ahcalhit->nscintillator;
					}
					if(islice==(*sii6).second) {  // cherenkov
						meancerHcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==2) {
					if( (islice==(*sii2).second)||(islice==(*sii4).second) ) { // either photo detector
						meanscinHcal+=ahcalhit->nscintillator;
					}
					if( (islice==(*sii5).second)||(islice==(*sii7).second)) {  // take light that hits photodetectors
						meancerHcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if( islice==(*sii3).second) { // PS (Polystyrene)
							meanscinHcal+=ahcalhit->energyDeposit;
							if(ievt<SCECOUNT) std::cout<<" meanscinHcal "<<meanscinHcal<<std::endl;
						}
						if( islice==(*sii6).second ) {  // quartz
							if(gendet==3) meancerHcal+=ahcalhit->edeprelativistic;
							if(gendet==4) meancerHcal+=ahcalhit->energyDeposit;
							if(ievt<SCECOUNT) std::cout<<" meancerHcal "<<meancerHcal<<std::endl;
						}
						for(size_t j=0;j<zxzz.size(); j++) {
							if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
						}
					}
				}
			}// end of sampling
		}  // end loop over hcal hits
	} // end of if dohcal
} // getMeanPhot

void getStuff(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input
		bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge, CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, //input
		float  &eesum,float &eesumcal,float &eesumem,float &eesumair,float &eesumdead, float &eesumcrystal,float &eesumPDe, //output
		float &eesumfiber1,float &eesumfiber2, //output
		float &eesumabs,float &eesumPDh,float &eesumedge,float &necertotecal,float &nescinttotecal,float &necertothcal, //output
		float &nescinttothcal, float &timecut, float &eecaltimecut, float &ehcaltimecut,//output
		TH1F* eecaltime, TH1F* ehcaltime, int &nine, int &ninh,
		TH1F *ecalpd1scint,TH1F *ecalpd1cer,TH1F *ecalpd2scint,TH1F *ecalpd2cer,
		TH1F *hcalpd1scint,TH1F *hcalpd1cer,TH1F *hcalpd2scint,TH1F *hcalpd2cer,
	      	TH1F *ecalpd1scintz,TH1F *ecalpd1cerz,TH1F *ecalpd2scintz,TH1F *ecalpd2cerz,
	        TH1F *hcalpd1scintz,TH1F *hcalpd1cerz,TH1F *hcalpd2scintz,TH1F *hcalpd2cerz,
	        TH2F *eecal2d, TH2F *ehcal2d)
        {
	cout<<"getStuff Func **************"<<endl;
	int nbyteecal(0), nbytehcal(0), nbyteedge(0);
	if(doecal) { // ecal hists
		nbyteecal = b_ecal->GetEntry(ievt);
		float eecaltimecut=0.;
		for(size_t i=0;i<ecalhits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
			long long int ihitchan=aecalhit->cellID;
			int idet = (ihitchan) & 0x07;
			int ix = (ihitchan >>3) & 0x3F ;  // is this right?
			if(ix>32) ix=ix-64;
			int iy =(ihitchan >>10) & 0x3F ; // is this right?
			if(iy>32) iy=iy-64;
			int  islice = (ihitchan >>17) & 0x07;
			int  ilayer = (ihitchan>> 20) & 0x07;
			if((ilayer!=0)&&(ilayer!=1)) cout<<"danger danger will robinson ilayer not zero"<<endl;
			if(islice>nsliceecal) {
				cout<<"  danger danger will robinson islice nsliceecal are "<<islice<<" "<<nsliceecal<<endl;
				cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<endl;
			} else {
				if(i<SCECOUNT&&ievt<SCECOUNT) cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<endl;
			}
			float ae=aecalhit->energyDeposit;
			nine+=aecalhit->n_inelastic; // what is this
			// check contribs
			Contributions zxzz=aecalhit->truth;
			float hacheck=0.;
			for(size_t j=0;j<zxzz.size(); j++) {
				hacheck+=(zxzz.at(j)).deposit;
				if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
				eecaltime->Fill((zxzz.at(j)).time);
			}
			if(ae>0.001) {
				if(hacheck/ae<0.99999) cout<<"missing contribs: ecal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ae is "<<ae<<" ratio "<<hacheck/ae<<endl;
			}
			eesum+=ae;
			if(islice==(*eii0).second){
				eesumair+=ae;
				eesumcal+=aecalhit->energyDeposit;
				eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii1).second){
				eesumPDe+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii2).second){
				eesumcrystal+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii3).second){
				eesumair+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii4).second){
				eesumdead+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii5).second){
				eesumair+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii6).second){
				eesumcrystal+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}
			if(islice==(*eii7).second){
				eesumPDe+=ae;
				eesumcal+=aecalhit->energyDeposit;
                                eesumem+=aecalhit->edeprelativistic;
			}

			if(islice==(*eii1).second) {  // pd on entrance to ecal
				for(int ijk=0;ijk<finenbin;ijk++){
					ecalpd1scint->Fill((ijk+0.5)*timebinsize,aecalhit->nscinttime[ijk]);
					ecalpd1cer->Fill((ijk+0.5)*timebinsize,aecalhit->ncertime[ijk]);
					ecalpd1scintz->Fill((ijk+0.5)*timebinsizez,aecalhit->nscinttimez[ijk]);
					ecalpd1cerz->Fill((ijk+0.5)*timebinsizez,aecalhit->ncertimez[ijk]);
				}
			}
			if(islice==(*eii7).second) {  // pd on exist of ecal
				for(int ijk=0;ijk<finenbin;ijk++){
					ecalpd2scint->Fill((ijk+0.5)*timebinsize,aecalhit->nscinttime[ijk]);
					ecalpd2cer->Fill((ijk+0.5)*timebinsize,aecalhit->ncertime[ijk]);
					ecalpd2scintz->Fill((ijk+0.5)*timebinsizez,aecalhit->nscinttimez[ijk]);
					ecalpd2cerz->Fill((ijk+0.5)*timebinsizez,aecalhit->ncertimez[ijk]);
				}
			}


			if(gendet==1) {   // use photons as generated in optical material
				if((islice==(*eii2).second)||(islice==(*eii6).second)) {  // crystal
					necertotecal+=aecalhit->ncerenkov;
					nescinttotecal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==2){
				if( (islice==(*eii1).second)||(islice==(*eii7).second) ) { // either photo detector
					necertotecal+=aecalhit->ncerenkov;
					nescinttotecal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==3||gendet==4){
				if(idet==5) {
					if((islice==(*eii2).second)||(islice==(*eii6).second) ) {  // crystal
						nescinttotecal+=aecalhit->energyDeposit;
						if(gendet==3) necertotecal+=aecalhit->edeprelativistic;
						if(gendet==4) necertotecal+=aecalhit->energyDeposit;
					}
				}
			}
			eecal2d->Fill(ix,iy,aecalhit->energyDeposit);
			//cout<<"idet=="<<idet<<", gendet=="<<gendet<<", islice=="<<islice<<endl;
		} //end loop over ecalhits
	} //end of doecal

	if(dohcal) {
		nbytehcal = b_hcal->GetEntry(ievt);
		// hcal hits
		float ehcaltimecut=0.;
		for(size_t i=0;i<hcalhits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
			float ah=ahcalhit->energyDeposit;
			ninh+=ahcalhit->n_inelastic;
			// check contribs
			Contributions zxzz=ahcalhit->truth;
			float hacheck=0.;
			for(size_t j=0;j<zxzz.size(); j++) {
				hacheck+=(zxzz.at(j)).deposit;
				if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
				ehcaltime->Fill((zxzz.at(j)).time);
			}
			if(ah>0.001) {
				if(hacheck/ah<0.99999) cout<<"missing contribs: hcal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ah is "<<ah<<" ratio "<<hacheck/ah<<endl;
			}
			eesum+=ah;
			eesumcal+=ahcalhit->energyDeposit;
			eesumem+=ahcalhit->edeprelativistic;

			long long int ihitchan=ahcalhit->cellID;
			if(hcaltype==0) { // fiber
				int idet = (ihitchan) & 0xFF;
				int ilayer = (ihitchan >>8) & 0xFFF;
				int itube = (ihitchan >>20) & 0xFFF;
				int iair = (ihitchan >>32) & 0x7;
				int itype = (ihitchan >>35) & 0x7;
				int ifiber=0; int iabs=0; int iphdet=0;  int ihole=0;
				int ix=0; int iy=0;
				if((itype==0)&&(iair==0)&&(itube!=0)) iabs=1;
				if(itype==1) ifiber=1; // scint
				if(itype==2) ifiber=2; // quartz
				if(itype==3) iphdet=1; //scint pt
				if(itype==4) iphdet=2; // quartz pt
				if(((iair==1)||(iair==2))&&(itype==0)) ihole=1;
				if(itube==0) ihole=1;
				ix=itube;
				iy=ilayer;
				if(gendet==1) {  // take light as generated in fiber
					if(ifiber==1) {  // scintillating fibers
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if(ifiber==2) {  // quartz fibers
						necertothcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==2) {
					if(iphdet==1) {  // take light that hits photodetectors
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if(iphdet==2) {  // take light that hits photodetectors
						necertothcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if(ifiber==1) {//scint fiber
							nescinttothcal+=ahcalhit->energyDeposit;
						}
						if(ifiber==2) {//quartz (cer) fiber
							if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
							if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
						}
					}
				}
				if(ifiber==1) {eesumfiber1+=ah;}
				if(ifiber==2) {eesumfiber2+=ah;}
				if(iabs==1) {eesumabs+=ah;}
				if(iphdet>1) {eesumPDh+=ah;}
			} //end of hcaltype=fiber
			else {  // sampling
				int idet = (ihitchan) & 0x07;
				int iy = (ihitchan >>3) & 0xFFF;
				int ix = (ihitchan >>15) & 0xFFF;
				int ilayer = (ihitchan >>27) & 0xFFF;
				int ibox2 = (ihitchan >> 39) & 0x03;
				int islice = (ihitchan >>41) & 0xF;
				cout<<" idet iy ix ilayer islice are "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<endl;
				cout<<"energy nscint ncer is "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<endl;
				cout<<" ps ii is "<<(*sii3).second<<endl;
				cout<<" quartz ii is "<<(*sii6).second<<endl;
				if(gendet==1) {  // take light as generated in media
					if(islice==(*sii3).second) {
						nescinttothcal+=ahcalhit->nscintillator;
						cout<<"add scint"<<endl;
					}
					if(islice==(*sii6).second) {  // chereknov
						necertothcal+=ahcalhit->ncerenkov;
						cout<<"add ceren"<<endl;
					}
				}
				else if(gendet==2) {
					if( (islice==(*sii2).second)||(islice==(*sii4).second) ) { // either photo detector
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if( (islice==(*sii5).second)||(islice==(*sii7).second)) {  // take light that hits photodetectors
						necertothcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if( islice==(*sii3).second) { //  ps
							nescinttothcal+=ahcalhit->energyDeposit;
						}
						if( islice==(*sii6).second ) {  // quartz
							if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
							if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
						}
					}
				}
				if( islice==(*sii3).second ) eesumfiber1+=ah; // scint
				if( islice==(*sii6).second ) eesumfiber2+=ah;  //cer
				if( (islice==(*sii1).second) || (islice==(*sii8).second) || (islice==(*sii9).second) ) eesumabs+=ah;
				if( (islice==(*sii2).second) || (islice==(*sii4).second) || (islice==(*sii5).second) || (islice==(*sii7).second)) eesumPDh+=ah;

				ehcal2d->Fill(ix,iy,ahcalhit->energyDeposit);
			} //end of hcal sampling
		}  // end loop over hcal hits
	}//end of if dohcal
	if(doedge) {
		nbyteedge = b_edge->GetEntry(ievt);
		if(ievt<SCECOUNT) cout<<" edgehits->size()= "<<edgehits->size()<<endl;
		for(size_t i=0;i<edgehits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* aedgehit =edgehits->at(i);
			float ae=aedgehit->energyDeposit;
			eesum+=ae;
			eesumedge+=ae;
		}  // end loop over escaping hits
	} // end doedge
} // end of getStuff

void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaEcal, float kappaHcal,
		float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype,
		TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, float &EEcal, float &EHcal,
		float &timecut, float &eecaltimecut, float &ehcaltimecut){

	float necertotecal(0),nescinttotecal(0),necertothcal(0),nescinttothcal(0);
	int nbyteecal, nbytehcal, nbyteedge;
	if(doecal) {
		nbyteecal = b_ecal->GetEntry(ievt);// ecal hits
		if(ievt<SCECOUNT) cout<<" ecalhits->size()="<<ecalhits->size()<<endl;
		eecaltimecut=0.;
		for(size_t i=0;i<ecalhits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
			long long int ihitchan=aecalhit->cellID;
			int idet = (ihitchan) & 0x07;
			int ix = (ihitchan >>3) & 0x3F ;  // is this right?
			if(ix>32) ix=ix-64;
			int iy =(ihitchan >>10) & 0x3F ; // is this right?
			if(iy>32) iy=iy-64;
			int  islice = (ihitchan >>17) & 0x07;
			int  ilayer = (ihitchan>> 20) & 0x07;
			Contributions zxzz=aecalhit->truth;
			if(gendet==1) {   // use photons as generated in optical material
				if((islice==(*eii2).second)||(islice==(*eii6).second)) {
					necertotecal+=aecalhit->ncerenkov;
					nescinttotecal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==2) {
				if( (islice==(*eii1).second)||(islice==(*eii7).second) ) {
					necertotecal+=aecalhit->ncerenkov;
					nescinttotecal+=aecalhit->nscintillator;
				}
			}
			else if(gendet==3||gendet==4){
				if(idet==5) {
					if((islice==(*eii2).second)||(islice==(*eii6).second) ) {  // crystal
						nescinttotecal+=aecalhit->energyDeposit;
						if(gendet==3) necertotecal+=aecalhit->edeprelativistic; //cerenkov light
						if(gendet==4) necertotecal+=aecalhit->energyDeposit;
						for(size_t j=0;j<zxzz.size(); j++) {if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit; }
					}
				}
			}
		}  // end loop over ecal hits
	}  // end do ecal

	if(dohcal) {
		nbytehcal = b_hcal->GetEntry(ievt);// hcal hits
		ehcaltimecut=0.;
		for(size_t i=0;i<hcalhits->size(); ++i) {
			CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
			Contributions zxzz=ahcalhit->truth;
			long long int ihitchan=ahcalhit->cellID;
			if(hcaltype==0) { // fiber
				int idet = (ihitchan) & 0xFF;
				int ilayer = (ihitchan >>8) & 0xFFF;
				int itube = (ihitchan >>20) & 0xFFF;
				int iair = (ihitchan >>32) & 0x7;
				int itype = (ihitchan >>35) & 0x7;
				int ifiber=0; int iabs=0; int iphdet=0;  int ihole=0;
				int ix=0; int iy=0;
				if((itype==0)&&(iair==0)&&(itube!=0)) iabs=1;
				if(itype==1) ifiber=1; // scint
				if(itype==2) ifiber=2; // quartz
				if(itype==3) iphdet=1; //scint pt
				if(itype==4) iphdet=2; // quartz pt
				if(((iair==1)||(iair==2))&&(itype==0)) ihole=1;
				if(itube==0) ihole=1;
				ix=itube;
				iy=ilayer;
				if(gendet==1) {  // take light as generated in fiber
					if(ifiber==1) {  // scintillating fibers
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if(ifiber==2) {  // quartz fibers
						necertothcal+=ahcalhit->ncerenkov;}
				}
				else if(gendet==2) {
					if(iphdet==1) {  // take light that hits photodetectors
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if(iphdet==2) {  // take light that hits photodetectors
						necertothcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if(ifiber==1) { nescinttothcal+=ahcalhit->energyDeposit;}
						if(ifiber==2) {
							if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
							if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
						}
						for(size_t j=0;j<zxzz.size(); j++) {
							if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
						}
					}
				}
			}//end of hcaltype=fiber
			else {  // sampling
				int idet = (ihitchan) & 0x07;
				int iy = (ihitchan >>3) & 0xFFF;
				int ix = (ihitchan >>15) & 0xFFF;
				int ilayer = (ihitchan >>27) & 0xFFF;
				int ibox2 = (ihitchan >> 39) & 0x03;
				int islice = (ihitchan >>41) & 0xF;
				if(gendet==1) {  // take light as generated in media
					if(islice==(*sii3).second) { nescinttothcal+=ahcalhit->nscintillator; }
					if(islice==(*sii6).second) { necertothcal+=ahcalhit->ncerenkov; } //add cerenkov
				}
				else if(gendet==2) {
					if( (islice==(*sii2).second)||(islice==(*sii4).second) ) { // either photo detector
						nescinttothcal+=ahcalhit->nscintillator;
					}
					if( (islice==(*sii5).second)||(islice==(*sii7).second)) {  // take light that hits photodetectors
						necertothcal+=ahcalhit->ncerenkov;
					}
				}
				else if(gendet==3||gendet==4) {
					if(idet==6) {
						if( islice==(*sii3).second ) { nescinttothcal+=ahcalhit->energyDeposit; } //ps
						if( islice==(*sii6).second ) {  // quartz
							if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
							if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
						}
						for(size_t j=0;j<zxzz.size(); j++) {if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;}
					}
				}
			} // end of hcal sampling
		}  // end for loop on hcalhits
	}  // end if dohcal
	if(doecal){
		float anecertotecal=necertotecal/meancerEcal;
		float anescinttotecal=nescinttotecal/meanscinEcal;
		EEcal=(anescinttotecal-kappaEcal*anecertotecal)/(1-kappaEcal);
                cout<<"ECAL:"<<endl;
		cout<<"		necertot="<<necertotecal   <<", meancer="  <<meancerEcal <<", necertot/meancer="   <<anecertotecal<<endl;
                cout<<"         nescintot="<<nescinttotecal<<", meanscint="<<meanscinEcal<<", nescintot/meanscint="<<anescinttotecal<<endl;
                cout<<"         kappaEcal="<<kappaEcal     <<", EEcal="    <<EEcal<<endl;
        }
	if(dohcal){
		float anecertothcal=necertothcal/meancerHcal;
		float anescinttothcal=nescinttothcal/meanscinHcal;
		EHcal=(anescinttothcal-kappaHcal*anecertothcal)/(1-kappaHcal);
		cout<<"HCAL: "<<endl;
		cout<<"		necertot="<<necertothcal   <<", meancer="  <<meancerHcal <<", necertot/meancer="   <<anecertothcal<<endl;
		cout<<"		nescintot="<<nescinttothcal<<", meanscint="<<meanscinHcal<<", nescintot/meanscint="<<anescinttothcal<<endl;
		cout<<"		kappaHcal="<<kappaHcal     <<", EHcal="    <<EHcal<<endl;
	}

} // end of getStuffDualCorr
