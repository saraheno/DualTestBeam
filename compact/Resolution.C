#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Factories.h"
#include "DDG4/Geant4Particle.h"
#include "DDG4/Geant4Data.h"
#include "../include/DualCrysCalorimeterHit.h"

#include <vector>
#include <functional>
#include <map>
#include <algorithm>

#include <iostream>
#include <sstream> // for ostringstream
#include <string>
using namespace std;


int SCECOUNT=5;
int SCECOUNT2=20;
int icountyuck=0;
int SCECOUNT3=10;


float timecut=10;
float betacut=1/1.5;
const int finenbin=40;
const float timemin=0.;
const float timemax=400.;
const float timemaxz=40.;
const float timebinsize=(timemax-timemin)/float(finenbin);
const float timebinsizez=(timemaxz-timemin)/float(finenbin);

float kappaEcal(1.),kappaHcal(1.);
float biggesttime=0.;
map<string, int> mapsampcalslice;

typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


void SCEDraw1 (TCanvas* canv, const char* name, TH1F* h1, const char* outfile, bool logy);
void SCEDraw1tp (TCanvas* canv, const char* name, TProfile* h1, const char* outfile);
void SCEDraw1_2D (TCanvas* canv, const char* name, TH2F* h1, const char* outfile,bool dline, float eohS,float eohC);
void SCEDraw2_2D (TCanvas* canv, const char* name, TH2F* h1, TH2F* h2, const char* outfile,bool doline, float eohS,float eohC);
void SCEDraw2 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile,bool logy);
void SCEDraw3 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile, bool logy);


void getStuff(map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, float &timecut, bool &fillhists,
	      float  &eesum,float &eesumcal,float &eesumem, float &eesumair,float &eesumdead, float &eesumcrystal,float &eesumPDe,float &eesumfiber1,float &eesumfiber2,float &eesumabs,float &eesumPDh,float &eesumairem, float &eesumdeadem, float &eesumcrystalem,float &eesumPDeem,float &eesumfiber1em, float &eesumfiber2em,float &eesumabsem,float &eesumPDhem,float &eesumedge,float &eesumedgerel, float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal,float &eecaltimecut, float &ehcaltimecut,float &erelecaltimecut, float &erelhcaltimecut,int &nine,int &ninh);
void FillTime(map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, float &timecut,
	      TH1F* eecaltime, TH1F* ehcaltime, TH1F *ecalpd1scint,TH1F *ecalpd1cer,TH1F *ecalpd2scint,TH1F *ecalpd2cer,TH1F *hcalpd1scint,TH1F *hcalpd1cer,TH1F *hcalpd2scint,TH1F *hcalpd2cer,TH1F *ecalpd1scintz,TH1F *ecalpd1cerz,TH1F *ecalpd2scintz,TH1F *ecalpd2cerz,TH1F *hcalpd1scintz,TH1F *hcalpd1cerz,TH1F *hcalpd2scintz,TH1F *hcalpd2cerz);

void getStuffDualCorr(bool domissCorr, float beamE, map<string, int> mapsampcalslice, int gendet, float kappaecal, float kappahcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype, bool doedge,float &eesumedge, float &eesumedgerel, TBranch* &b_ecal,TBranch* &b_hcal, TBranch* &b_edge,CalHits* &ecalhits, CalHits* &hcalhits,CalHits* &edgehits,float &EEcal, float &EHcal,float &timecut, float &eecaltimecut, float &ehcaltimecut, float &erelecaltimecut, float &erelhcaltimecut);

void getMeanPhot(map<string, int> mapsampcalslice,  int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,CalHits* &ecalhits, CalHits* &hcalhits,float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,float &timecut, float &eecaltimecut, float &ehcaltimecut, float &erelecaltimecut, float &erelhcaltimecut);

void CalibRefine(map<string, int> mapsampcalslice,  int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,CalHits* &ecalhits, CalHits* &hcalhits, float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,TH1F *CalEcalncer, TH1F *CalEcalnscint, TH1F *CalHcalncer, TH1F *CalHcalnscint
);


// hcal type 0=fiber, 1 = sampling
// gendet 1=active media photons, 2 = photodetector, 3=energy deposit 4 is a debug gendet
// ECALleaf is

//Resolution(2,"./output/out_DualTestBeam_20GeV_e-.root","./output/out_DualTestBeam_20GeV_pi-.root","./output/out_FSCEPonly_20GeV_e-.root","./output/out_FSCEPonly_20GeV_pi-.root",20,1,1,0,1,0,0,0,3,"hists_20GeV.root","DRCNoSegment","DRFNoSegment",1,0,1,1)



void Resolution(int num_evtsmax, const char* einputfilename, const char* piinputfilename,
		const char* hcalonlyefilename, const char* hcalonlypifilename,
		const float beamEE, bool doecal, bool dohcal, int hcaltype, bool doedge, bool domissCorr,bool doedgecut, float edgecut,int gendet, const char* outputfilename,const char* ECALleaf, const char* HCALleaf,bool doplots, bool dotimingplots,bool dodualcorr,bool twocalecalcorr) {

  // these must correspond to the "slice" physvolid used in DRCrys_geo
  // these correspond to slices in scepcal_drcrystal.xml
  // surely there is a better way to do this
  mapsampcalslice["air"]=0;
  mapsampcalslice["Iron"]=1;
  mapsampcalslice["PD1"]=2;
  mapsampcalslice["PS"]=3;
  mapsampcalslice["PD2"]=4;
  mapsampcalslice["PD3"]=5;
  mapsampcalslice["Quartz"]=6;
  mapsampcalslice["PD4"]=7;
  mapsampcalslice["Sep1"]=8;
  mapsampcalslice["Sep2"]=9;

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

  int ihaha,num_evt;
  float arms,amean;
  bool fillfill;  // filling some timing histograms in getstuff
  TBranch* b_mc;
  TBranch* b_ecal;
  TBranch* b_hcal;
  TBranch* b_edge;

  // define histograms

  TH1F *ehchan = new TH1F("ehchan","channel ID number",1028,0.,1028);
  TH1F *phchan = new TH1F("phchan","channel ID number",1028,0.,1028);
  TH1F *ehcEcalE = new TH1F("ehcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);
  TH1F *phcEcalE = new TH1F("phcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);
  TH1F *ehcHcalE = new TH1F("ehcHcalE","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE = new TH1F("phcHcalE","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalE1 = new TH1F("ehcHcalE1","fiber 1 hcal energy / cal E",200,0.,0.2);
  TH1F *phcHcalE1 = new TH1F("phcHcalE1","fiber 1 hcal energy / cal E",200,0.,0.2);
  TH1F *ehcHcalE2 = new TH1F("ehcHcalE2","fiber 2 hcal energy / cal E",200,0.,0.2);
  TH1F *phcHcalE2 = new TH1F("phcHcalE2","fiber 2 hcal energy / cal E",200,0.,0.2);
  TH2F *phcHcalvfE1 = new TH2F("phcHcalvfE1","fiber 1 hcal energy / cal E",200,0.,1.2,200,0.,0.2);
  TH2F *phcHcalvfE2 = new TH2F("phcHcalvfE2","fiber 2 hcal energy / cal E",200,0.,1.2,200,0.,0.2);
  TH1F *ehcEdgeRelf = new TH1F("ehcEdgeRelf","rel frac edge e",200,0.,1.);
  TH1F *phcEdgeRelf = new TH1F("phcEdgeRelf","rel frac edge e",200,0.,1.);
  TH1F *ehcEdgeE = new TH1F("ehcEdgeE","sum edge / beam E",200,0.,1.0);
  TH1F *phcEdgeE = new TH1F("phcEdgeE","sum edge / beamE-Erel/fnorm",200,0.,1.0);
  TH1F *ehcnonconsE = new TH1F("ehcnonconsE","non cons / beam E-Erel/fnorm",100,0.,1.5);
  TH1F *phcnonconsE = new TH1F("phcnonconsE","non cons ",100,0.,1.5);
  TH1F *phcnonconsE2 = new TH1F("phcnonconsE2","non cons / Ecal",100,0.,1.5);
  TH1F *phcnonconsE3 = new TH1F("phcnonconsE3","non cons / Ecal-Erel",100,0.,1.5);
  TH1F *phcnonconsE4 = new TH1F("phcnonconsE4","non cons / eCal-Erel/fnorm",100,0.,1.5);
  TH1F *ehcEdgeR = new TH1F("ehcEdgeR","beam - sum edge / beam E",100,0.,1.5);
  TH1F *phcEdgeR = new TH1F("phcEdgeR","beam - sum edge / beam E",100,0.,1.5);
  TH1F *ehcEcalncer = new TH1F("ehcEcalncer","total number of ecal cerenkov",  500,0.,1.5);
  TH1F *phcEcalncer = new TH1F("phcEcalncer","total number of ecal cerenkov", 500,0.,1.5);
  TH1F *phcEcalncer2 = new TH1F("phcEcalncer2","total number of ecal cerenkov", 500,0.,1.5);
  TH1F *phcEcalncer3 = new TH1F("phcEcalncer3","total number of ecal cerenkov", 500,0.,1.5);
  TH1F *ehcHcalncer = new TH1F("ehcHcalncer","total number of hcal cerenkov",   500,0.,1.5);
  TH1F *phcHcalncer = new TH1F("phcHcalncer","total number of hcal cerenkov",  500,0.,1.5);
  TH1F *phcHcalncer2 = new TH1F("phcHcalncer2","total number of hcal cerenkov",  500,0.,1.5);
  TH1F *phcHcalncer3 = new TH1F("phcHcalncer3","total number of hcal cerenkov",  500,0.,1.5);
  TH1F *phcEandHcalncer = new TH1F("phcEandHcalncer","total number of ecal+hcal cerenkov", 500,0.,1.5);
  TH1F *phcEandHcalnscint = new TH1F("phcEandHcalnscint","total number of ecal+hcal scintillation", 500,0.,1.5);
  TH1F *ehcEcalnscint = new TH1F("ehcEcalnscint","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalnscint = new TH1F("phcEcalnscint","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalnscint2 = new TH1F("phcEcalnscint2","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalnscint3 = new TH1F("phcEcalnscint3","total number of ecal scintillation", 500,0.,1.5);
  TH1F *ehcHcalnscint = new TH1F("ehcHcalnscint","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalnscint = new TH1F("phcHcalnscint","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalnscint2 = new TH1F("phcHcalnscint2","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalnscint3 = new TH1F("phcHcalnscint3","total number of hcal scintillation", 500,0.,1.5);
  TH1F *ehcEcalcorr = new TH1F("ehcEcalcorr","ecal dual corr", 500,0.,1.5);
  TH1F *phcEcalcorr = new TH1F("phcEcalcorr","ecal dual corr", 500,0.,1.5);
  TH1F *ehcHcalcorr = new TH1F("ehcHcalcorr","e hcal dual", 500,0.,1.5);
  TH1F *phcHcalcorr = new TH1F("phcHcalcorr","pi hcal dual", 500,0.,1.5);
  TH1F *phcEandHcalcorr = new TH1F("phcEandHcalcorr","pi hcal dual", 500,0.,1.5);
  TH2F *ehcEcalNsNc = new TH2F("ehcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcEcalNsNc = new TH2F("phcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalNsNc = new TH2F("ehcHcalNsNc","hcal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcHcalNsNc = new TH2F("phcHcalNsNc","hcal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *enonconsEcalcer = new TH2F("enonconsEcalcer","Ecal ncer versus non conservation of E",500,0.,1.0,500,0.,1.4);
  TH2F *enonconsEcalscin = new TH2F("enonconsEcalscin","Ecal nscin versus non conservation of E",500,0.,1.0,500,0.,1.4);
  TH2F *enonconsHcalcer = new TH2F("enonconsHcalcer","Hcal ncer versus non conservation of E",500,0.,1.0,500,0.,1.4);
  TH2F *enonconsHcalscin = new TH2F("enonconsHcalscin","Hcal nscin versus non conservation of E",500,0.,1.0,500,0.,1.4);
  TH2F *hcnonconsvesc = new TH2F("hcnonconsvesc","pion escaping versus noncons",500,0.,1.,500,0.,1.);
  TH2F *enonconsvnni  = new TH2F("enonconsvnni","non conservation of E versus num nucl",500,0.,1000,500,0.,1.0);
  TH2F *enonconsvf  = new TH2F("enonconsvf","non conservation of E versus f",500,0.,1.1,500,0.,1.0);
  TH2F *mes1Ecal = new TH2F("mes1Ecal","Ecal edge energy versus em frac",500,0.,1.5,500,0.,0.3);
  TH2F *mes2Ecal = new TH2F("mes2Ecal","Ecal edge energy versus num nucl int",500,0.,1000,500,0.,0.3);
  TH2F *mes3Ecal = new TH2F("mes3Ecal","Ecal edge energy versus ncer",500,0.,0.3,500,0.,1.5);
  TH2F *mes4Ecal = new TH2F("mes4Ecal","Ecal edge energy versus nscint",500,0.,0.3,500,0.,1.5);
  TH2F *mes1Hcal = new TH2F("mes1Hcal","Hcal edge energy versus em frac",500,0.,1.5,500,0.,1.0);
  TH2F *mes2Hcal = new TH2F("mes2Hcal","Hcal edge energy versus num nucl int",500,0.,1000,500,0.,1.0);
  TH2F *mes3Hcal = new TH2F("mes3Hcal","Hcal edge energy versus ncer",500,0.,0.3,500,0.,1.5);
  TH2F *mes4Hcal = new TH2F("mes4Hcal","Hcal edge energy versus nscint",500,0.,0.3,500,0.,1.5);
  TH2F *hfnscinEcal = new TH2F("hfnscinEcal","scin versus f ECAL",500,0.,1.5,500,0.,1.5);
  TH2F *hfncerEcal = new TH2F("hfncerEcal","cer versus f ECAL",500,0.,1.5,500,0.,1.5);
  TH2F *hfnscinHcal = new TH2F("hfnscinHcal","scin versus f ECAL",500,0.,1.5,500,0.,1.5);
  TH2F *hfncerHcal = new TH2F("hfncerHcal","cer versus f ECAL",500,0.,1.5,500,0.,1.5);
  TH2F *ehcEcalNsNctc = new TH2F("ehcEcalNsNctc","ecal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);
  TH2F *phcEcalNsNctc = new TH2F("phcEcalNsNctc","ecal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalNsNctc = new TH2F("ehcHcalNsNctc","hcal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);
  TH2F *phcHcalNsNctc = new TH2F("phcHcalNsNctc","hcal ncer versus nscint time cut",500,0.,1.5,500,0.,1.5);
  TH2F *ehcEcalMarco = new TH2F("ehcEcalMarco","ecal c v c/s",500,0.,1.5,500,0.,1.5);
  TH2F *phcEcalMarco = new TH2F("phcEcalMarco","ecal c v c/s",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalMarco = new TH2F("ehcHcalMarco","hcal c v c/s",500,0.,1.5,500,0.,1.5);
  TH2F *phcHcalMarco = new TH2F("phcHcalMarco","hcal c v c/s",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalf1f2 = new TH2F("ehcHcalf1f2","hcal f1e versus f2e",500,0.,2.5,500,0.,2.5);
  TH2F *phcHcalf1f2 = new TH2F("phcHcalf1f2","hcal f1e versus f2e",500,0.,2.5,500,0.,2.5);
  TH2F *ehecal2d = new TH2F("ehecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *phecal2d = new TH2F("phecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *ehhcal2d = new TH2F("ehhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *phhcal2d = new TH2F("phhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH1F *ehaphcal = new TH1F("ehaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);
  TH1F *phaphcal = new TH1F("phaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);
  TH1F *eheest = new TH1F("eheest","ratio estimated to true energy",500,0.,1.);
  TH1F *pheest = new TH1F("pheest","ratio estimated to true energy",500,0.,1.);
  TH1F *ehetrue = new TH1F("ehetrue","ratio deposited to incident energy",500,0.,1.1);
  TH1F *phetrue = new TH1F("phetrue","ratio deposited to incident energy",500,0.,1.1);
  TH1F *hedepcal = new TH1F("hedepcal","elec all deposited energies",500,0.,1.1);
  TH1F *hpdepcal = new TH1F("hpdepcal","pi all deposited energies",500,0.,1.1);
  TH1F *ehnecalcon = new TH1F("ehnecalcon","number contribs to ecal hit",1010,-10.,1000.);
  TH1F *phnecalcon = new TH1F("phnecalcon","number contribs to ecal hit",1010,-10.,1000.);
  TH2F *ehzvst = new TH2F("ehzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.);
  TH2F *phzvst = new TH2F("phzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.);
  TH1F *eecaltime = new TH1F("eecaltime","time of ecal const",100,0.,40.);
  TH1F *ehcaltime = new TH1F("ehcaltime","time of hcal const",100,0.,40.);
  TH1F *piecaltime = new TH1F("piecaltime","time of ecal const",100,0.,40.);
  TH1F *pihcaltime = new TH1F("pihcaltime","time of hcal const",100,0.,40.);
  TH2F *enscvni = new TH2F("enscvni","electron scint E versus number inelastic",100,0.,1000,100,0.,1.4);
  TH2F *pinalvni = new TH2F("pinalvni","pion all E versus number inelastic",100,0.,1000,100,0.,1.4);
  TH2F *pinscvni = new TH2F("pinscvni","pion scint  versus number inelastic",100,0.,1000,100,0.,1.4);
  TH2F *pincevni = new TH2F("pincevni","pion cherenkov  versus number inelastic",100,0.,1000,100,0.,1.4);
  TH1F *heesumcal = new TH1F("heesumcal","electron energy in calorimeter",400,0.,1.5);
  TH1F *heesumemcal = new TH1F("heesumemcal","electron relativistic energy in calorimeter",400,0.,1.5);
  TH1F *hefff = new TH1F("hefff","electron relativistic fraction in calorimeter",400,0.,1.5);
  TH1F *hpesumcal = new TH1F("hpesumcal","pion energy in calorimeter",400,0.,1.5);
  TH1F *hpesumemcal = new TH1F("hpesumemcal","pion relativistic energy in calorimeter",400,0.,1.5);
  TH1F *hpfff = new TH1F("hpfff","gamma fraction in calorimeter",400,0.,1.5);
  TH1F *hpfff2 = new TH1F("hpfff2","pion relativistic energy over beam energy",400,0.,1.5);
  TH1F *hpfff3 = new TH1F("hpfff3","pion relativistic energy over cal energy",400,0.,1.5);
  TH1F *hpfffabs = new TH1F("hpfffabs","pion relativistic fraction in absorber",400,0.,1.5);
  TH1F *hpffffib = new TH1F("hpffffib","pion relativistic fraction in fiber",400,0.,1.5);
  TH1F *CalEcalncer = new TH1F("CalEcalncer","total number of ecal cerenkov",  500,0.,1.5);
  TH1F *CalEcalnscint = new TH1F("CalEcalnscint","total number of ecal scint",  500,0.,1.5);
  TH1F *CalHcalncer = new TH1F("CalHcalncer","total number of hcal cerenkov",  500,0.,1.5);
  TH1F *CalHcalnscint = new TH1F("CalHcalnscint","total number of hcal scint",  500,0.,1.5);
  TH1F *acovECAL = new TH1F("acovECAL","covariance ecal",500,0.,1.);
  TH1F *acovHCAL = new TH1F("acovHCAL","covariance hcal",500,0.,1.);
  TH1F *ffscinffcer = new TH1F("ffscinffcer","f from scin and cer",100,0.,1.2);

  //********************************************************
  //  DANGER DANGER WILL ROBINSON
  //  THIS MUST ALIGN WITH TIMEMIN TIMEMAX AND FINENBIN IN ../src/DualCrysCalorimeterHit.h
  // if you change these, you MUST remake the samples

  std::cout<<"warning warning if you change the timing histograms, please read the comment in the code"<<std::endl;
  TH1F *eecalpd1scint = new TH1F("eecalpd1scint","electron scint photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
  TH1F *eecalpd1cer = new TH1F("eecalpd1cer","electron cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
  TH1F *pecalpd1scint = new TH1F("pecalpd1scint","pion scint photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
  TH1F *pecalpd1cer = new TH1F("pecalpd1cer","pion cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemax);
  TH1F *eecalpd2scint = new TH1F("eecalpd2scint","electron scint photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
  TH1F *eecalpd2cer = new TH1F("eecalpd2cer","electron cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
  TH1F *pecalpd2scint = new TH1F("pecalpd2scint","pion scint photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
  TH1F *pecalpd2cer = new TH1F("pecalpd2cer","pion cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemax);
  TH1F *ehcalpd1scint = new TH1F("ehcalpd1scint","elec scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
  TH1F *ehcalpd1cer = new TH1F("ehcalpd1cer","elec cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
  TH1F *phcalpd1scint = new TH1F("phcalpd1scint","pion scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
  TH1F *phcalpd1cer = new TH1F("phcalpd1cer","pion cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemax);
  TH1F *ehcalpd2scint = new TH1F("ehcalpd2scint","elec scint photon arrival time ns HCAL quartz fiber",finenbin,timemin,timemax);
  TH1F *ehcalpd2cer = new TH1F("ehcalpd2cer","elec cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemax);
  TH1F *phcalpd2scint = new TH1F("phcalpd2scint","pion scint photon arrival time ns quartz fiber",finenbin,timemin,timemax);
  TH1F *phcalpd2cer = new TH1F("phcalpd2cer","pion cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemax);

  TH1F *eecalpd1scintz = new TH1F("eecalpd1scintz","electron scint photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
  TH1F *eecalpd1cerz = new TH1F("eecalpd1cerz","electron cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
  TH1F *pecalpd1scintz = new TH1F("pecalpd1scintz","pion scint photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
  TH1F *pecalpd1cerz = new TH1F("pecalpd1cerz","pion cerenov photon arrival time ns ECAL PD1",finenbin,timemin,timemaxz);
  TH1F *eecalpd2scintz = new TH1F("eecalpd2scintz","electron scint photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
  TH1F *eecalpd2cerz = new TH1F("eecalpd2cerz","electron cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
  TH1F *pecalpd2scintz = new TH1F("pecalpd2scintz","pion scint photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
  TH1F *pecalpd2cerz = new TH1F("pecalpd2cerz","pion cerenov photon arrival time ns ECAL PD2",finenbin,timemin,timemaxz);
  TH1F *ehcalpd1scintz = new TH1F("ehcalpd1scintz","elec scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
  TH1F *ehcalpd1cerz = new TH1F("ehcalpd1cerz","elec cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
  TH1F *phcalpd1scintz = new TH1F("phcalpd1scintz","pion scint photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
  TH1F *phcalpd1cerz = new TH1F("phcalpd1cerz","pion cerenov photon arrival time ns HCAL scint fiber",finenbin,timemin,timemaxz);
  TH1F *ehcalpd2scintz = new TH1F("ehcalpd2scintz","elec scint photon arrival time ns HCAL quartz fiber",finenbin,timemin,timemaxz);
  TH1F *ehcalpd2cerz = new TH1F("ehcalpd2cerz","elec cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);
  TH1F *phcalpd2scintz = new TH1F("phcalpd2scintz","pion scint photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);
  TH1F *phcalpd2cerz = new TH1F("phcalpd2cerz","pion cerenov photon arrival time ns quartz fiber",finenbin,timemin,timemaxz);


  //****************************************************************************************************************************
  // process electrons

  TFile* ef = TFile::Open(einputfilename);
  TTree* et = (TTree*)ef->Get("EVENT;1");
  if(doecal) b_ecal = et->GetBranch(ECALleaf);
  if(dohcal) b_hcal = et->GetBranch(HCALleaf);
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");

  // get number of events in the file
  b_mc= et->GetBranch("MCParticles");
  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<std::endl<<std::endl<<"num_evt for electron file is  "<<num_evt<<std::endl;

  float meanscinEcal(0.),meanscinHcal(0.),meancerEcal(0),meancerHcal(0),egEcal(0.),egHcal(0.);
  float meaneecaltimecut(0),meanehcaltimecut(0),meanerelecaltimecut(0),meanerelhcaltimecut(0);
  float meanSEcal(0.),meanSHcal(0.),meanCEcal(0.),meanCHcal(0.);
  if(num_evt>0) {

    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);

    // first pass through file for rough calibration

    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<std::endl<<"  event number rough calibration "<<ievt<<std::endl;
      getMeanPhot(mapsampcalslice, gendet, ievt, doecal, dohcal, hcaltype, b_ecal,b_hcal, ecalhits, hcalhits, meanscinEcal, meanscinHcal, meancerEcal, meancerHcal, timecut, meaneecaltimecut, meanehcaltimecut, meanerelecaltimecut, meanerelhcaltimecut);
    }

    std::cout<<std::endl<<"done with getMeanPhot"<<std::endl<<std::endl;

    meanscinEcal=meanscinEcal/num_evt;
    meanscinHcal=meanscinHcal/num_evt;
    meancerEcal=meancerEcal/num_evt;
    meancerHcal=meancerHcal/num_evt;
    meaneecaltimecut=meaneecaltimecut/num_evt;
    meanehcaltimecut=meanehcaltimecut/num_evt;
    meanerelecaltimecut=meanerelecaltimecut/num_evt;
    meanerelhcaltimecut=meanerelhcaltimecut/num_evt;
    std::cout<<"mean scint ecal is "<<meanscinEcal/1000.<<std::endl;
    std::cout<<"mean scint hcal is "<<meanscinHcal/1000.<<std::endl;
    std::cout<<"mean cer ecal is "<<meancerEcal/1000.<<std::endl;
    std::cout<<"mean cer hcal is "<<meancerHcal/1000.<<std::endl<<std::endl;

    std::cout<<"mean e ecal timecut is "<<meaneecaltimecut/1000.<<std::endl;
    std::cout<<"mean e hcal timecut is "<<meanehcaltimecut/1000.<<std::endl;
    std::cout<<"mean rel ecal timecut is "<<meanerelecaltimecut/1000.<<std::endl;
    std::cout<<"mean rel hcal timecut is "<<meanerelhcaltimecut/1000.<<std::endl;


    // refine calibration

    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<std::endl<<"event number calibration refinement is "<<ievt<<std::endl;
      CalibRefine(mapsampcalslice, gendet, ievt, doecal, dohcal, hcaltype, b_ecal,b_hcal, ecalhits, hcalhits, meanscinEcal, meanscinHcal, meancerEcal, meancerHcal, CalEcalncer, CalEcalnscint, CalHcalncer, CalHcalnscint);
    }


    if(doecal ) {

      TF1 *gs = new TF1("gs", "gaus", 0, 1.5);

      arms = TMath::Max(CalEcalncer->GetRMS(),0.02); // keep at least +-0.02 (x1.5) width for fit for stability
      //arms = CalEcalncer->GetRMS();
      amean = CalEcalncer->GetMean();
      gs->SetParameter(1, amean);
      gs->SetParLimits(1, amean-0.02, amean+0.02); // keep mu within a "reasonable" range
      CalEcalncer->Fit("gs","R0L","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitEcalncer = (TF1*)CalEcalncer->GetListOfFunctions()->FindObject("gs");
      Double_t Ecalncer_p0= fitEcalncer->GetParameter(0);
      Double_t Ecalncer_p1= fitEcalncer->GetParameter(1);
      Double_t Ecalncer_p2= fitEcalncer->GetParameter(2);
      std::cout<<std::endl;
      std::cout<<"Ecal cer refine calib fit params "<<Ecalncer_p0<<" "<<Ecalncer_p1<<" "<<Ecalncer_p2<<std::endl;
      std::cout<<std::endl;
      meancerEcal=meancerEcal/Ecalncer_p1;

      arms = TMath::Max(CalEcalncer->GetRMS(),0.02); // keep at least +-0.02 (x1.5) width for fit for stability
      //arms = CalEcalnscint->GetRMS();
      amean = CalEcalnscint->GetMean();
      gs->SetParameter(1, amean);
      gs->SetParLimits(1, amean-0.02, amean+0.02); // keep mu within a "reasonable" range
      CalEcalnscint->Fit("gs","R0L","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitEcalnscint = (TF1*)CalEcalnscint->GetListOfFunctions()->FindObject("gs");
      Double_t Ecalnscint_p0= fitEcalnscint->GetParameter(0);
      Double_t Ecalnscint_p1= fitEcalnscint->GetParameter(1);
      Double_t Ecalnscint_p2= fitEcalnscint->GetParameter(2);
      meanscinEcal=meanscinEcal/Ecalnscint_p1;
      std::cout<<std::endl;
      std::cout<<"Ecal scint refine calib fit params "<<Ecalnscint_p0<<" "<<Ecalnscint_p1<<" "<<Ecalnscint_p2<<std::endl;
      std::cout<<std::endl;

    }

    if(dohcal) {
      arms = CalHcalncer->GetRMS();
      amean = CalHcalncer->GetMean();
      CalHcalncer->Fit("gaus","R0","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitHcalncer = (TF1*)CalHcalncer->GetListOfFunctions()->FindObject("gaus");
      Double_t Hcalncer_p0= fitHcalncer->GetParameter(0);
      Double_t Hcalncer_p1= fitHcalncer->GetParameter(1);
      Double_t Hcalncer_p2= fitHcalncer->GetParameter(2);
      meancerHcal=meancerHcal/Hcalncer_p1;
      std::cout<<std::endl;
      std::cout<<"Hcal cer refine calib fit params "<<Hcalncer_p0<<" "<<Hcalncer_p1<<" "<<Hcalncer_p2<<std::endl;
      std::cout<<std::endl;

      arms = CalHcalnscint->GetRMS();
      amean = CalHcalnscint->GetMean();
      CalHcalnscint->Fit("gaus","R0","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitHcalnscint = (TF1*)CalHcalnscint->GetListOfFunctions()->FindObject("gaus");
      Double_t Hcalnscint_p0= fitHcalnscint->GetParameter(0);
      Double_t Hcalnscint_p1= fitHcalnscint->GetParameter(1);
      Double_t Hcalnscint_p2= fitHcalnscint->GetParameter(2);
      meanscinHcal=meanscinHcal/Hcalnscint_p1;
      std::cout<<std::endl;
      std::cout<<"Hcal scint refine calib fit params "<<Hcalnscint_p0<<" "<<Hcalnscint_p1<<" "<<Hcalnscint_p2<<std::endl;
      std::cout<<std::endl;

    }

    if(doplots) {
      TCanvas* z1 = new TCanvas();
      SCEDraw1(z1,"z1",CalEcalncer,"junkz1.png",0);
      TCanvas* z2 = new TCanvas();
      SCEDraw1(z2,"z2",CalEcalnscint,"junkz2.png",0);
      TCanvas* z3 = new TCanvas();
      SCEDraw1(z3,"z3",CalHcalncer,"junkz3.png",0);
      TCanvas* z4 = new TCanvas();
      SCEDraw1(z4,"z4",CalHcalnscint,"junkz4.png",0);

    }


    // now that have calibration, do the real work for the electrons
    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<"event number second is "<<ievt<<std::endl;

      float eesum(0.),eesumcal(0.),eesumem(0.),eesumair(0.),eesumdead(0.),eesumcrystal(0.),eesumPDe(0.),eesumfiber1(0.),eesumfiber2(0.),eesumabs(0.),eesumPDh(0.),eesumedge(0.),eesumedgerel(0.),necertotecal(0.),nescinttotecal(0.),necertothcal(0.),nescinttothcal(0.),eecaltimecut(0.),ehcaltimecut(0.),erelecaltimecut(0.),erelhcaltimecut(0.),eesumairem(0.),eesumdeadem(0.),eesumcrystalem(0.),eesumPDeem(0.),eesumfiber1em(0.),eesumfiber2em(0.),eesumabsem(0.),eesumPDhem(0.);
      int nine(0),ninh(0);

      fillfill=1;
      getStuff(mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits, timecut,fillfill,eesum,eesumcal,eesumem,eesumair,eesumdead,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,eesumPDh,eesumairem,eesumdeadem,eesumcrystalem,eesumPDeem,eesumfiber1em,eesumfiber2em,eesumabsem,eesumPDhem,eesumedge,eesumedgerel,necertotecal,nescinttotecal,necertothcal,nescinttothcal,eecaltimecut, ehcaltimecut,erelecaltimecut,erelhcaltimecut,  nine,ninh);
      if(fillfill==1) FillTime(mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits, timecut,eecaltime,ehcaltime,eecalpd1scint,eecalpd1cer,eecalpd2scint,eecalpd2cer,ehcalpd1scint,ehcalpd1cer,ehcalpd2scint,ehcalpd2cer,eecalpd1scintz,eecalpd1cerz,eecalpd2scintz,eecalpd2cerz,ehcalpd1scintz,ehcalpd1cerz,ehcalpd2scintz,ehcalpd2cerz);

      float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge+eesumdead;
      float eedepcal=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumdead;
      float nonconse=(beamE-eachecks)/beamE;

      egEcal+=(eesumcrystal)/beamE;
      egHcal+=(eesumfiber1+eesumfiber2)/beamE;

      heesumcal->Fill(eesumcal/beamE);
      heesumemcal->Fill(eesumem/beamE);
      if(eesumem>0) {
	hefff->Fill(eesumem/eesumcal);

      } else {
	std::cout<<"esumemcal is zero "<<std::endl;
      }
      ehcEcalE->Fill(eesumcrystal/beamE);
      ehcHcalE->Fill((eesumfiber1+eesumfiber2)/eedepcal);
      ehcHcalE1->Fill(eesumfiber1/eedepcal);
      ehcHcalE2->Fill(eesumfiber2/eedepcal);
      ehcEdgeE->Fill(eesumedge/beamE);
      if(eesumedge>0) ehcEdgeRelf->Fill(eesumedgerel/eesumedge);
      ehcEdgeR->Fill((beamE-eesumedge)/beamE);
      ehcnonconsE->Fill(nonconse);
      ehcEcalncer->Fill(necertotecal/meancerEcal);
      ehcHcalncer->Fill(necertothcal/meancerHcal);
      ehcEcalnscint->Fill(nescinttotecal/meanscinEcal);
      ehcHcalnscint->Fill(nescinttothcal/meanscinHcal);
      ehcHcalf1f2->Fill(eesumfiber1/1000.,eesumfiber2/1000.);
      if((eesumfiber1+eesumfiber2)>0) ehaphcal->Fill((eesumfiber1+eesumfiber2)/(eesumabs+eesumfiber1+eesumfiber2));
      eheest->Fill((eesumcrystal+(eesumfiber1+eesumfiber2))/beamE);
      float ttt=nescinttotecal/meanscinEcal;
      float ttt2=necertotecal/meancerEcal;
      float tty=nescinttothcal/meanscinHcal;
      float tty2=necertothcal/meancerHcal;
      if( (ttt>0.1)&& (ttt2>0.2) )
      ehcEcalNsNc->Fill(ttt,ttt2);
      if( (tty>0.1)&& (tty2>0.2) )
      ehcHcalNsNc->Fill(tty,tty2);
      ehcEcalMarco->Fill(ttt2/ttt,ttt2);
      ehcHcalMarco->Fill(tty2/tty,tty2);
      ehetrue->Fill(eachecks/beamE);
      hedepcal->Fill(eedepcal/beamE);
      enscvni->Fill(nine+ninh,eedepcal/beamE);

      if(ievt<SCECOUNT) {
	std::cout<<std::endl<<std::endl<<"GETSTUFF electrons full calorimeter"<<std::endl;
	std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
	std::cout<<std::endl<<std::endl<<"total energy deposit "<<eesum/1000.<<std::endl;
	std::cout<<"       cal total energy deposit "<<eesumcal/1000.<<std::endl;
	std::cout<<"       cal EM total energy deposit "<<eesumem/1000.<<std::endl;
	std::cout<<"       in air "<<eesumair/1000.<<std::endl;
	std::cout<<"       in ecal dead material "<<eesumdead/1000.<<std::endl;
	std::cout<<"       in photodetector ecal "<<eesumPDe/1000.<<std::endl;
	std::cout<<"       in crystal "<<eesumcrystal/1000.<<std::endl;
	std::cout<<"       in fiber1 "<<eesumfiber1/1000.<<std::endl;
	std::cout<<"       in fiber2 "<<eesumfiber2/1000.<<std::endl;
	std::cout<<"       in absorber "<<eesumabs/1000.<<std::endl;
	std::cout<<"       in photodetect hcal "<<eesumPDh/1000.<<std::endl;
	std::cout<<"       edge detector "<<eesumedge/1000.<<std::endl;
	std::cout<<"       sum individual "<<eachecks/1000.<<std::endl;
	std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	std::cout<<"   ratio to incident energy "<<eachecks/beamE<<std::endl;
	std::cout<<"total number of cherenkov ecal is "<<necertotecal<<std::endl;
	std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl;
	std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
	std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl;
	std::cout<<"number inelastic is "<<nine+ninh<<std::endl;
	std::cout<<std::endl;
      }

    }  //end loop over events
  }  // end if no events
  ef->Close();
  std::cout<<"done with getstuff electrons"<<std::endl<<std::endl;

  egEcal/=num_evt;
  egHcal/=num_evt;
  std::cout<<std::endl<<"!! sampling fractions for electrons"<<std::endl;
  std::cout<<"   ecal sampling fraction "<<egEcal<<std::endl;
  std::cout<<"   hcal sampling fraction "<<egHcal<<std::endl;

  // get normalization for f, which is the shower em fraction, from electron plots
  std::cout<<std::endl<<" getting the em shower fraction "<<std::endl;
  arms = heesumemcal->GetRMS();
  amean = heesumemcal->GetMean();
  heesumemcal->Fit("gaus","R0","",amean-1.5*arms,amean+1.5*arms);
  float fnorm = heesumemcal->GetMean();
  TF1 *fitheesumemcal = (TF1*)heesumemcal->GetListOfFunctions()->FindObject("gaus");
  fnorm= fitheesumemcal->GetParameter(1);
  std::cout<<"!!  fnorm is "<<fnorm<<std::endl;

  //****************************************************************************************************************************
  // for calibration of hcal,  if have both an ecal and hcal


  if(doecal&&dohcal ) {
    TFile* efa = TFile::Open(hcalonlyefilename);
    TTree* eta = (TTree*)efa->Get("EVENT;1");
    b_mc= eta->GetBranch("MCParticles");
    b_hcal = eta->GetBranch(HCALleaf);
    std::cout<<std::endl<<std::endl<<"b_mc b_hcal are "<<b_mc<<" "<<b_hcal<<std::endl;
    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<std::endl<<std::endl<<"num_evt for electron hcal calibration file is  "<<num_evt<<std::endl;

  // loop over events
    if(num_evt>0) {
      CalHits* ecalhitsa = new CalHits();
      CalHits* hcalhitsa = new CalHits();
      CalHits* edgehitsa = new CalHits();
      if(dohcal) b_hcal->SetAddress(&hcalhitsa);
      if(doedge) b_edge = eta->GetBranch("EdgeDetNoSegment");
      if(doedge) b_edge->SetAddress(&edgehitsa);
            std::cout<<"7"<<std::endl;

      std::cout<<" branches set"<<std::endl;
      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<std::endl<<"event number first pass is "<<ievt<<std::endl;
	getMeanPhot(mapsampcalslice, gendet, ievt, 0, dohcal, hcaltype, b_ecal,b_hcal, ecalhitsa, hcalhitsa, meanscinEcal, meanscinHcal, meancerEcal, meancerHcal,timecut,meaneecaltimecut,meanehcaltimecut,meanerelecaltimecut,meanerelhcaltimecut);
      }
      std::cout<<"done with getMeanPhot for hcal calibration file"<<std::endl;
      meanscinHcal=meanscinHcal/num_evt;
      meancerHcal=meancerHcal/num_evt;
      std::cout<<"mean scint hcal is "<<meanscinHcal/1000.<<std::endl;
      std::cout<<"mean cer hcal is "<<meancerHcal/1000.<<std::endl;

    // now that have calibration, do the real work for the electrons
      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<"event number second is "<<ievt<<std::endl;

	float eesum(0.),eesumcal(0.),eesumem(0.),eesumair(0.),eesumdead(0.),eesumcrystal(0.),eesumPDe(0.),eesumfiber1(0.),eesumfiber2(0.),eesumabs(0.),eesumPDh(0.),eesumedge(0.),eesumedgerel(0.),necertotecal(0.),nescinttotecal(0.),necertothcal(0.),nescinttothcal(0.),eecaltimecut(0.),ehcaltimecut(0.),erelecaltimecut(0.),erelhcaltimecut(0.),eesumairem(0.),eesumdeadem(0.),eesumcrystalem(0.),eesumPDeem(0.),eesumfiber1em(0.),eesumfiber2em(0.),eesumabsem(0.),eesumPDhem(0.);
	int nine(0),ninh(0);

	fillfill=1;
	getStuff(mapsampcalslice,  gendet, ievt, 0, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhitsa,hcalhitsa,edgehitsa, timecut,fillfill,eesum,eesumcal,eesumem,eesumair,eesumdead,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,eesumPDh,eesumairem,eesumdeadem,eesumcrystalem,eesumPDeem,eesumfiber1em,eesumfiber2em,eesumabsem,eesumPDhem,eesumedge,eesumedgerel,necertotecal,nescinttotecal,necertothcal,nescinttothcal,eecaltimecut, ehcaltimecut,erelecaltimecut,erelhcaltimecut,  nine,ninh);


	float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge+eesumdead;
	float eedepcal=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumdead;
	float nonconse=(beamE-eachecks)/beamE;

	egHcal+=(eesumfiber1+eesumfiber2)/beamE;

	if(ievt<SCECOUNT) {
	  std::cout<<"GETSTUFF electrons hcal calorimeter"<<std::endl;
	  std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
	  std::cout<<std::endl<<std::endl<<"total energy deposit "<<eesum/1000.<<std::endl;
	  std::cout<<"       cal total energy deposit "<<eesumcal/1000.<<std::endl;
	  std::cout<<"       cal EM total energy deposit "<<eesumem/1000.<<std::endl;
	  std::cout<<"       in air "<<eesumair/1000.<<std::endl;
	  std::cout<<"       in ecal dead material "<<eesumdead/1000.<<std::endl;
	  std::cout<<"       in photodetector ecal "<<eesumPDe/1000.<<std::endl;
	  std::cout<<"       in crystal "<<eesumcrystal/1000.<<std::endl;
	  std::cout<<"       in fiber1 "<<eesumfiber1/1000.<<std::endl;
	  std::cout<<"       in fiber2 "<<eesumfiber2/1000.<<std::endl;
	  std::cout<<"       in absorber "<<eesumabs/1000.<<std::endl;
	  std::cout<<"       in photodetect hcal "<<eesumPDh/1000.<<std::endl;
	  std::cout<<"       edge detector "<<eesumedge/1000.<<std::endl;
	  std::cout<<"       sum individual "<<eachecks/1000.<<std::endl;
	  std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	  std::cout<<"   ratio to incident energy "<<eachecks/beamE<<std::endl;
	  std::cout<<"total number of cherenkov ecal is "<<necertotecal<<std::endl;
	  std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl;
	  std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
	  std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl;
	  std::cout<<"number inelastic is "<<nine+ninh<<std::endl;
	  std::cout<<std::endl;
	}
      }  //end event loop
      egHcal/=num_evt;
      std::cout<<std::endl<<"!! sampling fractions for electrons"<<std::endl;
      std::cout<<"   ecal sampling fraction "<<egEcal<<std::endl;
      std::cout<<"   hcal sampling fraction "<<egHcal<<std::endl;


    }
    efa->Close();
    std::cout<<"done with with electron calibration of hcal"<<std::endl;
  }  // end doecal&&dohcal

  //****************************************************************************************************************************
  // process pions for full calorimeter
  std::cout<<std::endl<<std::endl<<"!! processing pions for full calorimeter "<<std::endl;

  TFile* pif = TFile::Open(piinputfilename);
  TTree* pit = (TTree*)pif->Get("EVENT;1");
  if(pif==0) std::cout<<" no file "<<std::endl;
  if(pit==0) std::cout<<" no event "<<std::endl;
  std::cout<<"pion file open"<<std::endl;

  b_mc= pit->GetBranch("MCParticles");
  if(doecal) b_ecal = pit->GetBranch(ECALleaf);
  if(dohcal) b_hcal = pit->GetBranch(HCALleaf);
  if(doedge) b_edge = pit->GetBranch("EdgeDetNoSegment");
  std::cout<<"pion branches found"<<std::endl;
  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<"num_evt for pion file is  "<<num_evt<<std::endl;

  // loop over events

  float avedepecal(0.),avedephcal(0.);
  if(num_evt>0) {
    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);


    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<"event number pion is "<<ievt<<std::endl;

      float pesum(0.),pesumcal(0.),pesumem(0.),pesumair(0.),pesumdead(0.),pesumcrystal(0.),pesumPDe(0.),pesumfiber1(0.),pesumfiber2(0.),pesumabs(0.),pesumPDh(0.),pesumedge(0.),pesumedgerel(0.),npcertotecal(0.),npscinttotecal(0.),npcertothcal(0.),npscinttothcal(0.),pecaltimecut(0.),phcaltimecut(0.),prelecaltimecut(0.),prelhcaltimecut(0.),pesumairem(0.),pesumdeadem(0.),pesumcrystalem(0.),pesumPDeem(0.),pesumfiber1em(0.),pesumfiber2em(0.),pesumabsem(0.),pesumPDhem(0.);
      int nine(0),ninh(0);
      getStuff(mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,timecut,fillfill,pesum,pesumcal,pesumem,pesumair,pesumdead,pesumcrystal,pesumPDe,pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumairem,pesumdeadem,pesumcrystalem,pesumPDeem,pesumfiber1em,pesumfiber2em,pesumabsem,pesumPDhem,pesumedge,pesumedgerel,npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,pecaltimecut, phcaltimecut,prelecaltimecut,prelhcaltimecut,nine,ninh);
      fillfill=1;
      if(fillfill==1) FillTime(mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,timecut,piecaltime,pihcaltime,pecalpd1scint,pecalpd1cer,pecalpd2scint,pecalpd2cer,phcalpd1scint,phcalpd1cer,phcalpd2scint,phcalpd2cer,pecalpd1scintz,pecalpd1cerz,pecalpd2scintz,pecalpd2cerz,phcalpd1scintz,phcalpd1cerz,phcalpd2scintz,phcalpd2cerz);

      // gamma fraction from em fraction
      float pfff=0.;
      pfff=(pesumem/beamE/fnorm);

      //edge fraction
      float pedgeff=0.;
      pedgeff=pesumedge/(beamE-(pesumem/fnorm));

      // non cons fraction
      float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumedge+pesumdead;
      float pedepcal=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumdead;
      float edepEcal = pesumcrystal;
      avedepecal+=pesumcrystal;
      avedephcal+=pesumfiber1+pesumfiber2+pesumabs;
      float nonconsE=(beamE-pachecks);
      float noncons=0.;
      noncons=nonconsE/(beamE-(pesumem/fnorm));

      // calculate average energy deposited in each section

      if(domissCorr) {
	float ContainedEnergy=beamE-pesumedge;
	float ContainedFrac=ContainedEnergy/beamE;
	npcertotecal=npcertotecal/ContainedFrac;
	npscinttotecal=npscinttotecal/ContainedFrac;
	npcertothcal=npcertothcal/ContainedFrac;
	npscinttothcal=npscinttothcal/ContainedFrac;
      }
      float ppp = pesumedge/beamE;
      // std::cout<<"getstuff pion ppp edge cut are "<<ppp<<" "<<edgecut<<std::endl;

      if((!doedgecut)||(ppp<edgecut) ) {
	//	std::cout<<"passed edge cut "<<std::endl;

	phcnonconsE->Fill(noncons);
	phcnonconsE2->Fill(nonconsE/pesumcal);
	phcnonconsE3->Fill(nonconsE/(pesumcal-pesumem));
	phcnonconsE4->Fill(nonconsE/(pesumcal-(pesumem/fnorm)));
	mes1Ecal->Fill(pfff,pedgeff);
	mes2Ecal->Fill(nine,pedgeff);
	mes3Ecal->Fill(pedgeff,npcertotecal/meancerEcal);
	mes4Ecal->Fill(pedgeff,npscinttotecal/meanscinEcal);
	mes1Hcal->Fill(pfff,pedgeff);
	mes2Hcal->Fill(ninh,pedgeff);
	mes3Hcal->Fill(pedgeff,npcertothcal/meancerHcal);
	mes4Hcal->Fill(pedgeff,npscinttothcal/meanscinHcal);
	hpfff->Fill(pfff);
	hpfff2->Fill(pesumem/beamE);
	hpfff3->Fill(pesumem/pesumcal);
	if(pesumabs>0) hpfffabs->Fill(pesumabsem/pesumabs);
	if((pesumfiber1+pesumfiber2)>0) hpffffib->Fill((pesumfiber1em+pesumfiber2em)/(pesumfiber1+pesumfiber2));
	hpesumcal->Fill(pesumcal/beamE);
	hpesumemcal->Fill(pesumem/beamE);
	hfnscinEcal->Fill(pfff,npscinttotecal/meanscinEcal);
	hfncerEcal->Fill(pfff,npcertotecal/meancerEcal);
	hfnscinHcal->Fill(pfff,npscinttothcal/meanscinHcal);
	hfncerHcal->Fill(pfff,npcertothcal/meancerHcal);
	phcEcalE->Fill(pesumcrystal/beamE);
	phcHcalE->Fill((pesumfiber1+pesumfiber2)/pedepcal);
	phcHcalE1->Fill(pesumfiber1/pedepcal);
	phcHcalE2->Fill(pesumfiber2/pedepcal);
	phcHcalvfE1->Fill(pfff,pesumfiber1/pedepcal);
	phcHcalvfE2->Fill(pfff,pesumfiber2/pedepcal);
	phcEdgeE->Fill(pedgeff);
	if(pesumedge>0) phcEdgeRelf->Fill(pesumedgerel/pesumedge);
	phcEdgeR->Fill((beamE-pesumedge)/beamE);
	phcEandHcalncer->Fill((npcertotecal/meancerEcal)+(npcertothcal/meancerHcal));
	phcEandHcalnscint->Fill((npscinttotecal/meanscinEcal)+(npscinttothcal/meanscinHcal));
	phcEcalncer->Fill(npcertotecal/meancerEcal);
	phcEcalncer2->Fill(npcertotecal/meancerEcal);
	if(pesumcrystal>0.1*20000) phcEcalncer3->Fill(npcertotecal/(meancerEcal*edepEcal/beamE));
	if(ievt<SCECOUNT) std::cout<<" pesumcrystal npcertotecal meancerecal edepEcal beamE fill "<<pesumcrystal<<" "<<npcertotecal<<" "<<meancerEcal<<" "<<edepEcal<<" "<<beamE<<" "<<(npcertotecal/(meancerEcal*edepEcal/beamE))<<std::endl;
	phcHcalncer->Fill(npcertothcal/meancerHcal);
	phcHcalncer2->Fill(npcertothcal/meancerHcal);
	phcEcalnscint->Fill(npscinttotecal/meanscinEcal);
	phcEcalnscint2->Fill(npscinttotecal/meanscinEcal);
	if(pesumcrystal>0.1*20000) phcEcalnscint3->Fill(npscinttotecal/(meanscinEcal*edepEcal/beamE));
	if(ievt<SCECOUNT) std::cout<<" pesumcrystal npscinttotecal meanscinecal edepEcal beamE fill "<<pesumcrystal<<" "<<npscinttotecal<<" "<<meanscinEcal<<" "<<edepEcal<<" "<<beamE<<" "<<(npscinttotecal/(meanscinEcal*edepEcal/beamE))<<std::endl;
	phcHcalnscint->Fill(npscinttothcal/meanscinHcal);
	phcHcalnscint2->Fill(npscinttothcal/meanscinHcal);
	enonconsEcalcer->Fill(noncons,npcertotecal/meancerEcal);
	enonconsEcalscin->Fill(noncons,npscinttotecal/meanscinEcal);
	enonconsHcalcer->Fill(noncons,npcertothcal/meancerHcal);
	enonconsHcalscin->Fill(noncons,npscinttothcal/meanscinHcal);
	enonconsvnni->Fill(nine+ninh,noncons);
	enonconsvf->Fill(pfff,noncons);
	hcnonconsvesc->Fill(noncons,pedgeff);
	phcHcalf1f2->Fill(pesumfiber1/1000.,pesumfiber2/1000.);
	if((pesumfiber1+pesumfiber2)>0) phaphcal->Fill((pesumfiber1+pesumfiber2)/(pesumabs+pesumfiber1+pesumfiber2));
	pheest->Fill((pesumcrystal+(pesumfiber1+pesumfiber2))/beamE);

	float rrr=0.;float rrr2=0.;float rrx=0.;float rrx2=0.;
	if(meanscinEcal>0) rrr=npscinttotecal/meanscinEcal;
	if(meancerEcal>0) rrr2=npcertotecal/meancerEcal;
	if(meanscinHcal>0) rrx=npscinttothcal/meanscinHcal;
	if(meancerHcal>0) rrx2=npcertothcal/meancerHcal;
	phcEcalNsNc->Fill(rrr,rrr2);
	phcHcalNsNc->Fill(rrx,rrx2);
	phcEcalMarco->Fill(rrr2/rrr,rrr2);
	phcHcalMarco->Fill(rrx2/rrx,rrx2);

	float rrrtc=0.;float rrr2tc=0.;float rrxtc=0.;float rrx2tc=0.;
	if(meaneecaltimecut>0) rrrtc=pecaltimecut/meaneecaltimecut;
	if(meanerelecaltimecut>0) rrr2tc=prelecaltimecut/meanerelecaltimecut;
	if(meanehcaltimecut>0) rrxtc=phcaltimecut/meanehcaltimecut;
	if(meanerelhcaltimecut>0) rrx2tc=prelhcaltimecut/meanerelhcaltimecut;
	phcEcalNsNctc->Fill(rrrtc,rrr2tc);
	phcHcalNsNctc->Fill(rrxtc,rrx2tc);

	hpdepcal->Fill(pedepcal/beamE);
	phetrue->Fill(pachecks/beamE);
	pinalvni->Fill(nine+ninh,pedepcal/beamE);
	pinscvni->Fill(nine+ninh,rrr+rrx);
	pincevni->Fill(nine+ninh,rrr2+rrx2);


	if(ievt<SCECOUNT) {
	  std::cout<<std::endl<<std::endl;
	  std::cout<<"GETSTUFF pions"<<std::endl;
	  std::cout<<" phcaltimecut is "<<phcaltimecut/1000.<<std::endl;
	  std::cout<<" prelhcaltimecut is "<<prelhcaltimecut/1000.<<std::endl;
	  std::cout<<"total energy deposit "<<pesum/1000.<<std::endl;
	  std::cout<<"       cal total energy deposit "<<pesumcal/1000.<<std::endl;
	  std::cout<<"       cal EM total energy deposit "<<pesumem/1000.<<std::endl;
	  std::cout<<"       in air "<<pesumair/1000.<<std::endl;
	  std::cout<<"       in ecal dead "<<pesumdead/1000.<<std::endl;
	  std::cout<<"       in photodetector ecal "<<pesumPDe/1000.<<std::endl;
	  std::cout<<"       in crystal "<<pesumcrystal/1000.<<std::endl;
	  std::cout<<"       in fiber1 "<<pesumfiber1/1000.<<std::endl;
	  std::cout<<"       in fiber2 "<<pesumfiber2/1000.<<std::endl;
	  std::cout<<"       sum of fibers "<<(pesumfiber1+pesumfiber2)/1000.<<std::endl;
	  std::cout<<"       in absorber "<<pesumabs/1000.<<std::endl;
	  std::cout<<"       in photodetect hcal "<<pesumPDh/1000.<<std::endl;
	  std::cout<<"       edge detector "<<pesumedge/1000.<<std::endl;
	  std::cout<<"       sum individual "<<pachecks/1000.<<std::endl;
	  std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	  std::cout<<"   ratio to incident energy "<<pachecks/beamE<<std::endl;
	  std::cout<<"total number of cherenkov ecal is "<<npcertotecal<<std::endl;
	  std::cout<<"total number of scintillator ecal is "<<npscinttotecal<<std::endl;
	  std::cout<<"total number of cherenkov hcal is "<<npcertothcal<<std::endl;
	  std::cout<<"total number of scintillator hcal is "<<npscinttothcal<<std::endl<<std::endl;
	  std::cout<<" number of inelastic is "<<nine+ninh<<std::endl;
	}
      } //edge cut
    }  //end loop over events
  }  // end if no events
    // close pion file
  pif->Close();


  avedepecal/=num_evt;
  avedephcal/=num_evt;
  std::cout<<"average deposit in calorimeter is "<<avedepecal<<" in the ECAL and "<<avedephcal<<" in the HCAL"<<std::endl<<std::endl;

    // get mean of S and C for covariance calculation
  meanSEcal=phcEcalnscint->GetMean();
  meanSHcal=phcHcalnscint->GetMean();
  meanCEcal=phcEcalncer->GetMean();
  meanCHcal=phcHcalncer->GetMean();
  std::cout<<"meanSEcal is "<<meanSEcal<<std::endl;
  std::cout<<"meanCEcal is "<<meanCEcal<<std::endl;
  std::cout<<"meanSHcal is "<<meanSHcal<<std::endl;
  std::cout<<"meanCHcal is "<<meanCHcal<<std::endl;



  //****************************************************************************************************************************
  // process pions for hcal only calorimeter
  if(doecal&&dohcal) {
    std::cout<<std::endl<<std::endl<<"!! processing pions for hcal only calorimeter "<<std::endl;
    phcHcalnscint2->Reset();

    TFile* pifa = TFile::Open(hcalonlypifilename);
    TTree* pita = (TTree*)pifa->Get("EVENT;1");
    if(pifa==0) std::cout<<" no file "<<std::endl;
    if(pita==0) std::cout<<" no event "<<std::endl;
    std::cout<<"hcal only pion file open"<<std::endl;

    b_mc= pita->GetBranch("MCParticles");
    //if(doecal) b_ecal = pita->GetBranch(ECALleaf);
    if(dohcal) b_hcal = pita->GetBranch(HCALleaf);
    if(doedge) b_edge = pita->GetBranch("EdgeDetNoSegment");
    std::cout<<"pion branches found"<<std::endl;
    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<"num_evt for hcal only pion file is  "<<num_evt<<std::endl;

  // loop over events

    if(num_evt>0) {
      CalHits* ecalhitsa = new CalHits();
      //if(doecal) b_ecal->SetAddress(&ecalhitsa);
      CalHits* hcalhitsa = new CalHits();
      if(dohcal) b_hcal->SetAddress(&hcalhitsa);
      CalHits* edgehitsa = new CalHits();
      if(doedge) b_edge->SetAddress(&edgehitsa);


      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<"event number hcal only pion is "<<ievt<<std::endl;

	float pesum(0.),pesumcal(0.),pesumem(0.),pesumair(0.),pesumdead(0.),pesumcrystal(0.),pesumPDe(0.),pesumfiber1(0.),pesumfiber2(0.),pesumabs(0.),pesumPDh(0.),pesumedge(0.),pesumedgerel(0.),npcertotecal(0.),npscinttotecal(0.),npcertothcal(0.),npscinttothcal(0.),pecaltimecut(0.),phcaltimecut(0.),prelecaltimecut(0.),prelhcaltimecut(0.),pesumairem(0.),pesumdeadem(0.),pesumcrystalem(0.),pesumPDeem(0.),pesumfiber1em(0.),pesumfiber2em(0.),pesumabsem(0.),pesumPDhem(0.);
	int nine(0),ninh(0);
	// doecal is forced to zero here since the analysis is done for hcalonlypifile
	getStuff(mapsampcalslice,  gendet, ievt, 0, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhitsa,hcalhitsa,edgehitsa,timecut,fillfill,pesum,pesumcal,pesumem,pesumair,pesumdead,pesumcrystal,pesumPDe,pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumairem,pesumdeadem,pesumcrystalem,pesumPDeem,pesumfiber1em,pesumfiber2em,pesumabsem,pesumPDhem,pesumedge,pesumedgerel,npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,pecaltimecut, phcaltimecut,prelecaltimecut,prelhcaltimecut,nine,ninh);

      // gamma fraction from em fraction
	float pfff=0.;
	pfff=(pesumem/beamE/fnorm);
      //edge fraction
	float pedgeff=0.;
	pedgeff=pesumedge/(beamE-(pesumem/fnorm));
      // non cons fraction
	float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumedge+pesumdead;
	float pedepcal=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumdead;
	float nonconsE=(beamE-pachecks);
	float noncons=0.;
	noncons=nonconsE/(beamE-(pesumem/fnorm));

	if(domissCorr) {
	  float ContainedEnergy=beamE-pesumedge;
	  float ContainedFrac=ContainedEnergy/beamE;
	  npcertotecal=npcertotecal/ContainedFrac;
	  npscinttotecal=npscinttotecal/ContainedFrac;
	  npcertothcal=npcertothcal/ContainedFrac;
	  npscinttothcal=npscinttothcal/ContainedFrac;
	}
	float ppp = pesumedge/beamE;
      // std::cout<<"getstuff pion ppp edge cut are "<<ppp<<" "<<edgecut<<std::endl;

	if((!doedgecut)||(ppp<edgecut) ) {
	//	std::cout<<"passed edge cut "<<std::endl;
	  phcHcalnscint2->Fill(npscinttothcal/meanscinHcal);

	  if(ievt<SCECOUNT) {
	    std::cout<<std::endl<<std::endl;
	    std::cout<<"GETSTUFF pions hcal only"<<std::endl;
	    std::cout<<" phcaltimecut is "<<phcaltimecut/1000.<<std::endl;
	    std::cout<<" prelhcaltimecut is "<<prelhcaltimecut/1000.<<std::endl;
	    std::cout<<"total energy deposit "<<pesum/1000.<<std::endl;
	    std::cout<<"       cal total energy deposit "<<pesumcal/1000.<<std::endl;
	    std::cout<<"       cal EM total energy deposit "<<pesumem/1000.<<std::endl;
	    std::cout<<"       in air "<<pesumair/1000.<<std::endl;
	    std::cout<<"       in ecal dead "<<pesumdead/1000.<<std::endl;
	    std::cout<<"       in photodetector ecal "<<pesumPDe/1000.<<std::endl;
	    std::cout<<"       in crystal "<<pesumcrystal/1000.<<std::endl;
	    std::cout<<"       in fiber1 "<<pesumfiber1/1000.<<std::endl;
	    std::cout<<"       in fiber2 "<<pesumfiber2/1000.<<std::endl;
	    std::cout<<"       sum of fibers "<<(pesumfiber1+pesumfiber2)/1000.<<std::endl;
	    std::cout<<"       in absorber "<<pesumabs/1000.<<std::endl;
	    std::cout<<"       in photodetect hcal "<<pesumPDh/1000.<<std::endl;
	    std::cout<<"       edge detector "<<pesumedge/1000.<<std::endl;
	    std::cout<<"       sum individual "<<pachecks/1000.<<std::endl;
	    std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	    std::cout<<"   ratio to incident energy "<<pachecks/beamE<<std::endl;
	    std::cout<<"total number of cherenkov ecal is "<<npcertotecal<<std::endl;
	    std::cout<<"total number of scintillator ecal is "<<npscinttotecal<<std::endl;
	    std::cout<<"total number of cherenkov hcal is "<<npcertothcal<<std::endl;
	    std::cout<<"total number of scintillator hcal is "<<npscinttothcal<<std::endl<<std::endl;
	    std::cout<<" number of inelastic is "<<nine+ninh<<std::endl;
	  }
	} //edge cut
      }  //end loop over events

    }  // num evt >0
    // close pion file
    pifa->Close();

  } // do ecal && do hcal


    // get kappa for calorimeter in preparation for doing dual readout correction
  float hovereecalcer(1.),hovereecalscint(1.),hoverehcalcer(1.),hoverehcalscint(1.);
  if(dodualcorr) {
    if(doecal) {
      std::cout<<"!! calculating kappa for ecal "<<std::endl;
      if(doecal&&dohcal) {
	arms = phcEcalnscint3->GetRMS();
	amean = phcEcalnscint3->GetMean();
      } else {
	arms = phcEcalnscint2->GetRMS();
	amean = phcEcalnscint2->GetMean();
      }
      std::cout<<" hist mean rms are "<<amean<<" "<<arms<<std::endl;
      if(doecal&&dohcal) {
	if (arms==0.) arms=0.05; // protection against rare crash
	phcEcalnscint3->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      } else {
	phcEcalnscint2->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      }
      TF1 *fitphcEcalnscint;
      if(doecal&&dohcal) {
	fitphcEcalnscint= (TF1*)phcEcalnscint3->GetListOfFunctions()->FindObject("gaus");
      } else {
	fitphcEcalnscint= (TF1*)phcEcalnscint2->GetListOfFunctions()->FindObject("gaus");
      }
      Double_t fitphcEcalnscint_p0= fitphcEcalnscint->GetParameter(0);
      Double_t fitphcEcalnscint_p1= fitphcEcalnscint->GetParameter(1);
      Double_t fitphcEcalnscint_p2= fitphcEcalnscint->GetParameter(2);
      hovereecalscint=fitphcEcalnscint_p1;
      std::cout<<"hovereecalscint before correct is "<<hovereecalscint<<std::endl;
      if(twocalecalcorr) hovereecalscint*=beamE/avedepecal;
      std::cout<<"hovereecalscint after correct is "<<hovereecalscint<<std::endl;

      TCanvas* xx3 = new TCanvas();;
      if(doecal&&dohcal) {
	SCEDraw1(xx3,"xx3",phcEcalnscint3,"junkxx3.png",0);
      } else {
	SCEDraw1(xx3,"xx3",phcEcalnscint2,"junkxx3.png",0);
      }

      fitphcEcalnscint->Draw("same");
      std::cout<<std::endl;
      std::cout<<" Ecal nscint fit params "<<fitphcEcalnscint_p0<<" "<<fitphcEcalnscint_p1<<" "<<fitphcEcalnscint_p2<<std::endl;
      std::cout<<std::endl;

      if(doecal&&dohcal ) {
	arms = phcEcalncer3->GetRMS();
	amean = phcEcalncer3->GetMean();
      } else {
	arms = phcEcalncer2->GetRMS();
	amean = phcEcalncer2->GetMean();
      }
      std::cout<<" hist mean rms are "<<amean<<" "<<arms<<std::endl;
      if(doecal&&dohcal) {
	phcEcalncer3->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      } else {
	phcEcalncer2->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      }
      TF1 *fitphcEcalncer;
      if(doecal&&dohcal) {
	fitphcEcalncer= (TF1*)phcEcalncer3->GetListOfFunctions()->FindObject("gaus");
      } else {
	fitphcEcalncer= (TF1*)phcEcalncer2->GetListOfFunctions()->FindObject("gaus");
      }
      Double_t fitphcEcalncer_p0= fitphcEcalncer->GetParameter(0);
      Double_t fitphcEcalncer_p1= fitphcEcalncer->GetParameter(1);
      Double_t fitphcEcalncer_p2= fitphcEcalncer->GetParameter(2);
      hovereecalcer=fitphcEcalncer_p1;
      std::cout<<"hovereecalcer before correct is "<<hovereecalcer<<std::endl;
      if(twocalecalcorr) hovereecalcer*=beamE/avedepecal;
      std::cout<<"hovereecalcer after correct is "<<hovereecalcer<<std::endl;

      TCanvas* xx4 = new TCanvas();;
      if(doecal&&dohcal) {
	SCEDraw1(xx4,"xx4",phcEcalncer3,"junkxx4.png",0);
      } else {
	SCEDraw1(xx4,"xx4",phcEcalncer2,"junkxx4.png",0);
      }
      fitphcEcalncer->Draw("same");
      std::cout<<std::endl;
      std::cout<<" Ecal ncer fit params "<<fitphcEcalncer_p0<<" "<<fitphcEcalncer_p1<<" "<<fitphcEcalncer_p2<<std::endl;
      std::cout<<std::endl;
      kappaEcal= (1-hovereecalscint)/(1.-hovereecalcer);
      std::cout<<" hovereecalcer hovereecalscint "<<hovereecalcer<<" "<<hovereecalscint<<std::endl;
      std::cout<<" kappa ecal is "<<kappaEcal<<std::endl;
    }  // end doecal


    if(dohcal) {
      std::cout<<"!! calculating kappa for hcal "<<std::endl;
      float hoverehcalscint,hoverehcalcer;
      float arms,amean;
      arms = phcHcalnscint2->GetRMS();
      amean = phcHcalnscint2->GetMean();
      std::cout<<"arms amean are "<<arms<<" "<<amean<<std::endl;
      phcHcalnscint2->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitphcHcalnscint = (TF1*)phcHcalnscint2->GetListOfFunctions()->FindObject("gaus");
      Double_t fitphcHcalnscint_p0= fitphcHcalnscint->GetParameter(0);
      Double_t fitphcHcalnscint_p1= fitphcHcalnscint->GetParameter(1);
      Double_t fitphcHcalnscint_p2= fitphcHcalnscint->GetParameter(2);
      std::cout<<" p0 p1 p2 "<<fitphcHcalnscint_p0<<" "<<fitphcHcalnscint_p1<<" "<<fitphcHcalnscint_p2<<std::endl;
      hoverehcalscint=fitphcHcalnscint_p1;

      TCanvas* xx5 = new TCanvas();
      SCEDraw1(xx5,"xx5",phcHcalnscint2,"junkxx5.png",0);
      fitphcHcalnscint->Draw("same");
      std::cout<<std::endl;
      std::cout<<" Hcal nscint fit params "<<fitphcHcalnscint_p0<<" "<<fitphcHcalnscint_p1<<" "<<fitphcHcalnscint_p2<<std::endl;
      std::cout<<std::endl;

      arms = phcHcalncer2->GetRMS();
      amean = phcHcalncer2->GetMean();
      std::cout<<"arms amean are "<<arms<<" "<<amean<<std::endl;
      phcHcalncer2->Fit("gaus","R","",amean-1.5*arms,amean+1.5*arms);
      TF1 *fitphcHcalncer = (TF1*)phcHcalncer2->GetListOfFunctions()->FindObject("gaus");
      Double_t fitphcHcalncer_p0= fitphcHcalncer->GetParameter(0);
      Double_t fitphcHcalncer_p1= fitphcHcalncer->GetParameter(1);
      Double_t fitphcHcalncer_p2= fitphcHcalncer->GetParameter(2);
      std::cout<<" p0 p1 p2 "<<fitphcHcalnscint_p0<<" "<<fitphcHcalnscint_p1<<" "<<fitphcHcalnscint_p2<<std::endl;
      hoverehcalcer=fitphcHcalncer_p1;

      TCanvas* xx6 = new TCanvas();
      SCEDraw1(xx6,"xx6",phcHcalncer2,"junkxx6.png",0);
      fitphcHcalncer->Draw("same");
      std::cout<<std::endl;
      std::cout<<" Hcal ncer fit params "<<fitphcHcalncer_p0<<" "<<fitphcHcalncer_p1<<" "<<fitphcHcalncer_p2<<std::endl;
      std::cout<<std::endl;

      std::cout<<hoverehcalscint<<" "<<hoverehcalcer<<std::endl;
      kappaHcal= (1-hoverehcalscint)/(1.-hoverehcalcer);

      std::cout<<" hoverehcalscint hoverehcalcer kappa hcal are "<<hoverehcalscint<<" "<<hoverehcalcer<<" "<<kappaHcal<<std::endl;
    }  //end dohcal


  // now calculate with dual readout correction
    std::cout<<std::endl<<std::endl<<"!! staring dual readout loop "<<std::endl;
    TFile* pifc = TFile::Open(piinputfilename);
    TTree* pitc = (TTree*)pifc->Get("EVENT;1");
    if(pifc==0) std::cout<<" no file "<<std::endl;
    if(pitc==0) std::cout<<" no event "<<std::endl;
    std::cout<<"pion file open"<<std::endl;

    b_mc= pitc->GetBranch("MCParticles");
    if(doecal) b_ecal = pitc->GetBranch(ECALleaf);
    if(dohcal) b_hcal = pitc->GetBranch(HCALleaf);
    if(doedge) b_edge = pitc->GetBranch("EdgeDetNoSegment");
    std::cout<<"pion branches found"<<std::endl;
    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<"num_evt for pion file is  "<<num_evt<<std::endl;

    if(num_evt>0) {
      CalHits* ecalhitsc = new CalHits();
      if(doecal) b_ecal->SetAddress(&ecalhitsc);
      CalHits* hcalhitsc = new CalHits();
      if(dohcal) b_hcal->SetAddress(&hcalhitsc);
      CalHits* edgehitsc = new CalHits();
      if(doedge) b_edge->SetAddress(&edgehitsc);

      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||((ievt%SCECOUNT2)==0)) std::cout<<"event number pion is "<<ievt<<std::endl;
	float EcorEcal(0),EcorHcal(0),ecaltimecutcor(0),hcaltimecutcor(0),relecaltimecutcor(0),relhcaltimecutcor(0),desumedge(0.),desumedgerel(0.);
	getStuffDualCorr(domissCorr,beamE,mapsampcalslice, gendet, kappaEcal, kappaHcal, meanscinEcal, meancerEcal, meanscinHcal, meancerHcal,ievt,doecal,dohcal, hcaltype, doedge, desumedge,desumedge,b_ecal,b_hcal,b_edge,ecalhitsc,hcalhitsc, edgehitsc, EcorEcal, EcorHcal,timecut, ecaltimecutcor, hcaltimecutcor,relecaltimecutcor,relhcaltimecutcor);
	float pp2 = desumedge/beamE;
	if((!doedgecut)||(pp2<edgecut) ) {
	  phcEcalcorr->Fill(EcorEcal);
	  phcHcalcorr->Fill(EcorHcal);
	  phcEandHcalcorr->Fill(EcorEcal+EcorHcal);
	}
      }  //end loop over events
    }  // end if num_evt>0
  // close pion file
    pifc->Close();

  }  // end dodualcor

  float roughfit=0.4;
  float rmsscale=1.2;

  float rms1,ffm;
  // fits for testing the formula
  ffm = hpfff->GetBinCenter(hpfff->GetMaximumBin());
  hpfff->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
  TF1 *fitffn1 = (TF1*)hpfff->GetListOfFunctions()->FindObject("gaus");
  ffm = fitffn1->GetParameter(1);
  rms1 = rmsscale*fitffn1->GetParameter(2);
  hpfff->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
  TF1 *fitffn2 = (TF1*)hpfff->GetListOfFunctions()->FindObject("gaus");
  float ffn2_p0= fitffn2->GetParameter(0);
  float ffn2_p1= fitffn2->GetParameter(1);
  float ffn2_p2= fitffn2->GetParameter(2);
  if(doplots) {
    TCanvas* x1 = new TCanvas();
    SCEDraw1(x1,"x1",hpfff,"junkx1.png",0);
    fitffn2->Draw("same");
  }
  std::cout<<std::endl;
  std::cout<<" rel fract fit params "<<ffn2_p0<<" "<<ffn2_p1<<" "<<ffn2_p2<<std::endl;
  std::cout<<std::endl;

  if(doecal) {
    ffm = phcEcalncer->GetBinCenter(phcEcalncer->GetMaximumBin());
    phcEcalncer->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
    TF1 *fitEcalncer1 = (TF1*)phcEcalncer->GetListOfFunctions()->FindObject("gaus");
    ffm = fitEcalncer1->GetParameter(1);
    rms1 = rmsscale*fitEcalncer1->GetParameter(2);
    phcEcalncer->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
    TF1 *fitEcalncer2 = (TF1*)phcEcalncer->GetListOfFunctions()->FindObject("gaus");
    float Ecalncer2_p0= fitEcalncer2->GetParameter(0);
    float Ecalncer2_p1= fitEcalncer2->GetParameter(1);
    float Ecalncer2_p2= fitEcalncer2->GetParameter(2);
    if(doplots) {
      TCanvas* x2 = new TCanvas();
      SCEDraw1(x2,"x2",phcEcalncer,"junkx2.png",0);
      fitEcalncer2->Draw("same");
    }
    std::cout<<std::endl;
    std::cout<<" Ecal ncer fit params "<<Ecalncer2_p0<<" "<<Ecalncer2_p1<<" "<<Ecalncer2_p2<<std::endl;
    std::cout<<std::endl;

    ffm = phcEcalnscint->GetBinCenter(phcEcalnscint->GetMaximumBin());
    phcEcalnscint->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
    TF1 *fitEcalnscint1 = (TF1*)phcEcalnscint->GetListOfFunctions()->FindObject("gaus");
    ffm = fitEcalnscint1->GetParameter(1);
    rms1 = rmsscale*fitEcalnscint1->GetParameter(2);
    phcEcalnscint->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
    TF1 *fitEcalnscint2 = (TF1*)phcEcalnscint->GetListOfFunctions()->FindObject("gaus");
    float Ecalnscint2_p0= fitEcalnscint2->GetParameter(0);
    float Ecalnscint2_p1= fitEcalnscint2->GetParameter(1);
    float Ecalnscint2_p2= fitEcalnscint2->GetParameter(2);
    if(doplots) {
      TCanvas* x3 = new TCanvas();
      SCEDraw1(x3,"x3",phcEcalnscint,"junkx3.png",0);
      fitEcalnscint2->Draw("same");
    }
    std::cout<<std::endl;
    std::cout<<" Ecal nscint fit params "<<Ecalnscint2_p0<<" "<<Ecalnscint2_p1<<" "<<Ecalnscint2_p2<<std::endl;
    std::cout<<std::endl;

    if(dodualcorr) {
      ffm = phcEcalcorr->GetBinCenter(phcEcalcorr->GetMaximumBin());
      phcEcalcorr->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
      TF1 *fitEcalcorr1 = (TF1*)phcEcalcorr->GetListOfFunctions()->FindObject("gaus");
      ffm = fitEcalcorr1->GetParameter(1);
      rms1 = rmsscale*fitEcalcorr1->GetParameter(2);
      phcEcalcorr->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
      TF1 *fitEcalcorr2 = (TF1*)phcEcalcorr->GetListOfFunctions()->FindObject("gaus");
      float Ecalcorr2_p0= fitEcalcorr2->GetParameter(0);
      float Ecalcorr2_p1= fitEcalcorr2->GetParameter(1);
      float Ecalcorr2_p2= fitEcalcorr2->GetParameter(2);
      if(doplots) {
	TCanvas* x4 = new TCanvas();
	SCEDraw1(x4,"x4",phcEcalcorr,"junkx4.png",0);
	fitEcalcorr2->Draw("same");
      }
      std::cout<<std::endl;
      std::cout<<" Ecal corr fit params "<<Ecalcorr2_p0<<" "<<Ecalcorr2_p1<<" "<<Ecalcorr2_p2<<std::endl;
      std::cout<<std::endl;
    }

  }  // end doecal

  if(dohcal) {
    ffm = phcHcalncer->GetBinCenter(phcHcalncer->GetMaximumBin());
    phcHcalncer->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
    TF1 *fitHcalncer1 = (TF1*)phcHcalncer->GetListOfFunctions()->FindObject("gaus");
    ffm = fitHcalncer1->GetParameter(1);
    rms1 = rmsscale*fitHcalncer1->GetParameter(2);
    phcHcalncer->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
    TF1 *fitHcalncer2 = (TF1*)phcHcalncer->GetListOfFunctions()->FindObject("gaus");
    float Hcalncer2_p0= fitHcalncer2->GetParameter(0);
    float Hcalncer2_p1= fitHcalncer2->GetParameter(1);
    float Hcalncer2_p2= fitHcalncer2->GetParameter(2);
    if(doplots) {
      TCanvas* x5 = new TCanvas();
      SCEDraw1(x5,"x5",phcHcalncer,"junkx5.png",0);
      fitHcalncer2->Draw("same");
    }
    std::cout<<std::endl;
    std::cout<<" Hcal ncer fit params "<<Hcalncer2_p0<<" "<<Hcalncer2_p1<<" "<<Hcalncer2_p2<<std::endl;
    std::cout<<std::endl;

    ffm = phcHcalnscint->GetBinCenter(phcHcalnscint->GetMaximumBin());
    phcHcalnscint->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
    TF1 *fitHcalnscint1 = (TF1*)phcHcalnscint->GetListOfFunctions()->FindObject("gaus");
    ffm = fitHcalnscint1->GetParameter(1);
    rms1 = rmsscale*fitHcalnscint1->GetParameter(2);
    phcHcalnscint->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
    TF1 *fitHcalnscint2 = (TF1*)phcHcalnscint->GetListOfFunctions()->FindObject("gaus");
    float Hcalnscint2_p0= fitHcalnscint2->GetParameter(0);
    float Hcalnscint2_p1= fitHcalnscint2->GetParameter(1);
    float Hcalnscint2_p2= fitHcalnscint2->GetParameter(2);
    if(doplots) {
      TCanvas* x6 = new TCanvas();
      SCEDraw1(x6,"x6",phcHcalnscint,"junkx6.png",0);
      fitHcalnscint2->Draw("same");
    }
    std::cout<<std::endl;
    std::cout<<" Hcal nscint fit params "<<Hcalnscint2_p0<<" "<<Hcalnscint2_p1<<" "<<Hcalnscint2_p2<<std::endl;
    std::cout<<std::endl;

    if(dodualcorr) {
      ffm = phcHcalcorr->GetBinCenter(phcHcalcorr->GetMaximumBin());
      phcHcalcorr->Fit("gaus","R0","",ffm-roughfit,ffm+roughfit);
      TF1 *fitHcalcorr1 = (TF1*)phcHcalcorr->GetListOfFunctions()->FindObject("gaus");
      ffm = fitHcalcorr1->GetParameter(1);
      rms1 = rmsscale*fitHcalcorr1->GetParameter(2);
      phcHcalcorr->Fit("gaus","R0","",ffm-rms1,ffm+rms1);
      TF1 *fitHcalcorr2 = (TF1*)phcHcalcorr->GetListOfFunctions()->FindObject("gaus");
      float Hcalcorr2_p0= fitHcalcorr2->GetParameter(0);
      float Hcalcorr2_p1= fitHcalcorr2->GetParameter(1);
      float Hcalcorr2_p2= fitHcalcorr2->GetParameter(2);
      if(doplots) {
	TCanvas* x7 = new TCanvas();
	SCEDraw1(x7,"x7",phcHcalcorr,"junkx7.png",0);
	fitHcalcorr2->Draw("same");
      }
      std::cout<<std::endl;
      std::cout<<" Hcal corr fit params "<<Hcalcorr2_p0<<" "<<Hcalcorr2_p1<<" "<<Hcalcorr2_p2<<std::endl;
      std::cout<<std::endl;
    }  //end dodual cor

  }  //end do hcal


  //***********************************************************************************************************************
  if(doplots) {

  TCanvas* c1 = new TCanvas();
  SCEDraw2(c1,"c1",ehetrue,phetrue,"junk1.png",0);

  TCanvas* c1b = new TCanvas();
  SCEDraw2(c1b,"c1b",ehcEdgeE,phcEdgeE,"junk1b.png",0);
  TCanvas* c1bqq = new TCanvas();
  SCEDraw2(c1bqq,"c1bqq",ehcEdgeRelf,phcEdgeRelf,"junk1bqq.png",0);
  TCanvas* c1b2a = new TCanvas();
  SCEDraw1(c1b2a,"c1b2a",phcEdgeE,"junk1b2a.png",0);
  TCanvas* c1b2 = new TCanvas();
  SCEDraw1(c1b2,"c1b2",phcnonconsE,"junk1b2.png",0);
  //TCanvas* c1b21;
  //SCEDraw1(c1b21,"c1b21",phcnonconsE2,"junk1b21.png",0);
  //TCanvas* c1b22;
  //SCEDraw1(c1b22,"c1b22",phcnonconsE3,"junk1b22.png",0);
  //TCanvas* c1b23;
  //SCEDraw1(c1b23,"c1b23",phcnonconsE4,"junk1b23.png",0);

  TCanvas* c1c = new TCanvas();
  SCEDraw2(c1c,"c1c",ehcEdgeR,phcEdgeR,"junk1c.png",0);


  TCanvas* c1d = new TCanvas();
  SCEDraw1(c1d,"c1d",hedepcal,"junk1d.png",0);
  TCanvas* c1e = new TCanvas();
  SCEDraw1(c1e,"c1e",hpdepcal,"junk1e.png",0);



  if(doecal) {
    TCanvas* ce2 = new TCanvas();
    SCEDraw2(ce2,"ce2",ehcEcalE,phcEcalE,"junke2.png",0);
    TCanvas* ce3 = new TCanvas();
    SCEDraw2(ce3,"ce3",ehcEcalncer,ehcEcalnscint,"junke3.png",0);
    TCanvas* ce4 = new TCanvas();
    SCEDraw3(ce4,"ce4",phcEcalncer,phcEcalnscint,phcEcalcorr,"junke4.png",0);
    TCanvas* ce4dd = new TCanvas();
    SCEDraw3(ce4dd,"ce4dd",phcEandHcalncer,phcEandHcalnscint,phcEandHcalcorr,"junke4dd.png",0);
    TCanvas* ce5 = new TCanvas();
    SCEDraw1_2D(ce5,"ce5",ehcEcalNsNc,"junke5.png",0,0.,0.);
    //TCanvas* ce5b;
    //SCEDraw1_2D(ce5b,"ce5b",ehcEcalMarco,"junke5b.png",0,0.,0.);
    TCanvas* ce6 = new TCanvas();
    SCEDraw1_2D(ce6,"ce6",phcEcalNsNc,"junke6.png",1,0.5,(0.5/kappaEcal)+1-(1/kappaEcal));
    //TCanvas* ce6b;
    //SCEDraw1_2D(ce6b,"ce6b",phcEcalMarco,"junke6b.png",0,0.,0.);
    TCanvas* ce7 = new TCanvas();
    SCEDraw2(ce7,"ce7",eecaltime,piecaltime,"junke7.png",1);

    /*
    TCanvas* ll1;
    SCEDraw1(ll1,"ll1",phcEcalncer,"junkll1.png",0);
    TCanvas* ll2;
    SCEDraw1(ll2,"ll2",phcEcalnscint,"junkll2.png",0);
    TCanvas* ll3;
    SCEDraw1(ll3,"ll3",phcEcalcorr,"junkll3.png",0);
    */


    TCanvas* ecmes1 = new TCanvas();
    SCEDraw1_2D(ecmes1,"ecmes1",mes1Ecal,"junkemes1.png",0,0.,0.);
    TCanvas* ecmes2 = new TCanvas();
    SCEDraw1_2D(ecmes2,"ecmes2",mes2Ecal,"junkemes2.png",0,0.,0.);
    //TCanvas* ecmes3;
    //SCEDraw1_2D(ecmes3,"ecmes3",mes3Ecal,"junkemes3.png",0,0.,0.);
    //TCanvas* ecmes4;
    //SCEDraw1_2D(ecmes4,"ecmes4",mes4Ecal,"junkemes4.png",0,0.,0.);
    TCanvas* ecmes5 = new TCanvas();
    SCEDraw2_2D(ecmes5,"ecmes5",mes3Ecal,mes4Ecal,"junkemes5.png",0,0.,0.);

    //TCanvas* ecms1;
    //SCEDraw1_2D(ecms1,"ecms1",hfnscinEcal,"junkhes1.png",0,0.,0.);
    //TCanvas* ecms2;
    //SCEDraw1_2D(ecms2,"ecms2",hfncerEcal,"junkhes2.png",0,0.,0.);
    TCanvas* ecms3 = new TCanvas();
    SCEDraw2_2D(ecms3,"ecms3",hfncerEcal,hfnscinEcal,"junkhes3.png",0,0.,0.);


    //TCanvas* uu1;
    //SCEDraw1_2D(uu1,"uu1",enonconsEcalcer,"junkuu1.png",0,0.,0.);
    //TCanvas* uu2;
    //SCEDraw1_2D(uu2,"uu2",enonconsEcalscin,"junkuu2.png",0,0.,0.);
    TCanvas* kk2 = new TCanvas();
    SCEDraw2_2D(kk2,"kk2",enonconsEcalcer,enonconsEcalscin,"junkkk2.png",0,0.,0.);

  }



  if(dohcal) {
    TCanvas* ch2 = new TCanvas();
    SCEDraw2(ch2,"ch2",ehcHcalE,phcHcalE,"junkh2.png",0);
    TCanvas* ch2a = new TCanvas();
    SCEDraw2(ch2a,"ch2a",ehcHcalE1,phcHcalE1,"junkh2a.png",0);
    TCanvas* ch2b = new TCanvas();
    SCEDraw2(ch2b,"ch2b",ehcHcalE2,phcHcalE2,"junkh2b.png",0);
    TCanvas* ch3 = new TCanvas();
    SCEDraw2(ch3,"ch3",ehcHcalncer,ehcHcalnscint,"junkh3.png",0);
    TCanvas* ch4 = new TCanvas();
    SCEDraw3(ch4,"ch4",phcHcalncer,phcHcalnscint,phcHcalcorr,"junkh4.png",0);
    TCanvas* ch5 = new TCanvas();
    SCEDraw1_2D(ch5,"ch5",ehcHcalNsNc,"junkh5.png",0,0.,0.);
    TCanvas* ch5c = new TCanvas();
    SCEDraw1_2D(ch5c,"ch5c",ehcHcalNsNctc,"junkh5c.png",0,0.,0.);
    //TCanvas* ch5b;
    //SCEDraw1_2D(ch5b,"ch5b",ehcHcalMarco,"junkhb.png",0,0.,0.);
    TCanvas* ch6 = new TCanvas();
    SCEDraw1_2D(ch6,"ch6",phcHcalNsNc,"junkh6.png",0,0.5,(0.5/kappaHcal)+1-(1/kappaHcal));
    TCanvas* ch6c = new TCanvas();
    SCEDraw1_2D(ch6c,"ch6c",phcHcalNsNctc,"junkh6c.png",0,0.,0.);
    //TCanvas* ch6b;
    //SCEDraw1_2D(ch6b,"ch6b",phcHcalMarco,"junkh6b.png",0,0.,0.);
    TCanvas* ch7 = new TCanvas();
    SCEDraw2_2D(ch7,"ch7",ehcHcalf1f2,phcHcalf1f2,"junkh7.png",0,0.,0.);


    TCanvas* ch8 = new TCanvas();
    SCEDraw2(ch8,"ch8",ehcaltime,pihcaltime,"junkh8.png",1);


    TCanvas* hcmes1 = new TCanvas();
    SCEDraw1_2D(hcmes1,"hcmes1",mes1Hcal,"junkhmes1.png",0,0.,0.);
    TCanvas* hcmes2 = new TCanvas();
    SCEDraw1_2D(hcmes2,"hcmes2",mes2Hcal,"junkhmes2.png",0,0.,0.);
    //TCanvas* hcmes3;
    //SCEDraw1_2D(hcmes3,"hcmes3",mes3Hcal,"junkhmes3.png",0,0.,0.);
    //TCanvas* hcmes4;
    //SCEDraw1_2D(hcmes4,"hcmes4",mes4Hcal,"junkhmes4.png",0,0.,0.);
    TCanvas* hcmes5 = new TCanvas();
    SCEDraw2_2D(hcmes5,"hcmes5",mes3Hcal,mes4Hcal,"junkhmes5.png",0,0.,0.);

    //TCanvas* hcms1;
    //SCEDraw1_2D(hcms1,"hcms1",hfnscinHcal,"junkhms1.png",0,0.,0.);
    //TCanvas* hcms2;
    //SCEDraw1_2D(hcms2,"hcms2",hfncerHcal,"junkhms2.png",0,0.,0.);
    TCanvas* hcms3 = new TCanvas();
    SCEDraw2_2D(hcms3,"hcms3",hfncerHcal,hfnscinHcal,"junkhhes3.png",0,0.,0.);


    //TCanvas* uu3;
    //SCEDraw1_2D(uu3,"uu3",enonconsHcalcer,"junkuu3.png",0,0.,0.);
    //TCanvas* uu4;
    //SCEDraw1_2D(uu4,"uu2",enonconsHcalscin,"junkuu4.png",0,0.,0.);
    TCanvas* kk1 = new TCanvas();
    SCEDraw2_2D(kk1,"kk1",enonconsHcalcer,enonconsHcalscin,"junkkk1.png",0,0.,0.);

    TCanvas* ct4 = new TCanvas();
    SCEDraw3(ct4,"ct4",hpfff,hpfffabs,hpffffib,"junkct4.png",0);
  }


  //TCanvas* ch9a;
  // SCEDraw1_2D(ch9a,"ch9a",enscvni,"junkh9a.png",0,0.,0.);
    //TCanvas* ch9ba;
    //SCEDraw1_2D(ch9ba,"ch9ba",pinalvni,"junkh9bya.png",0,0.,0.);
    //TCanvas* ch9bb;
    //SCEDraw1_2D(ch9bb,"ch9bb",pinscvni,"junkh9byb.png",0,0.,0.);
    //TCanvas* ch9bc;
    //SCEDraw1_2D(ch9bc,"ch9bc",pincevni,"junkh9byc.png",0,0.,0.);
    TCanvas* ch9bd = new TCanvas();
    SCEDraw2_2D(ch9bd,"ch9bd",pincevni,pinscvni,"junkh9bd.png",0,0.,0.);


    TCanvas* nnn1 = new TCanvas();
    SCEDraw1_2D(nnn1,"nnn1",enonconsvnni,"junknnn1.png",0,0.,0.);
    TCanvas* nnn2 = new TCanvas();
    SCEDraw1_2D(nnn2,"nnn2",enonconsvf,"junknnn2.png",0,0.,0.);
    TCanvas* nnn3 = new TCanvas();
    SCEDraw1_2D(nnn3,"nnn3",hcnonconsvesc,"junknnn3.png",0,0.,0.);


  //TCanvas* c7;
  //SCEDrawp(c7,"c7",phcEcalNsNc_pfx,"junk7.png");

    TCanvas* cc1 = new TCanvas();
    SCEDraw2(cc1,"cc1",hefff,hpfff,"junkcc1.png",1);
    TCanvas* cd1 = new TCanvas();
    SCEDraw3(cd1,"cd1",hpfff,hpfff2,hpfff3,"junkcd1.png",1);


    TCanvas* bbv1 = new TCanvas();
    SCEDraw2(bbv1,"bbv1",acovECAL,acovHCAL,"junkbbv1.png",1);


    TCanvas* bbv2 = new TCanvas();
    SCEDraw1(bbv2,"bbv2",ffscinffcer,"junkbbv2.png",0);


  }






  if(dotimingplots) {



    TCanvas* cc2 = new TCanvas();
    SCEDraw2(cc2,"cc2",eecalpd1scint,eecalpd1cer,"junkcc2.png",1);
    TCanvas* cc3 = new TCanvas();
    SCEDraw2(cc3,"cc3",eecalpd2scint,eecalpd2cer,"junkcc3.png",1);
    TCanvas* cc4 = new TCanvas();
    SCEDraw2(cc4,"cc4",ehcalpd1scint,ehcalpd1cer,"junkcc4.png",1);
    TCanvas* cc5 = new TCanvas();
    SCEDraw2(cc5,"cc5",ehcalpd2scint,ehcalpd2cer,"junkcc5.png",1);


    TCanvas* cc6 = new TCanvas();
    SCEDraw2(cc6,"cc6",pecalpd1scint,pecalpd1cer,"junkcc6.png",1);
    TCanvas* cc7 = new TCanvas();
    SCEDraw2(cc7,"cc7",pecalpd2scint,pecalpd2cer,"junkcc7.png",1);
    TCanvas* cc8 = new TCanvas();
    SCEDraw2(cc8,"cc8",phcalpd1scint,phcalpd1cer,"junkcc8.png",1);
    TCanvas* cc9 = new TCanvas();
    SCEDraw2(cc9,"cc9",phcalpd2scint,phcalpd2cer,"junkcc9.png",1);

    TCanvas* cc2z = new TCanvas();
    SCEDraw2(cc2z,"cc2z",eecalpd1scintz,eecalpd1cerz,"junkcc2z.png",1);
    TCanvas* cc3z = new TCanvas();
    SCEDraw2(cc3z,"cc3z",eecalpd2scintz,eecalpd2cerz,"junkcc3z.png",1);
    TCanvas* cc4z = new TCanvas();
    SCEDraw2(cc4z,"cc4z",ehcalpd1scintz,ehcalpd1cerz,"junkcc4z.png",1);
    TCanvas* cc5z = new TCanvas();
    SCEDraw2(cc5z,"cc5z",ehcalpd2scintz,ehcalpd2cerz,"junkcc5z.png",1);


    TCanvas* cc6z = new TCanvas();
    SCEDraw2(cc6z,"cc6z",pecalpd1scintz,pecalpd1cerz,"junkcc6z.png",1);
    TCanvas* cc7z = new TCanvas();
    SCEDraw2(cc7z,"cc7z",pecalpd2scintz,pecalpd2cerz,"junkcc7z.png",1);
    TCanvas* cc8z = new TCanvas();
    SCEDraw2(cc8z,"cc8z",phcalpd1scintz,phcalpd1cerz,"junkcc8z.png",1);
    TCanvas* cc9z = new TCanvas();
    SCEDraw2(cc9z,"cc9z",phcalpd2scintz,phcalpd2cerz,"junkcc9z.png",1);



  }



  //TCanvas* ch6;
  //SCEDraw1_2D(ch6,"ch6",phcHcalNsNc,"junkh6.png",-b1Hcal/m1Hcal,0.);

    std::cout<<"biggest hit time was "<<biggesttime<<std::endl;

  //***********************************************************************************************************

  TFile * out = new TFile(outputfilename,"RECREATE");

  ffscinffcer->Write();

  acovECAL->Write();
  acovHCAL->Write();

  enonconsvnni->Write();
  enonconsvf->Write();
  hcnonconsvesc->Write();


  enonconsEcalcer->Write();
  enonconsEcalscin->Write();
  enonconsHcalcer->Write();
  enonconsHcalscin->Write();

  hfnscinEcal->Write();
  hfncerEcal->Write();
  hfnscinHcal->Write();
  hfncerHcal->Write();

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

  heesumcal->Write();
  heesumemcal->Write();
  hefff->Write();

  hpesumcal->Write();
  hpesumemcal->Write();
  hpfff->Write();
  hpfff2->Write();
  hpfff3->Write();
  hpfffabs->Write();
  hpffffib->Write();

  ehcEcalE->Write();
  phcEcalE->Write();

  ehcHcalE->Write();
  phcHcalE->Write();

  ehcHcalE1->Write();
  phcHcalE1->Write();
  ehcHcalE2->Write();
  phcHcalE2->Write();
  phcHcalvfE1->Write();
  phcHcalvfE2->Write();

  ehcEdgeRelf->Write();
  phcEdgeRelf->Write();
  ehcEdgeE->Write();
  phcEdgeE->Write();
  ehcnonconsE->Write();
  phcnonconsE->Write();
  phcnonconsE2->Write();
  phcnonconsE3->Write();
  phcnonconsE4->Write();



  ehcEdgeR->Write();
  phcEdgeR->Write();



   ehcEcalncer->Write();
   phcEcalncer->Write();
   phcEcalncer2->Write();
   phcEcalncer3->Write();

   ehcHcalncer->Write();
   phcHcalncer->Write();

   ehcEcalcorr->Write();
   phcEcalcorr->Write();
   phcEandHcalcorr->Write();

   ehcHcalcorr->Write();
   phcHcalcorr->Write();


  ehecal2d->Write();
  phecal2d->Write();


   ehcEcalnscint->Write();
   phcEcalnscint->Write();
   phcEcalnscint2->Write();
   phcEcalnscint3->Write();


   ehcHcalnscint->Write();
   phcHcalnscint->Write();


  ehhcal2d->Write();
  phhcal2d->Write();

  ehaphcal->Write();
  phaphcal->Write();

  eheest->Write();
  pheest->Write();


  ehetrue->Write();
  phetrue->Write();

  hedepcal->Write();
  hpdepcal->Write();

  ehnecalcon->Write();
  phnecalcon->Write();

  ehzvst->Write();
  phzvst->Write();

  ehchan->Write();
  phchan->Write();

  ehcEcalNsNc->Write();
  phcEcalNsNc->Write();
  ehcHcalNsNc->Write();
  phcHcalNsNc->Write();

  ehcEcalNsNctc->Write();
  phcEcalNsNctc->Write();
  ehcHcalNsNctc->Write();
  phcHcalNsNctc->Write();


  ehcEcalMarco->Write();
  phcEcalMarco->Write();
  ehcHcalMarco->Write();
  phcHcalMarco->Write();

  ehcHcalf1f2->Write();
  phcHcalf1f2->Write();

  //ehcEcalNsNc_pfx->Write();
  //ehcHcalNsNc_pfx->Write();
  //phcEcalNsNc_pfx->Write();
  //phcHcalNsNc_pfx->Write();

  eecaltime->Write();
  ehcaltime->Write();
  piecaltime->Write();
  pihcaltime->Write();

  enscvni->Write();
  pinalvni->Write();
  pinscvni->Write();
  pincevni->Write();

  CalEcalncer->Write();
  CalEcalnscint->Write();
  CalHcalncer->Write();
  CalHcalnscint->Write();


  out->Close();


}


void SCEDraw1 (TCanvas* canv,  const char* name,TH1F* h1, const char* outfile, bool logy) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
  gStyle->SetOptFit();
  gStyle->SetPalette(1,0);
  if(logy) canv->SetLogy();

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(1);
  h1->SetStats(111111);
  h1->Draw("HIST");



  canv->Print(outfile,".png");
  canv->Update();

  return;
}
void SCEDraw1tp (TCanvas* canv,  const char* name,TProfile* h1, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(kGreen);
  h1->SetStats(111111);
  h1->Draw("HIST");



  canv->Print(outfile,".png");
  canv->Update();

  return;
}

void SCEDrawp (TCanvas* canv,  const char* name,TProfile* h1, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(kGreen);
  h1->SetStats(111111);
  h1->SetMarkerSize(20);
  h1->SetMarkerStyle(4);
  h1->SetMarkerColor(6);
  gStyle->SetOptFit();
  h1->Draw("*");



  canv->Print(outfile,".png");
  canv->Update();

  return;
}


void SCEDraw1_2D (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile, bool doline, float eohS, float eohC) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(kGreen);
  h1->SetMarkerSize(0.2);
  h1->SetMarkerColor(kMagenta);
  h1->SetStats(111111);
  h1->Draw("colz");

  TLine line = TLine(eohS,eohC,1.,1.);
  line.SetLineColor(kBlue);
  line.SetLineWidth(2);
  if(doline) line.Draw("same");


  canv->Print(outfile,".png");
  canv->Update();

  return;
}

void SCEDraw2_2D (TCanvas* canv,  const char* name,TH2F* h1, TH2F* h2, const char* outfile, bool doline, float eohS, float eohC) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(kGreen);
  h1->SetLineWidth(kGreen);
  h1->SetMarkerColor(kGreen);
  h1->SetStats(111111);
  h1->Draw("colz");

  h2->SetLineColor(kBlue);
  h2->SetLineWidth(kBlue);
  h2->SetMarkerColor(kBlue);
  h2->SetStats(111111);
  h2->Draw("same");

  TLine line = TLine(eohS,eohC,1.,1.);
  line.SetLineColor(kCyan);
  line.SetLineWidth(2);
  if(doline) line.Draw("same");


  canv->Print(outfile,".png");
  canv->Update();

  return;
}


void SCEDraw2 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile, bool logy) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
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

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
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



  h3->SetLineColor(kBlue);
  h3->SetLineWidth(3);
  h3->SetStats(111111);
  h3->Draw("HIST same");


  canv->Print(outfile,".png");
  canv->Update();

  return;
}


void DecodeEcal (long long int ihitchan, int& idet, int& ix, int&iy, int& islice, int& ilayer, int&wc, int&type ) {
  idet = (ihitchan) & 0x07;
  ix = (ihitchan >>3) & 0x3F ;  // is this right?
  if(ix>32) ix=ix-64;
  iy =(ihitchan >>10) & 0x3F ; // is this right?
  if(iy>32) iy=iy-64;
  islice = (ihitchan >>17) & 0x07;
  ilayer = (ihitchan>> 20) & 0x07;
  wc=  (ihitchan>> 23) & 0x07;
  type=0;
  if((ilayer==0)&&(islice==1))  type=1;  //pd
  if((ilayer==0)&&(islice==2))  type=4;  //resin
  if((ilayer==0)&&(islice==3))  type=4;  //cookie
  if((ilayer==0)&&(islice==4))  type=2;  //crystal
  if((ilayer==0)&&(islice==5))  type=3;  //air




  if((ilayer==1)&&(islice==1))  type=2;  //crystal
  if((ilayer==1)&&(islice==2))  type=4;  //cookie
  if((ilayer==1)&&(islice==3))  type=4;  //resin
  if((ilayer==1)&&(islice==4))  type=1;  //pd
  if((ilayer==1)&&(islice==5))  type=3;  //air

  return;
}


void DecodeFiber (long long int ihitchan, int& idet, int& ilayer, int& itube, int& iair, int&itype, int& ifiber, int& iabs, int& iphdet, int& ihole, int& ix, int& iy) {
  idet = (ihitchan) & 0xFF;
  ilayer = (ihitchan >>8) & 0xFFF;
  itube = (ihitchan >>20) & 0xFFF;
  //int iair=0; int itype=0;
  iair = (ihitchan >>32) & 0x7;
  itype = (ihitchan >>35) & 0x7;
  ifiber=0; iabs=0; iphdet=0;  ihole=0;
  ix=0; iy=0;
  if((itype==0)&&(iair==0)&&(itube!=0)) iabs=1;
  if(itype==1) ifiber=1; // scint
  if(itype==2) ifiber=2; // quartz
  if(itype==3) iphdet=1; //scint pt
  if(itype==4) iphdet=2; // quartz pt
  if(((iair==1)||(iair==2))&&(itype==0)) ihole=1;
  if(itube==0) ihole=1;
  ix=itube;
  iy=ilayer;
  return;
}


void DecodeSampling(long long int ihitchan,int& idet, int& ix, int& iy, int& ilayer, int& ibox2, int& islice) {
  idet = (ihitchan) & 0x07;
  iy = (ihitchan >>3) & 0xFFF;
  ix = (ihitchan >>15) & 0xFFF;
  ilayer = (ihitchan >>27) & 0xFFF;
  ibox2 = (ihitchan >> 39) & 0x03;
  islice = (ihitchan >>41) & 0xF;
  return;
}

void CalibRefine(map<string, int> mapsampcalslice,  int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits,
	      float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,
		 TH1F *CalEcalncer, TH1F *CalEcalnscint, TH1F *CalHcalncer, TH1F *CalHcalnscint
		 ){

  int nbyteecal, nbytehcal, nbyteedge;
  float ameanscinEcal(0.);
  float ameancerEcal(0.);
  float ameanscinHcal(0.);
  float ameancerHcal(0.);


  if(doecal) {
    if(ievt<SCECOUNT) std::cout<<"CalibRefine phot ievt is "<<ievt<<std::endl;

    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      long long int ihitchan=aecalhit->cellID;
      int idet,ix,iy,islice,ilayer,wc,type;
      DecodeEcal(ihitchan,idet,ix,iy,islice,ilayer,wc,type );

      if(gendet==1) {   // use photons as generated in otical material
	if(type==2 ) {  // crystal
	  ameancerEcal+=aecalhit->ncerenkov;
	  ameanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( type==1 ) { // either photo detector
	  ameancerEcal+=aecalhit->ncerenkov;
	  ameanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==3||gendet==4){
	if(idet==5) {
	  if(type==2 ) {  // crystal
	    ameanscinEcal+=aecalhit->energyDeposit;
	    if(gendet==3) ameancerEcal+=aecalhit->edeprelativistic;
	    if(gendet==4) ameancerEcal+=aecalhit->energyDeposit;
	  }
	}
      }
    }
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);

      // hcal hits
    if(ievt<SCECOUNT) std::cout<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet ix iy ifiber iabs iphdet "<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet iy ix ilayer islice  "<<std::endl;

    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);

      long long int ihitchan=ahcalhit->cellID;

      if(hcaltype==0) { // fiber

	int idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy;
	DecodeFiber(ihitchan,idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy);

	//if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<ilayer<<" "<<itube<<" "<<iair<<" "<<itype<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>8)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>20)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>32)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>35)<<std::endl;


	if(gendet==1) {  // take light as generated in fiber
	  if(ifiber==1) {  // scintillating fibers
	    ameanscinHcal+=ahcalhit->nscintillator;
	  }
	  if(ifiber==2) {  // quartz fibers
	    ameancerHcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==2) {
	  if(iphdet==1) {  // take light that hits photodetectors
	    ameanscinHcal+=ahcalhit->nscintillator;
	  }
	  if(iphdet==2) {  // take light that hits photodetectors
	    ameancerHcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==3||gendet==4) {
	  if(idet==6) {
	    if(ifiber==1) {
	      ameanscinHcal+=ahcalhit->energyDeposit;
	    }
	    if(ifiber==2) {
	      if(gendet==3) ameancerHcal+=ahcalhit->edeprelativistic;
	      if(gendet==4) ameancerHcal+=ahcalhit->energyDeposit;
	    }
	  }
	}
      }
      else {  // sampling
	int idet,ix,iy,ilayer,ibox2,islice;
	DecodeSampling(ihitchan,idet,ix,iy,ilayer,ibox2,islice);

	if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;

	if(gendet==1) {  // take light as generated in media
	  if(islice==(*mapsampcalslice.find("PS")).second) {
	    ameanscinHcal+=ahcalhit->nscintillator;
	    //	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*mapsampcalslice.find("Quartz")).second) {  // cherenkov
	    ameancerHcal+=ahcalhit->ncerenkov;
	    //	    std::cout<<"add ceren"<<std::endl;
	  }
	}
	else if(gendet==2) {
	  if( (islice==(*mapsampcalslice.find("PD1")).second)||(islice==(*mapsampcalslice.find("PD2")).second) ) { // either photo detector
	    ameanscinHcal+=ahcalhit->nscintillator;
	  }
	  if( (islice==(*mapsampcalslice.find("PD3")).second)||(islice==(*mapsampcalslice.find("PD4")).second)) {  // take light that hits photodetectors
	    ameancerHcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==3||gendet==4) {
	  if(idet==6) {
	    if( islice==(*mapsampcalslice.find("PS")).second) { // PS
	      ameanscinHcal+=ahcalhit->energyDeposit;
	    }
	    if( islice==(*mapsampcalslice.find("Quartz")).second ) {  // quartz
	      if(gendet==3) ameancerHcal+=ahcalhit->edeprelativistic;
	      if(gendet==4) ameancerHcal+=ahcalhit->energyDeposit;
	    }
	  }
	}


      }


    }  // end loop over hcal hits
  }

  if(meancerEcal>0) CalEcalncer->Fill(ameancerEcal/meancerEcal);
  if(meanscinEcal>0) CalEcalnscint->Fill(ameanscinEcal/meanscinEcal);
  if(meancerHcal>0) CalHcalncer->Fill(ameancerHcal/meancerHcal);
  if(meanscinHcal>0) CalHcalnscint->Fill(ameanscinHcal/meanscinHcal);

}






void getMeanPhot(map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits,
		 float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,float &timecut, float &eecaltimecut, float &ehcaltimecut, float &erelecaltimecut, float &erelhcaltimecut){
  int nbyteecal, nbytehcal, nbyteedge;



  if(doecal) {
    if(ievt<SCECOUNT) std::cout<<"getMean phot ievt is "<<ievt<<std::endl;

    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      long long int ihitchan=aecalhit->cellID;
      int idet,ix,iy,islice,ilayer,wc,type;
      DecodeEcal(ihitchan,idet,ix,iy,islice,ilayer,wc,type );


      Contributions zxzz=aecalhit->truth;

      if(gendet==1) {   // use photons as generated in otical material
	if(type==2 ) {  // crystal
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( type==1 ) { // either photo detector
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==3||gendet==4){
	if(idet==5) {
	  if(type==2 ) {  // crystal
	    meanscinEcal+=aecalhit->energyDeposit;
	    if(gendet==3) meancerEcal+=aecalhit->edeprelativistic;
	    if(gendet==4) meancerEcal+=aecalhit->energyDeposit;

	    for(size_t j=0;j<zxzz.size(); j++) {
	      if( (zxzz.at(j)).time>biggesttime) biggesttime=(zxzz.at(j)).time;
	      if((zxzz.at(j)).time<timecut) {
		eecaltimecut+=(zxzz.at(j)).deposit;
		if(((aecalhit->contribBeta)[j])>betacut) erelecaltimecut+=(zxzz.at(j)).deposit;
	      }
	    }
	  }
	}
      }
    }
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);

      // hcal hits
    if(ievt<SCECOUNT) std::cout<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet ix iy ifiber iabs iphdet "<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet iy ix ilayer islice  "<<std::endl;

    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);

      long long int ihitchan=ahcalhit->cellID;

      Contributions zxzz=ahcalhit->truth;

      if(hcaltype==0) { // fiber

	int idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy;
	DecodeFiber(ihitchan,idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy);


	//if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<ilayer<<" "<<itube<<" "<<iair<<" "<<itype<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>8)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>20)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>32)<<std::endl;
	//std::cout<<std::hex<<(ihitchan>>35)<<std::endl;



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
	    if((ifiber==1)||(ifiber==2)) {
	      for(size_t j=0;j<zxzz.size(); j++) {
		if( (zxzz.at(j)).time>biggesttime) biggesttime=(zxzz.at(j)).time;
		if((zxzz.at(j)).time<timecut) {
		  if(ifiber==1) {
		  ehcaltimecut+=(zxzz.at(j)).deposit;
		  }
		  if(ifiber==2) {
		  if(((ahcalhit->contribBeta)[j])>betacut) erelhcaltimecut+=(zxzz.at(j)).deposit;
		  }
		}
	      }
	    }
	  }
	}
      }
      else {  // sampling

	int idet,ix,iy,ilayer,ibox2,islice;
	DecodeSampling(ihitchan,idet,ix,iy,ilayer,ibox2,islice);


	//	if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;

	if(gendet==1) {  // take light as generated in media
	  if(islice==(*mapsampcalslice.find("PS")).second) {
	    meanscinHcal+=ahcalhit->nscintillator;
	    //	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*mapsampcalslice.find("Quartz")).second) {  // cherenkov
	    meancerHcal+=ahcalhit->ncerenkov;
	    //	    std::cout<<"add ceren"<<std::endl;
	  }
	}
	else if(gendet==2) {
	  if( (islice==(*mapsampcalslice.find("PD1")).second)||(islice==(*mapsampcalslice.find("PD2")).second) ) { // either photo detector
	    meanscinHcal+=ahcalhit->nscintillator;
	  }
	  if( (islice==(*mapsampcalslice.find("PD3")).second)||(islice==(*mapsampcalslice.find("PD4")).second)) {  // take light that hits photodetectors
	    meancerHcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==3||gendet==4) {
	  if(idet==6) {
	    if( islice==(*mapsampcalslice.find("PS")).second) { // PS
	      meanscinHcal+=ahcalhit->energyDeposit;
	      if(ievt<SCECOUNT) std::cout<<" meanscinHcal "<<meanscinHcal<<std::endl;
	    }
	    if( islice==(*mapsampcalslice.find("Quartz")).second ) {  // quartz
	      if(gendet==3) meancerHcal+=ahcalhit->edeprelativistic;
	      if(gendet==4) meancerHcal+=ahcalhit->energyDeposit;
	      if(ievt<SCECOUNT) std::cout<<" meancerHcal "<<meancerHcal<<std::endl;
	    }
	    if(( islice==(*mapsampcalslice.find("PS")).second)||( islice==(*mapsampcalslice.find("Quartz")).second)) {
	      for(size_t j=0;j<zxzz.size(); j++) {
		if( (zxzz.at(j)).time>biggesttime) biggesttime=(zxzz.at(j)).time;
		if((zxzz.at(j)).time<timecut) {
		  ehcaltimecut+=(zxzz.at(j)).deposit;
		  if(((ahcalhit->contribBeta)[j])>betacut) erelhcaltimecut+=(zxzz.at(j)).deposit;
		}
	      }
	    }
	  }
	}


      }


    }  // end loop over hcal hits
  }


}





void getStuff(map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, float &timecut,bool &fillhists,
	      float  &eesum,float &eesumcal,float &eesumem, float &eesumair,float &eesumdead, float &eesumcrystal,float &eesumPDe,float &eesumfiber1,float &eesumfiber2,float &eesumabs,float &eesumPDh,
	      float &eesumairem, float &eesumdeadem, float &eesumcrystalem,float &eesumPDeem,float &eesumfiber1em, float &eesumfiber2em,float &eesumabsem,float &eesumPDhem,
	      float &eesumedge,float &eesumedgerel, float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal,
	      float &eecaltimecut, float &ehcaltimecut,float &erelecaltimecut, float &erelhcaltimecut,
	      int &nine,int &ninh
	      ){


  if(ievt<SCECOUNT) std::cout<<"getstuff phot ievt is "<<ievt<<std::endl;

  int nbyteecal, nbytehcal, nbyteedge;

  /*
  std::cout<<"entering getStuff"<<std::endl;
      std::cout<<"eesum now "<<eesum<<std::endl;
      std::cout<<"ehcaltimecut is "<<ehcaltimecut<<std::endl;
      std::cout<<eesumfiber1<<" "<<eesumfiber2<<" "<<eesumabs<<std::endl;
  */

  if(doecal) {
    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    eecaltimecut=0.;
    erelecaltimecut=0.;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);

      //      if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
      long long int ihitchan=aecalhit->cellID;
      int idet,ix,iy,islice,ilayer,wc,type;
      DecodeEcal(ihitchan,idet,ix,iy,islice,ilayer,wc,type );

      if((ilayer!=0)&&(ilayer!=1)) std::cout<<"danger danger will robinson ilayer not zero"<<std::endl;


      float ae=aecalhit->energyDeposit;
      nine+=aecalhit->n_inelastic;


      eesum+=ae;
      if(type==3) {eesumair+=ae;eesumcal+=aecalhit->energyDeposit;eesumem+=aecalhit->edeprelativistic;eesumairem+=aecalhit->edeprelativistic;};
      if(type==1) {eesumPDe+=ae;eesumcal+=aecalhit->energyDeposit;eesumem+=aecalhit->edeprelativistic;eesumPDeem+=aecalhit->edeprelativistic;};
      if(type==2) {eesumcrystal+=ae;eesumcal+=aecalhit->energyDeposit;eesumem+=aecalhit->edeprelativistic;eesumcrystalem+=aecalhit->edeprelativistic;};
      if(type==0||type==4) {eesumdead+=ae;eesumcal+=aecalhit->energyDeposit;eesumem+=aecalhit->edeprelativistic;eesumdeadem+=aecalhit->edeprelativistic;};





      if(gendet==1) {   // use photons as generated in otical material
	if(type==2) {  // crystal
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( type==1 ) { // either photo detector
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==3||gendet==4){
	if(idet==5) {
	  if( type==2 ) {  // crystal
	    nescinttotecal+=aecalhit->energyDeposit;
	    if(gendet==3) necertotecal+=aecalhit->edeprelativistic;
	    if(gendet==4) necertotecal+=aecalhit->energyDeposit;
	    Contributions zxzz=aecalhit->truth;
	    float hacheck=0.;
	    for(size_t j=0;j<zxzz.size(); j++) {
	      hacheck+=(zxzz.at(j)).deposit;
	      if((zxzz.at(j)).time<timecut) {
		eecaltimecut+=(zxzz.at(j)).deposit;
		if(((aecalhit->contribBeta)[j])>betacut) erelecaltimecut+=(zxzz.at(j)).deposit;
	      }
	      //if(fillhists) eecaltime->Fill((zxzz.at(j)).time);
	    }
	    if(ae>0.001) {
	      if(hacheck/ae<0.999) {

		if(icountyuck<SCECOUNT3) std::cout<<"missing contribs: ecal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ae is "<<ae<<" ratio "<<hacheck/ae<<std::endl;
				icountyuck+=1;
	      }
	    }
	  }
	}
      }

	//ehchan->Fill(aecalhit->cellID);
	//ehecal2d->Fill(ix,iy,aecalhit->energyDeposit);


    }  // end loop over ecal hits
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);


      // hcal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    ehcaltimecut=0.;
    erelhcaltimecut=0.;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      float ah=ahcalhit->energyDeposit;
      ninh+=ahcalhit->n_inelastic;

      float aarel = ahcalhit->edeprelativistic;
      eesum+=ah;
      eesumcal+=ahcalhit->energyDeposit;eesumem+=ahcalhit->edeprelativistic;
      /*
      std::cout<<"eesum now "<<eesum<<std::endl;
      std::cout<<"ehcaltimecut is "<<ehcaltimecut<<std::endl;
      std::cout<<eesumfiber1<<" "<<eesumfiber2<<" "<<eesumabs<<std::endl;
      */

      long long int ihitchan=ahcalhit->cellID;
      if(hcaltype==0) { // fiber

	int idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy;
	DecodeFiber(ihitchan,idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy);


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
	    if(ifiber==1) {
	      nescinttothcal+=ahcalhit->energyDeposit;
	    }
	    if(ifiber==2) {
	      if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
	      if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
	    }
	    if((ifiber==1)||(ifiber==2) ) {
      // check contribs
	      Contributions zxzz=ahcalhit->truth;
	      float hacheck=0.;
	      for(size_t j=0;j<zxzz.size(); j++) {
		hacheck+=(zxzz.at(j)).deposit;
		if((zxzz.at(j)).time<timecut) {
		  if(ifiber==1) {
		  ehcaltimecut+=(zxzz.at(j)).deposit;
		  }
		  if(ifiber==2) {
		  if(((ahcalhit->contribBeta)[j])>betacut) erelhcaltimecut+=(zxzz.at(j)).deposit;
		  }
		}
		//if(fillhists) ehcaltime->Fill((zxzz.at(j)).time);
	      }
	      if(ah>0.001) {
		if(hacheck/ah<0.99999) std::cout<<"missing contribs: hcal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ah is "<<ah<<" ratio "<<hacheck/ah<<std::endl;
	      }
	    }
	  }
	}

	if(ifiber==1) {eesumfiber1+=ah;eesumfiber1em+=aarel;}
	if(ifiber==2) {eesumfiber2+=ah;eesumfiber2em+=aarel;}
	if(iabs==1) {eesumabs+=ah;eesumabsem+=aarel;}
	if(iphdet>1) {eesumPDh+=ah;eesumPDhem+=aarel;}
	//std::cout<<"   "<<ihitchan<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<" "<<eesumfiber<<" "<<eesumabs<<" "<<eesumPDh<<std::endl;


      }
      else {  // sampling

	int idet,ix,iy,ilayer,ibox2,islice;
	DecodeSampling(ihitchan,idet,ix,iy,ilayer,ibox2,islice);



	//std::cout<<" idet iy ix ilayer islice are "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;

	//std::cout<<"energy nscint ncer is "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<std::endl;




	if(gendet==1) {  // take light as generated in media
	  if(islice==(*mapsampcalslice.find("PS")).second) {
	    nescinttothcal+=ahcalhit->nscintillator;
	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*mapsampcalslice.find("Quartz")).second) {  // chereknov
	    necertothcal+=ahcalhit->ncerenkov;
	    std::cout<<"add ceren"<<std::endl;
	  }
	}
	else if(gendet==2) {
	  if( (islice==(*mapsampcalslice.find("PD1")).second)||(islice==(*mapsampcalslice.find("PD2")).second) ) { // either photo detector
	    nescinttothcal+=ahcalhit->nscintillator;
	  }
	  if( (islice==(*mapsampcalslice.find("PD3")).second)||(islice==(*mapsampcalslice.find("PD4")).second)) {  // take light that hits photodetectors
	    necertothcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==3||gendet==4) {
	  //std::cout<<" here "<<std::endl;
	  if(idet==6) {
	  if( islice==(*mapsampcalslice.find("PS")).second) { //  ps
	    nescinttothcal+=ahcalhit->energyDeposit;
	    //	    if(i<10) std::cout<<" i nescinttothcal "<<i<<" "<<nescinttothcal<<std::endl;
	  }
	  if( islice==(*mapsampcalslice.find("Quartz")).second ) {  // quartz
	    if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
	    if(gendet==4) necertothcal+=ahcalhit->energyDeposit;

	  }
	  if((islice==(*mapsampcalslice.find("PS")).second)||(islice==(*mapsampcalslice.find("Quartz")).second) ){
      // check contribs
	      Contributions zxzz=ahcalhit->truth;
	      float hacheck=0.;
	      for(size_t j=0;j<zxzz.size(); j++) {
		hacheck+=(zxzz.at(j)).deposit;
		if((zxzz.at(j)).time<timecut) {
		  ehcaltimecut+=(zxzz.at(j)).deposit;
		  if(((ahcalhit->contribBeta)[j])>betacut) erelhcaltimecut+=(zxzz.at(j)).deposit;
		}
		//if(fillhists) ehcaltime->Fill((zxzz.at(j)).time);
	      }
	      if(ah>0.001) {
		if(hacheck/ah<0.99999) std::cout<<"missing contribs: hcal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ah is "<<ah<<" ratio "<<hacheck/ah<<std::endl;
	      }

	  }



	  //	  for(size_t j=0;j<zxzz.size(); j++) {
	  //  if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  //}
	  }
	}



	if( islice==(*mapsampcalslice.find("PS")).second ) {eesumfiber1+=ah;eesumfiber1em+=aarel;}; // scint
	if( islice==(*mapsampcalslice.find("Quartz")).second ) {eesumfiber2+=ah;eesumfiber2em+=aarel;};  //cer
	if( (islice==(*mapsampcalslice.find("Iron")).second)||(islice==(*mapsampcalslice.find("Sep1")).second)||(islice==(*mapsampcalslice.find("Sep2")).second) ) {eesumabs+=ah;eesumabsem+=aarel;};
	if(  (islice==(*mapsampcalslice.find("PD1")).second) || (islice==(*mapsampcalslice.find("PD2")).second) ||  (islice==(*mapsampcalslice.find("PD3")).second) || (islice==(*mapsampcalslice.find("PD4")).second)) {eesumPDh+=ah;eesumPDhem+=aarel;};


      }
	//	if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<islice<<" "<<ilayer<<std::endl;


	//ehchan->Fill(ahcalhit->cellID);
	//ehhcal2d->Fill(ix,iy,ahcalhit->energyDeposit);

    }  // end loop over hcal hits
  }

  if(doedge) {
    nbyteedge = b_edge->GetEntry(ievt);


    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of edge hits is "<<edgehits->size()<<std::endl;
    for(size_t i=0;i<edgehits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aedgehit =edgehits->at(i);

      float ae=aedgehit->energyDeposit;
      eesum+=ae;
      eesumedge+=ae;
      eesumedgerel=aedgehit->edepepgam;

    }  // end loop over escaping hits
  }

  //std::cout<<" yuck 2 ehcaltimecut is "<<ehcaltimecut<<std::endl;


}


void FillTime(map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, float &timecut,TH1F* eecaltime, TH1F* ehcaltime,TH1F *ecalpd1scint,TH1F *ecalpd1cer,TH1F *ecalpd2scint,TH1F *ecalpd2cer,TH1F *hcalpd1scint,TH1F *hcalpd1cer,TH1F *hcalpd2scint,TH1F *hcalpd2cer,TH1F *ecalpd1scintz,TH1F *ecalpd1cerz,TH1F *ecalpd2scintz,TH1F *ecalpd2cerz,TH1F *hcalpd1scintz,TH1F *hcalpd1cerz,TH1F *hcalpd2scintz,TH1F *hcalpd2cerz){

  if(ievt<SCECOUNT) std::cout<<"getstuff phot ievt is "<<ievt<<std::endl;
  int nbyteecal, nbytehcal, nbyteedge;


  if(doecal) {
    nbyteecal = b_ecal->GetEntry(ievt);
      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      long long int ihitchan=aecalhit->cellID;
      int idet,ix,iy,islice,ilayer,wc,type;
      DecodeEcal(ihitchan,idet,ix,iy,islice,ilayer,wc,type );
      float ae=aecalhit->energyDeposit;
      if((ilayer==0)&&(islice==1)) {  // pd on entrance to ecal
	for(int ijk=0;ijk<finenbin;ijk++){
	  ecalpd1scint->Fill((ijk+0.5)*timebinsize,aecalhit->nscinttime[ijk]);
	  ecalpd1cer->Fill((ijk+0.5)*timebinsize,aecalhit->ncertime[ijk]);
	  ecalpd1scintz->Fill((ijk+0.5)*timebinsizez,aecalhit->nscinttimez[ijk]);
	  ecalpd1cerz->Fill((ijk+0.5)*timebinsizez,aecalhit->ncertimez[ijk]);
	}
      }
      if((ilayer==1)&&(islice==2)) {  // pd on exist of ecal
	for(int ijk=0;ijk<finenbin;ijk++){
	  ecalpd2scint->Fill((ijk+0.5)*timebinsize,aecalhit->nscinttime[ijk]);
	  ecalpd2cer->Fill((ijk+0.5)*timebinsize,aecalhit->ncertime[ijk]);
	  ecalpd2scintz->Fill((ijk+0.5)*timebinsizez,aecalhit->nscinttimez[ijk]);
	  ecalpd2cerz->Fill((ijk+0.5)*timebinsizez,aecalhit->ncertimez[ijk]);
	}
      }
      if(gendet==3||gendet==4){
	if(idet==5) {
	  if( type==2 ) {  // crystal
	    Contributions zxzz=aecalhit->truth;
	    for(size_t j=0;j<zxzz.size(); j++) {
	      eecaltime->Fill((zxzz.at(j)).time);
	    }
	  }
	}
      }


    }  // end loop over ecal hits
  }  //end doecal

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);
      // hcal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      float ah=ahcalhit->energyDeposit;
      float aarel = ahcalhit->edeprelativistic;
      long long int ihitchan=ahcalhit->cellID;
      if(hcaltype==0) { // fiber
	int idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy;
	DecodeFiber(ihitchan,idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy);
	if(iphdet==1) {  // pd on scintillating fibers
	  for(int ijk=0;ijk<finenbin;ijk++){
	    hcalpd1scint->Fill((ijk+0.5)*timebinsize,ahcalhit->nscinttime[ijk]);
	    hcalpd1cer->Fill((ijk+0.5)*timebinsize,ahcalhit->ncertime[ijk]);
	    hcalpd1scintz->Fill((ijk+0.5)*timebinsizez,ahcalhit->nscinttimez[ijk]);
	    hcalpd1cerz->Fill((ijk+0.5)*timebinsizez,ahcalhit->ncertimez[ijk]);
	  }
	}
	if(iphdet==2) {  // pd on cherenkov fibers
	  for(int ijk=0;ijk<finenbin;ijk++){
	    hcalpd2scint->Fill((ijk+0.5)*timebinsize,ahcalhit->nscinttime[ijk]);
	    hcalpd2cer->Fill((ijk+0.5)*timebinsize,ahcalhit->ncertime[ijk]);
	    hcalpd2scintz->Fill((ijk+0.5)*timebinsizez,ahcalhit->nscinttimez[ijk]);
	    hcalpd2cerz->Fill((ijk+0.5)*timebinsizez,ahcalhit->ncertimez[ijk]);
	  }
	}
	if(gendet==3||gendet==4) {
	  if(idet==6) {
	    if((ifiber==1)||(ifiber==2) ) {
      // check contribs
	      Contributions zxzz=ahcalhit->truth;
	      for(size_t j=0;j<zxzz.size(); j++) {
		ehcaltime->Fill((zxzz.at(j)).time);
	      }
	    }
	  }
	}

      }  //end hcaltype==0
      else {  // sampling
	int idet,ix,iy,ilayer,ibox2,islice;
	DecodeSampling(ihitchan,idet,ix,iy,ilayer,ibox2,islice);
	if( (islice==(*mapsampcalslice.find("PD1")).second)||(islice==(*mapsampcalslice.find("PD2")).second)) {  // pd on scint?
	  for(int ijk=0;ijk<finenbin;ijk++){
	    hcalpd1scint->Fill((ijk+0.5)*timebinsize,ahcalhit->nscinttime[ijk]);
	    hcalpd1cer->Fill((ijk+0.5)*timebinsize,ahcalhit->ncertime[ijk]);
	    hcalpd1scintz->Fill((ijk+0.5)*timebinsizez,ahcalhit->nscinttimez[ijk]);
	    hcalpd1cerz->Fill((ijk+0.5)*timebinsizez,ahcalhit->ncertimez[ijk]);
	  }
	}
	if( (islice==(*mapsampcalslice.find("PD3")).second)||(islice==(*mapsampcalslice.find("PD4")).second)) {  // pd on quartz?
	  for(int ijk=0;ijk<finenbin;ijk++){
	    hcalpd2scint->Fill((ijk+0.5)*timebinsize,ahcalhit->nscinttime[ijk]);
	    hcalpd2cer->Fill((ijk+0.5)*timebinsize,ahcalhit->ncertime[ijk]);
	    hcalpd2scintz->Fill((ijk+0.5)*timebinsizez,ahcalhit->nscinttimez[ijk]);
	    hcalpd2cerz->Fill((ijk+0.5)*timebinsizez,ahcalhit->ncertimez[ijk]);
	  }
	}
	if(gendet==3||gendet==4) {
	  if((islice==(*mapsampcalslice.find("PS")).second)||(islice==(*mapsampcalslice.find("Quartz")).second) ){
      // check contribs
	    Contributions zxzz=ahcalhit->truth;
	    for(size_t j=0;j<zxzz.size(); j++) {
	      ehcaltime->Fill((zxzz.at(j)).time);
	    }
	  }
	}
      }  // end sampling

    }  // end loop over hcal hits
  }  //end dohcal

  return;

}


void getStuffDualCorr(bool domissCorr, float beamE, map<string, int> mapsampcalslice, int gendet, float kappaEcal, float kappaHcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype,
		      bool doedge, float &eesumedge, float&eesumedgerel, TBranch* &b_ecal,TBranch* &b_hcal, TBranch* &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
		      float &EEcal, float &EHcal, 	      float &timecut, float &eecaltimecut, float &ehcaltimecut,
		      float &erelecaltimecut, float &erelhcaltimecut)
{
  float necertotecal(0),nescinttotecal(0),necertothcal(0),nescinttothcal(0);
  int nbyteecal, nbytehcal, nbyteedge;


  if(ievt<SCECOUNT) std::cout<<"getstuffdualcorr phot ievt is "<<ievt<<std::endl;


  if(doecal) {
    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    eecaltimecut=0.;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      long long int ihitchan=aecalhit->cellID;
      int idet,ix,iy,islice,ilayer,wc,type;
      DecodeEcal(ihitchan,idet,ix,iy,islice,ilayer,wc,type );



      Contributions zxzz=aecalhit->truth;






      if(gendet==1) {   // use photons as generated in otical material
	if(type==2) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( type==1 ) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==3||gendet==4){
	if(idet==5) {
	  if(type==2) {  // crystal
	    nescinttotecal+=aecalhit->energyDeposit;
	    if(gendet==3) necertotecal+=aecalhit->edeprelativistic;
	    if(gendet==4) necertotecal+=aecalhit->energyDeposit;

	    for(size_t j=0;j<zxzz.size(); j++) {
	      if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
	    }
	  }
	}
      }

    }  // end loop over ecal hits
  }  // end do ecal

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);

      // hcal hits
    if(ievt<SCECOUNT) std::cout<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet ix iy ifiber iabs iphdet "<<std::endl;
    ehcaltimecut=0.;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      Contributions zxzz=ahcalhit->truth;
      long long int ihitchan=ahcalhit->cellID;
      if(hcaltype==0) { // fiber

	int idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy;
	DecodeFiber(ihitchan,idet,ilayer,itube,iair,itype,ifiber,iabs,iphdet,ihole,ix,iy);

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
	  if(ifiber==1) {
	    nescinttothcal+=ahcalhit->energyDeposit;
	  }
	  if(ifiber==2) {
	    if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
	    if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
	  }
	  for(size_t j=0;j<zxzz.size(); j++) {
	    if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  }
	  }
	}



      }
      else {  // sampling

	int idet,ix,iy,ilayer,ibox2,islice;
	DecodeSampling(ihitchan,idet,ix,iy,ilayer,ibox2,islice);


	//	if(ievt<SCECOUNT) std::cout<<"   "<<std::hex<<ihitchan<<std::dec<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<islice<<" "<<ilayer<<" "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<std::endl;



	if(gendet==1) {  // take light as generated in media
	  if(islice==(*mapsampcalslice.find("PS")).second) {
	    nescinttothcal+=ahcalhit->nscintillator;
	    //	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*mapsampcalslice.find("Quartz")).second) {  // cherenkov
	    necertothcal+=ahcalhit->ncerenkov;
	    //	    std::cout<<"add ceren"<<std::endl;
	  }
	}
	else if(gendet==2) {
	  if( (islice==(*mapsampcalslice.find("PD1")).second)||(islice==(*mapsampcalslice.find("PD2")).second) ) { // either photo detector
	    nescinttothcal+=ahcalhit->nscintillator;
	  }
	  if( (islice==(*mapsampcalslice.find("PD3")).second)||(islice==(*mapsampcalslice.find("PD4")).second)) {  // take light that hits photodetectors
	    necertothcal+=ahcalhit->ncerenkov;
	  }
	}
	else if(gendet==3||gendet==4) {
	  if(idet==6) {
	  if( islice==(*mapsampcalslice.find("PS")).second ) { // ps
	    nescinttothcal+=ahcalhit->energyDeposit;
	    //if(i<10) std::cout<<" i nescinttothcal "<<i<<" "<<nescinttothcal<<std::endl;
	  }
	  if( islice==(*mapsampcalslice.find("Quartz")).second ) {  // quartz
	    if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
	    if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
	    //if(i<10) std::cout<<" i necertothcal "<<i<<" "<<necertothcal<<std::endl;
	  }
	  for(size_t j=0;j<zxzz.size(); j++) {
	    if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  }
	  }
	}

      }

    }  // end hit loop
  }  // end do hcal



  float ContainedFrac(1.);
  if(doedge) {
    nbyteedge = b_edge->GetEntry(ievt);


    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of edge hits is "<<edgehits->size()<<std::endl;
    for(size_t i=0;i<edgehits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aedgehit =edgehits->at(i);

      float ae=aedgehit->energyDeposit;
      eesumedge+=ae;
      eesumedgerel=aedgehit->edeprelativistic;

    }  //f end loop over escaping hits

    float ContainedFrac=(beamE-eesumedge)/beamE;
  }




  //  std::cout<<"  getstuffdual ecal cer scint count "<<necertotecal<<" "<<nescinttotecal<<std::endl;
  float anecertotecal=necertotecal/meancerEcal/ContainedFrac;
  float anescinttotecal=nescinttotecal/meanscinEcal/ContainedFrac;
  //std::cout<<"  getstuffdual norm ecal cer scint count "<<anecertotecal<<" "<<anescinttotecal<<std::endl;
  EEcal=(anescinttotecal-kappaEcal*anecertotecal)/(1-kappaEcal);


  float anecertothcal=necertothcal/meancerHcal/ContainedFrac;
  float anescinttothcal=nescinttothcal/meanscinHcal/ContainedFrac;
  EHcal=(anescinttothcal-kappaHcal*anecertothcal)/(1-kappaHcal);
  //  std::cout<<"necertothcal nescinttothcal meancerhcal meanscinhcal anecertothcal anescinttothcal kappahcal "<<necertothcal<<" "<<nescinttothcal<<" "<<meancerHcal<<" "<<meanscinHcal<<" "<<anecertothcal<<" "<<anescinttothcal<<" "<<kappaHcal<<std::endl;

  //std::cout<<"getstuffdual outputing ecal hcal "<<EEcal<<" "<<EHcal<<std::endl;


}
