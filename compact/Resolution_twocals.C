//
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
#include <map>
#include <algorithm>





// IGNORE THIS!!!!!!!!!!!!!!
// LOOK AT THE DECLARATION JUST AT THE crystalana def!!!!!!!!!!!!!!!!!!!!
const int nsliceecal = 8;
std::string nameecalslice[nsliceecal] = {"air","PD1","crystal1","gap1","middlemat","gap2","crystal2","PD2"};
int SCECOUNT=1;
bool dodualcorr=1;
float timecut=1000;
double hoverehcalscint(1.),hoverehcalcer(1.);
double hovereecalscint(1.),hovereecalcer(1.);
double kappaHcal;
double kappaEcal;

// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
//const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
  typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
  typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


void SCEDraw1 (TCanvas* canv, const char* name, TH1F* h1, const char* outfile, bool logy);
void SCEDraw1tp (TCanvas* canv, const char* name, TProfile* h1, const char* outfile);
void SCEDraw1_2D (TCanvas* canv, const char* name, TH2F* h1, const char* outfile,float eohS,float eohC);
void SCEDraw2 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile,bool logy);
void SCEDraw3 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile, bool logy);


void getStuff(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype,  bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
	      float  &eesum,float &eesumair, float &eesumdead, float &eesumcrystal,float &eesumPDe,float &eesumfiber1, float &eesumfiber2,float &eesumabs,float &eesumPDh,float &eesumedge,float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal,
	      float &timecut, float &eecaltimecut, float &ehcaltimecut,
	      TH1F* eecaltime, TH1F* ehcaltime,
	      int &nin
);


void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaecal, float kappahcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal, 
		      CalHits* &ecalhits, CalHits* &hcalhits,
		      float &EEcal, float &EHcal,
	      float &timecut, float &eecaltimecut, float &ehcaltimecut
);

void getMeanPhot(map<string, int> mapecalslice, map<string, int> mapsampcalslice,  int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits, 
		 float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,
	      float &timecut, float &eecaltimecut, float &ehcaltimecut
);


map<string, int> mapecalslice; 
map<string,int>::iterator eii0;
map<string,int>::iterator eii1;
map<string,int>::iterator eii2;
map<string,int>::iterator eii3;
map<string,int>::iterator eii4;
map<string,int>::iterator eii5;
map<string,int>::iterator eii6;
map<string,int>::iterator eii7;

map<string, int> mapsampcalslice; 
map<string,int>::iterator sii1;
map<string,int>::iterator sii2;
map<string,int>::iterator sii3;
map<string,int>::iterator sii4;
map<string,int>::iterator sii5;
map<string,int>::iterator sii6;
map<string,int>::iterator sii7;
map<string,int>::iterator sii8;
map<string,int>::iterator sii9;


// hcal type 0=fiber, 1 = sampling
// gendet 1=active media photons, 2 = photodetector, 3=energy deposit 4 is a debug gendet
// ECALleaf is 
// Resolution(20,"./output/out_DualTestBeam_50GeV_e-.root","./output/out_DualTestBeam_50GeV_pi-.root",
// "./output/out_FSCEPonly_50GeV_e-.root","./output/out_FSCEPonly_50GeV_pi-.root",
// 50,1,1,0,1,3,"hists_50GeV.root","DRCNoSegment","DRFNoSegment",1)

void Resolution(int num_evtsmax, const char* einputfilename, const char* piinputfilename,
		const char* hcalonlyefilename, const char* hcalonlypifilename,
		const float beamEE, bool doecal, bool dohcal, int hcaltype, bool doedge, int gendet, const char* outputfilename,const char* ECALleaf, const char* HCALleaf,bool doplots) {


  // these must correspond to the "slice" physvolid used in DRCrys_geo
  // these correspond to slices in scepcal_drcrystal.xml
  // surely there is a better way to do this
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

  int ihaha,num_evt;


  TBranch* b_mc;
  TBranch*  b_ecal;
  TBranch* b_hcal;
  TBranch* b_edge;
   


  // define histograms

  // calorimeter infor
  TH1F *ehchan = new TH1F("ehchan","channel ID number",1028,0.,1028);
  TH1F *phchan = new TH1F("phchan","channel ID number",1028,0.,1028);

  TH1F *ehcEcalE = new TH1F("ehcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);
  TH1F *phcEcalE = new TH1F("phcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);

  TH1F *ehcHcalE = new TH1F("ehcHcalE","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalEHO = new TH1F("ehcHcalEHO","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE = new TH1F("phcHcalE","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalEHO = new TH1F("phcHcalEHO","sum fiber hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalE1 = new TH1F("ehcHcalE1","fiber 1 hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalE1HO = new TH1F("ehcHcalE1HO","fiber 1 hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE1 = new TH1F("phcHcalE1","fiber 1 hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE1HO = new TH1F("phcHcalE1HO","fiber 1 hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalE2 = new TH1F("ehcHcalE2","fiber 2 hcal energy / beam E",200,0.,0.5);
  TH1F *ehcHcalE2HO = new TH1F("ehcHcalE2HO","fiber 2 hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE2 = new TH1F("phcHcalE2","fiber 2 hcal energy / beam E",200,0.,0.5);
  TH1F *phcHcalE2HO = new TH1F("phcHcalE2HO","fiber 2 hcal energy / beam E",200,0.,0.5);

  TH1F *ehcEdgeE = new TH1F("ehcEdgeE","sum escaping / beam E",100,0.,1.5);
  TH1F *phcEdgeE = new TH1F("phcEdgeE","sum escaping / beam E",100,0.,1.5);


  TH1F *ehcEdgeR = new TH1F("ehcEdgeR","beam - sum escaping / beam E",100,0.,1.5);
  TH1F *phcEdgeR = new TH1F("phcEdgeR","beam - sum escaping / beam E",100,0.,1.5);


  TH1F *ehcEcalncer = new TH1F("ehcEcalncer","total number of ecal cerenkov",  500,0.,1.5);
  TH1F *phcEcalncer = new TH1F("phcEcalncer","total number of ecal cerenkov", 500,0.,1.5);

  TH1F *ehcHcalncer = new TH1F("ehcHcalncer","total number of hcal cerenkov",   500,0.,1.5);
  TH1F *phcHcalncer = new TH1F("phcHcalncer","total number of hcal cerenkov",  500,0.,1.5);
  TH1F *fullEcalCalncer = new TH1F("fullEcalCalncer","total number of cerenkov Ecal cal",  500,0.,1.5);


  TH1F *ehcEcalnscint = new TH1F("ehcEcalnscint","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalnscint = new TH1F("phcEcalnscint","total number of ecal scitillation", 500,0.,1.5);

  TH1F *ehcHcalnscint = new TH1F("ehcHcalnscint","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalnscint = new TH1F("phcHcalnscint","total number of hcal scitillation", 500,0.,1.5);
  TH1F *fullEcalCalnscint = new TH1F("fullEcalCalnscint","total number of scitillation Ecal cal", 500,0.,1.5);

  TH1F *ehcEcalcorr = new TH1F("ehcEcalcorr","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalcorr = new TH1F("phcEcalcorr","total number of ecal scitillation", 500,0.,1.5);

  TH1F *ehcHcalcorr = new TH1F("ehcHcalcorr","e hcal dual", 500,0.,1.5);
  TH1F *phcHcalcorr = new TH1F("phcHcalcorr","pi hcal dual", 500,0.,1.5);
  TH1F *phcFulldetcorr = new TH1F("phcFulldetcorr","pi full detector dual", 500,0.,1.5);

 
  TH2F *ehcEcalNsNc = new TH2F("ehcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcEcalNsNc = new TH2F("phcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalNsNc = new TH2F("ehcHcalNsNc","hcal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcHcalNsNc = new TH2F("phcHcalNsNc","hcal ncer versus nscint",500,0.,1.5,500,0.,1.5);


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

  TH2F *enscvni = new TH2F("enscvni","electron scint E versus number inelastic",100,0.,1.2,100,0.,1000);
  TH2F *pinscvni = new TH2F("pinscvni","pion scint E versus number inelastic",100,0.,1.2,100,0.,1000);




  //calibration constants
  float meanscinEcal(0),meanscinHcal(0),meancerEcal(0),meancerHcal(0);
  float meaneecaltimecut(0),meanehcaltimecut(0);


  
  //****************************************************************************************************************************
  // process electrons in full calorimeter first pass to get calibration of single detector or ECAL of combined detector

  TFile* ef = TFile::Open(einputfilename);
  TTree* et = (TTree*)ef->Get("EVENT;1");
  b_mc= et->GetBranch("MCParticles");
  if(doecal) b_ecal = et->GetBranch(ECALleaf);
  if(dohcal) b_hcal = et->GetBranch(HCALleaf);
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");
  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<std::endl<<std::endl<<"num_evt for electron file full calorimeter is  "<<num_evt<<std::endl<<std::endl;
  
  // loop over events 


  
  if(num_evt>0) {  

    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);

    float ajunk1,ajunk2;
    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<std::endl<<"event number first pass electron file full detector is "<<ievt<<std::endl;
      getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, doecal, dohcal, hcaltype, b_ecal,b_hcal, ecalhits, hcalhits, meanscinEcal, ajunk1, meancerEcal, ajunk2,timecut, meaneecaltimecut, meanehcaltimecut);
    }

    std::cout<<"done with getMeanPhot full detector"<<std::endl;
    meanscinEcal=meanscinEcal/num_evt;
    meancerEcal=meancerEcal/num_evt;
     std::cout<<"mean scint ecal is "<<meanscinEcal<<std::endl;
     std::cout<<"mean cer ecal is "<<meancerEcal<<std::endl;
     std::cout<<"done with electron calibration of ECAL"<<std::endl<<std::endl;
  }
  ef->Close();

  //****************************************************************************************************************************
  // for calibration of hcal,  if have both an ecal and hcal


  if(doecal&&dohcal ) {
    TFile* efa = TFile::Open(hcalonlyefilename);
    TTree* eta = (TTree*)efa->Get("EVENT;1");
    b_mc= eta->GetBranch("MCParticles");
    b_hcal = eta->GetBranch(HCALleaf);
    //std::cout<<"b_mc b_hcal are "<<b_mc<<" "<<b_hcal<<std::endl;
    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<std::endl<<std::endl<<"num_evt for electron hcal calibration file is  "<<num_evt<<std::endl;
  
  // loop over events 

    float meanscinEcala(0),meanscinHcalla(0),meancerEcala(0),meancerHcalla(0), meaneecaltimecutla(0),meanehcaltimecutla(0);

  
    if(num_evt>0) {  

      CalHits* ecalhitsa = new CalHits();
      CalHits* hcalhitsa = new CalHits();
      b_hcal->SetAddress(&hcalhitsa);

      std::cout<<" branches set"<<std::endl;

    // first pass through file
      float ajunk1,ajunk2;
      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<std::endl<<"event number first pass is "<<ievt<<std::endl;
	getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, 0, dohcal, hcaltype, b_ecal,b_hcal, ecalhitsa, hcalhitsa, ajunk1, meanscinHcalla, ajunk2, meancerHcalla,timecut,meaneecaltimecutla,meanehcaltimecutla);
      }

      std::cout<<"done with getMeanPhot for hcal calibration file"<<std::endl;
      meanscinHcal=meanscinHcalla/num_evt;
      meancerHcal=meancerHcalla/num_evt;
      std::cout<<"mean scint hcal is "<<meanscinHcal<<std::endl;
      std::cout<<"mean cer hcal is "<<meancerHcal<<std::endl;

    }

    efa->Close();
    std::cout<<"done with with electron calibration of hcal"<<std::endl<<std::endl;

  }



  //process electrons with calibration in full detector
 
  ef = TFile::Open(einputfilename);
  et = (TTree*)ef->Get("EVENT;1");

  b_mc= et->GetBranch("MCParticles");
  if(doecal) b_ecal = et->GetBranch(ECALleaf);
  if(dohcal) b_hcal = et->GetBranch(HCALleaf);
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");

  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<std::endl<<std::endl<<"num_evt for electron file full detector is  "<<num_evt<<std::endl;

  CalHits* ecalhits = new CalHits();
  if(doecal) b_ecal->SetAddress(&ecalhits);
  CalHits* hcalhits = new CalHits();
  if(dohcal) b_hcal->SetAddress(&hcalhits);
  CalHits* edgehits = new CalHits();
  if(doedge) b_edge->SetAddress(&edgehits);

   
  if(num_evt>0) {  

    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number second is "<<ievt<<std::endl;

      float eesum=0.;
      float eesumair=0;
      float eesumdead=0;
      float eesumcrystal=0;
      float eesumPDe=0;
      float eesumfiber1=0;
      float eesumfiber2=0.;
      float eesumabs=0;
      float eesumPDh=0;
      float eesumedge=0;
      float necertotecal=0;
      float nescinttotecal=0;
      float necertothcal=0;
      float nescinttothcal=0;
      float eecaltimecut=0.;
      float ehcaltimecut=0.;
      int nin=0;

      getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,eesum,eesumair,eesumdead,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal,timecut, eecaltimecut, ehcaltimecut,eecaltime,ehcaltime,nin);

    
      ehcEcalE->Fill(eesumcrystal/beamE);
      ehcHcalE->Fill((eesumfiber1+eesumfiber2)/beamE);
      ehcHcalE1->Fill(eesumfiber1/beamE);
      ehcHcalE2->Fill(eesumfiber2/beamE);
      ehcEdgeE->Fill(eesumedge/beamE);
      ehcEdgeR->Fill((beamE-eesumedge)/beamE);
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


      float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge+eesumdead;
      ehetrue->Fill(eachecks/beamE);

      float eedepcal=eesumair+eesumdead+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh;
      hedepcal->Fill(eedepcal/beamE);

      enscvni->Fill(eedepcal/beamE,nin);

      std::cout<<"GETSTUFF electrons full detector"<<std::endl;
      std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
      std::cout<<std::endl<<std::endl<<"total energy deposit "<<eesum/1000.<<std::endl;
      std::cout<<"       in air "<<eesumair/1000.<<std::endl;
      std::cout<<"       in ecal dead material "<<eesumdead/1000.<<std::endl;
      std::cout<<"       in photodetector ecal "<<eesumPDe/1000.<<std::endl;
      std::cout<<"       in crystal "<<eesumcrystal/1000.<<std::endl;
      std::cout<<"       in fiber1 "<<eesumfiber1/1000.<<std::endl;
      std::cout<<"       in fiber2 "<<eesumfiber2/1000.<<std::endl;
      std::cout<<"       in absorber "<<eesumabs/1000.<<std::endl;
      std::cout<<"       in photodetect hcal "<<eesumPDh/1000.<<std::endl;
      std::cout<<"       escaping detector "<<eesumedge/1000.<<std::endl;

      std::cout<<"       sum individual "<<eachecks/1000.<<std::endl;
      std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
      std::cout<<"   ratio to incident energy "<<eachecks/beamE<<std::endl;

      std::cout<<"total number of cherenkov ecal is "<<necertotecal<<std::endl;
      std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl;
      std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
      std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl;
      std::cout<<"number inelastic is "<<nin<<std::endl<<std::endl;


    }  //end loop over events
  }  // end if no events

  //  float amean = hceest->GetMean();

  ef->Close();
  std::cout<<"done with getstuff electrons full detector"<<std::endl<<std::endl;

  //process electrons with calibration in HCAL only

  if(doecal&&dohcal ) {
    TFile* efa = TFile::Open(hcalonlyefilename);
    TTree* eta = (TTree*)efa->Get("EVENT;1");

    b_mc= eta->GetBranch("MCParticles");
    if(dohcal) b_hcal = eta->GetBranch(HCALleaf);
    if(doedge) b_edge = eta->GetBranch("EdgeDetNoSegment");

    
    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<std::endl<<std::endl<<"num_evt for hcal only electron file is  "<<num_evt<<std::endl;


    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);

    if(num_evt>0) {  

      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number second is "<<ievt<<std::endl;
	float eesum=0.;
	float eesumair=0;
	float eesumdead=0;
	float eesumcrystal=0;
	float eesumPDe=0;
	float eesumfiber1=0;
	float eesumfiber2=0.;
	float eesumabs=0;
	float eesumPDh=0;
	float eesumedge=0;
	float necertotecal=0;
	float nescinttotecal=0;
	float necertothcal=0;
	float nescinttothcal=0;
	float eecaltimecut=0.;
	float ehcaltimecut=0.;
	
	int nin=0;

	getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, 0, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,eesum,eesumair,eesumdead,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal,timecut, eecaltimecut, ehcaltimecut,eecaltime,ehcaltime,nin);


	ehcHcalEHO->Fill((eesumfiber1+eesumfiber2)/beamE);
	ehcHcalE1HO->Fill(eesumfiber1/beamE);
	ehcHcalE2HO->Fill(eesumfiber2/beamE);

	float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge+eesumdead;
	float eedepcal=eesumair+eesumdead+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh;



	std::cout<<"GETSTUFF hcal only electrons"<<std::endl;
	std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
	std::cout<<std::endl<<std::endl<<"total energy deposit "<<eesum/1000.<<std::endl;
	std::cout<<"       in air "<<eesumair/1000.<<std::endl;
	std::cout<<"       in ecal dead material "<<eesumdead/1000.<<std::endl;
	std::cout<<"       in photodetector ecal "<<eesumPDe/1000.<<std::endl;
	std::cout<<"       in crystal "<<eesumcrystal/1000.<<std::endl;
	std::cout<<"       in fiber1 "<<eesumfiber1/1000.<<std::endl;
	std::cout<<"       in fiber2 "<<eesumfiber2/1000.<<std::endl;
	std::cout<<"       in absorber "<<eesumabs/1000.<<std::endl;
	std::cout<<"       in photodetect hcal "<<eesumPDh/1000.<<std::endl;
	std::cout<<"       escaping detector "<<eesumedge/1000.<<std::endl;
	
	std::cout<<"       sum individual "<<eachecks/1000.<<std::endl;
	std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	std::cout<<"   ratio to incident energy "<<eachecks/beamE<<std::endl;

	std::cout<<"total number of cherenkov ecal is "<<necertotecal<<std::endl;
	std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl;
	std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
	std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl;
	std::cout<<"number inelastic is "<<nin<<std::endl<<std::endl;

      }  //end loop over events
    }  // end if no events
    efa->Close();
  } // end both ecal and hcal



  std::cout<<"done with getstuff electrons hcal only"<<std::endl<<std::endl;

  
  //****************************************************************************************************************************
  // process pions for hcalonly to get fiber e/h

  if(doecal&&dohcal) {
    TFile* pifa = TFile::Open(hcalonlypifilename);
    TTree* pita = (TTree*)pifa->Get("EVENT;1");
    if(pifa==0) std::cout<<" no file "<<std::endl;
    if(pita==0) std::cout<<" no event "<<std::endl;
    std::cout<<"pion file hcal only open"<<std::endl;

    b_mc= pita->GetBranch("MCParticles");
    if(doecal) b_ecal = pita->GetBranch(ECALleaf);
    if(dohcal) b_hcal = pita->GetBranch(HCALleaf);
    if(doedge) b_edge = pita->GetBranch("EdgeDetNoSegment");
    std::cout<<"pion branches found"<<std::endl;

    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<"num_evt for pion file is  "<<num_evt<<std::endl;
  
  // loop over events 
  
    if(num_evt>0) {  

      CalHits* hcalhits = new CalHits();
      if(dohcal) b_hcal->SetAddress(&hcalhits);
      CalHits* edgehits = new CalHits();
      if(doedge) b_edge->SetAddress(&edgehits);

      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number pion hcanonly detector "<<ievt<<std::endl;

	float pesum(0.),pesumair(0.),pesumdead(0.),pesumcrystal(0.),pesumPDe(0.),pesumfiber1(0.),pesumfiber2(0.),pesumabs(0.),pesumPDh(0.),pesumedge(0.),npcertotecal(0.),npscinttotecal(0.),npcertothcal(0.),npscinttothcal(0.),eecaltimecut(0.),ehcaltimecut(0.);
	int nin=0;

	getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, 0, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,pesum,pesumair,pesumdead,pesumcrystal,pesumPDe,pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumedge,npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,timecut, eecaltimecut, ehcaltimecut,piecaltime,pihcaltime,nin);
	std::cout<<" yuck ehcaltimecut is "<<ehcaltimecut<<std::endl;

	phcEcalE->Fill(pesumcrystal/beamE);
	phcHcalE->Fill((pesumfiber1+pesumfiber2)/beamE);
	phcHcalE1->Fill(pesumfiber1/beamE);
	phcHcalE2->Fill(pesumfiber2/beamE);
	phcEdgeE->Fill(pesumedge/beamE);
	phcEdgeR->Fill((beamE-pesumedge)/beamE);
	phcEcalncer->Fill(npcertotecal/meancerEcal);
	phcHcalncer->Fill(npcertothcal/meancerHcal);
	phcEcalnscint->Fill(npscinttotecal/meanscinEcal);
	phcHcalnscint->Fill(npscinttothcal/meanscinHcal);
	phcHcalf1f2->Fill(pesumfiber1/1000.,pesumfiber2/1000.);
	if((pesumfiber1+pesumfiber2)>0) phaphcal->Fill((pesumfiber1+pesumfiber2)/(pesumabs+pesumfiber1+pesumfiber2));
	pheest->Fill((pesumcrystal+(pesumfiber1+pesumfiber2))/beamE);

	float rrr=npscinttotecal/meanscinEcal;
	float rrr2=npcertotecal/meancerEcal;
	float rrx=npscinttothcal/meanscinHcal;
	float rrx2=npcertothcal/meancerHcal;
	if( (rrr>0.1)&&(rrr2>0.2) )
	  phcEcalNsNc->Fill(rrr,rrr2);  
	if( (rrx>0.1)&&(rrx2>0.2) )
	  phcHcalNsNc->Fill(rrx,rrx2);  

	phcEcalMarco->Fill(rrr2/rrr,rrr2);
	phcHcalMarco->Fill(rrx2/rrx,rrx2);

	float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumedge+pesumdead;
	float pedepcal=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumdead;
	hpdepcal->Fill(pedepcal/beamE);
	phetrue->Fill(pachecks/beamE);
	pinscvni->Fill(pedepcal/beamE,nin);

	std::cout<<std::endl<<std::endl;
	std::cout<<"GETSTUFF pions HCALONLY detector"<<std::endl;
	std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
	std::cout<<"total energy deposit "<<pesum/1000.<<std::endl;
	std::cout<<"       in air "<<pesumair/1000.<<std::endl;
	std::cout<<"       in ecal dead "<<pesumdead/1000.<<std::endl;
	std::cout<<"       in photodetector ecal "<<pesumPDe/1000.<<std::endl;
	std::cout<<"       in crystal "<<pesumcrystal/1000.<<std::endl;
	std::cout<<"       in fiber1 "<<pesumfiber1/1000.<<std::endl;
	std::cout<<"       in fiber2 "<<pesumfiber2/1000.<<std::endl;
	std::cout<<"       in absorber "<<pesumabs/1000.<<std::endl;
	std::cout<<"       in photodetect hcal "<<pesumPDh/1000.<<std::endl;
	std::cout<<"       escaping detector "<<pesumedge/1000.<<std::endl;

	std::cout<<"       sum individual "<<pachecks/1000.<<std::endl;
	std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	std::cout<<"   ratio to incident energy "<<pachecks/beamE<<std::endl;

	std::cout<<"total number of cherenkov ecal is "<<npcertotecal<<std::endl;
	std::cout<<"total number of scintillator ecal is "<<npscinttotecal<<std::endl;
	std::cout<<"total number of cherenkov hcal is "<<npcertothcal<<std::endl;
	std::cout<<"total number of scintillator hcal is "<<npscinttothcal<<std::endl;
	std::cout<<" number of inelastic is "<<nin<<std::endl<<std::endl;


      }  //end loop over events
    }  // end if no events
    pifa->Close();
    hoverehcalscint=phcHcalnscint->GetMean();
    hoverehcalcer=phcHcalncer->GetMean();
    kappaHcal= (1-hoverehcalscint)/(1.-hoverehcalcer);
    std::cout<<" hoverehcalscint hoverehcalcer kappa hcal are "<<hoverehcalscint<<" "<<hoverehcalcer<<" "<<kappaHcal<<std::endl<<std::endl;

  
  }
  
      


  
  //****************************************************************************************************************************
  // process pions for full detector and try to get e/h for ecal

  if(doecal&&dohcal) {
    TFile* pif = TFile::Open(piinputfilename);
    TTree* pit = (TTree*)pif->Get("EVENT;1");
    if(pif==0) std::cout<<" no file "<<std::endl;
    if(pit==0) std::cout<<" no event "<<std::endl;
    std::cout<<"pion file full detector open"<<std::endl;

    b_mc= pit->GetBranch("MCParticles");
    if(doecal) b_ecal = pit->GetBranch(ECALleaf);
    if(dohcal) b_hcal = pit->GetBranch(HCALleaf);
    if(doedge) b_edge = pit->GetBranch("EdgeDetNoSegment");
    std::cout<<"pion branches found"<<std::endl;

    ihaha = b_mc->GetEntries();
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<"num_evt for pion file full is  "<<num_evt<<std::endl;
  
  // loop over events 
  
    if(num_evt>0) {  

      CalHits* ecalhits = new CalHits();
      if(doecal) b_ecal->SetAddress(&ecalhits);
      CalHits* hcalhits = new CalHits();
      if(dohcal) b_hcal->SetAddress(&hcalhits);
      CalHits* edgehits = new CalHits();
      if(doedge) b_edge->SetAddress(&edgehits);

      for(int ievt=0;ievt<num_evt; ++ievt) {
	if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number pion full detector "<<ievt<<std::endl;

	float pesum(0.),pesumair(0.),pesumdead(0.),pesumcrystal(0.),pesumPDe(0.),pesumfiber1(0.),pesumfiber2(0.),pesumabs(0.),pesumPDh(0.),pesumedge(0.),npcertotecal(0.),npscinttotecal(0.),npcertothcal(0.),npscinttothcal(0.),eecaltimecut(0.),ehcaltimecut(0.);
	int nin=0;

	getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,pesum,pesumair,pesumdead,pesumcrystal,pesumPDe,pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumedge,npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,timecut, eecaltimecut, ehcaltimecut,piecaltime,pihcaltime,nin);
	std::cout<<" yuck ehcaltimecut is "<<ehcaltimecut<<std::endl;


	fullEcalCalnscint->Fill((npscinttotecal/meanscinEcal)/(1.-npscinttothcal/meanscinHcal));
	fullEcalCalncer->Fill((npcertotecal/meancerEcal)/(1.-npcertothcal/meancerHcal));
	std::cout<<(npscinttotecal/meanscinEcal)<<" "<<(1.-npscinttothcal/meanscinHcal)<<std::endl;
	std::cout<<(npcertotecal/meancerEcal)<<" "<<(1.-npcertothcal/meancerHcal)<<std::endl;

	
	float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumedge+pesumdead;
	float pedepcal=pesumair+pesumPDe+pesumcrystal+pesumfiber1+pesumfiber2+pesumabs+pesumPDh+pesumdead;
	std::cout<<std::endl<<std::endl;
	std::cout<<"GETSTUFF pion full detector"<<std::endl;
	std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
	std::cout<<"total energy deposit "<<pesum/1000.<<std::endl;
	std::cout<<"       in air "<<pesumair/1000.<<std::endl;
	std::cout<<"       in ecal dead "<<pesumdead/1000.<<std::endl;
	std::cout<<"       in photodetector ecal "<<pesumPDe/1000.<<std::endl;
	std::cout<<"       in crystal "<<pesumcrystal/1000.<<std::endl;
	std::cout<<"       in fiber1 "<<pesumfiber1/1000.<<std::endl;
	std::cout<<"       in fiber2 "<<pesumfiber2/1000.<<std::endl;
	std::cout<<"       in absorber "<<pesumabs/1000.<<std::endl;
	std::cout<<"       in photodetect hcal "<<pesumPDh/1000.<<std::endl;
	std::cout<<"       escaping detector "<<pesumedge/1000.<<std::endl;

	std::cout<<"       sum individual "<<pachecks/1000.<<std::endl;
	std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
	std::cout<<"   ratio to incident energy "<<pachecks/beamE<<std::endl;

	std::cout<<"total number of cherenkov ecal is "<<npcertotecal<<std::endl;
	std::cout<<"total number of scintillator ecal is "<<npscinttotecal<<std::endl;
	std::cout<<"total number of cherenkov hcal is "<<npcertothcal<<std::endl;
	std::cout<<"total number of scintillator hcal is "<<npscinttothcal<<std::endl;
	std::cout<<" number of inelastic is "<<nin<<std::endl<<std::endl;


      }  //end loop over events
    }  // end if no events

    hovereecalscint=fullEcalCalnscint->GetMean();
    hovereecalcer=fullEcalCalncer->GetMean();
    kappaEcal= (1-hovereecalscint)/(1.-hovereecalcer);
    std::cout<<" hoverehcalscint hoverehcalcer kappa ecal are "<<hovereecalscint<<" "<<hovereecalcer<<" "<<kappaEcal<<std::endl<<std::endl;

      pif->Close();

  
  }
  
      





  
  // do dual corr  full detector

  if(dodualcorr) {

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
  
     if(num_evt>0) {  


       CalHits* ecalhits = new CalHits();
       if(doecal) b_ecal->SetAddress(&ecalhits);
       CalHits* hcalhits = new CalHits();
       if(dohcal) b_hcal->SetAddress(&hcalhits);
       CalHits* edgehits = new CalHits();
       if(doedge) b_edge->SetAddress(&edgehits);

       
       for(int ievt=0;ievt<num_evt; ++ievt) {
	 if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number pion is "<<ievt<<std::endl;

	 float EcorEcal(0),EcorHcal(0),eecaltimecutcor(0),ehcaltimecutcor(0);


	 getStuffDualCorr(mapecalslice, mapsampcalslice, gendet, kappaEcal, kappaHcal, meanscinEcal, meancerEcal, meanscinHcal, meancerHcal,  ievt,doecal,dohcal, hcaltype, b_ecal,b_hcal, ecalhits,hcalhits, EcorEcal, EcorHcal,timecut, eecaltimecutcor, ehcaltimecutcor);

	 phcEcalcorr->Fill(EcorEcal);
	 phcHcalcorr->Fill(EcorHcal);
	 phcFulldetcorr->Fill(EcorEcal+EcorHcal);



       }   //end loop over events


  // close pion file
	 pif->Close();
    

     }  // end if no events

  }


  
  //***********************************************************************************************************************
  if(doplots) {

  TCanvas* c1;
  SCEDraw2(c1,"c1",ehetrue,phetrue,"junk1.png",0);

  TCanvas* c1b;
  SCEDraw2(c1b,"c1b",ehcEdgeE,phcEdgeE,"junk1b.png",0);

  TCanvas* c1c;
  SCEDraw2(c1c,"c1c",ehcEdgeR,phcEdgeR,"junk1c.png",0);


  TCanvas* c1d;
  SCEDraw1(c1d,"c1d",hedepcal,"junk1d.png",0);
  TCanvas* c1e;
  SCEDraw1(c1e,"c1e",hpdepcal,"junk1e.png",0);

  

  if(doecal) {
    TCanvas* ce2;
    SCEDraw2(ce2,"ce2",ehcEcalE,phcEcalE,"junke2.png",0);
    TCanvas* ce3;
    SCEDraw2(ce3,"ce3",ehcEcalncer,ehcEcalnscint,"junke3.png",0);
    TCanvas* ce4;
    SCEDraw3(ce4,"ce4",phcEcalncer,phcEcalnscint,phcEcalcorr,"junke4.png",0);
    TCanvas* ce5;
    SCEDraw1_2D(ce5,"ce5",ehcEcalNsNc,"junke5.png",0.,0.);
    TCanvas* ce5b;
    SCEDraw1_2D(ce5b,"ce5b",ehcEcalMarco,"junke5b.png",0.,0.);
    TCanvas* ce6;
    SCEDraw1_2D(ce6,"ce6",phcEcalNsNc,"junke6.png",0.,0.);
    TCanvas* ce6b;
    SCEDraw1_2D(ce6b,"ce6b",phcEcalMarco,"junke6b.png",0.,0.);
    TCanvas* ce7;
    SCEDraw2(ce7,"ce7",eecaltime,piecaltime,"junke7.png",1);

  }


  
  if(dohcal) {
    TCanvas* ch2;
    SCEDraw2(ch2,"ch2",ehcHcalE,phcHcalE,"junkh2.png",0);
    TCanvas* ch2a;
    SCEDraw2(ch2a,"ch2a",ehcHcalE1,phcHcalE1,"junkh2a.png",0);
    TCanvas* ch2b;
    SCEDraw2(ch2b,"ch2b",ehcHcalE2,phcHcalE2,"junkh2b.png",0);
    TCanvas* ch2HO;
    SCEDraw2(ch2HO,"ch2HO",ehcHcalEHO,phcHcalEHO,"junkh2HO.png",0);
    TCanvas* ch2aHO;
    SCEDraw2(ch2aHO,"ch2aHO",ehcHcalE1HO,phcHcalE1HO,"junkh2aHO.png",0);
    TCanvas* ch2bHO;
    SCEDraw2(ch2bHO,"ch2bHO",ehcHcalE2HO,phcHcalE2HO,"junkh2bHO.png",0);
    TCanvas* ch4a;
    SCEDraw2(ch4a,"ch4a",fullEcalCalncer,fullEcalCalnscint,"junkh4a.png",0);

    
    TCanvas* ch3;
    SCEDraw2(ch3,"ch3",ehcHcalncer,ehcHcalnscint,"junkh3.png",0);
    TCanvas* ch4;
    SCEDraw3(ch4,"ch4",phcHcalncer,phcHcalnscint,phcFulldetcorr,"junkh4.png",0);
    TCanvas* ch5;
    SCEDraw1_2D(ch5,"ch5",ehcHcalNsNc,"junkh5.png",0.,0.);
    TCanvas* ch5b;
    SCEDraw1_2D(ch5b,"ch5b",ehcHcalMarco,"junkhb.png",0.,0.);
    TCanvas* ch7;
    SCEDraw1_2D(ch7,"ch7",ehcHcalf1f2,"junkh7.png",0.,0.);
    TCanvas* ch7b;
    SCEDraw1_2D(ch7b,"ch7b",phcHcalf1f2,"junkh7b.png",0.,0.);

    TCanvas* ch8;
    SCEDraw2(ch8,"ch8",ehcaltime,pihcaltime,"junkh8.png",1);

  }


    TCanvas* ch9a;
    SCEDraw1_2D(ch9a,"ch9a",enscvni,"junkh9a.png",0.,0.);
    TCanvas* ch9b;
    SCEDraw1_2D(ch9b,"ch9b",pinscvni,"junkh9by.png",0.,0.);



  //TCanvas* c7;
  //SCEDrawp(c7,"c7",phcEcalNsNc_pfx,"junk7.png");

  }

  //***********************************************************************************************************

  TFile * out = new TFile(outputfilename,"RECREATE");



  ehcEcalE->Write();
  phcEcalE->Write();

  ehcHcalE->Write();
  phcHcalE->Write();
  ehcHcalEHO->Write();
  phcHcalEHO->Write();

  ehcHcalE1->Write();
  phcHcalE1->Write();
  ehcHcalE1HO->Write();
  phcHcalE1HO->Write();


  ehcHcalE2->Write();
  phcHcalE2->Write();
  ehcHcalE2HO->Write();
  phcHcalE2HO->Write();


  ehcEdgeE->Write();
  phcEdgeE->Write();



  ehcEdgeR->Write();
  phcEdgeR->Write();



   ehcEcalncer->Write();
   phcEcalncer->Write();

   ehcHcalncer->Write();
   phcHcalncer->Write();
   fullEcalCalncer->Write();

   ehcEcalcorr->Write();
   phcEcalcorr->Write();

   ehcHcalcorr->Write();
   phcHcalcorr->Write();
   phcFulldetcorr->Write();


  ehecal2d->Write();
  phecal2d->Write();


   ehcEcalnscint->Write();
   phcEcalnscint->Write();


   ehcHcalnscint->Write();
   phcHcalnscint->Write();
   fullEcalCalnscint->Write();

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
   pinscvni->Write();

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
  h1->SetMarkerColor(3);
  gStyle->SetOptFit();
  h1->Draw("*");



  canv->Print(outfile,".png");
  canv->Update();

  return;
}


void SCEDraw1_2D (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile, float eohS, float eohC) {

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
  h1->Draw("");

  TLine line = TLine(eohS,eohC,1.,1.);
  line.SetLineColor(kYellow);
  line.SetLineWidth(2);
  line.Draw("same");


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


void getMeanPhot(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits, 
		 float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal, 	      float &timecut, float &eecaltimecut, float &ehcaltimecut

){
  int nbyteecal, nbytehcal, nbyteedge;


  

  if(doecal) {
    if(ievt<SCECOUNT) std::cout<<"getMean phot ievt is "<<ievt<<std::endl;

    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
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

      if(gendet==1) {   // use photons as generated in otical material
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
    }
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);
    
      // hcal hits
    if(ievt<SCECOUNT) std::cout<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet ix iy ifiber iabs iphdet "<<std::endl;
    //if(ievt<SCECOUNT) std::cout<<"    ihitchan idet iy ix ilayer islice  "<<std::endl;

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
	  for(size_t j=0;j<zxzz.size(); j++) {
	    if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  }
	  }
	}
      }
      else {  // sampling

	int idet = (ihitchan) & 0x07;
	int iy = (ihitchan >>3) & 0xFFF;  
	int ix = (ihitchan >>15) & 0xFFF;  
	int ilayer = (ihitchan >>27) & 0xFFF;  
	int ibox2 = (ihitchan >> 39) & 0x03;
	int islice = (ihitchan >>41) & 0xF;  
 
	if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;

	if(gendet==1) {  // take light as generated in media
	  if(islice==(*sii3).second) {
	    meanscinHcal+=ahcalhit->nscintillator;
	    //	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*sii6).second) {  // cherenkov
	    meancerHcal+=ahcalhit->ncerenkov;
	    //	    std::cout<<"add ceren"<<std::endl;
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
	  if( islice==(*sii3).second) { // PS
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


      }


    }  // end loop over hcal hits
  }


}





void getStuff(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
	      float  &eesum,float &eesumair,float &eesumdead, float &eesumcrystal,float &eesumPDe,float &eesumfiber1,float &eesumfiber2,float &eesumabs,float &eesumPDh,float &eesumedge,float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal,
	      float &timecut, float &eecaltimecut, float &ehcaltimecut,
	      TH1F* eecaltime, TH1F* ehcaltime, int &nin){


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
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);





      //      if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
      long long int ihitchan=aecalhit->cellID;
      int idet = (ihitchan) & 0x07;
      int ix = (ihitchan >>3) & 0x3F ;  // is this right?
      if(ix>32) ix=ix-64;
      int iy =(ihitchan >>10) & 0x3F ; // is this right?
      if(iy>32) iy=iy-64;
      int  islice = (ihitchan >>17) & 0x07;
      int  ilayer = (ihitchan>> 20) & 0x07;



      	
      if((ilayer!=0)&&(ilayer!=1)) std::cout<<"danger danger will robinson ilayer not zero"<<std::endl;
      if(islice>nsliceecal) {
	std::cout<<"  danger danger will robinson islice nsliceecal are "<<islice<<" "<<nsliceecal<<std::endl;
	std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<std::endl;
      } else {
	if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;
      }
      

      float ae=aecalhit->energyDeposit;
      nin+=aecalhit->n_inelastic;


      // check contribs
      Contributions zxzz=aecalhit->truth;
      float hacheck=0.;
      for(size_t j=0;j<zxzz.size(); j++) {
	hacheck+=(zxzz.at(j)).deposit;
	if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
	eecaltime->Fill((zxzz.at(j)).time);
      }
      if(ae>0.001) {
	if(hacheck/ae<0.99999) std::cout<<"missing contribs: ecal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ae is "<<ae<<" ratio "<<hacheck/ae<<std::endl;
      }




      eesum+=ae;
      if(islice==(*eii0).second)eesumair+=ae;
      if(islice==(*eii1).second)eesumPDe+=ae;
      if(islice==(*eii2).second)eesumcrystal+=ae;
      if(islice==(*eii3).second)eesumair+=ae;
      if(islice==(*eii4).second)eesumdead+=ae;
      if(islice==(*eii5).second)eesumair+=ae;
      if(islice==(*eii6).second)eesumcrystal+=ae;
      if(islice==(*eii7).second)eesumPDe+=ae;




      if(gendet==1) {   // use photons as generated in otical material
	if((islice==(*eii2).second)||(islice==(*eii6).second)) {  // crystal
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
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
	//for(size_t j=0;j<zxzz.size(); j++) {
	//  if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
	//}
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
    //float ehcaltimecut=0.;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      float ah=ahcalhit->energyDeposit;
      nin+=ahcalhit->n_inelastic;

      // check contribs
      Contributions zxzz=ahcalhit->truth;
      float hacheck=0.;
      for(size_t j=0;j<zxzz.size(); j++) {
	hacheck+=(zxzz.at(j)).deposit;
	if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	ehcaltime->Fill((zxzz.at(j)).time);
      }
      if(ah>0.001) {
	if(hacheck/ah<0.99999) std::cout<<"missing contribs: hcal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ah is "<<ah<<" ratio "<<hacheck/ah<<std::endl;
      }




      eesum+=ah;
      /*
      std::cout<<"eesum now "<<eesum<<std::endl;
      std::cout<<"ehcaltimecut is "<<ehcaltimecut<<std::endl;
      std::cout<<eesumfiber1<<" "<<eesumfiber2<<" "<<eesumabs<<std::endl;
      */

      long long int ihitchan=ahcalhit->cellID;
      if(hcaltype==0) { // fiber

	int idet = (ihitchan) & 0xFF;
	int ilayer = (ihitchan >>8) & 0xFFF;  
	int itube = (ihitchan >>20) & 0xFFF;  
	//int iair=0; int itype=0;
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




	//int ifiber  =(ihitchan >>21) & 0x03;
	//int iabs=(ihitchan >>23) & 0x03;
	//int iphdet=(ihitchan >>25) & 0x03;
	//if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<std::endl;


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
	  //	  for(size_t j=0;j<zxzz.size(); j++) { 
	  //  if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  //}
	  }
	}

	if(ifiber==1) {eesumfiber1+=ah;}
	if(ifiber==2) {eesumfiber2+=ah;}
	if(iabs==1) {eesumabs+=ah;}
	if(iphdet>1) {eesumPDh+=ah;}
	//std::cout<<"   "<<ihitchan<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<" "<<eesumfiber<<" "<<eesumabs<<" "<<eesumPDh<<std::endl;


      }
      else {  // sampling


	int idet = (ihitchan) & 0x07;
	int iy = (ihitchan >>3) & 0xFFF;  
	int ix = (ihitchan >>15) & 0xFFF;  
	int ilayer = (ihitchan >>27) & 0xFFF;  
	int ibox2 = (ihitchan >> 39) & 0x03;
	int islice = (ihitchan >>41) & 0xF;  



	//std::cout<<" idet iy ix ilayer islice are "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;

	//std::cout<<"energy nscint ncer is "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<std::endl;


	//std::cout<<" ps ii is "<<(*sii3).second<<std::endl;
	//std::cout<<" quartz ii is "<<(*sii6).second<<std::endl;


	if(gendet==1) {  // take light as generated in media
	  if(islice==(*sii3).second) {
	    nescinttothcal+=ahcalhit->nscintillator;
	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*sii6).second) {  // chereknov
	    necertothcal+=ahcalhit->ncerenkov;
	    std::cout<<"add ceren"<<std::endl;
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
	  //std::cout<<" here "<<std::endl;
	  if(idet==6) {
	  if( islice==(*sii3).second) { //  ps
	    nescinttothcal+=ahcalhit->energyDeposit;
	    //	    if(i<10) std::cout<<" i nescinttothcal "<<i<<" "<<nescinttothcal<<std::endl;
	  }
	  if( islice==(*sii6).second ) {  // quartz
	    if(gendet==3) necertothcal+=ahcalhit->edeprelativistic;
	    if(gendet==4) necertothcal+=ahcalhit->energyDeposit;
	    //if(i<10) std::cout<<" i necertothcal "<<i<<" "<<necertothcal<<std::endl;
	  }
	  //	  for(size_t j=0;j<zxzz.size(); j++) {
	  //  if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	  //}
	  }
	}



	if( islice==(*sii3).second ) eesumfiber1+=ah; // scint
	if( islice==(*sii6).second ) eesumfiber2+=ah;  //cer
	if( (islice==(*sii1).second)||(islice==(*sii8).second)||(islice==(*sii9).second) ) eesumabs+=ah;
	if(  (islice==(*sii2).second) || (islice==(*sii4).second) ||  (islice==(*sii5).second) || (islice==(*sii7).second)) eesumPDh+=ah;


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

    }  // end loop over escaping hits
  }

  //std::cout<<" yuck 2 ehcaltimecut is "<<ehcaltimecut<<std::endl;


}


void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaEcal, float kappaHcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal, 
	      CalHits* &ecalhits, CalHits* &hcalhits, 
		      float &EEcal, float &EHcal, 	      float &timecut, float &eecaltimecut, float &ehcaltimecut)
{
  float necertotecal(0),nescinttotecal(0),necertothcal(0),nescinttothcal(0);
  int nbyteecal, nbytehcal, nbyteedge;

  std::cout<<" getstuff meanscinEcal is "<<meanscinEcal<<std::endl;
  if(ievt<SCECOUNT) std::cout<<"getstuffdualcorr phot ievt is "<<ievt<<std::endl;


  if(doecal) {
    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
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

      //      if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"getstuffdualcorr  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;



      
      if(gendet==1) {   // use photons as generated in otical material
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


	int idet = (ihitchan) & 0xFF;
	int ilayer = (ihitchan >>8) & 0xFFF;  
	int itube = (ihitchan >>20) & 0xFFF; 
	//int iair=0; int itype=0; 
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

	//int ifiber  =(ihitchan >>21) & 0x03;
	//int iabs=(ihitchan >>23) & 0x03;
	//int iphdet=(ihitchan >>25) & 0x03;
	//	if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<ifiber<<" "<<iabs<<" "<<iphdet<<std::endl;
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

	int idet = (ihitchan) & 0x07;
	int iy = (ihitchan >>3) & 0xFFF;  
	int ix = (ihitchan >>15) & 0xFFF;  
	int ilayer = (ihitchan >>27) & 0xFFF;  
	int ibox2 = (ihitchan >> 39) & 0x03;
	int islice = (ihitchan >>41) & 0xF;  


	//	if(ievt<SCECOUNT) std::cout<<"   "<<std::hex<<ihitchan<<std::dec<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<islice<<" "<<ilayer<<" "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<std::endl;
 


	if(gendet==1) {  // take light as generated in media
	  if(islice==(*sii3).second) {
	    nescinttothcal+=ahcalhit->nscintillator;
	    //	    std::cout<<"add scint"<<std::endl;
	  }
	  if(islice==(*sii6).second) {  // cherenkov
	    necertothcal+=ahcalhit->ncerenkov;
	    //	    std::cout<<"add ceren"<<std::endl;
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
	  if( islice==(*sii3).second ) { // ps
	    nescinttothcal+=ahcalhit->energyDeposit;
	    //if(i<10) std::cout<<" i nescinttothcal "<<i<<" "<<nescinttothcal<<std::endl;
	  }
	  if( islice==(*sii6).second ) {  // quartz
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
  std::cout<<"  getstuffdual ecal cer scint count "<<necertotecal<<" "<<nescinttotecal<<std::endl;
  float anecertotecal=necertotecal/meancerEcal;
  float anescinttotecal=nescinttotecal/meanscinEcal;
  std::cout<<"  getstuffdual norm ecal cer scint count "<<anecertotecal<<" "<<anescinttotecal<<std::endl;
  EEcal=(anescinttotecal-kappaEcal*anecertotecal)/(1-kappaEcal);


  float anecertothcal=necertothcal/meancerHcal;
  float anescinttothcal=nescinttothcal/meanscinHcal;
  EHcal=(anescinttothcal-kappaHcal*anecertothcal)/(1-kappaHcal);
  std::cout<<"necertothcal nescinttothcal meancerhcal meanscinhcal anecertothcal anescinttothcal kappahcal "<<necertothcal<<" "<<nescinttothcal<<" "<<meancerHcal<<" "<<meanscinHcal<<" "<<anecertothcal<<" "<<anescinttothcal<<" "<<kappaHcal<<std::endl;
  
  std::cout<<"getstuffdual outputing ecal hcal "<<EEcal<<" "<<EHcal<<std::endl;


}


