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

#include <iostream>
#include <sstream> // for ostringstream
#include <string>
using namespace std;


// This file is modified to accound for two-segmented ECAL crystals
// to run in batch copy: 
//     root -l -b -q 'Resolution.C(nevt,"electron_rootfile.root", "pion_rootfile.root","hcalonly_rootfile.root",energy,doecal,dohcal,doedge,gendet,hcaltype,"output.root","ECALleaf","HCALleaf")'
//     nevt                         number of events to process, e.g. 50
//     electron_rootfile.root       ddsim output with electron gun, e.g. out_DualTestBeam-dial_e-10_.root
//     pion_rootfile.root           ddsim output with electron gun, e.g. out_DualTestBeam-dial_e-10_.root
//     hcalonly_rootfile.root       for calibration hcal if both ecal + hcal, e.g. " "
//     energy                       particle gun energy, e.g. 20
//     doecal & dohcal & doedge     include ECAL & HCAL & Edges
//     gendet                       photon in:
//                                  active media (e.g. ecal crystal) = 1;
//                                  photodetector                    = 2;
//                                  energy deposit                   = 3;
//                                  debug                            = 4; 
//     hcaltype                     fiber=0, sampling=1
//     output.root                  output root file with all histograms --> to be used in ResvE.C
//     ECALleaf                     ECAL ttree in electron_rootfile.root and pion_rootfile.root, e.g. DRCNoSegment
//     HCALleaf                     HCAL ttree in electron_rootfile.root and pion_rootfile.root, e.g. DRFNoSegment



// IGNORE THIS!!!!!!!!!!!!!!
// LOOK AT THE DECLARATION JUST AT THE crystalana def!!!!!!!!!!!!!!!!!!!!
const int nsliceecal = 6;
std::string nameecalslice[nsliceecal] = {"air","PD1","crystal1","gap","crystal2","PD2"};
int SCECOUNT=1;
bool dodualcorr=1;
bool doplots=1;
float timecut=1000;

// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
//const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry
typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


void SCEDraw1(TCanvas* canv, const char* name, TH1F* h1, const char* outfile, bool logy);
void SCEDraw1tp(TCanvas* canv, const char* name, TProfile* h1, const char* outfile);
void SCEDraw1_2D(TCanvas* canv, const char* name, TH2F* h1, const char* outfile,float eohS,float eohC);
void SCEDraw2(TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile,bool logy);
void SCEDraw3(TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile, bool logy);


void getStuff(map<string,int> mapecalslice, map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, 
   int hcaltype,  bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge, CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
   float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber1, float &eesumfiber2,float &eesumabs,float &eesumPDh,
   float &eesumedge,float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal, float &timecut, float &eecaltimecut, 
   float &ehcaltimecut, TH1F* eecaltime, TH1F* ehcaltime, int &nin);


void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaecal, float kappahcal, 
   float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype, 
   TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, float &EEcal, float &EHcal,
   float &timecut, float &eecaltimecut, float &ehcaltimecut); //just for pions (??)

void getMeanPhot(map<string, int> mapecalslice, map<string, int> mapsampcalslice,  //input
   int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input
   TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, //input
   float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal,float &timecut, float &eecaltimecut, float &ehcaltimecut); // output


map<string, int> mapecalslice;
map<string,int>::iterator eii0;
map<string,int>::iterator eii1;
map<string,int>::iterator eii2;
map<string,int>::iterator eii3;
map<string,int>::iterator eii4;
map<string,int>::iterator eii5;

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

void Resolution(int num_evtsmax, const char* einputfilename, const char* piinputfilename, const char* hcalonlyefilename,
                const float beamEE, bool doecal, bool dohcal, int hcaltype, bool doedge, int gendet, const char* outputfilename,
                const char* ECALleaf, const char* HCALleaf){

  mapecalslice["air"]=0;
  mapecalslice["PD1"]=1;
  mapecalslice["crystal1"]=2;
  mapecalslice["gap"]=3;
  mapecalslice["crystal2"]=4;
  mapecalslice["PD2"]=5;

  eii0 = mapecalslice.find("air");
  eii1 = mapecalslice.find("PD1");
  eii2 = mapecalslice.find("crystal1");
  eii3 = mapecalslice.find("gap");
  eii4 = mapecalslice.find("crystal2");
  eii5 = mapecalslice.find("PD2");

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
  TBranch* b_ecal;
  TBranch* b_hcal;
  TBranch* b_edge;      

  const int num_histograms = 6;
  vector<TH1F*> eenergy;
  vector<TH1F*> pienergy;
  vector<TH1F*> ephoton;
  vector<TH1F*> piphoton;
  vector<string> energy_title = {"sum(ecalcrys_E) / beamE", "sum(fiber_E) / beamE", "sum(scintfiber_E) / beamE", "sum(cerfiber_E) / beamE", "sum(edge_E) / beamE", "sum(beamE-edge_E) / beamE"};
  vector<string> nphoton_title = {"ntotcer_ecal/meancer_ecal", "ntotcer_hcal/meancer_hcal", "ntotscint_ecal/meanscint_ecal", "ntotscint_hcal/meanscint_hcal"};
  for (size_t i = 0; i < num_histograms; ++i) {
    eenergy.push_back(new TH1F(Form("electron_%zu", i), energy_title[i].c_str(), 100, 0.0, 1.5));
    pienergy.push_back(new TH1F(Form("pion_%zu", i), energy_title[i].c_str(), 100, 0.0, 1.5));
    if (i>3) { continue;}
    ephoton.push_back(new TH1F(Form("eph_%zu",i), nphoton_title[i].c_str(),   100, 0.0, 1.5));
    piphoton.push_back(new TH1F(Form("piph_%zu",i), nphoton_title[i].c_str(), 100, 0.0, 1.5));
  }

  TH1F *phcEcalcorr = new TH1F("phcEcalcorr","pi- total number of ecal scitillation", 50,0.,1.5);  
  TH1F *phcHcalcorr = new TH1F("phcHcalcorr","pi- hcal dual", 50,0.,1.5);
  
  TH2F *ehcEcalNsNc = new TH2F("ehcEcalNsNc","e- ecal ncer versus nscint",50,0.,1.5,50,0.,1.5);
  TH2F *phcEcalNsNc = new TH2F("phcEcalNsNc","pi- ecal ncer versus nscint",50,0.,1.5,50,0.,1.5);
  TH2F *ehcHcalNsNc = new TH2F("ehcHcalNsNc","e- hcal ncer versus nscint",50,0.,1.5,50,0.,1.5);
  TH2F *phcHcalNsNc = new TH2F("phcHcalNsNc","pi- hcal ncer versus nscint",50,0.,1.5,50,0.,1.5);
  
  TH2F *ehcEcalMarco = new TH2F("ehcEcalMarco","e- ecal c v c/s",50,0.,1.5,50,0.,1.5);
  TH2F *phcEcalMarco = new TH2F("phcEcalMarco","pi- ecal c v c/s",50,0.,1.5,50,0.,1.5);
  TH2F *ehcHcalMarco = new TH2F("ehcHcalMarco","e- hcal c v c/s",50,0.,1.5,50,0.,1.5);
  TH2F *phcHcalMarco = new TH2F("phcHcalMarco","pi- hcal c v c/s",50,0.,1.5,50,0.,1.5);
  
  TH2F *ehcHcalf1f2 = new TH2F("ehcHcalf1f2","e- hcal f1e versus f2e",50,0.,2.5,50,0.,2.5);
  TH2F *phcHcalf1f2 = new TH2F("phcHcalf1f2","pi- hcal f1e versus f2e",50,0.,2.5,50,0.,2.5);
  
  TH2F *ehecal2d = new TH2F("ehecal2d","e- lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *phecal2d = new TH2F("phecal2d","e- lego of ecal", 41,-20.,20.,41,-20.,20.);
  
  TH1F *ehaphcal = new TH1F("ehaphcal","e- (scint_fiber + quartz_fiber)/ (total hcal energy)" ,50,0.,0.2);
  TH1F *phaphcal = new TH1F("phaphcal","pi- (scint_fiber + quartz_fiber)/ (total hcal energy)",50,0.,0.2);
  
  TH1F *eheest = new TH1F("eheest","e- ratio estimated to true energy",50,0.,1.);
  TH1F *pheest = new TH1F("pheest","pi- ratio estimated to true energy",50,0.,1.);
  
  TH1F *ehetrue = new TH1F("ehetrue","e- ratio deposited to incident energy",50,0.,1.1);
  TH1F *phetrue = new TH1F("phetrue","pi- ratio deposited to incident energy",50,0.,1.1);
  
  TH1F *hedepcal = new TH1F("hedepcal","e- all deposited energies",50,0.,1.1);
  TH1F *hpdepcal = new TH1F("hpdepcal","pi- all deposited energies",50,0.,1.1);
  
  TH1F *ehnecalcon = new TH1F("ehnecalcon","e- number contribs to ecal hit",1010,-10.,1000.);
  TH1F *phnecalcon = new TH1F("phnecalcon","pi- number contribs to ecal hit",1010,-10.,1000.);
  
  TH2F *ehzvst = new TH2F("ehzvst","e- z position of hit versus time ",100,-50.,300.,100,0.,100.); 
  TH2F *phzvst = new TH2F("phzvst","pi- z position of hit versus time ",100,-50.,300.,100,0.,100.); 
  
  TH1F *eecaltime = new TH1F("eecaltime","e- time of ecal const",100,0.,40.);
  TH1F *ehcaltime = new TH1F("ehcaltime","e- time of hcal const",100,0.,40.);
  TH1F *piecaltime = new TH1F("piecaltime","pi- time of ecal const",100,0.,40.);
  TH1F *pihcaltime = new TH1F("pihcaltime","pi- time of hcal const",100,0.,40.);
  
  TH2F *enscvni = new TH2F("enscvni","e- scint E versus number inelastic",100,0.,1.2,100,0.,1000);
  TH2F *pinscvni = new TH2F("pinscvni","pi- scint E versus number inelastic",100,0.,1.2,100,0.,1000);
  
  
  //****************************************************************************************************************************
  // process electrons
  
  TFile* ef = TFile::Open(einputfilename);
  TTree* et = (TTree*)ef->Get("EVENT;1");
  
  b_mc= et->GetBranch("MCParticles");
  if(doecal) b_ecal = et->GetBranch(ECALleaf);
  if(dohcal) b_hcal = et->GetBranch(HCALleaf);
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment"); 
  
  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  cout<<"num_evt for electron file is  "<<num_evt<<endl;
    
  float meanscinEcal(0),meanscinHcal(0),meancerEcal(0),meancerHcal(0);  
  float meaneecaltimecut(0),meanehcaltimecut(0);
  if(num_evt>0) { // loop over events
    std::cout<<" looping over events: "<<std::endl;
    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);
  
    for(int ievt=0;ievt<num_evt; ++ievt) {// first pass through file get hcal & ecal: meanCer and meanScint + timing info
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) cout<<"event number first pass is "<<endl;
      getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, doecal, dohcal, hcaltype, b_ecal,b_hcal, ecalhits, hcalhits, //input
      meanscinEcal, meanscinHcal, meancerEcal, meancerHcal,timecut, meaneecaltimecut, meanehcaltimecut); // output: 
    }
  
    cout<<"done with getMeanPhot"<<endl;
    meanscinEcal=meanscinEcal/num_evt;
    meanscinHcal=meanscinHcal/num_evt;
    meancerEcal=meancerEcal/num_evt;
    meancerHcal=meancerHcal/num_evt;
    std::cout<<"mean scint ecal is "<<meanscinEcal<<std::endl;
    std::cout<<"mean scint hcal is "<<meanscinHcal<<std::endl;
    std::cout<<"mean cer ecal is "<<meancerEcal<<std::endl;
    std::cout<<"mean cer hcal is "<<meancerHcal<<std::endl;
  
    // second pass through file
    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) cout<<"event number second is "<<ievt<<endl;
      float eesum(0), eesumair(0), eesumcrystal(0), eesumPDe(0), eesumfiber1(0), eesumfiber2(0), eesumabs(0), eesumPDh(0), eesumedge(0);
      float necertotecal(0), nescinttotecal(0), necertothcal(0), nescinttothcal(0), eecaltimecut(0), ehcaltimecut(0);
      int nin(0);
      getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,eesum,eesumair,eesumcrystal,eesumPDe,eesumfiber1,eesumfiber2,eesumabs,eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal,timecut, eecaltimecut, ehcaltimecut,eecaltime,ehcaltime,nin);

      vector<float> eesums   = {eesumcrystal, eesumfiber1, eesumfiber2, eesumfiber1 + eesumfiber2, eesumedge};
      vector<float> nphotons = {necertotecal/meancerEcal, necertothcal/meancerHcal, nescinttotecal/meanscinEcal, nescinttothcal/meanscinHcal};

      auto eit = eesums.begin();
      for (auto hist : eenergy) {
        hist->Fill(*eit);
        ++eit;
      }
      auto nit = nphotons.begin();
      for (auto hist: ephoton) {
        hist->Fill(*nit);
	++nit;
      }

  
      ehcHcalf1f2->Fill(eesumfiber1/1000.,eesumfiber2/1000.);
  
      if((eesumfiber1+eesumfiber2)>0) ehaphcal->Fill((eesumfiber1+eesumfiber2)/(eesumabs+eesumfiber1+eesumfiber2));
      eheest->Fill((eesumcrystal+(eesumfiber1+eesumfiber2))/beamE);
      float ttt=nescinttotecal/meanscinEcal;
      float ttt2=necertotecal/meancerEcal;
      float tty=nescinttothcal/meanscinHcal;
      float tty2=necertothcal/meancerHcal;
      if( (ttt>0.1)&& (ttt2>0.2) ) ehcEcalNsNc->Fill(ttt,ttt2);  
      if( (tty>0.1)&& (tty2>0.2) ) ehcHcalNsNc->Fill(tty,tty2);  
  
      ehcEcalMarco->Fill(ttt2/ttt,ttt2);
      ehcHcalMarco->Fill(tty2/tty,tty2);
  
  
      float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh+eesumedge;
      ehetrue->Fill(eachecks/beamE);
  
      float eedepcal=eesumair+eesumPDe+eesumcrystal+eesumfiber1+eesumfiber2+eesumabs+eesumPDh;
      hedepcal->Fill(eedepcal/beamE);
  
      enscvni->Fill(eedepcal/beamE,nin);
  
      cout<<"GETSTUFF electrons"<<endl;
      cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<endl;
      cout<<"   total energy deposit "<<eesum/1000.<<endl;
      cout<<"       in air "<<eesumair/1000.<<endl;
      cout<<"       in photodetector ecal "<<eesumPDe/1000.<<endl;
      cout<<"       in crystal "<<eesumcrystal/1000.<<endl;
      cout<<"       in fiber1 "<<eesumfiber1/1000.<<endl;
      cout<<"       in fiber2 "<<eesumfiber2/1000.<<endl;
      cout<<"       in absorber "<<eesumabs/1000.<<endl;
      cout<<"       in photodetect hcal "<<eesumPDh/1000.<<endl;
      cout<<"       escaping detector "<<eesumedge/1000.<<endl;
  
      cout<<"       sum individual "<<eachecks/1000.<<endl;
      cout<<"   incident energy "<<beamE/1000.<<endl;
      cout<<"   ratio to incident energy "<<eachecks/beamE<<endl;
  
      cout<<"total number of cherenkov ecal is "<<necertotecal<<endl;
      cout<<"total number of scintillator ecal is "<<nescinttotecal<<endl;
      cout<<"total number of cherenkov hcal is "<<necertothcal<<endl;
      cout<<"total number of scintillator hcal is "<<nescinttothcal<<endl;
      cout<<"number inelastic is "<<nin<<endl; 
    }  //end loop over events
  }  // end process electron if no events
  //  float amean = hceest->GetMean();
  ef->Close();
  cout<<"done with getstuff electrons"<<endl;
  //****************************************************************************************************************************
  // for calibration of hcal,  if have both an ecal and hcal 
  if(doecal&&dohcal ) {
    TFile* efa = TFile::Open(hcalonlyefilename);
    TTree* eta = (TTree*)efa->Get("EVENT;1");
  
    b_mc= eta->GetBranch("MCParticles");
    b_hcal = eta->GetBranch(HCALleaf);
    std::cout<<"b_mc b_hcal are "<<b_mc<<" "<<b_hcal<<std::endl;
 
    ihaha = b_mc->GetEntries();
    cout<<"ihaha "<<ihaha<<endl;
    num_evt= std::min(ihaha,num_evtsmax);
    std::cout<<std::endl<<std::endl<<"num_evt for electron hcal calibration file is  "<<num_evt<<std::endl;
    // loop over events 
    float meanscinEcala(0),meanfscinHcala(0),meancerEcala(0),meancerHcala(0), meaneecaltimecut(0),meanehcaltimecut(0);
    
    if(num_evt>0) { 
      CalHits* ecalhitsa = new CalHits();
      CalHits* hcalhitsa = new CalHits();
      b_hcal->SetAddress(&hcalhitsa);
  
      std::cout<<" branches set"<<std::endl;
      // first pass through file
      for(int ievt=0;ievt<num_evt; ++ievt) {
         if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<std::endl<<"event number first pass is "<<ievt<<std::endl;
         getMeanPhot(mapecalslice, mapsampcalslice, gendet, ievt, 0, dohcal, hcaltype, b_ecal,b_hcal, ecalhitsa, hcalhitsa, 
        meanscinEcal, meanscinHcal, meancerEcal, meancerHcal,timecut,meaneecaltimecut,meanehcaltimecut);
      }
  
      std::cout<<"done with getMeanPhot for hcal calibration file"<<std::endl;
      meanscinHcal=meanscinHcal/num_evt;
      meancerHcal=meancerHcal/num_evt;
      std::cout<<"mean scint hcal is "<<meanscinHcal<<std::endl;
      std::cout<<"mean cer hcal is "<<meancerHcal<<std::endl;
  
    } 
    efa->Close();
    std::cout<<"done with with electron calibration of hcal"<<std::endl;
  } // end if(doecal&&dohcal)
  
  //****************************************************************************************************************************
  // process pions  
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
  
  float b1Ecal=0.;float m1Ecal=1.;
  float b1Hcal=0.;float m1Hcal=1.;
  
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
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) cout<<"event number pion is "<<ievt<<endl;
      float pesum(0), pesumair(0), pesumcrystal(0), pesumPDe(0), pesumfiber1(0), pesumfiber2(0), pesumabs(0), pesumPDh(0), pesumedge(0);
      float npcertotecal(0), npscinttotecal(0), npcertothcal(0), npscinttothcal(0), eecaltimecut(0), ehcaltimecut(0);
      int nin(0);
      getStuff(mapecalslice, mapsampcalslice,  gendet, ievt, doecal, dohcal, hcaltype, doedge, b_ecal,b_hcal,b_edge,
        ecalhits,hcalhits,edgehits,pesum,pesumair,pesumcrystal,pesumPDe,pesumfiber1,pesumfiber2,pesumabs,pesumPDh,pesumedge,
	npcertotecal,npscinttotecal,npcertothcal,npscinttothcal,timecut, eecaltimecut, ehcaltimecut,piecaltime,pihcaltime,nin);



      vector<float> piesums   = {pesumcrystal, pesumfiber1, pesumfiber2, pesumfiber1 + pesumfiber2, pesumedge};
      vector<float> npiphotons = {npcertotecal/meancerEcal, npcertothcal/meancerHcal, npscinttotecal/meanscinEcal, npscinttothcal/meanscinHcal};

      auto pit = piesums.begin();
      for (auto hist : pienergy) {
        hist->Fill(*pit);
        ++pit;
      }
      auto nit = npiphotons.begin();
      for (auto hist : piphoton) {
	hist->Fill(*nit);
	++nit;
      }
  
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
      pinscvni->Fill(pedepcal/beamE,nin);
  
      std::cout<<std::endl<<std::endl;
      std::cout<<"GETSTUFF pions"<<std::endl;
      std::cout<<" ehcaltimecut is "<<ehcaltimecut/1000.<<std::endl;
      std::cout<<"total energy deposit "<<pesum/1000.<<std::endl;
      std::cout<<"       in air "<<pesumair/1000.<<std::endl;
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
      std::cout<<"total number of scintillator hcal is "<<npscinttothcal<<std::endl<<std::endl;
      std::cout<<" number of inelastic is "<<nin<<std::endl;
    }  //end loop over events
    if(dodualcorr) {
      std::cout<<" starting fits"<<std::endl;
      //** fits
      TF1 *gEcale = new TF1("gEcale","[0]*(x-1.)+1",0.,1.);
      TF1 *gHcale = new TF1("gHcale","[0]*(x-1.)+1",0.,1.);
      TF1 *gEcalp = new TF1("gEcalp","[0]*(x-1.)+1",0.,1.);
      TF1 *gHcalp = new TF1("gHcalp","[0]*(x-1.)+1",0.,1.);
  
      // fit to get e/h
      if(doecal) {
        TProfile* ehcEcalNsNc_pfx = ehcEcalNsNc->ProfileX();
        ehcEcalNsNc_pfx->Fit("gEcale","W");
      }
      if(dohcal) {
        TProfile* ehcHcalNsNc_pfx = ehcHcalNsNc->ProfileX();
        ehcHcalNsNc_pfx->Fit("gHcale","W");
      }
      if(doecal) {
        TProfile* phcEcalNsNc_pfx = phcEcalNsNc->ProfileX();
        phcEcalNsNc_pfx->Fit("gEcalp","W");
        //b1Ecal=gEcalp->GetParameter(0);
        //m1Ecal=gEcalp->GetParameter(1);
        m1Ecal=gEcalp->GetParameter(0);
        b1Ecal=1-m1Ecal;
        std::cout<<"for ecal b m are "<<b1Ecal<<" "<<m1Ecal<<std::endl;
      }
      double kappaEcal = 1+(b1Ecal/m1Ecal); 
      double hovereecalscint=piphoton[2]->GetMean();;
      double hovereecalcer=piphoton[0]->GetMean();
      kappaEcal= (1-hovereecalscint)/(1.-hovereecalcer);
      std::cout<<" kappa ecal is "<<kappaEcal<<std::endl; 
  
      if(dohcal) {
        TProfile* phcHcalNsNc_pfx = phcHcalNsNc->ProfileX();
        phcHcalNsNc_pfx->Fit("gHcalp","W");
        //    b1Hcal=gHcalp->GetParameter(0);
        //    m1Hcal=gHcalp->GetParameter(1);
        m1Hcal=gHcalp->GetParameter(0);
        b1Hcal=1-m1Hcal;
  
        std::cout<<"for hcal b m are "<<b1Hcal<<" "<<m1Hcal<<std::endl;
        TCanvas* cyuck;
        SCEDraw1tp(cyuck,"cyuck",phcHcalNsNc_pfx,"junkyuck.png");
      }
      double kappaHcal = 1+(b1Hcal/m1Hcal); 
      double hoverehcalscint=piphoton[3]->GetMean();;
      double hoverehcalcer=piphoton[1]->GetMean();
      kappaHcal= (1-hoverehcalscint)/(1.-hoverehcalcer);
  
      std::cout<<" hoverehcalscint hoverehcalcer kappa hcal are "<<hoverehcalscint<<" "<<hoverehcalcer<<" "<<kappaHcal<<std::endl;
      
      for(int ievt=0;ievt<num_evt; ++ievt) {
        if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number pion is "<<ievt<<std::endl;
        float EcorEcal(0),EcorHcal(0),eecaltimecutcor(0),ehcaltimecutcor(0); 
        getStuffDualCorr(mapecalslice, mapsampcalslice, gendet, kappaEcal, kappaHcal, meanscinEcal, meancerEcal, meanscinHcal, meancerHcal,  
           ievt,doecal,dohcal, hcaltype, b_ecal,b_hcal, ecalhits,hcalhits, EcorEcal, EcorHcal,timecut, eecaltimecutcor, ehcaltimecutcor);
        phcEcalcorr->Fill(EcorEcal);
        phcHcalcorr->Fill(EcorHcal);
      }  //end loop over events
    }
  }  // end if no events
  // close pion file
  pif->Close();
  
  
  //***********************************************************************************************************
  TFile * out = new TFile(outputfilename,"RECREATE");
  for (auto ee : eenergy)     {ee->Write();}
  for (auto pie : pienergy)   {pie->Write();}
  for (auto enph : ephoton)   {enph->Write();}
  for (auto pinph : piphoton) {pinph->Write();}
  phcEcalcorr->Write();
  ehcHcalcorr->Write();
  phcHcalcorr->Write();
  ehecal2d->Write();
  phecal2d->Write();
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
  //ehchan->Write();
  //phchan->Write();   
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
  //hmass_p[0]->Write(); 
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
} // end of Resolution

void SCEDraw1(TCanvas* canv,const char* name,TH1F* h1,const char* outfile, bool logy)
  {
  canv= new TCanvas(name,name,200,10,700,500);
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


void getMeanPhot(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, //input
   int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input
   TBranch* &b_ecal,TBranch* &b_hcal, CalHits* &ecalhits, CalHits* &hcalhits, //input
   float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal, //output
   float &timecut, float &eecaltimecut, float &ehcaltimecut){//output

   int nbyteecal, nbytehcal, nbyteedge;
   if(doecal) {
    cout<<"getMean i events = "<<ievt<<" SCECOUNT="<<SCECOUNT<<endl;
    if(ievt<SCECOUNT) cout<<"getMean phot ievt is "<<ievt<<endl;
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
        if((islice==(*eii2).second)||(islice==(*eii4).second) ) {  // crystal
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
        if( (islice==(*eii1).second)||(islice==(*eii5).second) ) { // either photo detector
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
      } 
      else if(gendet==3||gendet==4){
	if(idet==5) {
	  meanscinEcal+=aecalhit->energyDeposit;
	  if(gendet==3) meancerEcal+=aecalhit->edeprelativistic;
	  if(gendet==4) meancerEcal+=aecalhit->energyDeposit;
	  for(size_t j=0;j<zxzz.size(); j++) {
	    if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
	  }
	}
      }
    } //end of ecalhitsize
  } // end of if doecal

  if(dohcal) { // hcal hits
    nbytehcal = b_hcal->GetEntry(ievt);
    if(ievt<SCECOUNT) std::cout<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    ehcaltimecut=0.;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      long long int ihitchan=ahcalhit->cellID;
      Contributions zxzz=ahcalhit->truth;
      if(hcaltype==0) { // fiber
	int idet = (ihitchan) & 0xFF;
	int ilayer = (ihitchan >>8) & 0xFFF;  
	int itube = (ihitchan >>20) & 0xFFF;  
	int iair = (ihitchan >>32) & 0x3;  
	int itype = (ihitchan >>35) & 0x3;  
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
	//if(ievt<SCECOUNT) std::cout<<"   "<<ihitchan<<" " <<std::hex<<ihitchan<<std::dec<<" "<<idet<<" "<<iy<<" "<<ix<<" "<<ilayer<<" "<<islice<<std::endl;
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
      }// end of sampling
    }  // end loop over hcal hits
  } // end of if dohcal
} // getMeanPhot

void getStuff(map<string, int> mapecalslice,  map<string, int> mapsampcalslice, int gendet, int ievt, bool doecal, bool dohcal, int hcaltype, //input 
  bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge, CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits, //input
  float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber1,float &eesumfiber2,float &eesumabs,float &eesumPDh,float &eesumedge,//output
  float &necertotecal,float &nescinttotecal,float &necertothcal,float &nescinttothcal, float &timecut, float &eecaltimecut, float &ehcaltimecut,//output
  TH1F* eecaltime, TH1F* ehcaltime, int &nin){

  if(ievt<SCECOUNT) cout<<"getstuff phot ievt is "<<ievt<<endl;
  int nbyteecal(0), nbytehcal(0), nbyteedge(0);
  
  std::cout<<"entering getStuff # ecalhits="<<ecalhits->size()<<std::endl;
  //std::cout<<"eesum = "<<eesum<<" ehcaltimecut= "<<ehcaltimecut<<" eesumfiber1="<<eesumfiber1<<" eesumfiber2= "<<eesumfiber2<<" eesumabs= "<<eesumabs<<std::endl;
  
  if(doecal) { // ecal hists
    nbyteecal = b_ecal->GetEntry(ievt);
    if(ievt<SCECOUNT) cout<<" number of ecal hits is "<<ecalhits->size()<<endl;
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
      nin+=aecalhit->n_inelastic; // what is this
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
      if(islice==(*eii0).second)eesumair+=ae;
      if(islice==(*eii1).second)eesumPDe+=ae;
      if(islice==(*eii5).second)eesumPDe+=ae;
      if(islice==(*eii2).second)eesumcrystal+=ae;
      if(islice==(*eii4).second)eesumcrystal+=ae;

      if(gendet==1) {   // use photons as generated in optical material
        if((islice==(*eii2).second)||(islice==(*eii4).second)) {  // crystal
          necertotecal+=aecalhit->ncerenkov;
          nescinttotecal+=aecalhit->nscintillator;
        }
      }
      else if(gendet==2) {
        if( (islice==(*eii1).second)||(islice==(*eii5).second) ) { // either photo detector
          necertotecal+=aecalhit->ncerenkov;
          nescinttotecal+=aecalhit->nscintillator;
        }
      } 
      else if(gendet==3||gendet==4){
	 if(idet==5) {
          nescinttotecal+=aecalhit->energyDeposit;
          if(gendet==3) necertotecal+=aecalhit->edeprelativistic;
          if(gendet==4) necertotecal+=aecalhit->energyDeposit;
          //for(size_t j=0;j<zxzz.size(); j++) {
          //   if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
          //}
        }
      }
    }  // end loop over ecal hits
  }// end if doecal

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);
    // hcal hits
    if(ievt<SCECOUNT) cout<<" number of hcal hits is "<<hcalhits->size()<<endl;
    float ehcaltimecut=0.;
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
	if(hacheck/ah<0.99999) cout<<"missing contribs: hcal check contributions Ncontrib is "<<zxzz.size()<<" hackec is  "<<hacheck<<" ah is "<<ah<<" ratio "<<hacheck/ah<<endl;
      }
      eesum+=ah;
      long long int ihitchan=ahcalhit->cellID;
      if(hcaltype==0) { // fiber
	int idet = (ihitchan) & 0xFF;
	int ilayer = (ihitchan >>8) & 0xFFF;  
	int itube = (ihitchan >>20) & 0xFFF;  
	//int iair=0; int itype=0;
	int iair = (ihitchan >>32) & 0x3;  
	int itype = (ihitchan >>35) & 0x3; 
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
           //for(size_t j=0;j<zxzz.size(); j++) { 
	   // if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	   //}
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
	    //for(size_t j=0;j<zxzz.size(); j++) {
	    //  if((zxzz.at(j)).time<timecut) ehcaltimecut+=(zxzz.at(j)).deposit;
	    //}
	  }
	}
	if( islice==(*sii3).second ) eesumfiber1+=ah; // scint
	if( islice==(*sii6).second ) eesumfiber2+=ah;  //cer
	if( (islice==(*sii1).second)||(islice==(*sii8).second)||(islice==(*sii9).second) ) eesumabs+=ah;
	if(  (islice==(*sii2).second) || (islice==(*sii4).second) ||  (islice==(*sii5).second) || (islice==(*sii7).second)) eesumPDh+=ah;
      } //end of hcal sampling
    }  // end loop over hcal hits
  }//end of if dohcal
  if(doedge) {
    nbyteedge = b_edge->GetEntry(ievt);
    if(ievt<SCECOUNT) cout<<" number of edge hits is "<<edgehits->size()<<endl;
    for(size_t i=0;i<edgehits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aedgehit =edgehits->at(i);
      float ae=aedgehit->energyDeposit;
      eesum+=ae;
      eesumedge+=ae;
    }  // end loop over escaping hits
  } // end doedge
} // end of getStuff 

void getStuffDualCorr(map<string, int> mapecalslice, map<string, int> mapsampcalslice, int gendet, float kappaEcal, float kappaHcal, 
  float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, int hcaltype, TBranch* &b_ecal,TBranch* &b_hcal,
  CalHits* &ecalhits, CalHits* &hcalhits, float &EEcal, float &EHcal, float &timecut, float &eecaltimecut, float &ehcaltimecut){

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
      if(gendet==1) {   // use photons as generated in otical material
        if((islice==(*eii2).second)||(islice==(*eii4).second)) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( (islice==(*eii1).second)||(islice==(*eii5).second) ) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      } 
      else if(gendet==3||gendet==4){
	if(idet==5) {
          nescinttotecal+=aecalhit->energyDeposit;
	  if(gendet==3) necertotecal+=aecalhit->edeprelativistic;
	  if(gendet==4) necertotecal+=aecalhit->energyDeposit;
          for(size_t j=0;j<zxzz.size(); j++) {
	    if((zxzz.at(j)).time<timecut) eecaltimecut+=(zxzz.at(j)).deposit;
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
	int iair = (ihitchan >>32) & 0x3;  
	int itype = (ihitchan >>35) & 0x3;  
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
      }//end of hcaltype=fiber
      else {  // sampling
	int idet = (ihitchan) & 0x07;
	int iy = (ihitchan >>3) & 0xFFF;  
	int ix = (ihitchan >>15) & 0xFFF;  
	int ilayer = (ihitchan >>27) & 0xFFF;  
	int ibox2 = (ihitchan >> 39) & 0x03;
	int islice = (ihitchan >>41) & 0xF;  
        //if(ievt<SCECOUNT) std::cout<<"   "<<std::hex<<ihitchan<<std::dec<<" " <<idet<<" "<<ix<<" "<<iy<<" "<<islice<<" "<<ilayer<<" "<<ahcalhit->energyDeposit<<" "<<ahcalhit->nscintillator<<" "<<ahcalhit->ncerenkov<<std::endl;
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
      } // end of if sampling 
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
} // end of getStuffDualCorr
