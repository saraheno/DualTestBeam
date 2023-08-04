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
const int nsliceecal = 4;
std::string nameecalslice[nsliceecal] = {"air","PD1","crystal","PD2"};
int SCECOUNT=5;


// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
//const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
  typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
  typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;


void SCEDraw1 (TCanvas* canv, const char* name, TH1F* h1, const char* outfile);
void SCEDrawp (TCanvas* canv, const char* name, TProfile* h1, const char* outfile);
void SCEDraw1_2D (TCanvas* canv, const char* name, TH2F* h1, const char* outfile,float eohS,float eohC);
void SCEDraw2 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile);
void SCEDraw3 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile);


void getStuff(map<string, int> mapecalslice, int gendet, int ievt, bool doecal, bool dohcal, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber,float &eesumabs,float &eesumPDh,float &eesumedge,int &necertotecal,int &nescinttotecal,int &necertothcal,int &nescinttothcal);


void getStuffDualCorr(map<string, int> mapecalslice, int gendet, float kappaecal, float kappahcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, TBranch* &b_ecal,TBranch* &b_hcal, 
		      CalHits* &ecalhits, CalHits* &hcalhits,
float &EEcal, float &EHcal);

void getMeanPhot(map<string, int> mapecalslice,  int gendet, int ievt, bool doecal, bool dohcal,TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits, 
float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal);

double EOVERH(double b, double m);




void crystalana(int num_evtsmax, const char* einputfilename, const char* piinputfilename, const float beamEE, bool doecal, bool dohcal, bool doedge, int gendet, const char* outputfilename) {


map<string, int> mapecalslice; 
mapecalslice["air"]=0;
mapecalslice["PD1"]=1;
mapecalslice["crystal"]=2;
mapecalslice["PD2"]=3;





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

  TH1F *ehcHcalE = new TH1F("ehcHcalE","sum fiber hcal energy / beam E",100,0.,1.5);
  TH1F *phcHcalE = new TH1F("phcHcalE","sum fiber hcal energy / beam E",100,0.,1.5);

  TH1F *ehcEdgeE = new TH1F("ehcEdgeE","sum escaping / beam E",100,0.,1.5);
  TH1F *phcEdgeE = new TH1F("phcEdgeE","sum escaping / beam E",100,0.,1.5);

  TH1F *ehcEcalncer = new TH1F("ehcEcalncer","total number of ecal cerenkov",  500,0.,1.5);
  TH1F *phcEcalncer = new TH1F("phcEcalncer","total number of ecal cerenkov", 500,0.,1.5);

  TH1F *ehcHcalncer = new TH1F("ehcHcalncer","total number of hcal cerenkov",   500,0.,1.5);
  TH1F *phcHcalncer = new TH1F("phcHcalncer","total number of hcal cerenkov",  500,0.,1.5);


  TH1F *ehcEcalnscint = new TH1F("ehcEcalnscint","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalnscint = new TH1F("phcEcalnscint","total number of ecal scitillation", 500,0.,1.5);

  TH1F *ehcHcalnscint = new TH1F("ehcHcalnscint","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalnscint = new TH1F("phcHcalnscint","total number of hcal scitillation", 500,0.,1.5);

  TH1F *ehcEcalcorr = new TH1F("ehcEcalcorr","total number of ecal scintillation", 500,0.,1.5);
  TH1F *phcEcalcorr = new TH1F("phcEcalcorr","total number of ecal scitillation", 500,0.,1.5);

  TH1F *ehcHcalcorr = new TH1F("ehcHcalcorr","total number of hcal scintillation", 500,0.,1.5);
  TH1F *phcHcalcorr = new TH1F("phcHcalcorr","total number of hcal scitillation", 500,0.,1.5);


  TH2F *ehcEcalNsNc = new TH2F("ehcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcEcalNsNc = new TH2F("phcEcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *ehcHcalNsNc = new TH2F("ehcHcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);
  TH2F *phcHcalNsNc = new TH2F("phcHcalNsNc","ecal ncer versus nscint",500,0.,1.5,500,0.,1.5);

  TH2F *ehecal2d = new TH2F("ehecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *phecal2d = new TH2F("phecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);

  TH2F *ehhcal2d = new TH2F("ehhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *phhcal2d = new TH2F("phhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);

  TH1F *ehaphcal = new TH1F("ehaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);
  TH1F *phaphcal = new TH1F("phaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);

  TH1F *eheest = new TH1F("eheest","ratio estimated to true energy",500,0.,1.);
  TH1F *pheest = new TH1F("pheest","ratio estimated to true energy",500,0.,1.);

  TH1F *ehetrue = new TH1F("ehetrue","ratio deposited to incident energy",500,0.,1.);
  TH1F *phetrue = new TH1F("phetrue","ratio deposited to incident energy",500,0.,1.);

  TH1F *ehnecalcon = new TH1F("ehnecalcon","number contribs to ecal hit",1010,-10.,1000.);
  TH1F *phnecalcon = new TH1F("phnecalcon","number contribs to ecal hit",1010,-10.,1000.);

  TH2F *ehzvst = new TH2F("ehzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.); 
  TH2F *phzvst = new TH2F("phzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.); 




  //****************************************************************************************************************************
  // process electrons

  TFile* ef = TFile::Open(einputfilename);
  TTree* et = (TTree*)ef->Get("EVENT;1");

  b_mc= et->GetBranch("MCParticles");
  if(doecal) b_ecal = et->GetBranch("DRCNoSegment");
  if(dohcal) b_hcal = et->GetBranch("DRFNoSegment");
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");





  ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<"num_evt for electron file is  "<<num_evt<<std::endl;
  
  // loop over events 

  float meanscinEcal(0),meanscinHcal(0),meancerEcal(0),meancerHcal(0);
  float hcalSampf;

  
  if(num_evt>0) {  

    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);

    // first pass through file

    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number first pass is "<<ievt<<std::endl;
      getMeanPhot(mapecalslice, gendet, ievt, doecal, dohcal, b_ecal,b_hcal, ecalhits, hcalhits, meanscinEcal, meanscinHcal, meancerEcal, meancerHcal);
    }
    meanscinEcal=meanscinEcal/num_evt;
    meanscinHcal=meanscinHcal/num_evt;
    meancerEcal=meancerEcal/num_evt;
    meancerHcal=meancerHcal/num_evt;
    std::cout<<"mean scint ecal is "<<meanscinEcal<<std::endl;
    std::cout<<"mean scint hcal is "<<meanscinHcal<<std::endl;
    std::cout<<"mean cer ecal is "<<meancerEcal<<std::endl;
    std::cout<<"mean cer hcal is "<<meancerHcal<<std::endl;
    hcalSampf=0.01;  // THIS IS WRONG!!!!!!!!
    std::cout<<"hcal sampling fraction is "<<hcalSampf<<std::endl;


    // second pass through file
    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number second is "<<ievt<<std::endl;


      float eesum=0.;
      float eesumair=0;
      float eesumcrystal=0;
      float eesumPDe=0;
      float eesumfiber=0;
      float eesumabs=0;
      float eesumPDh=0;
      float eesumedge=0;
      int necertotecal=0;
      int nescinttotecal=0;
      int necertothcal=0;
      int nescinttothcal=0;

      getStuff(mapecalslice,  gendet, ievt, doecal, dohcal, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,eesum,eesumair,eesumcrystal,eesumPDe,eesumfiber,eesumabs,eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal);

    
      ehcEcalE->Fill(eesumcrystal/beamE);
      ehcHcalE->Fill(eesumfiber*hcalSampf/beamE);
      ehcEdgeE->Fill(eesumedge/beamE);
      ehcEcalncer->Fill(necertotecal/meancerEcal);
      ehcHcalncer->Fill(necertothcal/meancerHcal);
      ehcEcalnscint->Fill(nescinttotecal/meanscinEcal);
      ehcHcalnscint->Fill(nescinttothcal/meanscinHcal);
      if(eesumfiber>0) ehaphcal->Fill(eesumfiber/(eesumabs+eesumfiber));
      eheest->Fill((eesumcrystal+eesumfiber*hcalSampf)/beamE);
      ehcEcalNsNc->Fill(nescinttotecal/meanscinEcal,necertotecal/meancerEcal);  
      ehcHcalNsNc->Fill(nescinttothcal/meanscinHcal,necertothcal/meancerHcal);  


      float eachecks=eesumair+eesumPDe+eesumcrystal+eesumfiber+eesumabs+eesumPDh+eesumedge;
      ehetrue->Fill(eachecks/beamE);


      std::cout<<std::endl<<std::endl<<"total energy deposit "<<eesum/1000.<<std::endl;
      std::cout<<"       in air "<<eesumair/1000.<<std::endl;
      std::cout<<"       in photodetector ecal "<<eesumPDe/1000.<<std::endl;
      std::cout<<"       in crystal "<<eesumcrystal/1000.<<std::endl;
      std::cout<<"       in fiber "<<eesumfiber/1000.<<std::endl;
      std::cout<<"       in absorber "<<eesumabs/1000.<<std::endl;
      std::cout<<"       in photodetect hcal "<<eesumPDh/1000.<<std::endl;
      std::cout<<"       escaping detector "<<eesumedge/1000.<<std::endl;

      std::cout<<"       sum individual "<<eachecks/1000.<<std::endl;
      std::cout<<"   incident energy "<<beamE/1000.<<std::endl;
      std::cout<<"   ratio to incident energy "<<eachecks/beamE<<std::endl;

      std::cout<<"total number of cherenkov ecal is "<<necertotecal<<std::endl;
      std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl;
      std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
      std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl<<std::endl;









    }  //end loop over events
  }  // end if no events

  //  float amean = hceest->GetMean();

  ef->Close();


  //****************************************************************************************************************************
  // process pions

  TFile* pif = TFile::Open(piinputfilename);
  TTree* pit = (TTree*)pif->Get("EVENT;1");

  b_mc= pit->GetBranch("MCParticles");
  if(doecal) b_ecal = et->GetBranch("DRCNoSegment");
  if(dohcal) b_hcal = et->GetBranch("DRFNoSegment");
  if(doedge) b_edge = et->GetBranch("EdgeDetNoSegment");


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


      float pesum=0.;
      float pesumair=0;
      float pesumcrystal=0;
      float pesumPDe=0;
      float pesumfiber=0;
      float pesumabs=0;
      float pesumPDh=0;
      float pesumedge=0;
      int npcertotecal=0;
      int npscinttotecal=0;
      int npcertothcal=0;
      int npscinttothcal=0;

      getStuff(mapecalslice,  gendet, ievt, doecal, dohcal, doedge, b_ecal,b_hcal,b_edge,ecalhits,hcalhits,edgehits,pesum,pesumair,pesumcrystal,pesumPDe,pesumfiber,pesumabs,pesumPDh,pesumedge,npcertotecal,npscinttotecal,npcertothcal,npscinttothcal);

    
      phcEcalE->Fill(pesumcrystal/beamE);
      phcHcalE->Fill(pesumfiber*hcalSampf/beamE);
      phcEdgeE->Fill(pesumedge/beamE);
      phcEcalncer->Fill(npcertotecal/meancerEcal);
      phcHcalncer->Fill(npcertothcal/meancerHcal);
      phcEcalnscint->Fill(npscinttotecal/meanscinEcal);
      phcHcalnscint->Fill(npscinttothcal/meanscinHcal);
      if(pesumfiber>0) phaphcal->Fill(pesumfiber/(pesumabs+pesumfiber));
      pheest->Fill((pesumcrystal+pesumfiber*hcalSampf)/beamE);
      phcEcalNsNc->Fill(npscinttotecal/meanscinEcal,npcertotecal/meancerEcal);  
      phcHcalNsNc->Fill(npscinttothcal/meanscinHcal,npcertothcal/meancerHcal);  




      float pachecks=pesumair+pesumPDe+pesumcrystal+pesumfiber+pesumabs+pesumPDh+pesumedge;
      phetrue->Fill(pachecks/beamE);


      std::cout<<std::endl<<std::endl<<"total energy deposit "<<pesum/1000.<<std::endl;
      std::cout<<"       in air "<<pesumair/1000.<<std::endl;
      std::cout<<"       in photodetector ecal "<<pesumPDe/1000.<<std::endl;
      std::cout<<"       in crystal "<<pesumcrystal/1000.<<std::endl;
      std::cout<<"       in fiber "<<pesumfiber/1000.<<std::endl;
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









    }  //end loop over events
  }  // end if no events

  //  float amean = hceest->GetMean();





  //** fits
  TF1 *gEcale = new TF1("gEcale","pol1",0.,1.);
  TF1 *gHcale = new TF1("gHcale","pol1",0.,1.);
  TF1 *gEcalp = new TF1("gEcalp","pol1",0.,1.);
  TF1 *gHcalp = new TF1("gHcalp","pol1",0.,1.);


      // fit to get e/h
  TProfile* ehcEcalNsNc_pfx = ehcEcalNsNc->ProfileX();
  ehcEcalNsNc_pfx->Fit("gEcale","W");
  TProfile* ehcHcalNsNc_pfx = ehcHcalNsNc->ProfileX();
  ehcHcalNsNc_pfx->Fit("gHcale","W");


  TProfile* phcEcalNsNc_pfx = phcEcalNsNc->ProfileX();
  phcEcalNsNc_pfx->Fit("gEcalp","W");
  float b1Ecal=gEcalp->GetParameter(0);
  float m1Ecal=gEcalp->GetParameter(1);
  std::cout<<"for ecal b m are "<<b1Ecal<<" "<<m1Ecal<<std::endl;
  //double eohSEcal = EOVERH(b1Ecal,m1Ecal);
  //double eohCEcal = 1/(1-(1-(1/eohSEcal))*m1Ecal);
  //std::cout<<"Ecal eohS eohC are "<<eohSEcal<<" "<<eohCEcal<<std::endl;
  double kappaEcal = 1+(b1Ecal/m1Ecal);
  std::cout<<" kappa ecal is "<<kappaEcal<<std::endl;


  TProfile* phcHcalNsNc_pfx = phcHcalNsNc->ProfileX();
  phcEcalNsNc_pfx->Fit("gHcalp","W");
  float b1Hcal=gHcalp->GetParameter(0);
  float m1Hcal=gHcalp->GetParameter(1);
  std::cout<<"for hcal b m are "<<b1Hcal<<" "<<m1Hcal<<std::endl;
  //double eohSHcal = EOVERH(b1Ecal,m1Ecal);
  //double eohCHcal = 1/(1-(1-(1/eohSHcal))*m1Hcal);
  //std::cout<<"Hcal eohS eohC are "<<eohSHcal<<" "<<eohCHcal<<std::endl;
  double kappaHcal = 1+(b1Hcal/m1Hcal);
  std::cout<<" kappa hcal is "<<kappaHcal<<std::endl;

      
  // no calculate with dual readout correction  




  if(num_evt>0) {  
    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);


    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number pion is "<<ievt<<std::endl;

      float EcorEcal(0),EcorHcal(0);

      std::cout<<" meanscinEcal is "<<meanscinEcal<<std::endl;
      getStuffDualCorr(mapecalslice, gendet, kappaEcal, kappaHcal, meanscinEcal, meancerEcal, meanscinHcal, meancerHcal,  ievt,doecal,dohcal, b_ecal,b_hcal, ecalhits,hcalhits, EcorEcal, EcorHcal);

      phcEcalcorr->Fill(EcorEcal);
      phcHcalcorr->Fill(EcorHcal);



    }  //end loop over events
  }  // end if no events



  // close pion file
  pif->Close();


  //***********************************************************************************************************************

  TCanvas* c1;
  SCEDraw2(c1,"c1",ehetrue,phetrue,"junk1.png");

  TCanvas* c2;
  SCEDraw2(c2,"c2",ehcEcalE,phcEcalE,"junk2.png");


  TCanvas* c3;
  SCEDraw2(c3,"c3",ehcEcalncer,ehcEcalnscint,"junk3.png");

  TCanvas* c4;
  SCEDraw3(c4,"c4",phcEcalncer,phcEcalnscint,phcEcalcorr,"junk4.png");


  TCanvas* c5;
  SCEDraw1_2D(c5,"c5",ehcEcalNsNc,"junk5.png",0.,0.);



  TCanvas* c6;
  SCEDraw1_2D(c6,"c6",phcEcalNsNc,"junk6.png",-b1Ecal/m1Ecal,0.);

  std::cout<<"haha0"<<std::endl;
  //TCanvas* c7;
  //SCEDrawp(c7,"c7",phcEcalNsNc_pfx,"junk7.png");

  std::cout<<"haha1"<<std::endl;

  //***********************************************************************************************************

  TFile * out = new TFile(outputfilename,"RECREATE");

  ehcEcalE->Write();
  phcEcalE->Write();

  ehcHcalE->Write();
  phcHcalE->Write();

  ehcEdgeE->Write();
  phcEdgeE->Write();

  std::cout<<"haha2"<<std::endl;;

   ehcEcalncer->Write();
   phcEcalncer->Write();

   ehcHcalncer->Write();
   phcHcalncer->Write();

   ehcEcalcorr->Write();
   phcEcalcorr->Write();

   ehcHcalcorr->Write();
   phcHcalcorr->Write();


  ehecal2d->Write();
  phecal2d->Write();


   ehcEcalnscint->Write();
   phcEcalnscint->Write();


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

   //ehcEcalNsNc_pfx->Write();
   //ehcHcalNsNc_pfx->Write();
   //phcEcalNsNc_pfx->Write();
   //phcHcalNsNc_pfx->Write();



  out->Close();

}

void SCEDraw1 (TCanvas* canv,  const char* name,TH1F* h1, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(111111);  
  h1->Draw("HIST");



  canv->Print(outfile,".png");

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

  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(111111);  
  h1->SetMarkerSize(20);
  h1->SetMarkerStyle(4);
  h1->SetMarkerColor(3);
  gStyle->SetOptFit();
  h1->Draw("*");



  canv->Print(outfile,".png");

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

  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(111111);  
  h1->Draw("");

  TLine line = TLine(eohS,eohC,1.,1.);
  line.SetLineColor(kYellow);
  line.SetLineWidth(2);
  line.Draw("same");


  canv->Print(outfile,".png");

  return;
}


void SCEDraw2 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);


  float max = std::max(h1->GetMaximum(),h2->GetMaximum());
  h1->SetMaximum(max*1.3);



  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(111111);  
  h1->Draw("HIST");



  h2->SetLineColor(2);
  h2->SetLineWidth(3);
  h2->SetStats(111111);  
  h2->Draw("HIST same");


  canv->Print(outfile,".png");

  return;
}


void SCEDraw3 (TCanvas* canv,  const char* name, TH1F* h1, TH1F* h2, TH1F* h3, const char* outfile) {

  canv= new TCanvas(name,name,200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);


  float max = std::max(h1->GetMaximum(),h2->GetMaximum());
  h1->SetMaximum(max*1.3);



  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(111111);  
  h1->Draw("HIST");



  h2->SetLineColor(2);
  h2->SetLineWidth(3);
  h2->SetStats(111111);  
  h2->Draw("HIST same");



  h3->SetLineColor(4);
  h3->SetLineWidth(3);
  h3->SetStats(111111);  
  h3->Draw("HIST same");


  canv->Print(outfile,".png");

  return;
}


void getMeanPhot(map<string, int> mapecalslice,  int gendet, int ievt, bool doecal, bool dohcal, TBranch* &b_ecal,TBranch* &b_hcal,
	      CalHits* &ecalhits, CalHits* &hcalhits, 
float &meanscinEcal, float &meanscinHcal, float &meancerEcal, float &meancerHcal){
  int nbyteecal, nbytehcal, nbyteedge;


  

  if(doecal) {
    if(ievt<SCECOUNT) std::cout<<"getMean phot ievt is "<<ievt<<std::endl;

    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      int ihitchan=aecalhit->cellID;
      int idet = (ihitchan) & 0x07;
      int ix = (ihitchan >>3) & 0x3F ;  // is this right?
      if(ix>32) ix=ix-64;
      int iy =(ihitchan >>10) & 0x3F ; // is this right?
      if(iy>32) iy=iy-64;
      int  islice = (ihitchan >>17) & 0x07;
      int  ilayer = (ihitchan>> 20) & 0x07;
      
      map<string,int>::iterator ii1 = mapecalslice.find("crystal");
      map<string,int>::iterator ii2 = mapecalslice.find("PD1");
      map<string,int>::iterator ii3 = mapecalslice.find("PD2");

      if(gendet==1) {   // use photons as generated in otical material
	if(islice==(*ii1).second) {
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( (islice==(*ii2).second)||(islice==(*ii3).second) ) {
	  meancerEcal+=aecalhit->ncerenkov;
	  meanscinEcal+=aecalhit->nscintillator;
	}
	else{   // use photons detected in photodetectors  
	  std::cout<<" illegal gendet"<<std::endl;
	  return;
	}



      }  // end loop over ecal hits
    }
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);
    

      // hcal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      int ihitchan=ahcalhit->cellID;
      int idet = (ihitchan) & 0x07;
      int ix = (ihitchan >>3) & 0xFF;   // is this right?
      if(ix>128) ix=ix-256;
      int iy =(ihitchan >>11) & 0xFF;   // is this right?
      if(iy>128) iy=iy-256;
      int ifiber  =(ihitchan >>21) & 0x03;
      int iabs=(ihitchan >>23) & 0x03;
      int iphdet=(ihitchan >>25) & 0x03;

      
      if(gendet==1) {   // use photons as generated in otical material
	if(ifiber==1) {
	  meancerHcal+=ahcalhit->ncerenkov;
	  meanscinHcal+=ahcalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( iabs==1) {
	  meancerHcal+=ahcalhit->ncerenkov;
	  meanscinHcal+=ahcalhit->nscintillator;
	}

      }
      else {   // use photons detected in photodetectors  
	std::cout<<" illegal gendet"<<std::endl;
	return;
      }

    }  // end loop over hcal hits
  }




}





void getStuff(map<string, int> mapecalslice,  int gendet, int ievt, bool doecal, bool dohcal, bool doedge,TBranch* &b_ecal,TBranch* &b_hcal,TBranch*  &b_edge,
	      CalHits* &ecalhits, CalHits* &hcalhits, CalHits* &edgehits,
float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber,float &eesumabs,float &eesumPDh,float &eesumedge,int &necertotecal,int &nescinttotecal,int &necertothcal,int &nescinttothcal){


    if(ievt<SCECOUNT) std::cout<<"getstuff phot ievt is "<<ievt<<std::endl;

  int nbyteecal, nbytehcal, nbyteedge;



      if(doecal) {
	nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
	if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
      for(size_t i=0;i<ecalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);

	necertotecal+=aecalhit->ncerenkov;
	nescinttotecal+=aecalhit->nscintillator;
	if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
        int ihitchan=aecalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0x3F ;  // is this right?
	if(ix>32) ix=ix-64;
	int iy =(ihitchan >>10) & 0x3F ; // is this right?
	if(iy>32) iy=iy-64;
        int  islice = (ihitchan >>17) & 0x07;
        int  ilayer = (ihitchan>> 20) & 0x07;


      map<string,int>::iterator ii0 = mapecalslice.find("air");
      map<string,int>::iterator ii1 = mapecalslice.find("crystal");
      map<string,int>::iterator ii2 = mapecalslice.find("PD1");
      map<string,int>::iterator ii3 = mapecalslice.find("PD2");


	
	if((ilayer!=0)&&(ilayer!=1)) std::cout<<"danger danger will robinson ilayer not zero"<<std::endl;
	if(islice>nsliceecal) {
	  std::cout<<"  danger danger will robinson islice nsliceecal are "<<islice<<" "<<nsliceecal<<std::endl;
	  std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<std::endl;
	} else {
	  if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;
	}

	float ae=aecalhit->energyDeposit;
	eesum+=ae;
	if(islice==(*ii0).second)eesumair+=ae;
	if(islice==(*ii1).second)eesumPDe+=ae;
	if(islice==(*ii2).second)eesumcrystal+=ae;
	if(islice==(*ii3).second)eesumPDe+=ae;

	//ehchan->Fill(aecalhit->cellID);
	//ehecal2d->Fill(ix,iy,aecalhit->energyDeposit);


      }  // end loop over ecal hits
      }


      if(dohcal) {
	nbytehcal = b_hcal->GetEntry(ievt);


      // hcal hits
	if(ievt<SCECOUNT) std::cout<<std::endl<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
      for(size_t i=0;i<hcalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);

	necertothcal+=ahcalhit->ncerenkov;
	nescinttothcal+=ahcalhit->nscintillator;
	if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<ahcalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<ahcalhit->energyDeposit<<","<<ahcalhit->ncerenkov<<","<<ahcalhit->nscintillator<<")"<<std::endl;
        int ihitchan=ahcalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0xFF;   // is this right?
	if(ix>128) ix=ix-256;
	int iy =(ihitchan >>11) & 0xFF;   // is this right?
	if(iy>128) iy=iy-256;
	int ifiber  =(ihitchan >>21) & 0x03;
	int iabs=(ihitchan >>23) & 0x03;
	int iphdet=(ihitchan >>25) & 0x03;
	if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"  idet,ix,iy is ("<<idet<<","<<ix<<","<<iy<<")"<<std::endl;
	if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"  ifiber,iabs,iphdet is ("<<ifiber<<","<<iabs<<","<<iphdet<<")"<<std::endl;


	float ah=ahcalhit->energyDeposit;
	eesum+=ah;

	if(ifiber==1) eesumfiber+=ah;
	if(iabs==1) eesumabs+=ah;
	if(iphdet==1) eesumPDh+=ah;

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




}


void getStuffDualCorr(map<string, int> mapecalslice, int gendet, float kappaEcal, float kappaHcal, float meanscinEcal, float meancerEcal, float meanscinHcal, float meancerHcal, int  ievt,bool doecal,bool dohcal, TBranch* &b_ecal,TBranch* &b_hcal, 
	      CalHits* &ecalhits, CalHits* &hcalhits, 
float &EEcal, float &EHcal)
{
  int necertotecal(0),nescinttotecal(0),necertothcal(0),nescinttothcal(0);
  int nbyteecal, nbytehcal, nbyteedge;

  std::cout<<" getstuff meanscinEcal is "<<meanscinEcal<<std::endl;
  if(ievt<SCECOUNT) std::cout<<"getstuffdualcorr phot ievt is "<<ievt<<std::endl;


  if(doecal) {
    nbyteecal = b_ecal->GetEntry(ievt);

      // ecal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of ecal hits is "<<ecalhits->size()<<std::endl;
    for(size_t i=0;i<ecalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
      int ihitchan=aecalhit->cellID;
      int idet = (ihitchan) & 0x07;
      int ix = (ihitchan >>3) & 0x3F ;  // is this right?
      if(ix>32) ix=ix-64;
      int iy =(ihitchan >>10) & 0x3F ; // is this right?
      if(iy>32) iy=iy-64;
      int  islice = (ihitchan >>17) & 0x07;
      int  ilayer = (ihitchan>> 20) & 0x07;
      

	  if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"getstuffdualcorr  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;


      map<string,int>::iterator ii1 = mapecalslice.find("crystal");
      map<string,int>::iterator ii2 = mapecalslice.find("PD1");
      map<string,int>::iterator ii3 = mapecalslice.find("PD2");

      if(gendet==1) {   // use photons as generated in otical material
	if(islice==(*ii1).second) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( (islice==(*ii2).second)||(islice==(*ii3).second) ) {
	  necertotecal+=aecalhit->ncerenkov;
	  nescinttotecal+=aecalhit->nscintillator;
	}
	else{   // use photons detected in photodetectors  
	  std::cout<<" illegal gendet"<<std::endl;
	  return;
	}



      }  // end loop over ecal hits
    }
  }

  if(dohcal) {
    nbytehcal = b_hcal->GetEntry(ievt);
    

      // hcal hits
    if(ievt<SCECOUNT) std::cout<<std::endl<<" number of hcal hits is "<<hcalhits->size()<<std::endl;
    for(size_t i=0;i<hcalhits->size(); ++i) {
      CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);
      int ihitchan=ahcalhit->cellID;
      int idet = (ihitchan) & 0x07;
      int ix = (ihitchan >>3) & 0xFF;   // is this right?
      if(ix>128) ix=ix-256;
      int iy =(ihitchan >>11) & 0xFF;   // is this right?
      if(iy>128) iy=iy-256;
      int ifiber  =(ihitchan >>21) & 0x03;
      int iabs=(ihitchan >>23) & 0x03;
      int iphdet=(ihitchan >>25) & 0x03;

      
      if(gendet==1) {   // use photons as generated in otical material
	if(ifiber==1) {
	  necertothcal+=ahcalhit->ncerenkov;
	  nescinttothcal+=ahcalhit->nscintillator;
	}
      }
      else if(gendet==2) {
	if( iabs==1) {
	  necertothcal+=ahcalhit->ncerenkov;
	  nescinttothcal+=ahcalhit->nscintillator;
	}

      }
      else {   // use photons detected in photodetectors  
	std::cout<<" illegal gendet"<<std::endl;
	return;
      }

    }  // end loop over hcal hits
  }




  std::cout<<"  getstuffdual cer scint count "<<necertotecal<<" "<<nescinttotecal<<std::endl;
  float anecertotecal=necertotecal/meancerEcal;
  float anescinttotecal=nescinttotecal/meanscinEcal;
  std::cout<<"  getstuffdual cer scint count "<<anecertotecal<<" "<<anescinttotecal<<std::endl;
  EEcal=(anescinttotecal-kappaEcal*anecertotecal)/(1-kappaEcal);


  float anecertothcal=necertothcal/meancerHcal;
  float anescinttothcal=nescinttothcal/meanscinHcal;
  EHcal=(anescinttothcal-kappaHcal*anecertothcal)/(1-kappaHcal);

  
  std::cout<<"getstuffdual outputing ecal hcal "<<EEcal<<" "<<EHcal<<std::endl;






}


double EOVERH(double b,double m) {
  std::cout<<"in EOVERH"<<std::endl;
 std:cout<<" in;puts are b m "<<b<<" "<<m<<std::endl;

  double eoh=1.;
  double A = m;
  double B=-2.*m;
  double C= m-1-(b/m);
  double rt=B*B-4*A*C;
  double s1(0),s2(0);
  if(rt>0) {
    rt=sqrt(rt);
    s1=(-B+rt)/2/A;
    s2=(-B-rt)/2/A;
  }
  else{std::cout<<"bad quad A B C rt are "<<A<<" "<<B<<" "<<C<<" "<<rt<<std::endl;}
  std::cout<<"s1 s2 are "<<s1<<" "<<s2<<std::endl;
  
  eoh=1/s1;
  std::cout<<"choosing larger one so invert is less than one "<<eoh<<std::endl;

  return eoh;
}
