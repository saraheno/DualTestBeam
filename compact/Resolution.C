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
#include <algorithm>


// information about hit channel IT numbers
//const int nchan = 4;
//const int ichan[nchan] = {64,73,74,75};  // channel 74 is the crystal, 73 and 75 the two kill media

const int nsliceecal = 4;
std::string nameecalslice[nsliceecal] = {"air","PD1","crystal","PD2"};
int SCECOUNT=1;


// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
//const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry
float hcalcalib=1./0.036;

void SCEDraw1 (TCanvas* canv, TH1F* h1, const char* outfile);
void SCEDraw2 (TCanvas* canv, TH1F* h1, TH1F* h2, const char* outfile);
void getStuff(int ievt, bool doecal, bool dohcal, bool doedge,TBranch* b_ecal,TBranch* b_hcal,TBranch*  b_edge,float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber,float &eesumabs,float &eesumPDh,float &eesumedge,int &necertotecal,int &nescinttotecal,int &necertothcal,int &nescinttothcal);



  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
  typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
  typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;



void crystalana(int num_evtsmax, const char* einputfilename, const char* piinputfilename, const float beamEE, bool doecal, bool dohcal, bool doedge, const char* outputfilename) {

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
  TH1F *pihchan = new TH1F("pihchan","channel ID number",1028,0.,1028);

  TH1F *ehcEcalE = new TH1F("ehcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);
  TH1F *pihcEcalE = new TH1F("pihcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);

  TH1F *ehcHcalE = new TH1F("ehcHcalE","sum fiber hcal energy / beam E",100,0.,1.5);
  TH1F *pihcHcalE = new TH1F("pihcHcalE","sum fiber hcal energy / beam E",100,0.,1.5);

  TH1F *ehcEdgeE = new TH1F("ehcEdgeE","sum escaping / beam E",100,0.,1.5);
  TH1F *pihcEdgeE = new TH1F("pihcEdgeE","sum escaping / beam E",100,0.,1.5);

  TH1F *ehcEcalncer = new TH1F("ehcEcalncer","total number of cerenkov",100000,0.,100000);
  TH1F *pihcEcalncer = new TH1F("pihcEcalncer","total number of cerenkov",100000,0.,100000);

  TH1F *ehcHcalncer = new TH1F("ehcHcalncer","total number of cerenkov",100000,0.,100000);
  TH1F *pihcHcalncer = new TH1F("pihcHcalncer","total number of cerenkov",100000,0.,100000);

  TH2F *ehecal2d = new TH2F("ehecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *pihecal2d = new TH2F("pihecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);

  TH2F *ehhcal2d = new TH2F("ehhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *pihhcal2d = new TH2F("pihhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);

  TH1F *ehaphcal = new TH1F("ehaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);
  TH1F *pihaphcal = new TH1F("pihaphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);

  TH1F *eheest = new TH1F("eheest","ratio estimated to true energy",500,0.,1.);
  TH1F *piheest = new TH1F("piheest","ratio estimated to true energy",500,0.,1.);

  TH1F *ehetrue = new TH1F("ehetrue","ratio deposited to incident energy",500,0.,1.);
  TH1F *pihetrue = new TH1F("pihetrue","ratio deposited to incident energy",500,0.,1.);

  TH1F *ehnecalcon = new TH1F("ehnecalcon","number contribs to ecal hit",1010,-10.,1000.);
  TH1F *pihnecalcon = new TH1F("pihnecalcon","number contribs to ecal hit",1010,-10.,1000.);

  TH2F *ehzvst = new TH2F("ehzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.); 
  TH2F *pihzvst = new TH2F("pihzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.); 


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
  
  if(num_evt>0) {  



    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number is "<<ievt<<std::endl;

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


      getStuff(ievt, doecal, dohcal, doedge, b_ecal,b_hcal,b_edge,eesum,eesumair,eesumcrystal,eesumPDe,eesumfiber,eesumabs,eesumPDh,eesumedge,necertotecal,nescinttotecal,necertothcal,nescinttothcal);

    
      ehcEcalE->Fill(eesumcrystal/beamE);
      ehcHcalE->Fill(eesumfiber*hcalcalib/beamE);
      ehcEdgeE->Fill(eesumedge/beamE);
      ehcEcalncer->Fill(necertotecal);
      ehcHcalncer->Fill(necertothcal);
      if(eesumfiber>0) ehaphcal->Fill(eesumfiber/(eesumabs+eesumfiber));
      eheest->Fill((eesumcrystal+eesumfiber*hcalcalib)/beamE);



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
      std::cout<<"total number of scintillator ecal is "<<nescinttotecal<<std::endl<<std::endl;
      std::cout<<"total number of cherenkov hcal is "<<necertothcal<<std::endl;
      std::cout<<"total number of scintillator hcal is "<<nescinttothcal<<std::endl<<std::endl;









    }  //end loop over events
  }  // end if no events

  //  float amean = hceest->GetMean();

  ef->Close();

  TCanvas* c1;
  SCEDraw1(c1,ehetrue,"junk1.png");

  TCanvas* c2;
  SCEDraw2(c2,ehcEcalE,ehcHcalE,"junk2.png");


  //***********************************************************************************************************

  TFile * out = new TFile(outputfilename,"RECREATE");

  ehcEcalE->Write();
  pihcEcalE->Write();

  ehcHcalE->Write();
  pihcHcalE->Write();

  ehcEdgeE->Write();
  pihcEdgeE->Write();

   ehcEcalncer->Write();
   pihcEcalncer->Write();

  ehecal2d->Write();
  pihecal2d->Write();

  ehhcal2d->Write();
  pihhcal2d->Write();

  ehaphcal->Write();
  pihaphcal->Write();

  eheest->Write();
  piheest->Write();


  ehetrue->Write();
  pihetrue->Write();

  ehnecalcon->Write();
  pihnecalcon->Write();

  ehzvst->Write();
  pihzvst->Write();

   ehchan->Write();
   pihchan->Write();



  out->Close();

}

void SCEDraw1 (TCanvas* canv, TH1F* h1, const char* outfile) {

  canv= new TCanvas("Canvas1","Canvas1",200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(0);  
  h1->Draw("HIST");



  canv->Print(outfile,".png");

  return;
}
void SCEDraw2 (TCanvas* canv, TH1F* h1, TH1F* h2, const char* outfile) {

  canv= new TCanvas("Canvas2","Canvas2",200,10,700,500);


  //canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);

  h1->SetLineColor(3);
  h1->SetLineWidth(3);
  h1->SetStats(0);  
  h1->Draw("HIST");



  h2->SetLineColor(2);
  h2->SetLineWidth(3);
  h2->SetStats(0);  
  h2->Draw("HIST same");


  canv->Print(outfile,".png");

  return;
}


void getStuff(int ievt, bool doecal, bool dohcal, bool doedge,TBranch* b_ecal,TBranch* b_hcal,TBranch*  b_edge,float  &eesum,float &eesumair,float &eesumcrystal,float &eesumPDe,float &eesumfiber,float &eesumabs,float &eesumPDh,float &eesumedge,int &necertotecal,int &nescinttotecal,int &necertothcal,int &nescinttothcal)
{




  int nbyteecal, nbytehcal, nbyteedge;

    CalHits* ecalhits = new CalHits();
    if(doecal) b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    if(dohcal) b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    if(doedge) b_edge->SetAddress(&edgehits);


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
	
	if((ilayer!=0)&&(ilayer!=1)) std::cout<<"danger danger will robinson ilayer not zero"<<std::endl;
	if(islice>nsliceecal) {
	  std::cout<<"  danger danger will robinson islice nsliceecal are "<<islice<<" "<<nsliceecal<<std::endl;
	  std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<std::endl;
	} else {
	  if(i<SCECOUNT&&ievt<SCECOUNT) std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;
	}

	float ae=aecalhit->energyDeposit;
	eesum+=ae;
	if(islice==0)eesumair+=ae;
	if(islice==1)eesumPDe+=ae;
	if(islice==2)eesumcrystal+=ae;
	if(islice==3)eesumPDe+=ae;

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
