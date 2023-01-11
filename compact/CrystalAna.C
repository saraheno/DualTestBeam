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
bool dogen=1;
bool doecal=1;
bool dohcal=1;
bool doedge=1;
int SCECOUNT=50;



// this is now hardwared in DualCrysCalorimeterHit.h
// need to figure out how to charge this
//const int HARDWIREDmax=1000;

// DANGER DANGER WILL ROBINSON!!!!!!!!!!!!!!!!!!!!!!!!
//  this must be changed whenever you change the hcalgeometry
float hcalcalib=1./0.036;


void crystalana(int num_evtsmax, const char* inputfilename, const float beamE) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;
  typedef dd4hep::sim::Geant4HitData::MonteCarloContrib Contribution;
  typedef std::vector<dd4hep::sim::Geant4HitData::MonteCarloContrib> Contributions;

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


  // define histograms

  //gen particles
  TH1F *hgenPsize = new TH1F("hgenPsize","number of generator particles",600,0.,40000);
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",1000,-500,500);
  TH1F *hgenfrstpid = new TH1F("hgenfrstpid","pdgID of first gen",1000,-500,500);
  TH1F *hgenfrstE = new TH1F("hgenfrstE","energy of first gen",500,0.,500.);


  // calorimeter infor
  TH1F *hchan = new TH1F("hchan","channel ID number",1028,0.,1028);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy / beam E",100,0.,1.5);
  TH1F *hcHcalE = new TH1F("hcHcalE","sum fiber hcal energy / beam E",100,0.,1.5);
  TH1F *hcEdgeE = new TH1F("hcEdgeE","sum escaping / beam E",100,0.,1.5);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100000,0.,100000);
  TH1F *hcHcalncer = new TH1F("hcHcalncer","total number of cerenkov",100000,0.,100000);


  TH2F *hecal2d = new TH2F("hecal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);
  TH2F *hhcal2d = new TH2F("hhcal2d","lego of ecal", 41,-20.,20.,41,-20.,20.);


  TH1F *haphcal = new TH1F("haphcal","ratio of fiber to total to  energy hcal",50,0.,0.2);
  TH1F *heest = new TH1F("heest","ratio estimated to true energy",500,0.,1.);
  TH1F *hetrue = new TH1F("hetrue","ratio deposited to incident energy",500,0.,1.);

  TH1F *hnecalcon = new TH1F("hnecalcon","number contribs to ecal hit",1010,-10.,1000.);

  TH2F *hzvst = new TH2F("hzvst","z position of hit versus time ",100,-50.,300.,100,0.,100.); 


  // open data and output file for histograms

  //  const char* inputfilename="/data/users/eno/dd4hep/DD4hep/DDDetectors/compact/testSid.root";
  const char* outputfilename="hist.root";

  // get Tree
  //  TFile *f = new TFile(inputfilename);
  //f->Print();

  TFile* f = TFile::Open(inputfilename);
  TTree* t = (TTree*)f->Get("EVENT;1");
  t->Print();

  

  TBranch* b_mc = t->GetBranch("MCParticles");
  TBranch* b_ecal = t->GetBranch("DRCNoSegment");
  TBranch* b_hcal = t->GetBranch("DRFNoSegment");
  TBranch* b_edge = t->GetBranch("EdgeDetNoSegment");


  int ihaha = b_mc->GetEntries();
  int num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<"num_evt is  "<<num_evt<<std::endl;
  
  
   

  // loop over events 
  
  if(num_evt>0) {

    GenParts* gens = new GenParts();
    b_mc->SetAddress(&gens);
    CalHits* ecalhits = new CalHits();
    b_ecal->SetAddress(&ecalhits);
    CalHits* hcalhits = new CalHits();
    b_hcal->SetAddress(&hcalhits);
    CalHits* edgehits = new CalHits();
    b_edge->SetAddress(&edgehits);


    for(int ievt=0;ievt<num_evt; ++ievt) {
      if((ievt<SCECOUNT)||(ievt%SCECOUNT)==0) std::cout<<"event number is "<<ievt<<std::endl;


      float mainee=-1.;
      if(dogen){
      // gen particles
	int nbytegen = b_mc->GetEntry(ievt);


	if(ievt<SCECOUNT) std::cout<<"  gen parts size "<<gens->size()<<std::endl;
	hgenPsize->Fill(gens->size());


	for(size_t i=0;i<gens->size(); ++i) {
	  dd4hep::sim::Geant4Particle* agen =gens->at(i);
	  float px=agen->psx;
	  float py=agen->psy;
	  float pz=agen->psz;
	  float mass=agen->mass;
	  float ee=sqrt(mass*mass+px*px+py*py+pz*pz);
	  if(ee>mainee) {
	    mainee=ee;
	    hgenfrstpid->Fill(agen->pdgID);
	    hgenfrstE->Fill(ee/1000);
	    if(i<SCECOUNT) std::cout<<"  gen pid "<<agen->pdgID<<" energy "<<ee<<std::endl;
	  }
	  float vsz=agen->vsz;
	  if(vsz<-100) std::cout<<" pid "<<agen->pdgID<<" "
				<<"vs ("<<agen->vsx<<","<<agen->vsy<<","<<agen->vsz<<") "
				<<"ve ("<<agen->vex<<","<<agen->vey<<","<<agen->vez<<") "
				<<"ps ("<<agen->psx<<","<<agen->psy<<","<<agen->psz<<") "
				<<"pe ("<<agen->pex<<","<<agen->pey<<","<<agen->pez<<") "
				<<std::endl;
	  hgenPdgID->Fill(agen->pdgID);
	}
      }

      // kludge for now
      mainee=beamE*1000;



      float esum=0.;
      float esumair=0;
      float esumcrystal=0;
      float esumPDe=0;
      float esumfiber=0;
      float esumabs=0;
      float esumPDh=0;
      float esumedge=0;
      int ncertotecal=0;
      int nscinttotecal=0;
      int ncertothcal=0;
      int nscinttothcal=0;


      if(doecal) {
      int nbyteecal = b_ecal->GetEntry(ievt);


      // ecal hits
      if(ievt<SCECOUNT) std::cout<<"   ecal size "<<ecalhits->size()<<std::endl;
      for(size_t i=0;i<ecalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);

	ncertotecal+=aecalhit->ncerenkov;
	nscinttotecal+=aecalhit->nscintillator;
	if(i<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
        int ihitchan=aecalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0x1F ;  // is this right?
	if(ix>16) ix=ix-32;
	int iy =(ihitchan >>8) & 0x1F ; // is this right?
	if(iy>16) iy=iy-32;
        int  islice = (ihitchan >>13) & 0x07;
        int  ilayer = (ihitchan>> 16) & 0x07;
	
	if((ilayer!=0)&&(ilayer!=1)) std::cout<<"danger danger will robinson ilayer not zero"<<std::endl;
	if(islice>nsliceecal) {
	  std::cout<<"  danger danger will robinson islice nsliceecal are "<<islice<<" "<<nsliceecal<<std::endl;
	  std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<std::endl;
	} else {
	  if(i<SCECOUNT) std::cout<<"  idet,ix,iy,ilayer, islice is ("<<idet<<","<<ix<<","<<iy<<","<<std::dec<<ilayer<<","<<islice<<")"<<" slice name is "<<nameecalslice[islice]<<std::endl;
	}

	float ae=aecalhit->energyDeposit;
	esum+=ae;
	if(islice==0)esumair+=ae;
	if(islice==1)esumPDe+=ae;
	if(islice==2)esumcrystal+=ae;
	if(islice==3)esumPDe+=ae;

	hchan->Fill(aecalhit->cellID);
	hecal2d->Fill(ix,iy,aecalhit->energyDeposit);

	// get MC truth information about individual contributions to this hit
	Contributions zxzz=aecalhit->truth;
	float hacheck=0.;
	if(i<SCECOUNT) std::cout<<"    number of contributes "<<zxzz.size()<<std::endl;
	hnecalcon->Fill(zxzz.size());
	//	if(zxzz.size()>HARDWIREDmax) std::cout<<" number of const "<<zxzz.size()<<" greater than hardwared limit in SDAction so that information is missing"<<std::endl;
	for(size_t j=0;j<zxzz.size(); j++) {
	  //	  std::cout<<"testing truth truth number "<<i<<" with pdgID "<<(zxzz.at(i)).pdgID<<std::endl;
	  // other member functions are trackID, deposit, time, length, x,y,z
	  // see DDG4/include/DDG4/Geant4Data.h
	  hacheck+=(zxzz.at(j)).deposit;
	  // right now can only save 1000.  lose info for more
	  //if(j<HARDWIREDmax) std::cout<<"     contrib charge ["<<j<<"] "<<aecalhit->contribCharge[j]<<" pid is "<<(zxzz.at(j)).pdgID<<" velo "<<aecalhit->contribBeta[j]<<std::endl;
	  if((i<SCECOUNT)&&(j<10)) std::cout<<"     contrib charge ["<<j<<"] "<<aecalhit->contribCharge[j]<<" pid is "<<(zxzz.at(j)).pdgID<<" velo "<<aecalhit->contribBeta[j]<<std::endl;
	  hzvst->Fill((zxzz.at(j)).z,(zxzz.at(j)).time);
	}
	if(i<SCECOUNT) {
	  std::cout<<"    difference between truth sum and total deposit is "<<hacheck-ae<<" where "<<ae<<" is the hit size."<<std::endl;
	  if(ae>0) std::cout<<"      percent diff is "<<(hacheck-ae)/ae<<std::endl;
	}


      }  // end loop over ecal hits
      }


      if(dohcal) {
      // hcal hits


      int nbytehcal = b_hcal->GetEntry(ievt);


      if(ievt<SCECOUNT) std::cout<<"   hcal size "<<hcalhits->size()<<std::endl;
      for(size_t i=0;i<hcalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* ahcalhit =hcalhits->at(i);

	ncertothcal+=ahcalhit->ncerenkov;
	nscinttothcal+=ahcalhit->nscintillator;
	if(i<SCECOUNT) std::cout<<std::endl<<" hit channel (hex) is "<< std::hex<<ahcalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<ahcalhit->energyDeposit<<","<<ahcalhit->ncerenkov<<","<<ahcalhit->nscintillator<<")"<<std::endl;
        int ihitchan=ahcalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0x7F;   // is this right?
	if(ix>64) ix=ix-128;
	int iy =(ihitchan >>10) & 0x7F;   // is this right?
	if(iy>64) iy=iy-128;
	int ifiber  =(ihitchan >>17) & 0x03;
	int iabs=(ihitchan >>19) & 0x03;
	int iphdet=(ihitchan >>21) & 0x03;
	if(i<SCECOUNT) std::cout<<"  idet,ix,iy is ("<<idet<<","<<ix<<","<<iy<<")"<<std::endl;
	if(i<SCECOUNT) std::cout<<"  ifiber,iabs,iphdet is ("<<ifiber<<","<<iabs<<","<<iphdet<<")"<<std::endl;


	float ah=ahcalhit->energyDeposit;
	esum+=ah;

	if(ifiber==1) esumfiber+=ah;
	if(iabs==1) esumabs+=ah;
	if(iphdet==1) esumPDh+=ah;




	hchan->Fill(ahcalhit->cellID);
	hhcal2d->Fill(ix,iy,ahcalhit->energyDeposit);
      }  // end loop over hcal hits
      }





      if(doedge) {
      int nbyteedge = b_edge->GetEntry(ievt);

      // energies escaping calroimeter
      std::cout<<"   edge size "<<edgehits->size()<<std::endl;
      for(size_t i=0;i<edgehits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* aedgehit =edgehits->at(i);


	float ae=aedgehit->energyDeposit;
	esum+=ae;
	esumedge+=ae;

      }  // end loop over escaping hits
      }





    
      hcEcalE->Fill(esumcrystal/1000./mainee);
      hcHcalE->Fill(esumfiber*hcalcalib/1000./mainee);
      hcEdgeE->Fill(esumedge/1000./mainee);
      hcEcalncer->Fill(ncertotecal);
      hcHcalncer->Fill(ncertothcal);
      if(esumfiber>0) haphcal->Fill(esumfiber/(esumabs+esumfiber));
      heest->Fill((esumcrystal+esumfiber*hcalcalib)/mainee);



      float achecks=esumair+esumPDe+esumcrystal+esumfiber+esumabs+esumPDh+esumedge;
      hetrue->Fill(achecks/mainee);


      std::cout<<std::endl<<std::endl<<"total energy deposit "<<esum/1000.<<std::endl;
      std::cout<<"       in air "<<esumair/1000.<<std::endl;
      std::cout<<"       in photodetector ecal "<<esumPDe/1000.<<std::endl;
      std::cout<<"       in crystal "<<esumcrystal/1000.<<std::endl;
      std::cout<<"       in fiber "<<esumfiber/1000.<<std::endl;
      std::cout<<"       in absorber "<<esumabs/1000.<<std::endl;
      std::cout<<"       in photodetect hcal "<<esumPDh/1000.<<std::endl;
      std::cout<<"       escaping detector "<<esumedge/1000.<<std::endl;

      std::cout<<"       sum individual "<<achecks/1000.<<std::endl;
      std::cout<<"   incident energy "<<mainee/1000.<<std::endl;
      std::cout<<"   ratio to incident energy "<<achecks/mainee<<std::endl;

      std::cout<<"total number of cherenkov ecal is "<<ncertotecal<<std::endl;
      std::cout<<"total number of scintillator ecal is "<<nscinttotecal<<std::endl<<std::endl;
      std::cout<<"total number of cherenkov hcal is "<<ncertothcal<<std::endl;
      std::cout<<"total number of scintillator hcal is "<<nscinttothcal<<std::endl<<std::endl;









    }  //end loop over events
  }  // end if no events

  //  float amean = hceest->GetMean();

    
  


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcHcalE->Write();
  hcEdgeE->Write();
  hcEcalncer->Write();
  hecal2d->Write();
  hhcal2d->Write();
  haphcal->Write();
  heest->Write();
  hgenfrstpid->Write();
  hgenfrstE->Write();
  hetrue->Write();
  hnecalcon->Write();
  hzvst->Write();
  out->Close();

}



void CrystalAna() {
  crystalana(5,"out.root",10.);
  return;
}
