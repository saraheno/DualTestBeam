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


void crystalana(int num_evtsmax, const char* inputfilename) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrysCalorimeterHit*> CalHits;


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
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",600,-200,200);


  // calorimeter infor
  TH1F *hchan = new TH1F("hchan","channel ID number",1028,0.,1028);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy",100,0.,100.);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000);

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


    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;


      // gen particles
      int nbytegen = b_mc->GetEntry(ievt);
      if( nbytegen>0) {
      std::cout<<" gen byte "<<nbytegen<<" bytes "<<std::endl;
      }

      std::cout<<"gen parts size "<<gens->size()<<std::endl;
     hgenPsize->Fill(gens->size());
      for(size_t i=0;i<gens->size(); ++i) {
	dd4hep::sim::Geant4Particle* agen =gens->at(i);
	hgenPdgID->Fill(agen->pdgID);
      }

      int nbyteecal = b_ecal->GetEntry(ievt);
      if( nbyteecal>0) {
      std::cout<<" Ecal Hits bytes "<<nbyteecal<<" bytes "<<std::endl;
      }
      float esum=0.;
      int ncertot=0;
      int nscinttot=0;

      // ecal hits
      std::cout<<"ecal size "<<ecalhits->size()<<std::endl;
      for(size_t i=0;i<ecalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* aecalhit =ecalhits->at(i);
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
	nscinttot+=aecalhit->nscintillator;
	std::cout<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
        int ihitchan=aecalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0x1F;
	int iy =(ihitchan >>6) & 0x1F;
        int  ilayer = (ihitchan >>13) & 0x07;
        int  islice = (ihitchan>> 16) & 0x07;
std::cout<<"idet,iy,ix,ilayer, islice is ("<<idet<<","<<std::hex<<iy<<","<<ix<<","<<std::dec<<ilayer<<","<<islice<<")"<<std::endl;
	hchan->Fill(aecalhit->cellID);
      }  // end loop over ecal hits
      // hcal hits


      int nbytehcal = b_hcal->GetEntry(ievt);
      if( nbytehcal>0) {
      std::cout<<" Ecal Hits  byte"<<nbytehcal<<" bytes "<<std::endl;
      }

      std::cout<<"hcal size "<<hcalhits->size()<<std::endl;
      for(size_t i=0;i<hcalhits->size(); ++i) {
	CalVision::DualCrysCalorimeterHit* aecalhit =hcalhits->at(i);
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
	nscinttot+=aecalhit->nscintillator;
	std::cout<<" hit channel (hex) is "<< std::hex<<aecalhit->cellID<<std::dec<<" (energy,nceren,nscin)=("<<aecalhit->energyDeposit<<","<<aecalhit->ncerenkov<<","<<aecalhit->nscintillator<<")"<<std::endl;
        int ihitchan=aecalhit->cellID;
        int idet = (ihitchan) & 0x07;
	int ix = (ihitchan >>3) & 0x1F;
	int iy =(ihitchan >>6) & 0x1F;
std::cout<<"idet,iy,ix is ("<<idet<<","<<std::hex<<iy<<","<<ix<<","<<std::dec<<")"<<std::endl;
	hchan->Fill(aecalhit->cellID);
      }  // end loop over hcal hits

    
      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);


      std::cout<<" total energy deposit "<<esum<<std::endl;
      std::cout<<" total number of cherenkov is "<<ncertot<<std::endl;
      std::cout<<" total number of scintillator is "<<nscinttot<<std::endl;


    }  //end loop over events
  }  // end if no events
    
  


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcEcalncer->Write();
  out->Close();

}

