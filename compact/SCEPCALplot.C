
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
#include "../src/DualCrystalCalorimeterHit.h"

#include <vector>
#include <algorithm>


const float sampling_fraction=100.;

void SCEPCALplot(int num_evtsmax, const char* inputfilename) {


  typedef std::vector<dd4hep::sim::Geant4Particle*> GenParts;
  typedef std::vector<CalVision::DualCrystalCalorimeterHit*> CalHits;

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
  result = gSystem->Load("libDDDualCrystal");
  result = gSystem->Load("libDDG4Plugins");


  // define histograms
  TH1F *hgenPsize = new TH1F("hgenPsize","number of generator particles",600,0.,40000);
  TH1F *hgenPdgID = new TH1F("hgenpdgID","pdgID of generator particles",600,-200,200);
  TH1F *hcEcalE = new TH1F("hcEcalE","sum crystal ecal energy",100,0.,100.);
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000000);

  // open data and output file for histograms

  //  const char* inputfilename="/data/users/eno/dd4hep/DD4hep/DDDetectors/compact/testSid.root";
  const char* outputfilename="hist.root";

  // get Tree
  //  TFile *f = new TFile(inputfilename);
  //f->Print();
  GenParts* pgenparts = new GenParts();
  CalHits* pcalhits = new CalHits();
  int num_evt;

  TFile* f = TFile::Open(inputfilename);
  TTree* t = (TTree*)f->Get("EVENT;1");
  t->Print();




  // loop over events in the gen loop
  TBranch* b_mc = t->GetBranch("MCParticles");
  int ihaha = b_mc->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_mc->GetName()<<std::endl;
  std::cout<<"num_evt gen loop is "<<num_evt<<std::endl;
  
  if(num_evt>0) {
    GenParts* gens = new GenParts();
    b_mc->SetAddress(&gens);
    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;
      int nbyte = b_mc->GetEntry(ievt);
      if( nbyte>0) {
	std::cout<<" Gen parts "<<nbyte<<" bytes "<<gens->size() <<std::endl;
      }
      //      for(size_t i=0;i<gens->size(); ++i) {
      hgenPsize->Fill(gens->size());
      for(size_t i=0;i<10; ++i) {
	dd4hep::sim::Geant4Particle* agen =gens->at(i);
	//	std::cout<<"       "<<i<<" pdgid "<<agen->pdgID<<std::endl;
	hgenPdgID->Fill(agen->pdgID);
      }
    }
  }





  // loop over events in the ecal loop
  TBranch* b_ecal = t->GetBranch("CrystalEcalBarrelHits");
  ihaha=b_ecal->GetEntries();
  num_evt= std::min(ihaha,num_evtsmax);
  std::cout<<" doing "<<b_ecal->GetName()<<std::endl;
  std::cout<<"num_evt ecal hit loop is "<<num_evt<<std::endl;
  
  if(num_evt>0) {
    CalHits* ecalhits = new CalHits();
    b_ecal->SetAddress(&ecalhits);
    for(int ievt=0;ievt<num_evt; ++ievt) {
      std::cout<<"event number is "<<ievt<<std::endl;
      int nbyte = b_ecal->GetEntry(ievt);
      if( nbyte>0) {
	std::cout<<" Ecal Hits "<<nbyte<<" bytes "<<ecalhits->size() <<std::endl;
      }
      float esum=0.;
      int ncertot=0;
      for(size_t i=0;i<10; ++i) {
	CalVision::DualCrystalCalorimeterHit* aecalhit =ecalhits->at(i);
	//	std::cout<<"       "<<i<<" energy "<<aecalhit->energyDeposit<<std::endl;
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
      }
      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);
    }
  }
    


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcEcalncer->Write();
  out->Close();

}


