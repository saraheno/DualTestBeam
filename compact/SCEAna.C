
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

const int nchan = 4;
const int ichan[nchan] = {198,166,230,6};


void SCEana(int num_evtsmax, const char* inputfilename) {


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
  TH1F *hcEcalncer = new TH1F("hcEcalncer","total number of cerenkov",100,0.,10000);
  TH1F *hcEcalncer0 = new TH1F("hcEcalncer0","total number of cerenkov chan 1",100,0.,10000);
  TH1F *hcEcalncer1 = new TH1F("hcEcalncer1","total number of cerenkov chan 2",100,0.,10000);
  TH1F *hcEcalncer2 = new TH1F("hcEcalncer2","total number of cerenkov chan 3",100,0.,10000);
  TH1F *hcEcalncer3 = new TH1F("hcEcalncer3","total number of cerenkov chan 4",100,0.,10000);

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
      hgenPsize->Fill(gens->size());
      for(size_t i=0;i<gens->size(); ++i) {
        dd4hep::sim::Geant4Particle* agen =gens->at(i);
        hgenPdgID->Fill(agen->pdgID);
      }
    }
  }
  


 
  

  // loop over events in the ecal loop
  TBranch* b_ecal = t->GetBranch("DRCNoSegment");
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
      float esumchan[nchan]={0.,0.,0.,0.};
      int ncerchan[nchan]={0,0,0,0};
      int nscintchan[nchan]={0,0,0,0};
      int ncertot=0;
      int nscinttot=0;
      for(size_t i=0;i<ecalhits->size(); ++i) {
	CalVision::DualCrystalCalorimeterHit* aecalhit =ecalhits->at(i);
	//	std::cout<<"       "<<i<<" energy "<<aecalhit->energyDeposit<<std::endl;
	esum+=aecalhit->energyDeposit;
	ncertot+=aecalhit->ncerenkov;
	nscinttot+=aecalhit->nscintillator;
	std::cout<<" hit channel is "<< aecalhit->cellID<<" "<<aecalhit->energyDeposit<<" "<<aecalhit->ncerenkov<<" "<<aecalhit->nscintillator<<std::endl;
	int ijchan=aecalhit->nbin;
	for (int j=0;j<ijchan;j++) {
	  std::cout<<"  ncerwave["<<j<<"]="<<(aecalhit->ncerwave)[j]<<std::endl;
	  std::cout<<"  nscintwave["<<j<<"]="<<(aecalhit->nscintwave)[j]<<std::endl;
	}

      // there is a better way to do this
	int jchan=aecalhit->cellID;
	int kchan=-1;
	for( int i=0;i<nchan;i++ ) {
	  if(ichan[i]==jchan) kchan=i;
	}
	if(kchan==-1) {
	  std::cout<<"unknown hit channel is "<< aecalhit->cellID<<std::endl;
	} else {
	  esumchan[kchan]+=aecalhit->energyDeposit;
	  ncerchan[kchan]+=aecalhit->ncerenkov;
	  nscintchan[kchan]+=aecalhit->nscintillator;
	}


      }  // end loop over hits

    
      hcEcalE->Fill(esum/1000.);
      hcEcalncer->Fill(ncertot);
      hcEcalncer0->Fill(ncerchan[0]);
      hcEcalncer1->Fill(ncerchan[1]);
      hcEcalncer2->Fill(ncerchan[2]);
      hcEcalncer3->Fill(ncerchan[3]);



      std::cout<<" total energy deposit "<<esum<<std::endl;
      float check=0.;
      for( int i=0;i<nchan;i++) {
	std::cout<<"esum ["<<ichan[i]<<"]="<<esumchan[i]<<std::endl;
	check+=esumchan[i];
      }
      std::cout<<" check total energy desposit "<<check<<std::endl;

      std::cout<<" total number of cherenkov is "<<ncertot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	std::cout<<"ncerenkov ["<<ichan[i]<<"]="<<ncerchan[i]<<std::endl;
	check+=ncerchan[i];
      }
      std::cout<<" check ncerenkov "<<check<<std::endl;


      std::cout<<" total number of scintillator is "<<nscinttot<<std::endl;
      check=0;
      for( int i=0;i<nchan;i++) {
	std::cout<<"nscintillator ["<<ichan[i]<<"]="<<nscintchan[i]<<std::endl;
	check+=nscintchan[i];
      }
      std::cout<<" check nscintillator "<<check<<std::endl;


    }  //end loop over events
  }  // end if no events
    
  


  
 
 
  f->Close();

  TFile * out = new TFile(outputfilename,"RECREATE");
  hgenPsize->Write();
  hgenPdgID->Write();
  hcEcalE->Write();
  hcEcalncer->Write();
  hcEcalncer0->Write();
  hcEcalncer1->Write();
  hcEcalncer2->Write();
  hcEcalncer3->Write();
  out->Close();

}


