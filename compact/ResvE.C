#include <iostream>
#include <iomanip>
#include <locale>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
using std::vector;

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
Int_t           fCurrent; //!current Tree number in a TChain                       





void resolution(const char* inputfilename,const char* histname,double* aamean,double* aarms) {

  std::string name1(histname);
  std::string name2(inputfilename);
  std::string name=name1+name2;
  TCanvas* canv= new TCanvas(name.c_str(),name.c_str(),200,10,700,500);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
  gStyle->SetOptFit();


  std::cout<<"file name is "<<inputfilename<<std::endl;
  TFile *f = new TFile(inputfilename);
  TH1F* nhist = static_cast<TH1F*>(f->Get(histname)->Clone());
  Int_t imax = nhist->GetMaximumBin();
  Double_t amax = nhist->GetBinCenter(imax);
  std::cout<<"hist max at "<<amax<<std::endl;
  Double_t arms = nhist->GetRMS();
  std::cout<<"hist rms is "<<arms<<std::endl;
  TF1 *f1 = new TF1("f1","gaus",amax-2*arms,amax+2*arms);
  nhist->Fit("f1","R");
  TF1 *fit=nhist->GetFunction("f1");
  Double_t p0= f1->GetParameter(0);
  Double_t p1= f1->GetParameter(1);
  Double_t p2= f1->GetParameter(2);
  std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
  *aamean=p1;
  *aarms=p2;

  nhist->Draw("");
  canv->Print(name.c_str(),".png");
  canv->Update();


}



void res() {

  const int npoints=2;
  const char* filenames[npoints];

  filenames[0]="hists_10.root"; 
  filenames[1]="hists_30.root"; 

  

  double aatruemean[npoints];
  aatruemean[0]=10;aatruemean[1]=30;

  const int nhst=3;
  double aaamean[npoints][nhst],aarms[npoints][nhst],rrres[npoints][nhst];
 
  vector<string> hnam(nhst);
  hnam[0]="phcHcalncer";
  hnam[1]="phcHcalnscint";
  hnam[2]="phcHcalcorr";




  double abc,dej;
  for(int k=0;k<nhst;k++) {
    for(int j=0;j<npoints;j++){
      std::cout<<"k j are "<<k<<" "<<j<<" fitting "<<hnam[k]<<" at energy "<<aatruemean[j]<<std::endl;
      resolution(filenames[j],hnam[k].c_str(),&abc,&dej);
      aaamean[j][k]=abc;
      aarms[j][k]=dej;
      rrres[j][k]=0;
      if(abc!=0) rrres[j][k]=aarms[j][k]/abc;
    }
  }

  


  TF1 *f2 = new TF1("f2","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
  TF1 *calice = new TF1("calice","sqrt(0.57*0.57/x+0.016*0.016)",5,100);  // arXiv:1507.05892 taken from pg 21 but needs to be checked

  double arrres[npoints];
  for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][0];}
  auto g1 = new TGraph(npoints,aatruemean,arrres);
  g1->Fit("f2");


  for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][1];}
  auto g2 = new TGraph(npoints,aatruemean,arrres);
  g2->Fit("f2");


  for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][2];}
  auto g3 = new TGraph(npoints,aatruemean,arrres);
  g3->Fit("f2");


  auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
  //  gStyle->SetOptStat(111111);
  //gStyle->SetOptFit();

  float x1_l = 0.9;
  float y1_l = 0.80;
  float dx_l = 0.60;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);


  TH1 *frame = new TH1F("frame","",1000,0,120);
  frame->SetMinimum(0.);
  frame->SetMaximum(1.0);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("true energy (GeV)");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetYaxis()->SetTitle("percent resolution");
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->Draw("");

  calice->Draw("same");
  lgd->AddEntry(calice, "calice detector resolution", "l");

  g1->SetMarkerColor(kBlue);
  g1->SetMarkerStyle(21);
  g1->SetMarkerSize(1.5);
  g1->Draw("P");
  lgd->AddEntry(g1, "ncer", "l");

  g2->SetMarkerColor(kGreen);
  g2->SetMarkerStyle(23);
  g2->SetMarkerSize(1.5);
  g2->Draw("P");
  lgd->AddEntry(g2, "nscint", "l");


  g3->SetMarkerColor(kRed);
  g3->SetMarkerStyle(23);
  g3->SetMarkerSize(1.5);
  g3->Draw("P");
  lgd->AddEntry(g3, "dual", "l");

  

  lgd->Draw();
  Canvas->Print("resolution.png",".png");


  return;


}