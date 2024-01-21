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
  Double_t amean = nhist->GetMean();
  std::cout<<"hist mean rms is "<<amean<<" "<<arms<<std::endl;
  TF1 *f1 = new TF1("f1","gaus",amean-1.5*arms,amean+1.5*arms);
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

  const int npoints=7;
  const char* filenames[npoints];

  filenames[0]="hists_10GeV.root"; 
  filenames[1]="hists_20GeV.root"; 
  filenames[2]="hists_25GeV.root"; 
  filenames[3]="hists_30GeV.root"; 
  filenames[4]="hists_35GeV.root"; 
  filenames[5]="hists_40GeV.root"; 
  filenames[6]="hists_45GeV.root"; 
  //filenames[7]="hists_50GeV_3.root"; 
  //  filenames[5]="hists_100GeV.root"; 

  

  double aatruemean[npoints];
  aatruemean[0]=10;
  aatruemean[1]=20;
  aatruemean[2]=25;
  aatruemean[3]=30;
  aatruemean[4]=35;
  aatruemean[5]=40;
  aatruemean[6]=45;
  //aatruemean[7]=50;
  //  aatruemean[5]=100;

  const int nhst=4;
  double aaamean[npoints][nhst],aarms[npoints][nhst],rrres[npoints][nhst];
 
  vector<string> hnam(nhst);
  hnam[0]="phcHcalncer";
  hnam[1]="phcHcalnscint";
  hnam[2]="phcHcalcorr";
  hnam[3]="phcEdgeR";




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


  // digitization of Marco's pure hcal resolution
  //double Marco_y[7]={0.11217,0.076845,0.064773,0.046040,0.030399,0.024674,0.019147};;
  //double Marco_x[7]={5.206,10.528,15.944,31.143,61.259,121.14,288.60};
  //auto marco_g = new TGraph(7,Marco_x,Marco_y);
  //marco_g->Fit("f2");
  TF1 *mC = new TF1("MC","sqrt(0.7*0.7/x+0.11*0.11)",10,110);
  TF1 *mS = new TF1("MS","sqrt(0.32*0.32/x+0.08*0.08)",10,110);
  TF1 *mD = new TF1("MD","sqrt(0.25*0.25/x+0.01*0.01)",10,110);




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


  for (int k=0;k<npoints;k++) {arrres[k]=rrres[k][3];}
  auto g4 = new TGraph(npoints,aatruemean,arrres);
  g4->Fit("f2");


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
  lgd->SetBorderSize(0); lgd->SetTextSize(0.03); lgd->SetTextFont(62); lgd->SetFillColor(0);


  TH1 *frame = new TH1F("frame","",1000,0,70);
  frame->SetMinimum(0.);
  frame->SetMaximum(0.25);
  frame->SetStats(0);
  frame->GetXaxis()->SetTitle("true energy (GeV)");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetYaxis()->SetTitle("percent resolution");
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->Draw("");

  //calice->Draw("same");
  //lgd->AddEntry(calice, "calice detector resolution", "l");


  mC->SetLineColor(kBlue);
  mC->Draw("same");
  lgd->AddEntry(mC, "Marco's C resolution", "l");
  mS->SetLineColor(kGreen);
  mS->Draw("same");
  lgd->AddEntry(mS, "Marco's S resolution", "l");
  mD->SetLineColor(kRed);
  mD->Draw("same");
  lgd->AddEntry(mD, "Marco's dual resolution", "l");


  g1->SetMarkerColor(kBlue);
  g1->SetMarkerStyle(21);
  g1->SetMarkerSize(1.0);
  g1->Draw("P");
  lgd->AddEntry(g1, "ncer", "P");

  g2->SetMarkerColor(kGreen);
  g2->SetMarkerStyle(23);
  g2->SetMarkerSize(1.0);
  g2->Draw("P");
  lgd->AddEntry(g2, "nscint", "P");


  g3->SetMarkerColor(kRed);
  g3->SetMarkerStyle(23);
  g3->SetMarkerSize(1.0);
  g3->Draw("P");
  lgd->AddEntry(g3, "dual", "P");

  g4->SetMarkerColor(kMagenta);
  g4->SetMarkerStyle(23);
  g4->SetMarkerSize(1.0);
  g4->Draw("P");
  lgd->AddEntry(g4, "escaping", "P");

  

  lgd->Draw();
  Canvas->Print("resolution.png",".png");


  return;


}
