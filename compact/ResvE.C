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
#include <TGraphErrors.h>
#include <cmath>

TTree          *fChain;   //!pointer to the analyzed TTree or TChain               
Int_t           fCurrent; //!current Tree number in a TChain                       

void resolution(const char* inputfilename,const char* histname,double* aamean,double* aarms, double* aamean_err, double* aarms_err) {

  std::string name1(histname);
  std::string name2(inputfilename);
  name2.erase(name2.size()-5,name2.size());
  name2.erase(0,9);
  std::string name=name1+"-"+name2+".png";
  TCanvas* canv= new TCanvas(name.c_str(),name.c_str(),200,10,700,500);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetTickx(0);
  canv->SetTicky(0);
  gStyle->SetOptFit(000000);


  std::cout<<"file name is "<<inputfilename<<std::endl;
  TFile *f = new TFile(inputfilename);
  TH1F* nhist = static_cast<TH1F*>(f->Get(histname)->Clone());
  cout<<"number of histo entries= "<<nhist->GetEntries()<<endl;
  nhist->GetXaxis()->SetLimits(0,1.5);
  Int_t imax = nhist->GetMaximumBin();
  Double_t amax = nhist->GetBinCenter(imax);
  std::cout<<"hist max at "<<amax<<std::endl;
  Double_t arms = nhist->GetRMS();
  Double_t amean = nhist->GetMean();
  std::cout<<"hist mean rms is "<<amean<<" "<<arms<<std::endl;
  TF1 *f1 = new TF1("f1","gaus",amean-1.5*arms,amean+1.5*arms);
  nhist->Fit("f1");
  TF1 *fit=nhist->GetFunction("f1");
  Double_t p0= f1->GetParameter(0);
  Double_t p1= f1->GetParameter(1);
  Double_t p2= f1->GetParameter(2);
  Double_t p0e = f1->GetParError(0);
  Double_t p1e = f1->GetParError(1);
  Double_t p2e = f1->GetParError(2);

  std::cout<<"fit parameters are "<<p0<<" "<<p1<<" "<<p2<<std::endl;
  std::cout<<"fit parameters errors: "<<p0e<<", "<<p1e<<", "<<p2e<<std::endl;
  std::cout<<"error bars = "<<p2/p1*(sqrt(pow(p1e/p1,2)+pow(p2e/p2,2)))<<std::endl;
  *aamean=p1;
  *aarms=p2;
  *aamean_err=p1e;
  *aarms_err=p2e;

  nhist->Draw("");
  canv->Print(name.c_str());
  canv->Update();
}


void ResvE() {

  const int npoints=8;
  //const int npoints=2;
  std::string filenames[npoints];
  double aatruemean[npoints];
  std::string aatype = "FSCEPonly";

  filenames[0]="./output/hists_"+aatype+"_10GeV.root";
  filenames[1]="./output/hists_"+aatype+"_15GeV.root";
  filenames[2]="./output/hists_"+aatype+"_20GeV.root";
  filenames[3]="./output/hists_"+aatype+"_25GeV.root";
  filenames[4]="./output/hists_"+aatype+"_30GeV.root";
  filenames[5]="./output/hists_"+aatype+"_35GeV.root";
  filenames[6]="./output/hists_"+aatype+"_40GeV.root";
  filenames[7]="./output/hists_"+aatype+"_45GeV.root";
  //filenames[8]="./output/hists_"+aatype+"_50GeV.root";
  //filenames[9]="./output/hists_"+aatype+"_100GeV.root";
  //filenames[0]="./output/hists_"+aatype+"_20GeV.root";
  //filenames[1]="./output/hists_"+aatype+"_50GeV.root";
  
  aatruemean[0]=10;
  aatruemean[1]=15;
  aatruemean[2]=20;
  aatruemean[3]=25;
  aatruemean[4]=30;
  aatruemean[5]=35;
  aatruemean[6]=40;
  aatruemean[7]=45;
  //aatruemean[8]=50;
  //aatruemean[9]=100;

  const int nhst=5;
  double aamean[npoints][nhst],aarms[npoints][nhst],rrres[npoints][nhst];
  double aamean_err[npoints][nhst],aarms_err[npoints][nhst],rrres_err[npoints][nhst];
 
  vector<string> hnam(nhst);
  hnam[0]="pinphoton_0";
  hnam[1]="pinphoton_2";
  hnam[2]="phcEcalcorr";
  hnam[3]="pienergy_5";
  hnam[4]="hpdepcal";


  //cout<<"number of histo entries= "<<hnam[0].GetEntries()<<endl;


  double abc,dej,feh,ikl;
  for(int k=0;k<nhst;k++) {
    for(int j=0;j<npoints;j++){
      std::cout<<"k j are "<<k<<" "<<j<<" fitting "<<hnam[k]<<" at energy "<<aatruemean[j]<<std::endl;
      resolution(filenames[j].c_str(),hnam[k].c_str(),&abc,&dej);
      aaamean[j][k]=abc;
      aarms[j][k]=dej;
      aamean_err[j][k]=feh;
      aarms_err[j][k]=ikl;
      rrres[j][k]=0;
      rrres_err[j][k]=0;
      if(abc!=0){ 
	      rrres[j][k]=aarms[j][k]/abc;
	      rrres_err[j][k]=sqrt(pow(ikl/dej,2)+pow(feh/abc,2))*rrres[j][k];
	      //cout<<"computed in the loop = "<<rrres_err[j][k]<<" error terms"<<pow(dej/ikl,2)+pow(abc/feh,2)<<" dej="<<dej<<ikl<<endl;
      }
    }
  }

  


  TF1 *f2 = new TF1("f2","sqrt([0]*[0]/x+[1]*[1]/x/x+[2]*[2])");
  TF1 *calice = new TF1("calice","sqrt(0.57*0.57/x+0.016*0.016)",5,100);  // arXiv:1507.05892 taken from pg 21 but needs to be checked


  // digitization of Marco's pure hcal resolution
  //double Marco_y[7]={0.11217,0.076845,0.064773,0.046040,0.030399,0.024674,0.019147};;
  //double Marco_x[7]={5.206,10.528,15.944,31.143,61.259,121.14,288.60};
  //auto marco_g = new TGraph(7,Marco_x,Marco_y);
  //marco_g->Fit("f2");
  TF1 *mC = new TF1("MC","sqrt(0.7*0.7/x+0.11*0.11)",6,110);
  TF1 *mS = new TF1("MS","sqrt(0.32*0.32/x+0.08*0.08)",6,110);
  TF1 *mD = new TF1("MD","sqrt(0.25*0.25/x+0.01*0.01)",6,110);





  // cerenkov
  double arrres[npoints];
  double arrres_err[npoints];
  for (int k=0;k<npoints;k++) {
    arrres[k]=rrres[k][0];
    arrres_err[k]=rrres_err[k][0];
    std::cout<<" cerenkov point "<<k<<" energy "<<aatruemean[k]<<" = "<<arrres[k]<<std::endl;
  }
  auto g1 = new TGraphErrors(npoints,aatruemean,arrres,0,arrres_err);
  g1->Fit("f2");

  // scintillator
  for (int k=0;k<npoints;k++) {
    arrres[k]=rrres[k][1];
    arrres_err[k]=rrres_err[k][1];
    std::cout<<" scint point "<<k<<" energy "<<aatruemean[k]<<" = "<<arrres[k]<<std::endl;
  }
  auto g2 = new TGraphErrors(npoints,aatruemean,arrres,0,arrres_err);
  g2->Fit("f2");

  // dual
  for (int k=0;k<npoints;k++) {
    arrres[k]=rrres[k][2];
    arrres_err[k]=rrres_err[k][2];
    std::cout<<" dual point "<<k<<" energy "<<aatruemean[k]<<" = "<<arrres[k]<<std::endl;
  }
  auto g3 = new TGraphErrors(npoints,aatruemean,arrres,0,arrres_err);
  g3->Fit("f2");

  // only escaping
  for (int k=0;k<npoints;k++) {
	  arrres[k]=rrres[k][3];
	  arrres_err[k]=rrres_err[k][3];
	  std::cout<<" escaping point "<<k<<" energy "<<aatruemean[k]<<" = "<<arrres[k]<<std::endl;
  }
  auto g4 = new TGraphErrors(npoints,aatruemean,arrres,0,arrres_err);
  g4->Fit("f2");

  // all deposited energies
  for (int k=0;k<npoints;k++) {
	  arrres[k]=rrres[k][4];
	  arrres_err[k]=rrres_err[k][4];
	  std::cout<<" alldepE point "<<k<<" energy "<<aatruemean[k]<<" = "<<arrres[k]<<std::endl;
  }
  auto g5 = new TGraphErrors(npoints,aatruemean,arrres,0,arrres_err);
  g5->Fit("f2");


  auto Canvas= new TCanvas("Canvas","Canvas",200,10,700,500);
  gStyle->SetOptStat(000000);
  //gStyle->SetOptFit();

  float x1_l = 0.85;
  float y1_l = 0.80;
  float dx_l = 0.35;
  float dy_l = 0.25;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);
  //gPad->SetLogy();


  TH1 *frame = new TH1F("frame","",1000,0,110);
  frame->SetMinimum(0.01);
  frame->SetMaximum(0.5);
  frame->SetStats(0);
  //frame->SetTitle("FiberOnly HCAL pi- Resolution vs Energy");
  frame->SetTitle("Sampling HCAL pi- Resolution vs Energy");
  frame->GetXaxis()->SetTitle("Beam Energy (GeV)");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.03);
  frame->GetYaxis()->SetTitle("resolution");
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->Draw("");

  //calice->Draw("same");
  //lgd->AddEntry(calice, "calice detector resolution", "l");


   
  mC->SetLineColor(kBlue);
  mC->SetLineStyle(2);
  //mC->Draw("same");
  //lgd->AddEntry(mC, "Marco's C resolution", "l");
  mS->SetLineColor(kGreen);
  mS->SetLineStyle(2);
  //mS->Draw("same");
  //lgd->AddEntry(mS, "Marco's S resolution", "l");
  mD->SetLineColor(kRed);
  mD->SetLineStyle(2);
  //mD->Draw("same");
  //lgd->AddEntry(mD, "Marco's dual resolution", "l");
  

  // digitization of Chekanov's 40L-PFQ
  const int npts =4;
  double C40LPFQ_s_y[npts]={0.32,0.23,0.195,0.17};;
  double C40LPFQ_s_x[npts]={5.,10.,20.,40.};
  double C40LPFQ_c_y[npts]={0.44,0.29,0.275,0.26};;
  double C40LPFQ_c_x[npts]={5.,10.,20.,40.};

  std::cout<<"Chekanov numbers "<<std::endl;
  std::cout<<"    energy  raw scint corr scint raw chek corrected chek "<<std::endl;
  // convert sigma90 to sigma
  for(int jjj=0;jjj<npts;jjj++) {
    std::cout<<"    "<<C40LPFQ_s_x[jjj]<<" "<<C40LPFQ_s_y[jjj]<<" "<<1.25*C40LPFQ_s_y[jjj]<<" "<<C40LPFQ_c_y[jjj]<<" "<<1.25*C40LPFQ_c_y[jjj]<<std::endl;
    C40LPFQ_s_y[jjj]=C40LPFQ_s_y[jjj]*1.25;
    C40LPFQ_c_y[jjj]=C40LPFQ_c_y[jjj]*1.25;

  }

  /*auto C40LPFQ_s_g = new TGraph(npts,C40LPFQ_s_x,C40LPFQ_s_y);
  f2->SetParameter(0,0.5);
  f2->SetParameter(1,0.5);
  f2->SetParameter(2,0.5);
  C40LPFQ_s_g->Fit("f2");
  f2 = C40LPFQ_s_g->GetFunction("f2");
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(2);
  f2->SetLineStyle(2);
  C40LPFQ_s_g->SetLineColor(kGreen);
  C40LPFQ_s_g->SetLineStyle(2);
  C40LPFQ_s_g->SetMarkerSize(1.0);
  C40LPFQ_s_g->SetMarkerStyle(4.0);
  C40LPFQ_s_g->SetMarkerColor(kGreen);
  //  C40LPFQ_s_g->Draw("P");
  //lgd->AddEntry(C40LPFQ_s_g, "Chekanov C40LPFQ's S resolution", "l");
  auto C40LPFQ_c_g = new TGraph(npts,C40LPFQ_c_x,C40LPFQ_c_y);
  f2->SetParameter(0,0.5);
  C40LPFQ_c_g->Fit("f2");
  f2 = C40LPFQ_c_g->GetFunction("f2");
  f2->SetLineColor(kBlue);
  f2->SetLineWidth(2);
  f2->SetLineStyle(2);
  C40LPFQ_c_g->SetLineColor(kBlue);
  C40LPFQ_c_g->SetLineStyle(2);
  C40LPFQ_c_g->SetMarkerSize(1.0);
  C40LPFQ_c_g->SetMarkerStyle(4.0);
  C40LPFQ_c_g->SetMarkerColor(kBlue);
  //C40LPFQ_c_g->Draw("P");
  //lgd->AddEntry(C40LPFQ_c_g, "Chekanov C40LPFQ's C resolution", "l");

  f2 = g1->GetFunction("f2");
  f2->SetLineColor(kBlue);
  f2->SetLineWidth(1);
  g1->SetMarkerColor(kBlue);
  g1->SetLineColor(kBlue);
  g1->SetMarkerStyle(21);
  g1->SetMarkerSize(1.0);
  g1->Draw("P");
  lgd->AddEntry(g1, "ncer", "P");

  f2 = g2->GetFunction("f2");
  f2->SetLineColor(kGreen);
  f2->SetLineWidth(1);
  g2->SetMarkerColor(kGreen);
  g2->SetLineColor(kGreen);
  g2->SetMarkerStyle(21);
  g2->SetMarkerSize(1.0);
  g2->Draw("P");
  lgd->AddEntry(g2, "nscint", "P");

  f2 = g3->GetFunction("f2");
  f2->SetLineColor(kRed);
  f2->SetLineWidth(1);
  g3->SetMarkerColor(kRed);
  g3->SetLineColor(kRed);
  g3->SetMarkerStyle(23);
  g3->SetMarkerSize(1.0);
  g3->Draw("P");
  lgd->AddEntry(g3, "dual", "P");


  f2 = g4->GetFunction("f2");
  f2->SetLineColor(kMagenta);
  f2->SetLineWidth(1);
  g4->SetMarkerColor(kMagenta);
  g4->SetLineColor(kMagenta);
  g4->SetMarkerStyle(23);
  g4->SetMarkerSize(1.0);
  g4->Draw("P");
  lgd->AddEntry(g4, "escaping", "P");

  
 
  f2 = g5->GetFunction("f2");
  f2->SetLineColor(kCyan);
  f2->SetLineWidth(1);
  g5->SetMarkerColor(kCyan);
  g5->SetLineColor(kCyan);
  g5->SetMarkerStyle(23);
  g5->SetMarkerSize(1.0);
  g5->Draw("P");
  lgd->AddEntry(g5, "all deposited energy", "P");
  

  lgd->Draw();
  Canvas->Print("resolution.png");


  return;


}
