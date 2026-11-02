#include "TH1.h"
#include "TH1F.h"

#include <string>



void rmsrms(const char* filename, const char* hname1) 
{ 
  TString canvName = "Fig_";
  canvName += "Crystal";

  std::string str1 = "shower rms (mm)";
  const char* atitle = str1.c_str();

  //std::string strn1="ervrms3";
  //const char* hname1 =strn1.c_str();

  TFile *f1 = new TFile(filename);



  
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  

  int W = 800;
  int H = 600;
  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  // references for T, B, L, R
  float T = 0.08*H;
  float B = 0.12*H; 
  float L = 0.12*W;
  float R = 0.04*W;
  

  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  //canv->SetLeftMargin( L/W );
  //canv->SetRightMargin( R/W );
  //canv->SetTopMargin( T/H );
  //canv->SetBottomMargin( B/H );
  //canv->SetTickx(0);
  //canv->SetTicky(0);
  


  
  int n_ = 2;
  
  float x1_l = 0.8;
  //  float x1_l = 0.75;
  float y1_l = 0.90;
  
  float dx_l = 0.30;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  


  std::cout<<"getting first"<<std::endl;
  TH2F *A_pt = static_cast<TH2F*>(f1->Get(hname1)->Clone());
  //A_pt->GetYaxis()->SetTitle("measu/true  ");  
  //A_pt->GetYaxis()->SetTitleSize(0.05);  
  //A_pt->GetXaxis()->SetTitle(atitle);  
  //A_pt->GetXaxis()->SetTitleSize(0.05);  
  //A_pt->SetMarkerColor(kGreen);
  //A_pt->SetMarkerSize(0.5);
  //A_pt->SetLineColor(1);
  //A_pt->SetLineWidth(3);
  //A_pt->SetStats(0);
  A_pt->Draw("SCAT");



  float t = canv->GetTopMargin();
  float r = canv->GetRightMargin();
  float Offset   = 0.2;
  TString alabel="simulation meas/true vs rms";
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  latex.SetTextFont(42);
  latex.SetTextAlign(31);
  latex.SetTextSize(0.75*t);
  latex.DrawLatex(1-r,1-t+Offset*t,alabel);
  
  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");

  
  Double_t intercept,slope,int_err,sl_err;

  TCanvas* canv2 = new TCanvas("yuck","yuck",50,50,W,H);
  std::cout<<"make profile"<<std::endl;
  TH2F *fat = (TH2F*)A_pt->Rebin2D(4,1);
  TProfile* A_pt_pfx = fat->ProfileX();
  int nnn = fat->GetNbinsX();
  float xxxmax = fat->GetXaxis()->GetXmax();
  float delta = xxxmax/float(nnn);
  for (int ii=0;ii<nnn;ii++ ){
    float aa = A_pt_pfx->GetBinContent(ii);
    std::cout<<" bin "<<ii<<" at x "<<delta*(ii+0.5)<<ii<<" content "<<aa<<std::endl;
  }




  std::cout<<"do fit "<<std::endl;
  float low=15.;
  float high=30.;
  std::cout<<"fitting over x range "<<low<<"-"<<high<<std::endl;
  A_pt_pfx->Fit("pol1","V","",low,high);
  std::cout<<" get fit function "<<std::endl;
  TF1 *fitFun = (TF1*)A_pt_pfx->GetListOfFunctions()->FindObject("pol1");
  intercept= fitFun->GetParameter(0);
  slope= fitFun->GetParameter(1);
  int_err=fitFun->GetParError(0);
  sl_err=fitFun->GetParError(1);


  std::cout<<" slope interecept are "<<slope<<"+-"<<sl_err<<" "<<intercept<<"+-"<<int_err<<std::endl;
  std::cout<<" value at f =1 is "<<slope+intercept<<"=-"<<sqrt(sl_err*sl_err+int_err*int_err)<<std::endl;
  
  
  return;
}



