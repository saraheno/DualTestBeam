#include "TH1.h"
#include "TH1F.h"

#include <string>



void figscvfPbWOS() 
{ 
  TString canvName = "Fig_";
  canvName += "scvfPbWOS";

  std::string str1 = "f";
  const char* atitle = str1.c_str();

  std::string strn1="hfncerEcal";
  const char* hname1 =strn1.c_str();
  std::string strn2="hfnscinEcal";
  const char* hname2 =strn2.c_str();

  TFile *f1 = new TFile("hists_20GeV_BigEcal2.root");

 
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
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  


  
  int n_ = 2;
  
  float x1_l = 0.8;
  //  float x1_l = 0.75;
  float y1_l = 0.90;
  
  float dx_l = 0.30;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  
  TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(32); lgd->SetFillColor(0);


  std::cout<<"getting first"<<std::endl;
  TH2F *A_pt = static_cast<TH2F*>(f1->Get(hname1)->Clone());
  A_pt->GetYaxis()->SetTitle("Signal (normalized to electron)  ");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  
  A_pt->GetXaxis()->SetTitle(atitle);  
  A_pt->GetXaxis()->SetTitleSize(0.05);  
  A_pt->SetMarkerColor(kGreen);
  A_pt->SetLineColor(1);
  A_pt->SetLineWidth(3);
  A_pt->SetStats(0);
  A_pt->Draw("");


  std::cout<<"getting second"<<std::endl;
  TH2F *B_pt = static_cast<TH2F*>(f1->Get(hname2)->Clone());
  B_pt->GetYaxis()->SetTitle("Signal (normalized to electron)  ");  
  B_pt->GetYaxis()->SetTitleSize(0.05);  
  B_pt->GetXaxis()->SetTitle(atitle);  
  B_pt->GetXaxis()->SetTitleSize(0.05);  
  B_pt->SetLineColor(1);
  B_pt->SetMarkerColor(kRed);
  B_pt->SetLineWidth(3);
  B_pt->SetStats(0);
  B_pt->Draw("same");
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  lgd->Draw();


  float t = canv->GetTopMargin();
  float r = canv->GetRightMargin();
  float Offset   = 0.2;
  TString alabel="20 GeV PbWO pion simulation";
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

 Double_t intercept,slope;

  TCanvas* canv2 = new TCanvas("yuck","yuck",50,50,W,H);
  TProfile* A_pt_pfx = A_pt->ProfileX();
  A_pt_pfx->Fit("pol1","WW","",0.5,0.9);
  TF1 *fitFun = (TF1*)A_pt_pfx->GetListOfFunctions()->FindObject("pol1");
  intercept= fitFun->GetParameter(0);
  slope= fitFun->GetParameter(1);


  TCanvas* canv3 = new TCanvas("yuck2","yuck2",50,50,W,H);
  TProfile* B_pt_pfx = B_pt->ProfileX();
  B_pt_pfx->Fit("pol1","WW","",0.5,0.9);
  TF1 *fitFun2 = (TF1*)B_pt_pfx->GetListOfFunctions()->FindObject("pol1");
  intercept= fitFun2->GetParameter(0);
  slope= fitFun2->GetParameter(1);

  
  return;
}



