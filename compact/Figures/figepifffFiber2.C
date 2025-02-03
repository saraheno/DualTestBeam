#include "TH1.h"
#include "TH1F.h"

#include <string>



void figepifffFiber2() 
{ 
  TString canvName = "Fig_";
  canvName += "epiffFiber2";

  std::string str1 = "Relativistic fraction E deposits, entire calorimeter";
  const char* atitle = str1.c_str();

  std::string strn1="hefff";
  const char* hname1 =strn1.c_str();

  std::string strn2="hpfff2";
  const char* hname2 =strn2.c_str();

  TFile *f1 = new TFile("hists_20GeV_FSCEPSAonly.root");

 
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
  

  TLatex latex;
  
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
  TH1F *A_pt = static_cast<TH1F*>(f1->Get(hname1)->Clone());
  double aaA = A_pt->Integral();
  std::cout<<" first entries is "<<aaA<<std::endl;
  A_pt->Scale(1./aaA);

  A_pt->GetYaxis()->SetTitle(" percent  ");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  
  A_pt->GetXaxis()->SetTitle(atitle);  
  A_pt->GetXaxis()->SetTitleSize(0.05);  
  A_pt->SetLineColor(1);
  A_pt->SetLineWidth(3);
  A_pt->SetStats(0);
  A_pt->Draw("HIST ");


  std::cout<<"getting first"<<std::endl;
  TH1F *B_pt = static_cast<TH1F*>(f1->Get(hname2)->Clone());
  double aaB = B_pt->Integral();
  std::cout<<" second entries is "<<aaB<<std::endl;
  B_pt->Scale(10./aaB);

  B_pt->GetYaxis()->SetTitle(" Cherenkov signal (arb. unit)  ");  
  B_pt->GetYaxis()->SetTitleSize(0.05);  
  B_pt->GetXaxis()->SetTitle(atitle);  
  B_pt->GetXaxis()->SetTitleSize(0.05);  
  B_pt->SetLineColor(2);
  B_pt->SetLineWidth(3);
  B_pt->SetStats(0);
  B_pt->Draw("HIST same");




    // Writing the lumi information and the CMS "logo"
   // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
   
  
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  lgd->Draw();


  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");

  return;
}



