#include "TH1.h"
#include "TH1F.h"

#include <string>


int dolog=0;
void figefff() 
{ 

  TString canvName = "Fig_";
  canvName += "efff";


  std::string str1 = "electron relativist fraction";
  const char* atitle = str1.c_str();

  std::string str2="hefff";
  const char* hname1 =str2.c_str();

  std::string str3="PbWO";
  const char* lgd1 = str3.c_str();
  std::string  str4="Fiber1";
  const char* lgd2 = str4.c_str();
  std::string  str5="SampL";
  const char* lgd3 = str5.c_str();
  std::string  str6="Fiber2";
  const char* lgd4 = str6.c_str();
  std::string  str7="SampS";
  const char* lgd5 = str7.c_str();

  std::string str8 = " ";
  const char* htitle = str8.c_str();


  TFile *f1 = new TFile("hists_20GeV_BigEcal2.root");
  TFile *f2 = new TFile("hists_20GeV_FSCEPonly.root");
  TFile *f3 = new TFile("hists_20GeV_SampOnly.root");
  TFile *f4 = new TFile("hists_20GeV_FSCEPSAonly.root");
  TFile *f5 = new TFile("hists_20GeV_SampOnly2.root");


 
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
  
  if (dolog) canv->SetLogy();


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


  //A_pt->Rebin(2);
  //  int nbin;
  //nbin=A_pt->FindBin(0.6);
  //    A_pt->GetXaxis()->SetRangeUser(0,nbin);
  A_pt->GetXaxis()->SetRangeUser(0.5,1.0);
  A_pt->SetDirectory(0);
  A_pt->SetTitle(htitle);
  double aaA = A_pt->Integral();
std::cout<<" first entries is "<<aaA<<std::endl;
  A_pt->Scale(1./aaA);

  
  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt = static_cast<TH1F*>(f2->Get(hname1)->Clone());
  //nbin=B_pt->FindBin(0.6);
  //B_pt->GetXaxis()->SetRangeUser(0,nbin);
  B_pt->GetXaxis()->SetRangeUser(0.5,1.0);
    //B_pt->Rebin(2);
  std::cout<<"ha"<<std::endl;
  B_pt->SetDirectory(0);
  double aaB = B_pt->Integral();
std::cout<<" second entries is "<<aaB<<std::endl;
  B_pt->Scale(1/aaB);

  
  std::cout<<"getting third"<<std::endl;
  TH1F *C_pt = static_cast<TH1F*>(f3->Get(hname1)->Clone());
  //nbin=C_pt->FindBin(0.6);
  //C_pt->GetXaxis()->SetRangeUser(0,nbin);
  C_pt->GetXaxis()->SetRangeUser(0.5,1.0);
    //C_pt->Rebin(2);
  std::cout<<"ha"<<std::endl;
  C_pt->SetDirectory(0);
  double aaC = C_pt->Integral();
std::cout<<" third entries is "<<aaC<<std::endl;
  C_pt->Scale(1/aaC);
  

  std::cout<<"getting fourth"<<std::endl;
  TH1F *D_pt = static_cast<TH1F*>(f4->Get(hname1)->Clone());
  //nbin=D_pt->FindBin(0.6);
  //D_pt->GetXaxis()->SetRangeUser(0,nbin);
  D_pt->GetXaxis()->SetRangeUser(0.5,1.0);
    //B_pt->Rebin(2);
  std::cout<<"ha"<<std::endl;
  D_pt->SetDirectory(0);
  double aaD = D_pt->Integral();
std::cout<<" fourth entries is "<<aaD<<std::endl;
  D_pt->Scale(1/aaD);

  
  std::cout<<"getting fifth"<<std::endl;
  TH1F *E_pt = static_cast<TH1F*>(f5->Get(hname1)->Clone());
  //nbin=E_pt->FindBin(0.6);
  //E_pt->GetXaxis()->SetRangeUser(0,nbin);
  E_pt->GetXaxis()->SetRangeUser(0.5,1.0);
    //E_pt->Rebin(2);
  std::cout<<"ha"<<std::endl;
  E_pt->SetDirectory(0);
  double aaE = E_pt->Integral();
std::cout<<" fift entries is "<<aaE<<std::endl;
  E_pt->Scale(1/aaE);
  

  double max = std::max(A_pt->GetMaximum(),B_pt->GetMaximum());
  //  max = std::max(max,C_pt->GetMaximum());
  A_pt->SetMaximum(max*1.3);

  A_pt->GetYaxis()->SetTitle(" percent  ");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  
  A_pt->GetXaxis()->SetTitle(atitle);  
  A_pt->GetXaxis()->SetTitleSize(0.05);  



  A_pt->SetLineColor(1);
  A_pt->SetLineWidth(3);
  A_pt->SetStats(0);
  A_pt->Draw("HIST ");

  

  B_pt->SetLineColor(2);
  B_pt->SetLineWidth(3);
  B_pt->SetStats(0);
  B_pt->Draw("HIST same");

  
  C_pt->SetLineColor(3);
  C_pt->SetLineWidth(3);
  C_pt->SetStats(0);
  C_pt->Draw("HIST same");
  
  
  D_pt->SetLineColor(4);
  D_pt->SetLineWidth(4);
  D_pt->SetStats(0);
  D_pt->Draw("HIST same");
  
  
  E_pt->SetLineColor(5);
  E_pt->SetLineWidth(5);
  E_pt->SetStats(0);
  E_pt->Draw("HIST same");
  
  


  lgd->AddEntry(A_pt, lgd1, "l");
  lgd->AddEntry(B_pt, lgd2, "l");
  lgd->AddEntry(C_pt, lgd3, "l");
  lgd->AddEntry(D_pt, lgd4, "l");
  lgd->AddEntry(E_pt, lgd5, "l");
  //lgd->AddEntry(C_pt, "ModelBx500", "l");

 lgd->Draw();
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

 
  if (dolog) {
    canv->Print(canvName+"_log.pdf",".pdf");
    canv->Print(canvName+"_log.png",".png");}
  else{ 
    canv->Print(canvName+".pdf",".pdf");
    canv->Print(canvName+".png",".png");}


  

  return;
}



