#include "TH1.h"
#include "TH1F.h"

#include <string>


int dolog=0;
void figthreeresPbWOS() 
{ 
  TString canvName = "Fig_";
  canvName += "threeresPbWOS";


  std::string str1 = "Signal calibrated to electrons";
  const char* atitle = str1.c_str();

  std::string strn1="phcEcalncer";
  const char* hname1 =strn1.c_str();
  std::string strn2="phcEcalnscint";
  const char* hname2 =strn2.c_str();
  std::string strn3="phcEcalcorr";
  const char* hname3 =strn3.c_str();

  std::string str3="Cherenkov proxy energy deposit";
  const char* lgd1 = str3.c_str();
  std::string  str4="Scintillator proxy energy deposit";
  const char* lgd2 = str4.c_str();
  std::string  str5="Dual corrected proxy energy deposit";
  const char* lgd3 = str5.c_str();

  std::string str6="Cherenkov light at creation";
  const char* lgd4 = str6.c_str();
  std::string  str7="Scintillator light at creation";
  const char* lgd5 = str7.c_str();
  std::string  str8="Dual corrected light at creation";
  const char* lgd6 = str8.c_str();

  std::string strt = " ";
  const char* htitle = strt.c_str();


  TFile *f1 = new TFile("hists_20GeV_BigEcal2.root");
  TFile *f2 = new TFile("hists_20GeV_BigEcal2_g1.root");



 
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



  
  int n_ = 2;
  
  float x1_l = 0.5;
  //  float x1_l = 0.75;
  float y1_l = 0.90;
  
  float dx_l = 0.30;
  float dy_l = 0.16;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  
 TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(32); lgd->SetFillColor(0);



  std::cout<<"getting first"<<std::endl;
  TH1F *A_pt = static_cast<TH1F*>(f1->Get(hname1)->Clone());
  A_pt->Rebin(4);
  //A_pt->GetXaxis()->SetRangeUser(0.,0.6);
  A_pt->SetDirectory(0);
  A_pt->SetTitle(htitle);
  double aaA = A_pt->Integral();
  std::cout<<" first entries is "<<aaA<<std::endl;
  A_pt->Scale(1./aaA);

  std::cout<<"getting first g1"<<std::endl;
  TH1F *Ag1_pt = static_cast<TH1F*>(f2->Get(hname1)->Clone());
  Ag1_pt->Rebin(4);
  //Ag1_pt->GetXaxis()->SetRangeUser(0.,0.6);
  Ag1_pt->SetDirectory(0);
  Ag1_pt->SetTitle(htitle);
  double aaAg1 = Ag1_pt->Integral();
  std::cout<<" first g1 entries is "<<aaAg1<<std::endl;
  Ag1_pt->Scale(1./aaAg1);


  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt = static_cast<TH1F*>(f1->Get(hname2)->Clone());
  //B_pt->GetXaxis()->SetRangeUser(0.,0.6);
  B_pt->Rebin(4);
  B_pt->SetDirectory(0);
  double aaB = B_pt->Integral();
  std::cout<<" second entries is "<<aaB<<std::endl;
  B_pt->Scale(1/aaB);

  std::cout<<"getting second g1"<<std::endl;
  TH1F *Bg1_pt = static_cast<TH1F*>(f2->Get(hname2)->Clone());
  //Bg1_pt->GetXaxis()->SetRangeUser(0.,0.6);
  Bg1_pt->Rebin(4);
  Bg1_pt->SetDirectory(0);
  double aaBg1 = Bg1_pt->Integral();
  std::cout<<" second g1 entries is "<<aaBg1<<std::endl;
  Bg1_pt->Scale(1/aaBg1);

  
  std::cout<<"getting third"<<std::endl;
  TH1F *C_pt = static_cast<TH1F*>(f1->Get(hname3)->Clone());
  //C_pt->GetXaxis()->SetRangeUser(0.,0.6);
  C_pt->Rebin(4);
  std::cout<<"ha"<<std::endl;
  C_pt->SetDirectory(0);
  double aaC = C_pt->Integral();
std::cout<<" third entries is "<<aaC<<std::endl;
  C_pt->Scale(1/aaC);
  
  
  std::cout<<"getting third g1"<<std::endl;
  TH1F *Cg1_pt = static_cast<TH1F*>(f2->Get(hname3)->Clone());
  //Cg1_pt->GetXaxis()->SetRangeUser(0.,0.6);
  Cg1_pt->Rebin(4);
  Cg1_pt->SetDirectory(0);
  double aaCg1 = Cg1_pt->Integral();
  std::cout<<" third entries g1 is "<<aaCg1<<std::endl;
  Cg1_pt->Scale(1/aaCg1);
  

  double max = std::max(A_pt->GetMaximum(),B_pt->GetMaximum());
  max = std::max(max,C_pt->GetMaximum());
  A_pt->SetMaximum(max*1.3);

  A_pt->GetYaxis()->SetTitle(" percent  ");  
  A_pt->GetYaxis()->SetTitleSize(0.05);  
  A_pt->GetXaxis()->SetTitle(atitle);  
  A_pt->GetXaxis()->SetTitleSize(0.05);  



  A_pt->SetLineColor(1);
  A_pt->SetLineWidth(6);
  A_pt->SetLineStyle(1);
  A_pt->SetStats(0);
  A_pt->Draw("HIST ");

  Ag1_pt->SetLineColor(1);
  Ag1_pt->SetLineWidth(2);
  Ag1_pt->SetLineStyle(1);
  Ag1_pt->SetStats(0);
  Ag1_pt->Draw("HIST same");

  

  B_pt->SetLineColor(2);
  B_pt->SetLineWidth(6);
  B_pt->SetLineStyle(2);
  B_pt->SetStats(0);
  B_pt->Draw("HIST same");

  Bg1_pt->SetLineColor(2);
  Bg1_pt->SetLineWidth(2);
  Bg1_pt->SetLineStyle(2);
  Bg1_pt->SetStats(0);
  Bg1_pt->Draw("HIST same");

  
  C_pt->SetLineColor(3);
  C_pt->SetLineWidth(6);
  C_pt->SetLineStyle(3);
  C_pt->SetStats(0);
  C_pt->Draw("HIST same");
  
  Cg1_pt->SetLineColor(3);
  Cg1_pt->SetLineWidth(2);
  Cg1_pt->SetLineStyle(3);
  Cg1_pt->SetStats(0);
  Cg1_pt->Draw("HIST same");
  
  


  lgd->AddEntry(A_pt, lgd1, "l");
  lgd->AddEntry(B_pt, lgd2, "l");
  lgd->AddEntry(C_pt, lgd3, "l");
  lgd->AddEntry(Ag1_pt, lgd4, "l");
  lgd->AddEntry(Bg1_pt, lgd5, "l");
  lgd->AddEntry(Cg1_pt, lgd6, "l");
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

 
  if (dolog) {
    canv->Print(canvName+"_log.pdf",".pdf");
    canv->Print(canvName+"_log.png",".png");}
  else{ 
    canv->Print(canvName+".pdf",".pdf");
    canv->Print(canvName+".png",".png");}
  return;
}



