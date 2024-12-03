#include "TH1.h"
#include "TH1F.h"

#include <string>


int dolog=0;
void figg1g2SampS() 
{ 
  TString canvName = "Fig_";
  canvName += "g1g2SampS";


  std::string strt = "Sampling Fraction";
  const char* atitle = strt.c_str();

  std::string strn1="ehcHcalE1";
  const char* hname1 =strn1.c_str();
  std::string strn2="phcHcalE1";
  const char* hname2 =strn2.c_str();
  std::string strn3="ehcHcalE2";
  const char* hname3 =strn3.c_str();
  std::string strn4="phcHcalE2";
  const char* hname4 =strn4.c_str();


  std::string str1="electron g1";
  const char* lgd1 = str1.c_str();
  std::string  str2="pion g1";
  const char* lgd2 = str2.c_str();

  std::string str3="electron g2";
  const char* lgd3 = str3.c_str();
  std::string  str4="pion g2";
  const char* lgd4 = str4.c_str();


  std::string str6 = " ";
  const char* htitle = str6.c_str();


  TFile *f1 = new TFile("hists_20GeV_SampOnly2.root");



 
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

  float arms,amean;

  std::cout<<"getting first"<<std::endl;
  TH1F *A_pt = static_cast<TH1F*>(f1->Get(hname1)->Clone());
  A_pt->SetDirectory(0);
  A_pt->SetTitle(htitle);
  double aaA = A_pt->Integral();
  std::cout<<" first entries is "<<aaA<<std::endl;
  A_pt->Scale(1./aaA);

  

  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt = static_cast<TH1F*>(f1->Get(hname2)->Clone());
  std::cout<<"ha"<<std::endl;
  B_pt->SetDirectory(0);
  double aaB = B_pt->Integral();
  std::cout<<" second entries is "<<aaB<<std::endl;
  B_pt->Scale(1/aaB);

  std::cout<<"getting third"<<std::endl;
  TH1F *C_pt = static_cast<TH1F*>(f1->Get(hname3)->Clone());
  std::cout<<"ha"<<std::endl;
  C_pt->SetDirectory(0);
  double aaC = C_pt->Integral();
  std::cout<<" third entries is "<<aaC<<std::endl;
  C_pt->Scale(1/aaC);

  std::cout<<"getting fourth"<<std::endl;
  TH1F *D_pt = static_cast<TH1F*>(f1->Get(hname4)->Clone());
  std::cout<<"ha"<<std::endl;
  D_pt->SetDirectory(0);
  double aaD = D_pt->Integral();
  std::cout<<" fourth entries is "<<aaD<<std::endl;
  D_pt->Scale(1/aaD);


  double max = std::max(A_pt->GetMaximum(),B_pt->GetMaximum());
  //  max = std::max(max,C_pt->GetMaximum());
  A_pt->SetMaximum(max*1.3);





 
  std::cout<<std::endl;
  std::cout<<"fitting first hist"<<std::endl;
  arms = A_pt->GetRMS();
  amean = A_pt->GetMean();
  std::cout<<"mean rms are "<<amean<<" "<<arms<<std::endl;
  A_pt->Fit("gaus","R0","",amean-5.5*arms,amean+5.5*arms);

  std::cout<<std::endl;
  std::cout<<"fitting second hist"<<std::endl;
  arms = B_pt->GetRMS();
  amean = B_pt->GetMean();
  std::cout<<"mean rms are "<<amean<<" "<<arms<<std::endl;
  B_pt->Fit("gaus","R0","",amean-5.5*arms,amean+5.5*arms);

  std::cout<<std::endl;
  std::cout<<"fitting third hist"<<std::endl;
  arms = C_pt->GetRMS();
  amean = C_pt->GetMean();
  std::cout<<"mean rms are "<<amean<<" "<<arms<<std::endl;
  C_pt->Fit("gaus","R0","",amean-5.5*arms,amean+5.5*arms);

  std::cout<<std::endl;
  std::cout<<"fitting fourth hist"<<std::endl;
  arms = D_pt->GetRMS();
  amean = D_pt->GetMean();
  std::cout<<"mean rms are "<<amean<<" "<<arms<<std::endl;
  D_pt->Fit("gaus","R0","",amean-5.5*arms,amean+5.5*arms);




  
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
  B_pt->Draw("HIST  same");


  C_pt->SetLineColor(3);
  C_pt->SetLineWidth(3);
  C_pt->SetStats(0);
  C_pt->Draw("HIST  same");


  D_pt->SetLineColor(4);
  D_pt->SetLineWidth(3);
  D_pt->SetStats(0);
  D_pt->Draw("HIST  same");

  


  lgd->AddEntry(A_pt, lgd1, "l");
  lgd->AddEntry(B_pt, lgd2, "l");
  lgd->AddEntry(C_pt, lgd3, "l");
  lgd->AddEntry(D_pt, lgd4, "l");

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



