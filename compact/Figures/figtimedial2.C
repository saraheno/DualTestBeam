#include "TH1.h"
#include "TH1F.h"

#include <string>


int dolog=0;
void figtimedial2(const char* filename1, const char* filename2,float nevents, float scales1, float scalec1, float scales2, float scalec2) 
{ 
  TString canvName = "Fig_";
  canvName += "timecomp";


  std::string str1 = "Calibrated signal";
  const char* atitle = str1.c_str();

  std::string strn1="eeecaltimes";
  const char* hname1 =strn1.c_str();
  std::string strn2="eeecaltimec";
  const char* hname2 =strn2.c_str();



  std::string str3="times electron scint";
  const char* lgd1 = str3.c_str();
  std::string  str4="times electron cer";
  const char* lgd2 = str4.c_str();


  
  std::string strt = " ";
  const char* htitle = strt.c_str();


  TFile *f1 = new TFile(filename1);
  
  std::cout<<"getting first"<<std::endl;
  TH1F *A_pt1 = static_cast<TH1F*>(f1->Get(hname1)->Clone());
  A_pt1->Rebin(4);
  //A_pt1->GetXaxis()->SetRangeUser(0.,0.6);
  A_pt1->SetDirectory(0);
  A_pt1->SetTitle(htitle);
  double aaA1 = A_pt1->Integral();
  std::cout<<" first entries is "<<aaA1<<std::endl;
  //A_pt1->Scale(1./aaA1);
  A_pt1->Scale(1./scales1/nevents);

  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt1 = static_cast<TH1F*>(f1->Get(hname2)->Clone());
  //B1_pt->GetXaxis()->SetRangeUser(0.,0.6);
  B_pt1->Rebin(4);
  B_pt1->SetDirectory(0);
  double aaB1 = B_pt1->Integral();
  std::cout<<" second entries is "<<aaB1<<std::endl;
  //B_pt1->Scale(1/aaB1);
  B_pt1->Scale(1./scalec1/nevents);
  
  f1->Close();


  TFile *f2 = new TFile(filename2);
  
  std::cout<<"getting first"<<std::endl;
  TH1F *A_pt2 = static_cast<TH1F*>(f2->Get(hname1)->Clone());
  A_pt2->Rebin(4);
  //A_pt2->GetXaxis()->SetRangeUser(0.,0.6);
  A_pt2->SetDirectory(0);
  A_pt2->SetTitle(htitle);
  double aaA2 = A_pt2->Integral();
  std::cout<<" first entries is "<<aaA2<<std::endl;
  //A_pt2->Scale(1./aaA2);
  A_pt2->Scale(1./scales2/nevents);

  std::cout<<"getting second"<<std::endl;
  TH1F *B_pt2 = static_cast<TH1F*>(f2->Get(hname2)->Clone());
  //B_pt->GetXaxis()->SetRangeUser(0.,0.6);
  B_pt2->Rebin(4);
  B_pt2->SetDirectory(0);
  double aaB2 = B_pt2->Integral();
  std::cout<<" second entries is "<<aaB2<<std::endl;
  //B_pt2->Scale(1/aaB2);
  B_pt2->Scale(1./scalec2/nevents);
  
  f2->Close();

 
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



  TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1.0);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();

  

  double max = std::max(A_pt1->GetMaximum(),A_pt2->GetMaximum());
  A_pt1->SetMaximum(max*1.3);
  A_pt1->GetYaxis()->SetTitle(" percent  ");  
  A_pt1->GetYaxis()->SetTitleSize(0.05);  
  A_pt1->GetXaxis()->SetTitle(atitle);  
  A_pt1->GetXaxis()->SetTitleSize(0.05);  
  A_pt1->SetLineColor(1);
  A_pt1->SetLineWidth(6);
  A_pt1->SetLineStyle(1);
  A_pt1->SetStats(0);
  A_pt1->Draw("HIST ");
  A_pt2->SetLineColor(2);
  A_pt2->SetLineWidth(6);
  A_pt2->SetLineStyle(2);
  A_pt2->SetStats(0);
  A_pt2->Draw("HIST same");
  //lgd->AddEntry(A_pt1, lgd1, "l");
  //lgd->AddEntry(A_pt2, lgd1, "l");
  //lgd->AddEntry(C_pt, "ModelBx500", "l");
  //lgd->Draw();

 
    // Writing the lumi information and the CMS "logo"
   // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

 canv->cd();
 TPad *pad2 = new TPad("pad2","pad2",0.,0.,1.,0.5);
 pad2->SetTopMargin(0.0);
 pad2->SetBottomMargin(0.3);
 pad2->Draw();
 pad2->cd();




  double max2 = std::max(B_pt1->GetMaximum(),B_pt2->GetMaximum());
  B_pt1->SetMaximum(max2*1.3);
  B_pt1->GetYaxis()->SetTitle(" percent  ");  
  B_pt1->GetYaxis()->SetTitleSize(0.05);  
  B_pt1->GetXaxis()->SetTitle(atitle);  
  B_pt1->GetXaxis()->SetTitleSize(0.05);  
  B_pt1->SetLineColor(1);
  B_pt1->SetLineWidth(6);
  B_pt1->SetLineStyle(1);
  B_pt1->SetStats(0);
  B_pt1->Draw("HIST ");
  B_pt2->SetLineColor(2);
  B_pt2->SetLineWidth(6);
  B_pt2->SetLineStyle(2);
  B_pt2->SetStats(0);
  B_pt2->Draw("HIST same");
  //lgd->AddEntry(A_pt1, lgd1, "l");
  //lgd->AddEntry(A_pt2, lgd1, "l");
  //lgd->AddEntry(C_pt, "ModelBx500", "l");
  //lgd->Draw();

 
  
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  lgd->Draw();

   float t = canv->GetTopMargin();
  float r = canv->GetRightMargin();
  float Offset   = 0.2;
  TString alabel="20 GeV simulation";
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



