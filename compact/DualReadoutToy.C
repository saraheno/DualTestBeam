// force update
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <TGaxis.h>

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;

/*
	this code is generating random number 
	it is using the initial values from the DualResolutionPaper
	to run the code do: root -l -b -q 'DualReadoutToy.C'

*/
void SCEDraw1		(TCanvas* canv, const char* name, TH1F* h1, const char* outfile, bool logy);
void SCEDraw1_2D	(TCanvas* canv, const char* name, TH2F* h1, const char* outfile);
void SCEDraw1_2D_2 	(TCanvas* canv, const char* name, TH2F* h1, const char* outfile);
void setCanvas  	(TCanvas* canv,const char* name, bool hist, bool tprofile, TH1F* h1, TProfile* t1);

void dotoy(bool doplot,double h_s,double h_c,double nscint,double ncer,double fmean,double frms, double &meanS,double &meanC,double &sigmaS, double &sigmaC, double &sigmaD, double &acov);
TRandom rrr;

int nshowers = 100000;
bool doplots = 0;
bool verbose = 1;

TFile * out = new TFile("dualtoy.root","RECREATE");

void DualReadoutToy() {
	// values from the paper
	double h_s	=  0.9;
	double h_c	=  0.6;
	double nscint	=  10000.0;
	double ncer	=  10000.0;
	double fmean	=  0.6;
	double frms	=  0.2;

	double meanS, meanC, sigmaS, sigmaC, sigmaD, acov;

	dotoy(1, h_s, h_c, nscint, ncer, fmean, frms, meanS, meanC, sigmaS, sigmaC, sigmaD, acov);

	double precov = (1 - h_s) * (1 - h_c) * pow(frms,2); // eq. 8 of DualResolutionPaper

	cout<<"meanS: fit = "	<< meanS  << " , formula = " << fmean + (1 - fmean) * h_s << endl;
	cout<<"meanC: fit = "	<< meanC  << " , formula = " << fmean + (1 - fmean) * h_c << endl;
	cout<<"cov:   fit = "   << acov   << " , formula = " << precov << endl;
	cout<<"sigmaS     = "	<< sigmaS << endl;
	cout<<"sigmaC     = "	<< sigmaC << endl;
	cout<<"sigmaD     = "   << sigmaD << endl;

	//computing the 3 terms of the resolution formula eq 9
	double term1	= pow((1 - h_c) * sigmaS, 2); // scint term
	double term2	= pow((1 - h_s) * sigmaC, 2); // ceren term
	double sum12    = term1 + term2;

	double term3_true	= 2 * (1 - h_s) * (1 - h_c) * acov;   // from dotoy 
	double term3_formula	= 2 * (1 - h_s) * (1 - h_c) * precov; // from frms=0.2

	double dualpred		= (1 / (h_s - h_c)) * sqrt(sum12 - term3_true);    // resolution with fres = 0.2
	double dualpred2	= (1 / (h_s - h_c)) * sqrt(sum12 - term3_formula); // resolution with fit cov

	cout<<"DualToy = " << dualpred << ", DualFormula = " << dualpred2 << ", SigmaD= " << sigmaD << endl;

	//TGaxis::SetMaxDigits(2); // Set axis to scientific notation
	gStyle->SetOptStat(0);

	TH2F *covcheck		= new TH2F("covcheck",	  "cov: comp vs formula;cov(formula);cov(fit)",  100, 0.,  0.002,	100,  0.,  0.002);
	TH2F *covcheckf1	= new TH2F("covcheckf1",  "#sigma_{D} - #sigma_{formula} vs fres",100, 1.6, 0.22,	100, -0.01,0.0);
	TH2F *scintdual		= new TH2F("scintdual",	  "#sigma_{D} vs #sigma_{S}",	   100, 0.,  0.1,	100,  0.,  0.1);
	TH2F *scintvscer	= new TH2F("scintvscer",  "#sigma_{S} vs #sigma_{C};#sigma_{C};#sigma_{S}", 100, 0, 0.1, 100, 0, 0.1);

	TH1F *cerRes		= new TH1F("cerRes",	  "#sigma_{C}",		100, 0,	0.1);
	TH1F *scintRes          = new TH1F("scintRes",    "#sigma_{S}",		100, 0, 0.1);
	TH1F *dualRes           = new TH1F("dualRes",     "#sigma_{D}",		100, 0, 0.02);
	TH1F *dualFormula	= new TH1F("dualFormula", "#sigma_{formula}",	100, 0, 0.02);
	TH1F *dualTrue		= new TH1F("dualTrue",    "#sigma_{true}",	100, 0, 0.02);

	TH1F *dualcheck		= new TH1F("dualcheck",	  "#sigma_{D}-dualtrue",	100,-0.02,-0.01);
	TH1F *dualcheckf	= new TH1F("dualcheckf",  "#sigma_{D}-dualformula",	100,-0.02, 0.0);
	TH1F *dualcheckab	= new TH1F("dualcheckab", "dualformul-dualtrue",100,-0.01,0.02);

	int jmax = 0;
	int npts = 200;
	double range = min(fmean ,1 - fmean) / 2.;

	if (verbose) cout << "sigmaS,         sigmaC,        sigmaD,         term1,        term2,         term1+term2,   term3_true,  term3_formula" << endl;

	for(int j=1; j < npts; j++) { // computing dualres with 200 variations of fres
		double frestry = j * (range / npts);
		dotoy(0, h_s, h_c, nscint, ncer, fmean, frestry, meanS, meanC, sigmaS, sigmaC, sigmaD, acov);
		
		precov = (1 - h_s) * (1 - h_c) * pow(frestry,2); // eq. 8 of the paper
		covcheck->Fill(acov, precov);
		
		term1		= pow((1 - h_c) * sigmaS, 2);
		term2		= pow((1 - h_s) * sigmaC, 2);
		term3_true	= 2 * (1 - h_s) * (1 - h_c) * acov;
		term3_formula	= 2 * (1 - h_s) * (1 - h_c) * precov;
		
		if(term3_formula>term1+term2) {
			cout << "invalid prediction frestry = " << frestry << endl;
			jmax=j-1;
		}

		if(jmax!=0) term3_formula = 2 * pow((1 - h_s) * (1 - h_c) * jmax * (1./40.), 2); // factor (1./40.) to catch odd values; jmax should not be == 0

		if (verbose){
			cout << scientific << setprecision(6);
			cout << sigmaS     << setw(15) << sigmaC << setw(15) << sigmaD << setw(15);
			cout << term1      << setw(15) << term2 << setw(15) << sum12 << setw(15);
			cout << term3_true << setw(15) << term3_formula << endl;
		}

		double dualpreda = (1 / (h_s - h_c)) * sqrt(term1 + term2 - term3_true);
		double dualpredb = (1 / (h_s - h_c)) * sqrt(term1 + term2 - term3_formula);

		cerRes->Fill(sigmaC);
		scintRes->Fill(sigmaS);
		dualRes->Fill(sigmaD);
		dualFormula->Fill(dualpredb);
		dualTrue->Fill(dualpreda);

		dualcheck->Fill(sigmaD-dualpreda);
		dualcheckf->Fill(sigmaD-dualpredb);
		dualcheckab->Fill(dualpreda-dualpredb);
		scintdual->Fill(sigmaS,sigmaD);
		scintvscer->Fill(sigmaC,sigmaS);
		covcheckf1->Fill(frestry,sigmaD-dualpredb);
	}

	int ncanv = 9;
	vector<TCanvas*> canv;
	for (int i = 0; i < ncanv; ++i) {canv.push_back(new TCanvas(Form("canvas%d", i), Form("Canvas %d", i), 800, 600));}

	if (doplots){
		SCEDraw1_2D	(canv[0],	"c7", 	covcheck,	"covcheck.png");
		SCEDraw1_2D	(canv[1],	"c7a",	covcheckf1,	"covcheckf1.png");
		SCEDraw1_2D_2	(canv[2],	"c10",	scintdual,	"scintdual.png");
		SCEDraw1	(canv[3],	"c8", 	dualcheck,	"dualcheck.png",  0);
		SCEDraw1	(canv[4],	"c9", 	dualcheckf,	"dualcheckf.png", 0);
		SCEDraw1	(canv[5],	"c91",	dualcheckab,	"dualcheckab.png",0);
	}

	TH2F *scintdual2	= new TH2F("scintdual2", "toys #sigma_{D} vs #sigma_{S}",	100,0.,0.1,100,0.,0.1);
	TH2F *dualcheckha	= new TH2F("dualcheckha", "toys dual: true vs formula",	100,0.,0.1,100,0.,0.1);

	for(int j=1; j < 500; j++) { // compute resolution by varying nscint=ncer ranged from 100 to 50000 (= 100 * 500)
		double nnn = 100. * j;
		dotoy(0, h_s, h_c, nnn, nnn, fmean, frms, meanS, meanC, sigmaS, sigmaC, sigmaD, acov);
		scintdual2->Fill(sigmaS,sigmaD);
		
		precov		= (1 - h_s) * (1 - h_c) * pow((frms * fmean) ,2)/ (pow(frms, 2) + pow(fmean, 2));
		term1		= (1 - h_c) * (1 - h_c) * pow(sigmaS, 2);
		term2		= (1 - h_s) * (1 - h_s) * pow(sigmaC, 2);
		term3_true	= 2 * pow((1 - h_s) * (1 - h_c) * acov, 2);
		term3_formula	= 2 * pow((1 - h_s) * (1 - h_c) * precov, 2);

		double dualpreda = (1/(h_s-h_c)) * sqrt(term1 + term2 - term3_true);
		double dualpredb = (1/(h_s-h_c)) * sqrt(term1 + term2 - term3_formula);
		dualcheckha->Fill(dualpreda,dualpredb);
	}
	cout<<"last dualres is "<<sigmaD<<endl;

	if (doplots){
		SCEDraw1_2D_2(canv[6],	"toys: scint vs dual res",	scintdual2,	"scintdual2.png");
		SCEDraw1_2D_2(canv[7],	"toys: dual true vs formula",	dualcheckha,	"dualcheckha.png");
	}
	
	TH2F *scintdual3 = new TH2F("scintdual3", "#sigma_{D} vs meanC", 100,0.,0.1,100,0.,0.1);

	for(int j=1;j<10;j++) {
		double fmeanaa = 0.2 + 0.05 * j;
		dotoy(0, h_s, h_c, nscint, ncer, fmeanaa, frms, meanS, meanC, sigmaS, sigmaC, sigmaD, acov);
		scintdual3->Fill(meanC, sigmaD);
	}
	if (doplots) SCEDraw1_2D(canv[8],"c12",scintdual3,"cer_vs_sigmaD.png");
	covcheck->Write();
        covcheckf1->Write();
        scintdual->Write();
        dualcheck->Write();
        dualcheckf->Write();
        dualcheckab->Write();
	scintdual2->Write();
        dualcheckha->Write();
	cerRes->Write();
        scintRes->Write();
        dualRes->Write();
	scintdual3->Write();
	dualFormula->Write();
	dualTrue->Write();
	scintvscer->Write();
	out->Close();
}

void dotoy(bool doplot, double h_s,double h_c,double nscint,double ncer,double fmean,double frms, double &meanS, double& meanC, double &sigmaS, double &sigmaC, double &sigmaD, double &acov) {
	TH1F *fff  = new TH1F("fff",	"shower em fraction",	300,	0.,2.0);
	TH1F *sss  = new TH1F("sss",	"shower scintillation",	300,	0.,2.0);
	TH1F *ccc  = new TH1F("ccc",	"shower cherenkov",	300,	0.,2.0);
	TH1F *ddd  = new TH1F("ddd",	"dual readout",		900,	0.,2.0);
	TH1F *cov  = new TH1F("cov",	"covariance",		3000,  -2.,2.0);
	TH2F *sscc = new TH2F("sscc",   "Cerenkov versus scint",  1000,   0.,2.0,1000,0.,2.0);

	fff->Reset();
	sss->Reset();
	ccc->Reset();
	sscc->Reset();
	ddd->Reset();
	cov->Reset();

	acov=0.;

	// two last terms to compute the cov(S,C) in formula 7 of DualResolutionPaper
	double pmeans = fmean + h_s * (1 - fmean); // == <S>
	double pmeanc = fmean + h_c * (1 - fmean); // == <C>

	for (int i=0; i < nshowers; i++) {
		double FFF = -1.; // pick a value for fraction EM in the shower
		while((FFF<0)||(FFF>1) ) {
			FFF = rrr.Gaus(fmean,frms); // generate single random number following the Gaus distr with mean=fmean and rms=frms
		}
		fff->Fill(FFF);

		double SSS = FFF + h_s * (1-FFF);
		double CCC = FFF + h_c * (1-FFF);

		SSS = SSS * (1 + rrr.Gaus(0.,1/sqrt(nscint)));
		CCC = CCC * (1 + rrr.Gaus(0.,1/sqrt(ncer)));

		double DDD = ((1 - h_c) * SSS - (1 - h_s) * CCC)/(h_s - h_c); // energy estimate eq. 3 of the paper

		sss->Fill(SSS);
		ccc->Fill(CCC);
		ddd->Fill(DDD);
		sscc->Fill(SSS,CCC);

		cov->Fill((SSS - pmeans) * (CCC - (fmean - pmeanc)));
		acov+= (SSS - pmeans) * (CCC - pmeanc);
	}

	acov	=	acov / (nshowers-1);
	sigmaS	=	sss->GetRMS();
	sigmaC	=	ccc->GetRMS();
	sigmaD	=	ddd->GetRMS();
	meanS	=	sss->GetMean();
	meanC	=	ccc->GetMean();

	
	if(doplot) {
		TCanvas* c1;
		SCEDraw1(c1,"EM Shower",fff,"emshower.png",0);
		TCanvas* c2;
		SCEDraw1(c2,"Scint Shower",sss,"scintshower.png",0);
		TCanvas* c3;
		SCEDraw1(c3,"Cer Shower",ccc,"cershower.png",0);
		TCanvas* c4;
		SCEDraw1_2D(c4,"Scint vs Cer",sscc,"scint_vs_cer.png");
		TCanvas* c5;
		SCEDraw1(c5,"Dual Readout",ddd,"dualreadout.png",0);
		TCanvas* c6;
		SCEDraw1(c6,"Covariance",cov,"cov.png",0);
	}
	delete fff;
	delete sss;
	delete ccc;
	delete sscc;
	delete ddd;
	delete cov;
	return;
}

void setCanvas(TCanvas* canv,const char* name) {
        canv= new TCanvas(name,name,200,10,700,500);
        canv->SetFillColor(0);
        canv->SetBorderMode(0);
        canv->SetFrameFillStyle(0);
        canv->SetFrameBorderMode(0);
        canv->SetTickx(0);
        canv->SetTicky(0);
        return;
}

void SCEDraw1 (TCanvas* canv,  const char* name,TH1F* h1, const char* outfile, bool logy) {
	setCanvas(canv, name);
	if(logy) canv->SetLogy();
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(1);
	h1->SetStats(111111);
	h1->Draw("HIST");	
	canv->Print(outfile,".png");
	return;
}

void SCEDraw1_2D (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile) {
	setCanvas(canv, name);
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(kGreen);
	h1->SetStats(111111);
	h1->Draw("");	
	canv->Print(outfile,".png");
	return;
}

void SCEDraw1_2D_2 (TCanvas* canv,  const char* name,TH2F* h1, const char* outfile) {
	setCanvas(canv, name);
	h1->SetLineColor(kGreen);
	h1->SetLineWidth(kGreen);
	h1->SetStats(111111);
	h1->Draw("");
	TLine line = TLine(0.,0.,1.,1.);
	line.SetLineColor(kYellow);
	line.SetLineWidth(2);
	line.Draw();
	canv->cd(0);
	canv->Modified();
	canv->Update();
	canv->Print(outfile,".png");
	return;
}
