#include <vector>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/bin.h"

using namespace std;
void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, double jetR, const char* applyCorr, bool applyJEC);
void getUnfoldingInput(const char* workDirName="DataFits", const char* rapRegion="", const char* DSTag="DATA", const char* fitType="ctauMass", bool wantPureSMC=0, double jetR=0.3, const char* applyCorr="AccEff", bool applyJEC=1, bool statErr=1);
void getUnfoldingInput_all(bool statErr=1);
void plot18012Comp(bool isPr);
double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax, int cntmin, int cntmax);
void plotMCMassPars(const char* workDirName,
		    const char* rapRegion,
		    const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
		    const char* fitType, // "mass", "ctau"...
		    bool wantPureSMC, // =false,
		    double jetR,
		    const char* applyCorr, // = "",
		    bool applyJEC // =false
		    )
{
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  //double zbins1624 [] = {0.16, 0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};

  int nzbins = 0;
  double zedmin, zedmax;
  if (strcmp(rapRegion,"1624")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else {
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }
  TCanvas* c = new TCanvas ("c","",1000,800);

  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/MCPars",workDirName, rapRegion, fitType, DSTag));
  ofstream fileOut (Form("Output/%s/DataFits_%s/%s/%s/MCPars/ParsEvol.txt",workDirName, rapRegion, fitType, DSTag));
  fileOut << "zmin  zmax  n  nerr  alpha  alphaerr"<<endl;

  TH1F* aevol = new TH1F ("aevol",";z(J/#psi);alpha",nzbins, (strcmp(rapRegion,"1624")?zbins016:zbins1624));
  TH1F* nevol = new TH1F ("nevol",";z(J/#psi);n", nzbins, (strcmp(rapRegion,"1624")?zbins016:zbins1624));
  //TH1F* aevolnw = new TH1F ("aevolnw",";z(J/#psi);alpha",10, 0, 1);
  //TH1F* nevolnw = new TH1F ("nevolnw",";z(J/#psi);n", 10, 0, 1);

  TH1F* atot = new TH1F ("atot",";z(J/#psi);alpha",1000, zedmin, zedmin+0.03);
  TH1F* ntot = new TH1F ("ntot",";z(J/#psi);n", 1000, zedmin, zedmin+0.03);
  //TH1F* atotnw = new TH1F ("atotnw",";z(J/#psi);alpha",10, 0, 0.00001);
  //TH1F* ntotnw = new TH1F ("ntotnw",";z(J/#psi);n", 10, 0, 0.00001);


  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, jetR, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  TString ShapeTag = "";
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float n, n_errL,n_errH;
  float alpha, alpha_errL,alpha_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("N_Jpsi_val",&val);
  tr->SetBranchAddress("N_Jpsi_errL",&errL);
  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("n_Jpsi_val",&n);
  tr->SetBranchAddress("n_Jpsi_errL",&n_errL);
  tr->SetBranchAddress("n_Jpsi_errH",&n_errH);
  tr->SetBranchAddress("alpha_Jpsi_val",&alpha);
  tr->SetBranchAddress("alpha_Jpsi_errL",&alpha_errL);
  tr->SetBranchAddress("alpha_Jpsi_errH",&alpha_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);

  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    if (zmin < zedmin+0.02 && zmax == 1){
      atot->SetBinContent(atot->FindBin(zmin+0.0001),alpha);
      atot->SetBinError(atot->FindBin(zmin+0.0001),alpha_errL);
      ntot->SetBinContent(ntot->FindBin(zmin+0.0001),n);
      ntot->SetBinError(ntot->FindBin(zmin+0.0001),n_errL);
      fileOut<< zmin << "  " << zmax << "  " <<n<<"  "<<n_errL<<"  "<<alpha<<"  "<<alpha_errL<<endl;
      ShapeTag = jpsiName;
    }
    else
      { 
	avr=avr+alpha;
	avrn=avrn+n;
	tot++;
	aevol->SetBinContent(aevol->FindBin(zmin+0.02),alpha);
	aevol->SetBinError(aevol->FindBin(zmin+0.02),alpha_errL);
	nevol->SetBinContent(nevol->FindBin(zmin+0.02),n);
	nevol->SetBinError(nevol->FindBin(zmin+0.02),n_errL);

	fileOut<< zmin << "  " << zmax << "  " <<n<<"  "<<n_errL<<"  "<<alpha<<"  "<<alpha_errL<<endl;
      }
  }
  fileOut.close();
  TLegend* leg = NULL;
  TPaveText* tbox = NULL;
  TLine * ave = NULL;
  //TLine * avenw = NULL;

  aevol->SetMinimum(0);
  aevol->SetMaximum(4);
  nevol->SetMinimum(0);
  nevol->SetMaximum(4);
  atot->SetMinimum(0);
  atot->SetMaximum(4);
  ntot->SetMinimum(0);
  ntot->SetMaximum(4);

  avr=avr/tot;
  ave= new TLine(zedmin,avr,1,avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);

  leg = new TLegend(0.2, 0.6, 0.4, 0.8);
  leg->AddEntry(aevol, "diff.", "lep");
  leg->AddEntry(atot, "int.", "lep");
  leg->AddEntry(ave, "average", "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  tbox = new TPaveText(0.2,0.8,0.4,0.9, "BRNDC");
  tbox->AddText(strcmp(DSTag,"MCJPSINOPR")?(strcmp(DSTag,"MCJPSIPR")?"Data":"Prompt J/#psi MC"):"Nonprompt J/#psi MC");
  tbox->AddText(Form("%.1f < |y| < %.1f", ymin, ymax));
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);


  c->cd();
  aevol->Draw();
  ave->Draw("same");
  atot->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));

  avrn=avrn/tot;
  ave = new TLine(zedmin,avrn,1, avrn);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                       

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);

  c->cd();
  nevol->Draw();
  ntot->Draw("same E1");
  ave->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  f->Close();
  delete f;delete atot; delete ntot; delete aevol; delete nevol; delete c; delete leg; delete tbox; delete ave;
}


void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, double jetR, const char* applyCorr, bool applyJEC) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};
  double zbins024 [] = {0.064, 0.220, 0.376, 0.532, 0.688, 0.844, 1.000};

  string binTag = workDirName;
  double jtptmin;
  double jtptmax;
  if (binTag.find("midJtPt")!=std::string::npos) {binTag = "midJtPt"; jtptmin =30; jtptmax=40;}
  else if (binTag.find("lowJtPt")!=std::string::npos) {binTag = "lowJtPt"; jtptmin =20; jtptmax=30;}
  else if (binTag.find("lowerJtPt")!=std::string::npos) {binTag = "lowerJtPt"; jtptmin =10; jtptmax=20;}
  else if (binTag.find("highJtPt")!=std::string::npos) {binTag = "highJtPt"; jtptmin =40; jtptmax=50;}
  else if (binTag.find("higherJtPt")!=std::string::npos) {binTag = "higherJtPt"; jtptmin =50; jtptmax=60;}

  int nzbins = 0;
  double zedmin, zedmax;
  if (!strcmp(rapRegion,"016")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else if (!strcmp(rapRegion,"1624")){
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }
  else {
    nzbins = sizeof(zbins024)/sizeof(double)-1;
    zedmin = zbins024[0];
    zedmax = zbins024[nzbins];
  }

  TCanvas* c = new TCanvas ("c","",900,900);
  TH1F* bkgOrd_pp = NULL;
  TH1F* bkgOrd_pbpb = NULL;
  //if (strcmp(rapRegion,"1624"))
  bkgOrd_pp = new TH1F ("bkgOrd_pp",";z(J/#psi);background order", nzbins, zedmin, zedmax);
  bkgOrd_pbpb = new TH1F ("bkgOrd_pbpb",";z(J/#psi);background order", nzbins, zedmin, zedmax);
  //else 
  //bkgOrd = new TH1F ("bkgOrd",";z(J/#psi);background order", 5, 0, 1);

  string bkgPol [] = {"Uniform", "Chebychev1", "Chebychev2", "Chebychev3", "Chebychev4", "Chebychev5", "Chebychev6"};
  string bkgExp [] = {"Uniform", "ExpChebychev1", "ExpChebychev2", "ExpChebychev3", "ExpChebychev4", "ExpChebychev5", "ExpChebychev6"};
  gSystem->mkdir(Form("Output/%s/%s/%s/fitsPars",workDirName, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/%s/%s/result/tree_allvars.root",workDirName, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(workDirName, DSTag,"", fitType, wantPureSMC, jetR, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  TString ShapeTag = "";
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float n, n_errL,n_errH;
  float alpha, alpha_errL,alpha_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("N_Jpsi_val",&val);
  tr->SetBranchAddress("N_Jpsi_errL",&errL);
  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("n_Jpsi_val",&n);
  tr->SetBranchAddress("n_Jpsi_errL",&n_errL);
  tr->SetBranchAddress("n_Jpsi_errH",&n_errH);
  tr->SetBranchAddress("alpha_Jpsi_val",&alpha);
  tr->SetBranchAddress("alpha_Jpsi_errL",&alpha_errL);
  tr->SetBranchAddress("alpha_Jpsi_errH",&alpha_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);
  int ord = 0;
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    for (int j = 0; j<7; j++){
      if (bkgName == bkgPol[j]) { ord = j; ShapeTag = "PolChebychev"; break;}
      else if (bkgName == bkgExp [j]) {ord = j; ShapeTag = "ExpChebychev"; break;}
    }
    cout<<"[INFO] z: "<<zmin<<"-"<<zmax<< "bkg: " << bkgName <<endl;
    //if (!(zmin < zedmin+0.02 && zmax == 1)){
    if (TString(collSystem)=="PP") {
      bkgOrd_pp->SetBinContent(bkgOrd_pp->FindBin(zmin+0.001),ord);
      bkgOrd_pp->SetBinError(bkgOrd_pp->FindBin(zmin+0.001), 0.0001);
    }
    else {
      bkgOrd_pbpb->SetBinContent(bkgOrd_pbpb->FindBin(zmin+0.001),ord);
      bkgOrd_pbpb->SetBinError(bkgOrd_pbpb->FindBin(zmin+0.001), 0.0001);
    }
      //}
  }

    TLatex *  text2 = new TLatex(0.15 ,0.8,strcmp(ShapeTag,"PolChebychev")?"Exp. Chebychev bkg.":"Pol. Chebychev bkg.");
    text2->SetNDC();
    text2->SetTextFont(42);
    text2->SetTextSize(0.05);
    text2->SetLineWidth(2);

    TLatex *  text3 = new TLatex(0.15 ,0.75, Form("%.0f < p_{T,jet} < %.0f GeV", jtptmin, jtptmax));
    text3->SetNDC();
    text3->SetTextFont(42);
    text3->SetTextSize(0.05);
    text3->SetLineWidth(2);

    TLatex *  text = new TLatex(0.73 ,0.83,"CMS");
    text->SetNDC();
    text->SetTextFont(42);
    text->SetTextSize(0.06708595);
    text->SetLineWidth(5);

    TLatex *  text1 = new TLatex(0.66 ,0.77,"Preliminary");
    text1->SetNDC();
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->SetLineWidth(2);

    TLatex *  text4 = new TLatex(0.51 ,0.91,"pp 311 pb^{-1} (5.02 TeV)");
    text4->SetNDC();
    text4->SetTextFont(42);
    text4->SetTextSize(0.04);
    text4->SetLineWidth(2);

    TLatex *  text5 = new TLatex(0.51,0.91,"pbpb 1.7 nb^{-1} (5.02 TeV)");
    text5->SetNDC();
    text5->SetTextFont(42);
    text5->SetTextSize(0.04);
    text5->SetLineWidth(2);


    bkgOrd_pp->GetYaxis()->SetRangeUser(0, 6);
    bkgOrd_pp->SetMarkerColor(kMagenta+3);
    bkgOrd_pp->SetMarkerStyle(33);
    bkgOrd_pp->SetMarkerSize(3);
    bkgOrd_pp->SetLineColor(kMagenta+2);
    //bkgOrd_pp->SetOption("E1");

    bkgOrd_pbpb->GetYaxis()->SetRangeUser(0, 6);
    bkgOrd_pbpb->SetMarkerColor(kGreen+3);
    bkgOrd_pbpb->SetMarkerStyle(33);
    bkgOrd_pbpb->SetMarkerSize(3);
    bkgOrd_pbpb->SetLineColor(kGreen+2);
    //bkgOrd_pp->SetOption("E1");

    c->cd();
    bkgOrd_pp->Draw("EP");
    text->Draw("same");
    text1->Draw("same");
    text2->Draw("same");
    text3->Draw("same");
    text4->Draw("same");
    //c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/%s/%s/fitsPars/bkgOrder_%s_%s_pp.pdf",workDirName, fitType, DSTag, ShapeTag.Data(), binTag.c_str()));
    //c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));

    c->cd();
    bkgOrd_pbpb->Draw("EP");
    text->Draw("same");
    text1->Draw("same");
    text2->Draw("same");
    text3->Draw("same");
    text5->Draw("same");
    //c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/%s/%s/fitsPars/bkgOrder_%s_%s_pbpb.pdf",workDirName, fitType, DSTag, ShapeTag.Data(), binTag.c_str()));
    //c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    f->Close();
    delete f; delete bkgOrd_pp; delete bkgOrd_pbpb;
}


void getUnfoldingInput(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, double jetR, const char* applyCorr, bool applyJEC, bool statErr) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};
  double zbins024 [] = {0.064, 0.220, 0.376, 0.532, 0.688, 0.844, 1.000};

  string binTag = workDirName;
  if (binTag.find("midJtPt")!=std::string::npos) binTag = "midJtPt";
  else if (binTag.find("lowJtPt")!=std::string::npos) binTag = "lowJtPt";
  else if (binTag.find("lowerJtPt")!=std::string::npos) binTag = "lowerJtPt";
  else if (binTag.find("lowestJtPt")!=std::string::npos) binTag = "lowestJtPt";
  else if (binTag.find("highJtPt")!=std::string::npos) binTag = "highJtPt";
  else if (binTag.find("higherJtPt")!=std::string::npos) binTag = "higherJtPt";

  double prSyst;
  double nprSyst;
  string systName [] = {"ctauBkg", "ctauErr", "ctauRes", "ctauTrue", "massBkg", "massSig","fullAccEff"};//, "AccEffMisMod", "tnpmuidSyst", "AccEffStat", "tnpstaStat", "tnpstaSyst", "tnptrgStat", "tnptrgSyst", "tnpbinned", "tnptrkSyst", "tnpmuidStat"};

  int nSyst = sizeof(systName)/sizeof(systName[0]);

  int nzbins = 0;
  double zedmin, zedmax;
  if (!strcmp(rapRegion,"016")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else if (!strcmp(rapRegion,"1624")){
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }
  else {
    nzbins = sizeof(zbins024)/sizeof(double)-1;
    zedmin = zbins024[0];
    zedmax = zbins024[nzbins];
  }

  TH1F* prNhist_pp = new TH1F ("prNhist_pp",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);
  TH1F* nprNhist_pp = new TH1F ("nprNhist_pp",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);

  TH1F* prNhist_PbPb = new TH1F ("prNhist_PbPb",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);
  TH1F* nprNhist_PbPb = new TH1F ("nprNhist_PbPb",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);

  TH1F* prNhist_PbPb_cent = new TH1F ("prNhist_PbPb_cent",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);
  TH1F* nprNhist_PbPb_cent = new TH1F ("nprNhist_PbPb_cent",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);

  TH1F* prNhist_PbPb_peri = new TH1F ("prNhist_PbPb_peri",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);
  TH1F* nprNhist_PbPb_peri = new TH1F ("nprNhist_PbPb_peri",";z(J/#psi);N(J/#psi)", nzbins, zedmin, zedmax);
  

  gSystem->mkdir(Form("Output/%s/%s/%s/fitsPars",workDirName, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/%s/%s/result/tree_allvars.root",workDirName, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl;
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s", workDirName), DSTag,"", fitType, wantPureSMC, jetR, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }

  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL=0, bfrac_errH=0;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  //tr->SetBranchAddress("N_Jpsi_val",&val);
  //tr->SetBranchAddress("N_Jpsi_errL",&errL);
  //tr->SetBranchAddress("N_Jpsi_errH",&errH);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);
  int ord = 0;
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    prSyst = 0;
    nprSyst = 0;
    if (!statErr) {
      for (int j=0; j<nSyst ; j++) {
	double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_prompt_%s_%s.csv", binTag.c_str(), collSystem, systName[j].c_str()), zmin, zmax, ymin, ymax, centmin, centmax);
	double v2 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_nonprompt_%s_%s.csv", binTag.c_str(), collSystem, systName[j].c_str()), zmin, zmax, ymin, ymax, centmin, centmax);
	//cout <<"reading from file "<<Form("../Fitter/Systematics/csv/syst_%s_NJpsi_prompt_%s_%s.csv", binTag.c_str(), collSystem, systName[j].c_str())<<endl;
	prSyst=sqrt(pow(prSyst,2)+pow(v1,2));
	nprSyst=sqrt(pow(nprSyst,2)+pow(v2,2));
      }
    }
    cout <<"collSystem: "<<collSystem<<",cent:["<<centmin<<"-"<<centmax<<"], z:["<<zmin<<"-"<<zmax<<"], N_Jpsi_parLoad_mass = "<<val<<", N_Jpsi_parLoad_mass_err = "<<errL<<", b_Jpsi_val = "<<bfrac<<", b_Jpsi_errL = "<<bfrac_errL<<endl;
    double err = 0;
    if (TString(collSystem)=="PP") {
      prNhist_pp->SetBinContent(prNhist_pp->FindBin(zmin+0.001),val*(1-bfrac));
      if (statErr)
	err = val*(1-bfrac)*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2));
      else 
	err = val*(1-bfrac)*prSyst;
      if (err > val*(1-bfrac)) {
	cout <<"[WARNING] very big error"<<endl;
	err = val*(1-bfrac);
      }
      prNhist_pp->SetBinError(prNhist_pp->FindBin(zmin+0.001), err);
      
      nprNhist_pp->SetBinContent(nprNhist_pp->FindBin(zmin+0.001),val*bfrac);
      if (statErr)
	err = val*bfrac*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2));
      else 
	err = val*bfrac*nprSyst;
      if (err > val*bfrac) {
        cout <<"[WARNING] very big error"<<endl;
        err = val*bfrac;
      }
      nprNhist_pp->SetBinError(nprNhist_pp->FindBin(zmin+0.001), err);
    }
    else {
      if (centmin==0 && centmax==180) prNhist_PbPb->SetBinContent(prNhist_PbPb->FindBin(zmin+0.001),val*(1-bfrac));
      else if (centmin==0 && centmax==40) prNhist_PbPb_cent->SetBinContent(prNhist_PbPb_cent->FindBin(zmin+0.001),val*(1-bfrac));
      else if (centmin==40 && centmax==180) prNhist_PbPb_peri->SetBinContent(prNhist_PbPb_peri->FindBin(zmin+0.001),val*(1-bfrac));
      else cout<<"[WARNING] check centrality selection"<<endl;
      if (statErr)
	err = val*(1-bfrac)*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2));
      else 
	err = val*(1-bfrac)*prSyst;
      
      if (err > val*(1-bfrac)) {
        cout <<"[WARNING] very big error: "<<err/(val*(1-bfrac))<<endl;
        err = val*(1-bfrac);
      }

      if (centmin==0 && centmax==180) prNhist_PbPb->SetBinError(prNhist_PbPb->FindBin(zmin+0.001), err);
      else if (centmin==0 && centmax==40) prNhist_PbPb_cent->SetBinError(prNhist_PbPb_cent->FindBin(zmin+0.001), err);
      else if (centmin==40 && centmax==180) prNhist_PbPb_peri->SetBinError(prNhist_PbPb_peri->FindBin(zmin+0.001), err);


      if (centmin==0 && centmax==180) nprNhist_PbPb->SetBinContent(nprNhist_PbPb->FindBin(zmin+0.001),val*bfrac);
      else if (centmin==0 && centmax==40) nprNhist_PbPb_cent->SetBinContent(nprNhist_PbPb_cent->FindBin(zmin+0.001),val*bfrac);
      else if (centmin==40 && centmax==180) nprNhist_PbPb_peri->SetBinContent(nprNhist_PbPb_peri->FindBin(zmin+0.001),val*bfrac);
      if (statErr)
	err = val*bfrac*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2));
      else err = val*bfrac*nprSyst;
      if (err > val*bfrac) {
        cout <<"[WARNING] very big error: "<<err/(val*bfrac)<<endl;
        err = val*bfrac;
      }
      if (centmin==0 && centmax==180) nprNhist_PbPb->SetBinError(nprNhist_PbPb->FindBin(zmin+0.001), err);
      else if (centmin==0 && centmax==40) nprNhist_PbPb_cent->SetBinError(nprNhist_PbPb_cent->FindBin(zmin+0.001), err);
      else if (centmin==40 && centmax==180) nprNhist_PbPb_peri->SetBinError(nprNhist_PbPb_peri->FindBin(zmin+0.001), err);
    }
    
  }
  TFile* fsave = new TFile (Form("Output/%s/%s/%s/fitsPars/unfoldingInput_%s_rap%s_jetR%d_%s.root", workDirName, fitType, DSTag, binTag.c_str(), rapRegion, (int)(jetR*10), statErr?"statErr":"systErr"),"RECREATE");
  fsave->ls();
  prNhist_pp->Write(Form("prHist_PP_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist_pp->Write(Form("nprHist_PP_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  prNhist_PbPb->Write(Form("prHist_PbPb_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist_PbPb->Write(Form("nprHist_PbPb_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  prNhist_PbPb_cent->Write(Form("prHist_PbPb_%s_rap%s_centBin_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist_PbPb_cent->Write(Form("nprHist_PbPb_%s_rap%s_centBin_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  prNhist_PbPb_peri->Write(Form("prHist_PbPb_%s_rap%s_periBin_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist_PbPb_peri->Write(Form("nprHist_PbPb_%s_rap%s_periBin_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));

  fsave->Close();
  delete prNhist_pp; delete nprNhist_pp; delete prNhist_PbPb; delete nprNhist_PbPb; delete fsave; delete f;
}
void getUnfoldingInput_all(bool statErr) {  
  getUnfoldingInput("DataFits/DataFits_midJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
  getUnfoldingInput("DataFits/DataFits_lowJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
  getUnfoldingInput("DataFits/DataFits_lowerJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
  getUnfoldingInput("DataFits/DataFits_lowestJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
  getUnfoldingInput("DataFits/DataFits_highJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
  getUnfoldingInput("DataFits/DataFits_higherJtPt", "", "DATA", "ctauMass", 0, 0.3, "AccEff", 1, statErr);
}

void plot18012Comp(bool isPr) {
  TFile* file19 = TFile::Open("Output/DataFits_HIN18012Comp/ctauMass/DATA/fitsPars/unfoldingInput_DataFits_HIN18012Comp_rap016_jetR4_statErr.root");
  TFile* file18 = TFile::Open("~/DimuonCADIs/HIN-16-004/Fitter/Output/DataFits/DataFits_midJtPt/DataFits_016/ctauMass/DATA/fitsPars/unfoldingInput_midJtPt_rap016_statErr.root");
  TFile* file18Syst = TFile::Open("~/DimuonCADIs/HIN-16-004/Fitter/Output/DataFits/DataFits_midJtPt/DataFits_016/ctauMass/DATA/fitsPars/unfoldingInput_midJtPt_rap016_systErr.root");
  TH1F* hist19 = (TH1F*) file19->Get(Form("%sprHist_PP_DataFits_HIN18012Comp_rap016_statErr",isPr?"":"n"));
  TH1F* hist18 = (TH1F*) file18->Get(Form("%sprHist_midJtPt_rap016_statErr",isPr?"":"n"));
  TH1F* hist18Syst = (TH1F*) file18Syst->Get(Form("%sprHist_midJtPt_rap016_systErr",isPr?"":"n"));

  int nBin = hist19->GetNbinsX()+1;
  for (int i=0; i<=nBin; i++) {
    if (hist19->GetBinCenter(i)<0.44 || hist19->GetBinCenter(i)>1) {
      hist19->SetBinContent(i,0);
      hist19->SetBinError(i,0);
    }
  }
  nBin = hist18->GetNbinsX()+1;
  for (int i=0; i<=nBin; i++) {
    if (hist18->GetBinCenter(i)<0.44 || hist18->GetBinCenter(i)>1) {
      hist18->SetBinContent(i,0);
      hist18->SetBinError(i,0);
      hist18Syst->SetBinContent(i,0);
      hist18Syst->SetBinError(i,0);
    }
  }

  hist19->Scale(1./hist19->Integral("width"));
  hist18->Scale(1./hist18->Integral("width"));
  hist18Syst->Scale(1./hist18Syst->Integral("width"));

  hist19->SetLineColor(8);
  hist19->SetMarkerColor(8);
  hist19->SetMarkerStyle(kFullCircle);

  hist18->SetLineColor(9);
  hist18->SetMarkerColor(9);
  hist18->SetMarkerStyle(kFullSquare);

  hist18Syst->SetLineColor(9);
  hist18Syst->SetMarkerColor(9);
  hist18Syst->SetMarkerStyle(kFullSquare);
  hist18Syst->SetFillColorAlpha(38,0.5);

  hist19->GetYaxis()->SetTitle("1/N dN/dz");
  hist19->GetXaxis()->SetTitle("z");
  hist19->GetXaxis()->SetRangeUser(0,1);
  if (hist18->GetMaximum()>hist19->GetMaximum()) hist19->GetYaxis()->SetRangeUser(0,1.5*hist18->GetMaximum());
  else hist19->GetYaxis()->SetRangeUser(0,1.5*hist19->GetMaximum());

  TLegend* leg = new TLegend(0.2,0.3,0.5,0.4);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist19,"2017 pp, R=0.4","lp");
  leg->AddEntry(hist18,"2015 pp, R=0.4","lp");
  TCanvas* c = new TCanvas("c","",900,900);
  hist19->Draw("e1");
  hist18Syst->Draw("e2 same");
  hist18->Draw("e1 same");
  leg->Draw("same");
  c->SaveAs(Form("Output/MeasureDistComp_18012vs19007_ppData%s_midRap_noUnderflow.pdf",isPr?"":"npr"));
  c->SaveAs(Form("Output/MeasureDistComp_18012vs19007_ppData%s_midRap_noUnderflow.png",isPr?"":"npr"));
}

void compareRaa() {
  gStyle->SetOptStat(0);
  double glb19007  = sqrt(0.011*0.011+0.015*0.015);
  double glb16025 = 0.023;
  double Raa16025 [] = {0.356, 0.351,0.345,0.335,0.324,0.318};
  double stat16025 [] = {0.008,0.008,0.008,0.008,0.008,0.013};
  double syst16025 [] = {0.024,0.024,0.024,0.027,0.026,0.024};

  int nBin = sizeof(Raa16025)/sizeof(double);
  double totPP = 0.0460229;
  double statPP = 0.000276909;
  double systPP = 0.000844524;

  double totPbPb = 0.0118678;
  double statPbPb = 0.000586712;
  double systPbPb = 0.00133272;

  double Raa19007 = (totPbPb/totPP);
  double stat19007 = Raa19007*sqrt((statPbPb/totPbPb)*(statPbPb/totPbPb)+(statPP/totPP)*(statPP/totPP));
  double syst19007 = Raa19007*sqrt((systPbPb/totPbPb)*(systPbPb/totPbPb)+(systPP/totPP)*(systPP/totPP));

  TH1D* hist16025 = new TH1D("hist16025",";|y|;R_{AA}",nBin,0,2.4);
  TH1D* hSyst16025 = new TH1D("hSyst16025",";|y|;R_{AA}",nBin,0,2.4);
  TH1D* hist19007 = new TH1D("hist19007",";|y|;R_{AA}",1,0,2.4);
  TH1D* hSyst19007 = new TH1D("hSyst19007",";|y|;R_{AA}",1,0,2.4);
  for (int i=0;i<nBin;i++) {
    hist16025->SetBinContent(i+1,Raa16025[i]);
    hist16025->SetBinError(i+1,stat16025[i]);
    hSyst16025->SetBinContent(i+1,Raa16025[i]);
    hSyst16025->SetBinError(i+1,syst16025[i]+(Raa16025[i]*glb16025));
    }
  hist19007->SetBinContent(1,Raa19007);
  hist19007->SetBinError(1,stat19007);
  hSyst19007->SetBinContent(1,Raa19007);
  hSyst19007->SetBinError(1,syst19007+(Raa19007*glb19007));

  hist16025->SetLineColor(kRed+2);
  hist16025->SetLineWidth(2);
  hist16025->SetMarkerColor(kRed+2);
  hist16025->SetMarkerStyle(kFullCircle);
  hist16025->SetFillColorAlpha(kRed-5,0.35);
  hSyst16025->SetLineColor(kRed+2);
  hSyst16025->SetLineWidth(2);
  hSyst16025->SetMarkerColor(kRed+2);
  hSyst16025->SetMarkerStyle(kFullCircle);
  hSyst16025->SetFillColorAlpha(kRed-5,0.35);

  hist19007->SetLineColor(kBlue+2);
  hist19007->SetLineWidth(2);
  hist19007->SetMarkerColor(kBlue+2);
  hist19007->SetMarkerStyle(kFullSquare);
  hist19007->SetFillColorAlpha(kBlue-5,0.35);
  hSyst19007->SetLineColor(kBlue+2);
  hSyst19007->SetLineWidth(2);
  hSyst19007->SetMarkerColor(kBlue+2);
  hSyst19007->SetMarkerStyle(kFullSquare);
  hSyst19007->SetFillColorAlpha(kBlue-5,0.35);

  TLegend* leg = new TLegend(0.6,0.6,0.8,0.7);
  leg->SetBorderSize(0);
  leg->AddEntry(hSyst16025,"HIN-16-025","lp");
  leg->AddEntry(hSyst19007,"HIN-19-007","lp");

  TCanvas* c = new TCanvas ("c","",800,800);
  hSyst16025->GetYaxis()->SetRangeUser(0,1);
  hSyst16025->Draw("E2");
  hSyst19007->Draw("same E2");
  hist19007->Draw("same E1");
  hSyst16025->Draw("same E2");
  hist16025->Draw("same E1");
  leg->Draw("same");
  c->SaveAs("Output/RaaComp.pdf");
  c->SaveAs("Output/RaaComp.png");
}

double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax, int cntmin,int cntmax) {
  double ans;
  ans = 0;
  ifstream file(systfile);
  if (!(file.good())) return ans;

  string systname; getline(file,systname);

  string line;
  double zmin=0, zmax=0, ymin=0, ymax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;

  while (file.good()) {
    getline(file,line);
    if (line.size()==0) break;
    TString tline(line.c_str());
    TString t; Int_t from = 0, cnt=0;
    while (tline.Tokenize(t, from , ",")) {
      t.Strip(TString::kBoth,' ');
      value = atof(t.Data());
      if (cnt==0) ymin = atof(t.Data());
      else if (cnt==1) ymax = value;
      else if (cnt==2) ptmin = value;
      else if (cnt==3) ptmax = value;
      else if (cnt==4) zmin = value;
      else if (cnt==5) zmax = value;
      else if (cnt==6) centmin = value;
      else if (cnt==7) centmax = value;
      else if (cnt>8) {
	cout << "Warning, too many fields, I'll take the last one." << endl;
	continue;
      }
      cnt++;
    }
    //if (!(zmin == 0.4 && zmax == 1.0) && !(zmin == 0.2 && zmax == 1.0)) ///// to not take the integrated results
    if (zmin<zedmin+0.001 && zmin>zedmin-0.001 && zmax<zedmax+0.001 && zmax>zedmax-0.001 && ymin<rapmin+0.001 && ymin>rapmin-0.001 && ymax<rapmax+0.001 && ymax>rapmax-0.001 && centmin<cntmin+0.001 && centmin>cntmin-0.001 && centmax<cntmax+0.001 && centmax>cntmax-0.001)
      //ans.push_back(value);
      ans = value;
  }
  file.close();
  if (ans == 0) cout <<"[WARNING] systematic = 0;"<<endl;
  if (ans >1) cout <<"[WARNING] huge systematic in "<<systfile<<", syst = "<<ans<<endl;
  return ans;
}

