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
void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC);
void getUnfoldingInput(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC, bool statErr);
double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax);
void plotMCMassPars(const char* workDirName,
		    const char* rapRegion,
		    const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
		    const char* fitType, // "mass", "ctau"...
		    bool wantPureSMC, // =false,
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
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
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


void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
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
  TH1F* bkgOrd = NULL;
  if (strcmp(rapRegion,"1624"))
    bkgOrd = new TH1F ("bkgOrd",";z(J/#psi);background order", 7, 0.02, 1);
  else 
    bkgOrd = new TH1F ("bkgOrd",";z(J/#psi);background order", 5, 0, 1);

  string bkgPol [] = {"Uniform", "Chebychev1", "Chebychev2", "Chebychev3", "Chebychev4", "Chebychev5", "Chebychev6"};
  string bkgExp [] = {"Uniform", "ExpChebychev1", "ExpChebychev2", "ExpChebychev3", "ExpChebychev4", "ExpChebychev5", "ExpChebychev6"};
  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/fitsPars",workDirName, rapRegion, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
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
      bkgOrd->SetBinContent(bkgOrd->FindBin(zmin+0.001),ord);
      bkgOrd->SetBinError(bkgOrd->FindBin(zmin+0.001), 0.0001);
      //}
  }

    TLatex *  text2 = new TLatex(0.175 ,0.8,strcmp(ShapeTag,"PolChebychev")?"Exp. Chebychev bkg.":"Pol. Chebychev bkg.");
    text2->SetNDC();
    text2->SetTextFont(42);
    text2->SetTextSize(0.05);
    text2->SetLineWidth(2);

    TLatex *  text3 = new TLatex(0.21 ,0.75, Form("%.1f < |y| < %.1f", ymin, ymax));
    text3->SetNDC();
    text3->SetTextFont(42);
    text3->SetTextSize(0.05);
    text3->SetLineWidth(2);

    TLatex *  text = new TLatex(0.75 ,0.8,"CMS");
    text->SetNDC();
    text->SetTextFont(42);
    text->SetTextSize(0.06708595);
    text->SetLineWidth(5);

    TLatex *  text1 = new TLatex(0.7 ,0.72,"Preliminary");
    text1->SetNDC();
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->SetLineWidth(2);

    TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
    text4->SetNDC();
    text4->SetTextFont(42);
    text4->SetTextSize(0.05);
    text4->SetLineWidth(2);


    bkgOrd->GetYaxis()->SetRangeUser(0, 3);
    bkgOrd->SetMarkerColor(kMagenta+3);
    bkgOrd->SetMarkerStyle(33);
    bkgOrd->SetMarkerSize(3);
    bkgOrd->SetLineColor(kMagenta+2);
    //bkgOrd->SetOption("E1");

    c->cd();
    bkgOrd->Draw("EP");
    text->Draw("same");
    text1->Draw("same");
    text2->Draw("same");
    text3->Draw("same");
    text4->Draw("same");
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    f->Close();
    delete f; delete bkgOrd;
}


void getUnfoldingInput(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC, bool statErr) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  //double zbins1624 [] = {0.16, 0.3, 0.44, 0.58, 0.72, 0.86, 1.0};

  //double zbins016 [] = {0.2, 0.4, 0.6, 0.8, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};


  string binTag = workDirName;
  if (binTag.find("midJtPt")!=std::string::npos) binTag = "midJtPt";
  else if (binTag.find("lowJtPt")!=std::string::npos) binTag = "lowJtPt";
  else if (binTag.find("highJtPt")!=std::string::npos) binTag = "highJtPt";

  double prSyst;
  double nprSyst;
  string systName [] = {"ctauBkg", "ctauErr", "ctauRes", "ctauTrue", "massBkg", "massSig", "AccEffMisMod", "tnpmuidSyst", "AccEffStat", "tnpstaStat", "tnpstaSyst", "tnptrgStat", "tnptrgSyst", "tnpbinned", "tnptrkSyst", "tnpmuidStat"};

  int nSyst = sizeof(systName)/sizeof(systName[0]);

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

  TH1F* prNhist = NULL;//new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
  TH1F* nprNhist = NULL;//new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);

  if (!strcmp(rapRegion,"016")){
    prNhist = new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 7, 0.02, 1);
    nprNhist = new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 7, 0.02, 1);
  }
  else {
    prNhist = new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
    nprNhist = new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
  }
  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/fitsPars",workDirName, rapRegion, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl;
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
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
    //if (zmax < zmin+0.22){
    prSyst = 0;
    nprSyst = 0;
    for (int j=0; j<nSyst ; j++) {
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_%s_NJpsi_prompt_PP_%s.csv", binTag.c_str(), rapRegion, systName[i].c_str()), zmin, zmax, ymin, ymax);
      double v2 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_%s_NJpsi_nonprompt_PP_%s.csv", binTag.c_str(), rapRegion, systName[i].c_str()), zmin, zmax, ymin, ymax);
      prSyst=sqrt(pow(prSyst,2)+pow(v1,2));
      nprSyst=sqrt(pow(nprSyst,2)+pow(v2,2));
    }
      prNhist->SetBinContent(prNhist->FindBin(zmin+0.001),val*(1-bfrac));
      if (statErr)
	prNhist->SetBinError(prNhist->FindBin(zmin+0.001), val*(1-bfrac)*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
      else
	prNhist->SetBinError(prNhist->FindBin(zmin+0.001), val*(1-bfrac)*prSyst);

      nprNhist->SetBinContent(nprNhist->FindBin(zmin+0.001),val*bfrac);
      if (statErr)
	nprNhist->SetBinError(nprNhist->FindBin(zmin+0.001), val*bfrac*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
      else 
	nprNhist->SetBinError(nprNhist->FindBin(zmin+0.001),val*bfrac*nprSyst);
      //}
  }
  TFile* fsave = new TFile (Form("Output/%s/DataFits_%s/%s/%s/fitsPars/unfoldingInput_%s_rap%s_%s.root", workDirName, rapRegion, fitType, DSTag, binTag.c_str(), rapRegion, statErr?"statErr":"systErr"),"RECREATE");
  fsave->ls();
  prNhist->Write(Form("prHist_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist->Write(Form("nprHist_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  fsave->Close();
  delete prNhist; delete nprNhist; delete fsave; delete f;
}

double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax) {
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
    if (zmin<zedmin+0.001 && zmin>zedmin-0.001 && zmax<zedmax+0.001 && zmax>zedmax-0.001 && ymin<rapmin+0.001 && ymin>rapmin-0.001 && ymax<rapmax+0.001 && ymax>rapmax-0.001)
      //ans.push_back(value);
      ans = value;
  }
  file.close();
  if (ans == 0) cout <<"[WARNING] systematic = 0;"<<endl;
  if (ans >10) cout <<"[WARNING] huge systematic in "<<systfile<<endl;
  return ans;
}

