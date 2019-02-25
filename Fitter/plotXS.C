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
#include "plotFitParams.C"
#include "Macros/CMS/tdrstyle.C"

using namespace std;

void plotXS(const char* workDirName, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool statErr) {
  gStyle->SetOptStat(0);
  double ptbins [] = {0, 1.6, 2.4};//{6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50}; //{0, 1.6, 2.4};
  double nbins = sizeof(ptbins)/sizeof(double)-1;

  string systName [] = {"ctauBkg", "ctauErr", "ctauRes", "ctauTrue", "massBkg", "massSig", "AccEffMisMod", "tnpmuidSyst", "AccEffStat", "tnpstaStat", "tnpstaSyst", "tnptrgStat", "tnptrgSyst", "tnpbinned", "tnptrkSyst", "tnpmuidStat"};
  int nSyst = sizeof(systName)/sizeof(systName[0]);

  string lumiTag = workDirName;

  gSystem->mkdir(Form("Output/%s/%s/%s/fitsPars",workDirName, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/%s/%s/result/tree_allvars.root",workDirName, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl;
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(workDirName, DSTag,"", fitType, wantPureSMC, applyCorr, 1);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  ofstream prOut (Form("Output/%s/%s/%s/fitsPars/prXC_vsPt.csv",workDirName, fitType, DSTag));
  ofstream nprOut (Form("Output/%s/%s/%s/fitsPars/nprXC_vsPt.csv",workDirName, fitType, DSTag));

  TH1F* prXS = new TH1F ("prXS", ";pt;XS", nbins, ptbins); prXS->Sumw2();
  TH1F* nprXS = new TH1F ("nprXS", ";pt;XS", nbins, ptbins); nprXS->Sumw2();

  TH1F* prTot = new TH1F ("prTot", "tot XS for pt[6.5-50];whatever;XS", 1, 0, 10); prTot->Sumw2();
  TH1F* nprTot = new TH1F ("nprTot", "tot XS for pt[6.5-50];whatever;XS", 1, 0, 10); nprTot->Sumw2();

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
  double deltaPt = 1.0;
  double deltaRap = 1.0;
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);

    if (lumiTag.find("muonJSON")!=std::string::npos || lumiTag.find("18025")!=std::string::npos)
      lumi = 28.0e6;
    else if (lumiTag.find("oldLumi")!=std::string::npos) 
      lumi = 25.8e6;
    cout<<"[INFO] lumi = "<<lumi<<endl;
    deltaPt = 1.0;//ptmax - ptmin;
    deltaRap = ymax - ymin;
    double normfact = 1.0;//(1.0/(lumi*deltaPt*deltaRap*1e-3));
    int xBin = prXS->FindBin((ymin+ymax)/2);
    //if ((ptmax-ptmin)<20)
    //{

    double prSyst = 0;
    double nprSyst = 0;
    for (int j=0; j<nSyst ; j++) {
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_NoJets_total_NJpsi_prompt_PP_%s.csv", systName[i].c_str()), zmin, zmax, ymin, ymax);
      double v2 = readSyst(Form("../Fitter/Systematics/csv/syst_NoJets_total_NJpsi_nonprompt_PP_%s.csv", systName[i].c_str()), zmin, zmax, ymin, ymax);
      prSyst=sqrt(pow(prSyst,2)+pow(v1,2));
      nprSyst=sqrt(pow(nprSyst,2)+pow(v2,2));
    }
    
    prXS->SetBinContent(xBin, val*(1-bfrac)*normfact);
    if (statErr)
      prXS->SetBinError(xBin, val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    else 
      prXS->SetBinError(xBin, val*(1-bfrac)*prSyst);

    nprXS->SetBinContent(xBin, val*bfrac*normfact);
    if (statErr)
      nprXS->SetBinError(xBin, val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    else nprXS->SetBinError(xBin, val*bfrac*prSyst);

    normfact = 1.0;
    cout<<"pr XS = " << prXS->GetBinContent(xBin)/normfact<<" +- "<<prXS->GetBinError(xBin)/normfact<< " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    cout<<"npr XS = " << nprXS->GetBinContent(xBin)/normfact<<" +- "<<nprXS->GetBinError(xBin)/normfact << " for "<<ptmin<<" < pt < "<<ptmax <<endl;
    prOut<<ptmin<<"  "<<ptmax<<"   "<<prXS->GetBinContent(xBin)/normfact<<"  "<<prXS->GetBinError(xBin)/normfact<<endl;
    nprOut<<ptmin<<"  "<<ptmax<<"   "<<nprXS->GetBinContent(xBin)/normfact<<"  "<<nprXS->GetBinError(xBin)/normfact<<endl;
    //}
    //else if (ptmax>35)
    //{
    //prTot->SetBinContent(prTot->FindBin(5), val*(1-bfrac)*normfact);
    //prTot->SetBinError(prTot->FindBin(5), val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    //nprTot->SetBinContent(nprTot->FindBin(5), val*bfrac*normfact);
    //nprTot->SetBinError(nprTot->FindBin(5), val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    
    //cout<<"pr XS = " << val*(1-bfrac)*normfact <<" +- "<< val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) << " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    //cout<<"npr XS = " << val*bfrac*normfact <<" +- "<< val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) << " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    //prOut<<ptmin<<"  "<<ptmax<<"  "<< val*(1-bfrac)*normfact <<"  "<< val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) <<endl;
    //nprOut<<ptmin<<"  "<<ptmax<<"  "<<val*bfrac*normfact<<"  "<< val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2))<<endl;
    //}
  }
  prOut.close();
  nprOut.close();
  prXS->Draw();
  TFile* fsave = new TFile (Form("Output/%s/%s/%s/fitsPars/XSPlot_%s.root", workDirName, fitType, DSTag, statErr?"statErr":"systErr"), "RECREATE");
  fsave->ls();
  prXS->Write("prtotXS");
  nprXS->Write("nprtotXS");
  prTot->Write("prtotIntXS");
  nprTot->Write("nprtotIntXS");
  fsave->Close();
  //delete prNhist; delete nprNhist; delete fsave; delete f;
}//end of plotXS function



void compXSPt(bool isPr, bool isDiff) {
  gStyle->SetOptStat(0);

  double ptbins [] ={6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50};
  int nptbins = sizeof(ptbins)/sizeof(double);
  double centbins [] = {0,10};

  double bins [15];
  int nbins;
  if (isDiff){
    nbins = sizeof(ptbins)/sizeof(double);
    for (int i =0; i<nbins; i++)
      bins[i] = ptbins [i];
  }
  else {
    nbins = sizeof(centbins)/sizeof(double);
    for (int i =0; i<nbins; i++)
      bins[i] = centbins [i];
  }
  nbins = nbins-1;

  TH1F* corrHist = new TH1F("corrHist","", nbins, bins); corrHist->Sumw2();
  TH1F* refHist = new TH1F("refHist","", nbins, bins); refHist->Sumw2();
  TH1F* goldHist = NULL;

  ifstream accEff;
  accEff.open(Form("../MyEfficiency/RatioStudy/AccEff16025_%s.csv", isDiff?"pt_rap0024":"cent_rap0024"));

  Float_t x,y,z,w;
  Int_t nlines = 0;

  while (1)
    {
      accEff >> x >> y >> z >> w;
      if (!accEff.good()) break;
      if (isPr) corrHist->SetBinContent(corrHist->FindBin((x+y)*0.5), z);
      else corrHist->SetBinContent(corrHist->FindBin((x+y)*0.5), w);
      corrHist->SetBinError(corrHist->FindBin((x+y)*0.5), 0.001);
      nlines++;
    }
  accEff.close();

  double refVal = 0;
  double refErr = 0;

  ifstream ref;
  ref.open(Form("Input/%sXS16025_pt.csv", isPr?"pr":"npr"));
  nlines = 0;
  while (1)
    {
      ref >> x >> y >> z >> w;
      if (!ref.good()) break;
      if (isDiff)
	{
	  refHist->SetBinContent(refHist->FindBin((x+y)*0.5), z);
	  refHist->SetBinError(refHist->FindBin((x+y)*0.5), w);
	}
      else { 
	refVal = refVal+z; 
	refErr = refErr+(w/z)*(w/z);
      }
      nlines++;
    }
  ref.close();
  //cout << "[INFO] intVal = "<<refVal<<" +- "<<refVal*sqrt(refErr)<<" intVal/(dpt) = "<<refVal/(50-6.5)<<" +- "<<refVal*sqrt(refErr)<<endl;
  refVal = refVal/(50-6.5);
  refErr = refVal*sqrt(refErr);
  if (!isDiff)
    refHist->SetBinContent(refHist->FindBin(5), refVal);

  TFile* goldFile = TFile::Open("Output/DataFits_totalN_y024/ctauMass/DATA/fitsPars/XSPlot.root");

  goldHist = (TH1F*) goldFile->Get(Form("%stot%sXS",isPr?"pr":"npr", isDiff?"":"Int"));

  goldHist->Divide(corrHist);

  refHist->SetLineColor(kRed+2);
  goldHist->SetLineColor(kBlue+2);

  refHist->SetMarkerColor(kRed);
  goldHist->SetMarkerColor(kBlue);

  refHist->SetMarkerStyle(4);
  goldHist->SetMarkerStyle(4);

  refHist->SetMarkerSize(1);
  goldHist->SetMarkerSize(1);

  TLegend* leg = new TLegend (0.5, 0.7, 0.9, 0.8);
  leg->AddEntry(refHist, "HIN-16-025 results", "lep");
  leg->AddEntry(goldHist, "HIN-18-012 (evt-by-evt)", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TH1F* axisHist = new TH1F("axisHist","", nbins, bins);
  axisHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));
  axisHist->GetYaxis()->SetLimits(0,(goldHist->GetMaximum()<1)?1:(1.2*goldHist->GetMaximum()));
  goldHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));
  refHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));

  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  if (!isDiff) axisHist->Draw();
  else refHist->Draw();
  refHist->Draw("same");
  goldHist->Draw("same");
  leg->Draw("same");

  gSystem->mkdir("Output/XSComparison");
  c->SaveAs(Form("Output/XSComparison/%s%scomparison.pdf", isPr?"pr":"npr", isDiff?"Diff":"Int"));

  if (!isDiff) cout<<"[INFO] xs(tot18-012) = "<<goldHist->GetBinContent(1)<<endl;//" xs(18012) = "<<muonHist->GetBinContent(muonHist->FindBin(6.5))<<endl;
  goldFile->Close();

}//end of compXS function

void plotXStot(bool isPr, bool fonllCorr) {
  gStyle->SetOptStat(0);

  double ybins [] = {0, 1.6, 2.4};

  TFile* totFile = TFile::Open("Output/DataFits/DataFits_NoJets/DataFits_total/ctauMass/DATA/fitsPars/XSPlot_statErr.root");
  TFile* totSystFile = TFile::Open("Output/DataFits/DataFits_NoJets/DataFits_total/ctauMass/DATA/fitsPars/XSPlot_systErr.root");
  TFile* mcFile = TFile::Open(Form("Output/MCResults%s/mcResult_%s_underflowOff.root", fonllCorr?"/fonllCorr":"", isPr?"prompt":"nonprompt"));

  TH1F* totHist = (TH1F*) totFile->Get(Form("%stotXS",isPr?"pr":"npr"));
  TH1F* totSystHist = (TH1F*) totSystFile->Get(Form("%stotXS",isPr?"pr":"npr"));

  TH1F* midMCHist = (TH1F*) mcFile->Get("zDist_mid");
  TH1F* midMCTotHist =(TH1F*) mcFile->Get("Ntot_mid");
  TH1F* fwdMCHist = (TH1F*) mcFile->Get("zDist_fwd");
  TH1F* fwdMCTotHist =(TH1F*) mcFile->Get("Ntot_fwd");
  TH1F* mcHist = new TH1F("mcHist", "", 2, ybins);

  TFile* jetMidFile = TFile::Open(Form("Output/unfData/results/unfResult_%s_mid_wStats.root",isPr?"prompt":"nonprompt"));
  TFile* jetMidSystFile = TFile::Open(Form("Output/unfData/results/systErrs_%s_Mid_Total.root",isPr?"Prompt":"NonPrompt"));
  TFile* jetFwdFile = TFile::Open(Form("Output/unfData/results/unfResult_%s_fwd_wStats.root",isPr?"prompt":"nonprompt"));
  TFile* jetFwdSystFile =TFile::Open(Form("Output/unfData/results/systErrs_%s_Fwd_Total.root",isPr?"Prompt":"NonPrompt"));

  TH1F* midHist = (TH1F*) jetMidFile->Get("zUnf"); 
  TH1F* fwdHist = (TH1F*)jetFwdFile->Get("zUnf");
  TH1F* midSystHist = (TH1F*) jetMidSystFile->Get("quarkoniaSyst");
  TH1F* fwdSystHist = (TH1F*) jetFwdSystFile->Get("quarkoniaSyst");

  TH1F* midTrMHist = (TH1F*) jetMidSystFile->Get("unfTrMatrixSyst");
  TH1F* fwdTrMHist = (TH1F*) jetFwdSystFile->Get("unfTrMatrixSyst");

  double valMid =0;
  double errMid =0;
  double systMid =0;
  double valFwd =0;
  double errFwd =0;
  double systFwd =0;
  double valMCMid =0;
  double valMCFwd =0;
  double errMCMid =0;
  double errMCFwd =0;

  for (int j =0; j<4 ; j++){
    valFwd =valFwd + fwdHist->GetBinContent(fwdHist->FindBin(0.2+j*(0.2)));
    errFwd = sqrt(errFwd*errFwd + fwdHist->GetBinError(fwdHist->FindBin(0.2+j*0.2))*fwdHist->GetBinError(fwdHist->FindBin(0.2+j*0.2)));
    systFwd = sqrt(systFwd*systFwd + fwdSystHist->GetBinError(fwdSystHist->FindBin(0.2+j*0.2))*fwdSystHist->GetBinError(fwdSystHist->FindBin(0.2+j*0.2)));
    systFwd = sqrt(systFwd*systFwd + fwdTrMHist->GetBinError(fwdTrMHist->FindBin(0.2+j*0.2))*fwdTrMHist->GetBinError(fwdTrMHist->FindBin(0.2+j*0.2)));
    valMCFwd= valMCFwd + fwdMCHist->GetBinContent(fwdMCHist->FindBin(0.2+j*(0.2)));
    errMCFwd= sqrt(errMCFwd*errMCFwd + fwdMCHist->GetBinError(fwdMCHist->FindBin(0.2+j*0.2))*fwdMCHist->GetBinError(fwdMCHist->FindBin(0.2+j*0.2)));
  }
  if (isPr){
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0313)*(valFwd*0.0313));//JES
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0199)*(valFwd*0.0199));//SI
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0745)*(valFwd*0.0745));//JES prompt
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0267)*(valFwd*0.0267));//JER
  }
  else {
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0313)*(valFwd*0.0313));//JES
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0126)*(valFwd*0.0126));//SI 
    systFwd = sqrt(systFwd*systFwd + (valFwd*0.0166)*(valFwd*0.0166));//JER
  }
  for (int i =0; i<5 ; i++){
    valMid = valMid+midHist->GetBinContent(midHist->FindBin(0.44+(i*0.14)));
    errMid = sqrt(errMid*errMid + midHist->GetBinError(midHist->FindBin(0.44+(i*0.16)))*midHist->GetBinError(midHist->FindBin(0.46+(i*0.16))));
    systMid = sqrt(systMid*systMid + midSystHist->GetBinError(midSystHist->FindBin(0.44+(i*0.16)))*midSystHist->GetBinError(midSystHist->FindBin(0.46+(i*0.16))));
    systMid = sqrt(systMid*systMid + midTrMHist->GetBinError(midTrMHist->FindBin(0.44+(i*0.16)))*midTrMHist->GetBinError(midTrMHist->FindBin(0.46+(i*0.16))));
    valMCMid= valMCMid + midMCHist->GetBinContent(midMCHist->FindBin(0.44+(i*0.14)));
    errMCMid= sqrt(errMCMid*errMCMid + midMCHist->GetBinError(midMCHist->FindBin(0.44+(i*0.16)))*midMCHist->GetBinError(midMCHist->FindBin(0.44+(i*0.14))));
  }
  if (isPr){
    systMid = sqrt(systMid*systMid + (valMid*0.0106)*(valMid*0.0106));//JES
    systMid = sqrt(systMid*systMid + (valMid*0.0052)*(valMid*0.0052));//SI
    systMid = sqrt(systMid*systMid + (valMid*0.0461)*(valMid*0.0461));//JES prompt
    systMid = sqrt(systMid*systMid + (valMid*0.0130)*(valMid*0.0130));//JER
  }
  else {
    systMid = sqrt(systMid*systMid + (valMid*0.0106)*(valMid*0.0106));//JES
    systMid = sqrt(systMid*systMid + (valMid*0.0093)*(valMid*0.0093));//SI 
    systMid = sqrt(systMid*systMid + (valMid*0.0069)*(valMid*0.0069));//JER
  }

  int bin =0;
  bin = totHist->FindBin(1);
  cout<<"totXS in mid rapidity = "<<totHist->GetBinContent(bin)/(27.39e3*1.6)<<", statErr =  "<<totHist->GetBinError(bin)/(27.39e3*1.6)<<", systErr "<<totSystHist->GetBinError(bin)/(27.39e3*1.6)<<endl;
  errMid = errMid/valMid;
  systMid = systMid/valMid;
  valMid = valMid/totHist->GetBinContent(bin);
  errMid = sqrt(errMid*errMid + totHist->GetBinError(bin)*totHist->GetBinError(bin)/(totHist->GetBinContent(bin)*totHist->GetBinContent(bin)));
  errMid = errMid*valMid;
  systMid = sqrt(systMid*systMid + totSystHist->GetBinError(bin)*totSystHist->GetBinError(bin)/(totSystHist->GetBinContent(bin)*totSystHist->GetBinContent(bin)));
  systMid = systMid*valMid;
  totHist->SetBinContent(bin, valMid);
  totHist->SetBinError(bin, errMid);
  totSystHist->SetBinContent(bin, valMid);
  totSystHist->SetBinError(bin, systMid);

  errMCMid = errMCMid/valMCMid;
  valMCMid = valMCMid/midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001));
  errMCMid = sqrt(errMCMid*errMCMid + midMCTotHist->GetBinError(midMCTotHist->FindBin(0.0001))*midMCTotHist->GetBinError(midMCTotHist->FindBin(0.0001))/(midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001))*midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001))));
  errMCMid=errMCMid*valMCMid;
  mcHist->SetBinContent(bin, valMCMid);
  mcHist->SetBinError(bin, errMCMid);

  bin = totHist->FindBin(2);
  cout<<"totXS in fwd rapidity = "<<totHist->GetBinContent(bin)/(27.39e3*1.6)<<", statErr =  "<<totHist->GetBinError(bin)/(27.39e3*1.6)<<", systErr "<<totSystHist->GetBinError(bin)/(27.39e3*1.6)<<endl;
  errFwd = errFwd/valFwd;
  systFwd = systFwd/valFwd;
  valFwd = valFwd/totHist->GetBinContent(bin);
  errFwd = sqrt(errFwd*errFwd + totHist->GetBinError(bin)*totHist->GetBinError(bin)/(totHist->GetBinContent(bin)*totHist->GetBinContent(bin)));
  errFwd = errFwd*valFwd;
  systFwd = sqrt(systFwd*systFwd + totSystHist->GetBinError(bin)*totSystHist->GetBinError(bin)/(totSystHist->GetBinContent(bin)*totSystHist->GetBinContent(bin)));
  systFwd = systFwd*valFwd;
  totHist->SetBinContent(bin, valFwd);
  totHist->SetBinError(bin, errFwd);
  totSystHist->SetBinContent(bin, valFwd);
  totSystHist->SetBinError(bin, systFwd);

  errMCFwd = errMCFwd/valMCFwd;
  valMCFwd = valMCFwd/fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001));
  errMCFwd = sqrt(errMCFwd*errMCFwd + fwdMCTotHist->GetBinError(fwdMCTotHist->FindBin(0.0001))*fwdMCTotHist->GetBinError(fwdMCTotHist->FindBin(0.0001))/(fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001))*fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001))));
  errMCFwd=errMCFwd*valMCFwd;
  mcHist->SetBinContent(bin, valMCFwd);
  mcHist->SetBinError(bin, errMCFwd);

  gSystem->mkdir("Output/XSComparison");
  TFile *fsave = new TFile(Form("Output/XSComparison/%sXSRatio%s.root", isPr?"pr":"npr", fonllCorr?"_fonllCorr":""), "RECREATE");
  totHist->Write("XS");
  totSystHist->Write("XS_syst");
  mcHist->Write("XS_mc");
  fsave->Close();
}//end of compXStot

void drawXStot(bool fonllCorr) {
  gStyle->SetOptStat(0);
  setTDRStyle();

  TFile *prf = TFile::Open("Output/XSComparison/prXSRatio.root"); 
  TFile *nprf = TFile::Open(Form("Output/XSComparison/nprXSRatio%s.root", fonllCorr?"_fonllCorr":""));
  TH1F *prHist = (TH1F*) prf->Get("XS");
  TH1F *nprHist = (TH1F*) nprf->Get("XS");
  TH1F *prSystHist = (TH1F*) prf->Get("XS_syst");
  TH1F *nprSystHist = (TH1F*) nprf->Get("XS_syst");
  TH1F *prMCHist = (TH1F*) prf->Get("XS_mc");
  TH1F *nprMCHist = (TH1F*) nprf->Get("XS_mc");

  prHist->Scale(100);
  nprHist->Scale(100);
  prSystHist->Scale(100);
  nprSystHist->Scale(100);
  prMCHist->Scale(100);
  nprMCHist->Scale(100);

  prHist->SetMarkerColor(kMagenta+2);
  prHist->SetMarkerStyle(kFullCircle);
  prHist->SetMarkerSize(2);
  prHist->SetLineColor(kMagenta+2);
  prHist->SetLineWidth(2);

  nprHist->SetMarkerColor(kMagenta+2);
  nprHist->SetMarkerStyle(kFullSquare);
  nprHist->SetMarkerSize(2);
  nprHist->SetLineColor(kMagenta+2);
  nprHist->SetLineWidth(2);

  prSystHist->SetMarkerColor(kMagenta+2);
  prSystHist->SetMarkerStyle(kFullCircle);
  prSystHist->SetMarkerSize(0);
  prSystHist->SetLineColor(kMagenta+2);
  prSystHist->SetFillColorAlpha(kMagenta-5, 0.5);
  prSystHist->SetLineWidth(2);

  nprSystHist->SetMarkerColor(kMagenta+2);
  nprSystHist->SetMarkerStyle(kFullSquare);
  nprSystHist->SetMarkerSize(0);
  nprSystHist->SetLineColor(kMagenta+2);
  nprSystHist->SetFillColorAlpha(kMagenta-5, 0.5);
  nprSystHist->SetLineWidth(2);

  prMCHist->SetMarkerColor(kCyan+2);
  prMCHist->SetMarkerStyle(kFullCircle);
  prMCHist->SetMarkerSize(2);
  prMCHist->SetLineColor(kCyan+2);
  prMCHist->SetLineWidth(2);

  nprMCHist->SetMarkerColor(kCyan+2);
  nprMCHist->SetMarkerStyle(kFullSquare);
  nprMCHist->SetMarkerSize(2);
  nprMCHist->SetLineColor(kCyan+2);
  nprMCHist->SetLineWidth(2);

  double ybins [] = {0, 1.6, 2.4};
  TH1F *axisHist = new TH1F("axisHist","", 2, ybins);
  axisHist->GetYaxis()->SetRangeUser(0.005, 500);
  axisHist->GetYaxis()->SetTitle("J/#psi-in-jet fraction (%)");
  axisHist->GetXaxis()->SetTitle("|y_{J/#psi}|");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->SetTitleOffset(0.8);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  TLegend* leg = new TLegend(0.47,0.70,0.88,0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  leg->AddEntry(prHist, "Prompt data", "lp");
  leg->AddEntry(nprHist, "Nonprompt data", "lp");
  leg->AddEntry(prMCHist, "Prompt PYTHIA 8", "lp");
  leg->AddEntry(nprMCHist, Form("Nonprompt PYTHIA 8%s",fonllCorr?" (fonll)":""), "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *  text2 = new TLatex(0.2 , 0.22, "6.5 < p_{T,J/#psi} < 35 GeV for |y_{J/#psi}| < 1.6");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.044);
  text2->SetLineWidth(2);
  
  TLatex *  text3 = new TLatex(0.2 , 0.16, "3 < p_{T,J/#psi} < 35 GeV for 1.6 < |y_{J/#psi}| < 2.4");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.044);
  text3->SetLineWidth(2);
  
  TLatex *  text = new TLatex(0.2 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  
  TLatex *  text1 = new TLatex(0.18 ,0.77,"Preliminary");
  text1->SetNDC();
  text1->SetTextFont(52);
  text1->SetTextSize(0.055);
  text1->SetLineWidth(2);
  
  TLatex *  text4 = new TLatex(0.57 ,0.93,"pp 27.39 pb^{-1} (5.02 TeV)");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.04);
  text4->SetLineWidth(2);
  

  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  axisHist->Draw();
  prSystHist->Draw("E2 same");
  prHist->Draw("E1 same");
  prMCHist->Draw("E1 same");
  nprSystHist->Draw("E2 same");
  nprHist->Draw("E1 same");
  nprMCHist->Draw("E1 same");
  leg->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  //c->SaveAs("Output/XSComparison/totXSPlot.pdf");
  c->SetLogy();
  c->SaveAs(Form("Output/XSComparison/totXSPlot_logScale%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/totXSPlot_logScale%s.C", fonllCorr?"_fonllCorr":""));

  //////draw prompt and nonprompt separately
  leg = new TLegend(0.68,0.70,0.88,0.80);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(prHist, "Data", "lp");
  leg->AddEntry(prMCHist, "PYTHIA 8", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TLatex * text0 = new TLatex(0.68 ,0.82,"Prompt J/#psi");
  text0->SetNDC();
  text0->SetTextFont(62);
  text0->SetTextSize(0.044);
  text0->SetLineWidth(2);

  axisHist->GetYaxis()->SetRangeUser(0.005, 200);
  axisHist->Draw();
  prSystHist->Draw("E2 same");
  prHist->Draw("E1 same");
  prMCHist->Draw("E1 same");
  leg->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  text0->Draw("same");
  c->SetLogy();
  c->SaveAs(Form("Output/XSComparison/prXSPlot_logScale%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/prXSPlot_logScale%s.C", fonllCorr?"_fonllCorr":""));

  leg = new TLegend(0.68,0.70,0.88,0.80);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(nprHist, "Data", "lp");
  leg->AddEntry(nprMCHist, Form("PYTHIA 8%s",fonllCorr?" (fonll)":""), "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  text0 = new TLatex(0.58 ,0.82,"Nonprompt J/#psi");
  text0->SetNDC();
  text0->SetTextFont(62);
  text0->SetTextSize(0.044);
  text0->SetLineWidth(2);

  axisHist->GetYaxis()->SetRangeUser(0.05, 500);
  axisHist->Draw();
  nprSystHist->Draw("E2 same");
  nprHist->Draw("E1 same");
  nprMCHist->Draw("E1 same");
  leg->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  text0->Draw("same");
  c->SetLogy();
  c->SaveAs(Form("Output/XSComparison/nprXSPlot_logScale%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/nprXSPlot_logScale%s.C", fonllCorr?"_fonllCorr":""));

  c = new TCanvas("cr", "", 1000, 1000);
  prHist->Divide(prMCHist);
  prSystHist->Divide(prMCHist);
  nprHist->Divide(nprMCHist);
  nprSystHist->Divide(nprMCHist);
  axisHist->GetYaxis()->SetRangeUser(0, 5);
  axisHist->GetYaxis()->SetTitle("J/#psi-in-jet fraction data/MC");
  axisHist->GetYaxis()->SetTitleOffset(1.0);
  axisHist->GetXaxis()->SetTitle("|y_{J/#psi}|");
  axisHist->GetXaxis()->SetTitleOffset(0.8);
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetYaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  prHist->SetMarkerColor(kGreen+2);
  prHist->SetMarkerStyle(kFullCircle);
  prHist->SetMarkerSize(2);
  prHist->SetLineColor(kGreen+2);
  prHist->SetLineWidth(2);

  nprHist->SetMarkerColor(9);
  nprHist->SetMarkerStyle(kFullSquare);
  nprHist->SetMarkerSize(2);
  nprHist->SetLineColor(9);
  nprHist->SetLineWidth(2);

  prSystHist->SetMarkerColor(29);
  prSystHist->SetMarkerStyle(kFullCircle);
  prSystHist->SetMarkerSize(0);
  prSystHist->SetLineColor(29);
  prSystHist->SetFillColorAlpha(29, 0.5);
  prSystHist->SetLineWidth(2);

  nprSystHist->SetMarkerColor(40);
  nprSystHist->SetMarkerStyle(kFullSquare);
  nprSystHist->SetMarkerSize(0);
  nprSystHist->SetLineColor(40);
  nprSystHist->SetFillColorAlpha(40, 0.5);
  nprSystHist->SetLineWidth(2);

  leg = new TLegend(0.62,0.72,0.88,0.85);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(prHist, "Prompt", "lp");
  leg->AddEntry(nprHist, "Nonprompt", "lp");
  //leg->AddEntry(prMCHist, "prompt MC", "lp");
  //leg->AddEntry(nprMCHist, "nonprompt MC", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  axisHist->Draw();
  prSystHist->Draw("E2 same");
  prHist->Draw("E1 same");
  nprSystHist->Draw("E2 same");
  nprHist->Draw("E1 same");
  leg->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  c->SaveAs(Form("Output/XSComparison/totXSPlot_DataMCFraction%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/totXSPlot_DataMCFraction%s.C", fonllCorr?"_fonllCorr":""));

  //////draw prompt and nonprompt separately
  text0 = new TLatex(0.68 ,0.82,"Prompt J/#psi");
  text0->SetNDC();
  text0->SetTextFont(62);
  text0->SetTextSize(0.044);
  text0->SetLineWidth(2);

  axisHist->Draw();
  prSystHist->Draw("E2 same");
  prHist->Draw("E1 same");
  text0->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  c->SaveAs(Form("Output/XSComparison/prXSPlot_DataMCFraction%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/prXSPlot_DataMCFraction%s.C", fonllCorr?"_fonllCorr":""));

  text0 = new TLatex(0.58 ,0.82,"Nonprompt J/#psi");
  text0->SetNDC();
  text0->SetTextFont(62);
  text0->SetTextSize(0.044);
  text0->SetLineWidth(2);

  axisHist->Draw();
  nprSystHist->Draw("E2 same");
  nprHist->Draw("E1 same");
  text0->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  c->SaveAs(Form("Output/XSComparison/nprXSPlot_DataMCFraction%s.pdf", fonllCorr?"_fonllCorr":""));
  c->SaveAs(Form("Output/XSComparison/nprXSPlot_DataMCFraction%s.C", fonllCorr?"_fonllCorr":""));
}
