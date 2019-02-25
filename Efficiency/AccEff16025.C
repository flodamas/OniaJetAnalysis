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

void AccEff16025 (bool getDiff)
{
  TFile* praccFile = TFile::Open("/home/llr/cms/blanco/Analysis/Onia/RAA/Results/DimuonCADIs/HIN-16-004/Efficiency/files_acc/nominal/histos_jpsi_pp.root");
  TFile* preffFile = TFile::Open("/home/llr/cms/blanco/Analysis/Onia/RAA/Results/DimuonCADIs/HIN-16-004/Efficiency/files_eff/nominal/histos_jpsi_pp.root");
  TFile* npraccFile = TFile::Open("/home/llr/cms/blanco/Analysis/Onia/RAA/Results/DimuonCADIs/HIN-16-004/Efficiency/files_acc/nominal/histos_npjpsi_pp.root");
  TFile* npreffFile = TFile::Open("/home/llr/cms/blanco/Analysis/Onia/RAA/Results/DimuonCADIs/HIN-16-004/Efficiency/files_eff/nominal/histos_npjpsi_pp.root");

  double ptbins [] ={6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50};
  double centbins [] = {0,10};

  double bins [15];
  int nbins;
  if (getDiff){
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

  string tag = "";
  if (getDiff) tag = "pt_rap0024";
  else tag = "cent_rap0024";

  TH1F* praccNum = (TH1F*) praccFile->Get(Form("hnum_%s", tag.c_str()));
  TH1F* praccDen = (TH1F*) praccFile->Get(Form("hden_%s", tag.c_str()));
  TH1F* preffNum = (TH1F*) preffFile->Get(Form("hnum_%s", tag.c_str()));
  TH1F* preffDen = (TH1F*) preffFile->Get(Form("hden_%s", tag.c_str()));

  TH1F* npraccNum = (TH1F*) npraccFile->Get(Form("hnum_%s", tag.c_str()));
  TH1F* npraccDen = (TH1F*) npraccFile->Get(Form("hden_%s", tag.c_str()));
  TH1F* npreffNum = (TH1F*) npreffFile->Get(Form("hnum_%s", tag.c_str()));
  TH1F* npreffDen = (TH1F*) npreffFile->Get(Form("hden_%s", tag.c_str()));
  ofstream fileOut (Form("RatioStudy/AccEff16025_%s.csv", tag.c_str()));

  double val = 1.0;
  if (getDiff) fileOut<<"ptmin, ptmax, ";
  else fileOut<<"centmin, centmax, ";
  fileOut<<"prAccxEff, nprAccxEff"<<endl;
  for (int i=0; i<nbins; i++)
    {
      fileOut<<bins[i]<<", "<<bins[i+1]<<", ";
      val = 1.0;
      val = praccNum->GetBinContent(praccNum->FindBin(bins[i]));
      val = val/(praccDen->GetBinContent(praccDen->FindBin(bins[i])));
      val = val*preffNum->GetBinContent(preffNum->FindBin(bins[i]));
      val = val/(preffDen->GetBinContent(preffDen->FindBin(bins[i])));
      fileOut<<val<<", ";
      val = 1.0;
      val = npraccNum->GetBinContent(npraccNum->FindBin(bins[i]));
      val = val/(npraccDen->GetBinContent(npraccDen->FindBin(bins[i])));
      val = val*npreffNum->GetBinContent(npreffNum->FindBin(bins[i]));
      val = val/(npreffDen->GetBinContent(npreffDen->FindBin(bins[i])));
      fileOut<<val<<endl;
    }

  fileOut.close();
  npreffFile->Close();
  npraccFile->Close();
  preffFile->Close();
  praccFile->Close();
}
