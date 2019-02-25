#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>

void Eff2DRatios(){

  Double_t ptbins []={3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35};
  int nptbins = ((sizeof(ptbins)/sizeof(double))-1);

  gStyle->SetOptStat(0);
  TEfficiency* prEff = NULL;
  TEfficiency* nprEff = NULL;

  TFile* prFile = TFile::Open("../Fitter/Input/pr_correction_AccEff.root", "READ");
  TFile* nprFile = TFile::Open("../Fitter/Input/npr_correction_AccEff.root", "READ");
  prEff = (TEfficiency*) prFile->Get("hcorr_Jpsi_PP");
  nprEff = (TEfficiency*) nprFile->Get("hcorr_Jpsi_PP");

  TH2F* temp = (TH2F*) prFile->Get("hcorr_his_num");

  gSystem->mkdir("RatioStudy");
  TFile* fsave = new TFile ("RatioStudy/accEff_prNpr_Ratios.root","RECREATE");
  TH2F* effRatio = (TH2F*) temp->Clone("EffRatio"); effRatio->Sumw2(); effRatio->SetTitle("pr/Npr AccxEff ratio");
  TH2F* errRatioH = (TH2F*) temp->Clone("ErrRatioH"); errRatioH->Sumw2(); errRatioH->SetTitle("pr/Npr error high ratio");
  TH2F* errRatioL = (TH2F*) temp->Clone("ErrRatioL"); errRatioL->Sumw2(); errRatioL->SetTitle("pr/Npr error low ratio");

  int nBinsX = effRatio->GetNbinsX();
  int nBinsY = effRatio->GetNbinsY();

  for (int i=1; i<=nBinsX; i++){
    for (int j=1; j<=nBinsY; j++) {

      int effBin = effRatio->GetBin(i,j);

      effRatio->SetBinContent(effBin, prEff->GetEfficiency(effBin)/nprEff->GetEfficiency(effBin));
      effRatio->SetBinError(effBin, sqrt(prEff->GetEfficiencyErrorUp(effBin)*prEff->GetEfficiencyErrorUp(effBin)+nprEff->GetEfficiencyErrorUp(effBin)*nprEff->GetEfficiencyErrorUp(effBin)));
      errRatioH->SetBinContent(effBin, prEff->GetEfficiencyErrorUp(effBin)/nprEff->GetEfficiencyErrorUp(effBin));
      errRatioL->SetBinContent(effBin, prEff->GetEfficiencyErrorLow(effBin)/nprEff->GetEfficiencyErrorLow(effBin));
    }//end of y-bins loop
  }//end of x-bins loop
  prFile->Close();
  nprFile->Close();

  cout<<"[INFO] done with the efficiencies and start to save"<<endl;
  effRatio->Write("accEffRatio");
  errRatioH->Write("errorHRatio");
  errRatioL->Write("errorLRatio");

  //draw
  TCanvas* c = new TCanvas("c","",1000,1000);
  TH2F *pl = new TH2F ("pl",";y;pt",10, -2.4, 2.4, nptbins, ptbins);
  TLine* pt65 = new TLine(-1.6,6.5,1.6,6.5);
  pt65->SetLineColor(kRed);
  pt65->SetLineStyle(2);
  TLine* y1 = new TLine(-1.6,3,-1.6,6.5);
  y1->SetLineColor(kRed);
  y1->SetLineStyle(2);
  TLine* y2 = new TLine(1.6,3,1.6,6.5);
  y2->SetLineColor(kRed);
  y2->SetLineStyle(2);

  c->cd();
  pl->Draw();
  pl->SetTitle("pr/Npr AccxEff ratio");
  effRatio->Draw("colz same");
  pt65->Draw("same");
  y1->Draw("same");
  y2->Draw("same");
  c->SaveAs("RatioStudy/AccEff2DRatio.pdf");

  pl->Draw();
  pl->SetTitle("pr/Npr error high ratio");
  errRatioH->Draw("colz same");
  pt65->Draw("same");
  y1->Draw("same");
  y2->Draw("same");
  c->SaveAs("RatioStudy/errH2DRatio.pdf");

  pl->Draw();
  pl->SetTitle("pr/Npr error low ratio");
  errRatioL->Draw("colz same");
  pt65->Draw("same");
  y1->Draw("same");
  y2->Draw("same");
  c->SaveAs("RatioStudy/errL2DRatio.pdf");

  //fsave->Close();
  //delete fsave; delete prFile; delete nprFile; delete temp; delete effRatio; delete errRatioH; delete errRatioH; delete prEff; delete nprEff;
}//end of Eff2DRatios function

void Eff1DRatios(){
  Double_t ptbins []={3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35};
  int nptbins = ((sizeof(ptbins)/sizeof(double))-1);
  Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
  int netabins = ((sizeof(etabins)/sizeof(double))-1);
  Double_t ybins []={0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
  int nybins = ((sizeof(ybins)/sizeof(double))-1);

  string hisNames [] = {"AccEff_pt_016","AccEff_pt_1624","Acc_pt_016","Acc_pt_1624", "Eff_pt_016", "Eff_pt_1624"};
  string labels [] = {"pr-npr AccxEff in mid rapidity", "pr-npr AccxEff in fwd rapidity", "pr-npr Acc in mid rapidity", "pr-npr Acc in fwd rapidity", "pr-npr Eff in mid rapidity","pr-npr Eff in fwd rapidity"};

  gStyle->SetOptStat(0);
  TFile* prFile = TFile::Open("ANPlots/prAccEffCorr.root");
  TFile* nprFile = TFile::Open("ANPlots/nprAccEffCorr.root");
  TH1F* prhis = NULL;
  TH1F* nprhis = NULL;
  TCanvas* c = new TCanvas("c","",1000,1000);
  TH1F *hispt = new TH1F ("histpt",";pt", nptbins, ptbins);
  hispt->SetMinimum(0);
  hispt->SetMaximum(1.2);

  TLine* pt1 = new TLine(3,1,35,1);
  pt1->SetLineColor(kRed);
  pt1->SetLineStyle(2);
  TLegend* leg = NULL;

  for (int i =0; i<sizeof(hisNames)/sizeof(hisNames[0]); i++){
    prhis = (TH1F*) prFile->Get(hisNames[i].c_str());
    nprhis =(TH1F*) nprFile->Get(hisNames[i].c_str());
    
    prhis->SetMarkerStyle(33);
    prhis->SetMarkerColor(kMagenta+2);
    prhis->SetMarkerSize(2);
    //prhis->SetLineColor(39);
    
    nprhis->SetMarkerStyle(34);
    nprhis->SetMarkerColor(kBlue);
    nprhis->SetMarkerSize(2);
    //nprhis->SetLineColor(39);
    
    
    leg = new TLegend(0.6, 0.775, 0.8, 0.875);
    leg->AddEntry(prhis, "Prompt MC", "lep");
    leg->AddEntry(nprhis, "Nonprompt MC", "lep");
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    
    hispt->SetTitle(labels[i].c_str());
    c->cd();
    hispt->Draw();
    nprhis->Draw("same");
    prhis->Draw("same");
    pt1->Draw("same");
    leg->Draw("same");
    c->SaveAs(Form("RatioStudy/%s.pdf",hisNames[i].c_str()));
  }//end of the histograms loop
  prFile->Close();
  nprFile->Close();
  //delete prFile; delete nprFile; delete prhis; delete nprhis; delete c; delete leg; delete pt1; delete hispt;
}//end Eff1DRatios
