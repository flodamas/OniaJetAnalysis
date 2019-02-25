#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TBranchElement.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TProfile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;

void plot(bool doPrompt = true, bool doMid = false){

  string filename_nominal = "";
  string filename_up = "";
  string filename_down = "";

  if(doPrompt && doMid){
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_down.root";
  }
  
  if(doPrompt && !doMid){
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_down.root";
  }
  
  if(!doPrompt && doMid){
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_down.root";
  }
  
  if(!doPrompt && !doMid){
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_down.root";
  }

  TFile *file_nominal = new TFile(filename_nominal.c_str());
  TFile *file_up = new TFile(filename_up.c_str());
  TFile *file_down = new TFile(filename_down.c_str());

  
  TH2D *h2Unf_SI3_nominal = (TH2D*)file_nominal->Get("hReco_Iter3;1");
  TH2D *h2Unf_SI3_up = (TH2D*)file_up->Get("hReco_Iter3;1");
  TH2D *h2Unf_SI3_down = (TH2D*)file_down->Get("hReco_Iter3;1");
  
  TH1D *hZUnf_SI3_nominal;
  TH1D *hZUnf_SI3_up;
  TH1D *hZUnf_SI3_down;

  // take it from 2D output and project

  hZUnf_SI3_nominal=(TH1D*)h2Unf_SI3_nominal->ProjectionX("hZUnf_SI3_nominal",6,10);
  hZUnf_SI3_up=(TH1D*)h2Unf_SI3_up->ProjectionX("hZUnf_SI3_up",6,10);
  hZUnf_SI3_down=(TH1D*)h2Unf_SI3_down->ProjectionX("hZUnf_SI3_down",6,10);
  
  //rebin
  if(doMid){
    hZUnf_SI3_nominal->Rebin(7);
    hZUnf_SI3_up->Rebin(7);
    hZUnf_SI3_down->Rebin(7);
  }
  else{
    hZUnf_SI3_nominal->Rebin(10);
    hZUnf_SI3_up->Rebin(10);
    hZUnf_SI3_down->Rebin(10);
  }
  
  int binMin = 2;
  int binMax = 5;
  if(doMid) binMin = 4;
  if(doMid) binMax = 7;
  
  hZUnf_SI3_nominal->Scale(1/hZUnf_SI3_nominal->Integral(binMin,binMax));
  hZUnf_SI3_up->Scale(1/hZUnf_SI3_up->Integral(binMin,binMax));
  hZUnf_SI3_down->Scale(1/hZUnf_SI3_down->Integral(binMin,binMax));
  
  int zBins = 5;
  if(doMid) zBins = 7;

  double nominal_content = 0.;
  double up_content = 0.;
  double down_content = 0.;

  TH1D *hZUnf_SI3_up_ratioNominal;
  TH1D *hZUnf_SI3_down_ratioNominal;

  hZUnf_SI3_up_ratioNominal = (TH1D*)hZUnf_SI3_up->Clone();
  hZUnf_SI3_up_ratioNominal->Reset();
  //hZUnf_SI3_up_ratioNominal->Divide(hZUnf_SI3_nominal);

  hZUnf_SI3_down_ratioNominal = (TH1D*)hZUnf_SI3_down->Clone();
  hZUnf_SI3_down_ratioNominal->Reset();
  //hZUnf_SI3_down_ratioNominal->Divide(hZUnf_SI3_nominal);
    
  for(int i = 0; i < zBins; i++){

    nominal_content = hZUnf_SI3_nominal->GetBinContent(i+1);
    up_content = hZUnf_SI3_up->GetBinContent(i+1);
    down_content = hZUnf_SI3_down->GetBinContent(i+1);
    
    cout << "bin = " << i+1 << endl;
    if(nominal_content > 0) cout << "rel unc (in %) up = " << (nominal_content-up_content)*100/nominal_content  << " , down = " << (down_content-nominal_content)*100/nominal_content  << endl;

    double relEff_up = 0;
    double relEff_down = 0;

    if(nominal_content>0) relEff_up = (nominal_content-up_content)*100/nominal_content;
    if(nominal_content>0) relEff_down = (nominal_content-down_content)*100/nominal_content;
    
    hZUnf_SI3_up_ratioNominal->SetBinContent(i+1,relEff_up);
    hZUnf_SI3_down_ratioNominal->SetBinContent(i+1,relEff_down);
    
  }
  
  /*
  TH1D *hZUnf_SI3_up_ratioNominal;
  TH1D *hZUnf_SI3_down_ratioNominal;
  
  hZUnf_SI3_up_ratioNominal = (TH1D*)hZUnf_SI3_up->Clone();
  hZUnf_SI3_up_ratioNominal->Divide(hZUnf_SI3_nominal);

  hZUnf_SI3_down_ratioNominal = (TH1D*)hZUnf_SI3_down->Clone();
  hZUnf_SI3_down_ratioNominal->Divide(hZUnf_SI3_nominal);
  */
  
  TCanvas * mycan1 = new TCanvas("mycan1","mycan1",900,500);
  mycan1->Divide(2,1);
  
  TLegend *legend = new TLegend(0.7,0.6,.85,0.9,"","brNDC");
  legend->SetHeader("");

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);

  mycan1->cd(1);

  //gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.1);
  
  hZUnf_SI3_nominal->SetStats(0);

  hZUnf_SI3_nominal->SetMaximum(hZUnf_SI3_nominal->GetMaximum()*1.1);
  
  hZUnf_SI3_nominal->GetXaxis()->SetTitle("z");
  hZUnf_SI3_nominal->GetYaxis()->SetTitleOffset(1.9);
  hZUnf_SI3_nominal->GetYaxis()->SetTitle("1/N dN/dz");

  hZUnf_SI3_nominal->SetLineColor(kBlack);
  hZUnf_SI3_nominal->SetMarkerColor(kBlack);
  hZUnf_SI3_nominal->SetMarkerStyle(28);
  hZUnf_SI3_nominal->SetMarkerSize(1.5);
  
  hZUnf_SI3_up->SetLineColor(kgreen);
  hZUnf_SI3_down->SetLineColor(kred);

  hZUnf_SI3_up->SetMarkerColor(kgreen);
  hZUnf_SI3_down->SetMarkerColor(kred);
  
  hZUnf_SI3_up->SetMarkerStyle(24);
  hZUnf_SI3_down->SetMarkerStyle(22);
  
  hZUnf_SI3_up->SetMarkerSize(1.5);
  hZUnf_SI3_down->SetMarkerSize(1.5);

  hZUnf_SI3_down->SetLineWidth(2);
  hZUnf_SI3_up->SetLineWidth(2);
  hZUnf_SI3_nominal->SetLineWidth(2);
  
  hZUnf_SI3_nominal->Draw("HIST");
  hZUnf_SI3_up->Draw("HISTsame");
  hZUnf_SI3_down->Draw("HISTsame");
  
  legend->AddEntry(hZUnf_SI3_nominal, "nominal","l");
  legend->AddEntry(hZUnf_SI3_up, "up","l");
  legend->AddEntry(hZUnf_SI3_down, "down","l");
  
  legend->Draw("same");

  mycan1->Update();

  float xCoord;
  if(doMid) xCoord = 0.44;
  else xCoord = 0.2;
    
  TLine *line0 = new TLine(xCoord,0,xCoord,gPad->GetUymax());
  line0->SetLineColor(kred);
  line0->SetLineStyle(2);
  line0->SetLineWidth(2);
  line0->Draw("same");
    
  mycan1->cd(2);
  
  //gPad->SetRightMargin(0.2);

  hZUnf_SI3_up_ratioNominal->SetStats(0);
  hZUnf_SI3_up_ratioNominal->GetYaxis()->SetTitle("relative error in %");
  hZUnf_SI3_up_ratioNominal->GetYaxis()->SetRangeUser(-10.,10);
  hZUnf_SI3_up_ratioNominal->GetXaxis()->SetTitle("z");

  hZUnf_SI3_up_ratioNominal->SetLineColor(kgreen);
  hZUnf_SI3_down_ratioNominal->SetLineColor(kred);

  hZUnf_SI3_up_ratioNominal->SetMarkerColor(kgreen);
  hZUnf_SI3_down_ratioNominal->SetMarkerColor(kred);
  
  hZUnf_SI3_down_ratioNominal->SetLineWidth(2);
  hZUnf_SI3_up_ratioNominal->SetLineWidth(2);
      
  hZUnf_SI3_up_ratioNominal->Draw("HIST");
  hZUnf_SI3_down_ratioNominal->Draw("HISTsame");
  
  mycan1->Update();

  TLine *line1 = new TLine(xCoord,0,xCoord,gPad->GetUymax());
  line1->SetLineColor(kred);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");

  mycan1->Update();
  
  if(doPrompt && doMid){
    mycan1->SaveAs("../../plots/unf_data_prompt_mid_JER_Scaled.pdf");
  }

  if(doPrompt && !doMid){
    mycan1->SaveAs("../../plots/unf_data_prompt_fwd_JER_Scaled.pdf");
  }

  if(!doPrompt && doMid){
    mycan1->SaveAs("../../plots/unf_data_nonprompt_mid_JER_Scaled.pdf");
  }

  if(!doPrompt && !doMid){
    mycan1->SaveAs("../../plots/unf_data_nonprompt_fwd_JER_Scaled.pdf");

  }
  
  
}

void PlotRatios_JER(){

  plot(true,true);

  //plot(true,false);
  
  //plot(false,true);
  //plot(false,false);
  
}
