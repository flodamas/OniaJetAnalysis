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

void plot(bool doPrompt = false, bool doMid = true){

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";

  if(doPrompt && !doMid){
    filename1 ="../unfOutput/step1/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename2 ="../unfOutput/step2/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename3 ="../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename4 ="../unfOutput/step4/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBinsDiag.root";
  }

  if(!doPrompt && !doMid){
    filename1 ="../unfOutput/step1/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename2 ="../unfOutput/step2/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename3 ="../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBinsDiag.root";
    filename4 ="../unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBinsDiag.root";
  }
  
  
  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  /*
  TH1D *hZMeas;
  hZMeas=(TH1D*)file1->Get("hMSme_1;1");
  */
  
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)file3->Get("mcMeas;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",2,2);
  
  TH2D *h2ZTrue;
  h2ZTrue=(TH2D*)file3->Get("mcPrior;1");
  TH1D *hZTrue = h2ZTrue->ProjectionX("hZTrue",6,10);
  hZTrue->Rebin(10);
  
  /*
  TH1D *hZTrue;
  hZTrue=(TH1D*)file1->Get("hMTru_1;1");
  */

  TH2D *h2ZUnf_SI1 = (TH2D*)file1->Get("hReco_Iter1;1");
  TH2D *h2ZUnf_SI2 = (TH2D*)file2->Get("hReco_Iter1;1");
  TH2D *h2ZUnf_SI3 = (TH2D*)file3->Get("hReco_Iter1;1");
  TH2D *h2ZUnf_SI4 = (TH2D*)file4->Get("hReco_Iter1;1");
    
  TH1D *hZUnf_SI1 = (TH1D*)h2ZUnf_SI1->ProjectionX("hZUnf_SI1",2,2);
  TH1D *hZUnf_SI2 = (TH1D*)h2ZUnf_SI2->ProjectionX("hZUnf_SI2",2,2);
  TH1D *hZUnf_SI3 = (TH1D*)h2ZUnf_SI3->ProjectionX("hZUnf_SI3",2,2);
  TH1D *hZUnf_SI4 = (TH1D*)h2ZUnf_SI4->ProjectionX("hZUnf_SI4",2,2);

  /*
  hZUnf_SI1=(TH1D*)file1->Get("hMUnf_1_Iter3;1");
  hZUnf_SI2=(TH1D*)file2->Get("hMUnf_1_Iter3;1");
  hZUnf_SI3=(TH1D*)file3->Get("hMUnf_1_Iter3;1");
  hZUnf_SI4=(TH1D*)file4->Get("hMUnf_1_Iter3;1");
  */
  
  /*
  hZUnf_SI1->Scale(1/hZUnf_SI1->Integral());
  hZUnf_SI2->Scale(1/hZUnf_SI2->Integral());
  hZUnf_SI3->Scale(1/hZUnf_SI3->Integral());
  */

  TH1D *hZTrue_ratioMeas;
  TH1D *hZUnf_SI1_ratioMeas;
  TH1D *hZUnf_SI2_ratioMeas;
  TH1D *hZUnf_SI3_ratioMeas;
  TH1D *hZUnf_SI4_ratioMeas;
  
  hZTrue_ratioMeas = (TH1D*)hZTrue->Clone();
  hZTrue_ratioMeas->Divide(hZMeas);  
  
  hZUnf_SI1_ratioMeas = (TH1D*)hZUnf_SI1->Clone();
  hZUnf_SI1_ratioMeas->Divide(hZMeas);

  hZUnf_SI2_ratioMeas = (TH1D*)hZUnf_SI2->Clone();
  hZUnf_SI2_ratioMeas->Divide(hZMeas);

  hZUnf_SI3_ratioMeas = (TH1D*)hZUnf_SI3->Clone();
  hZUnf_SI3_ratioMeas->Divide(hZMeas);

  hZUnf_SI4_ratioMeas = (TH1D*)hZUnf_SI4->Clone();
  hZUnf_SI4_ratioMeas->Divide(hZMeas);
  
  
  TCanvas * mycan1 = new TCanvas("mycan1","mycan1",900,500);
  mycan1->Divide(2,1);
  
  TLegend *legend = new TLegend(0.8,0.15,.95,0.9,"","brNDC");
  legend->SetHeader("");

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);

  mycan1->cd(1);

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.2);
  
  hZMeas->SetStats(0);
  if(doPrompt & doMid){
    hZMeas->SetTitle("prompt mid-rapidity");
  }
  if(doPrompt & !doMid){
    hZMeas->SetTitle("prompt fow-rapidity");
  }
  if(!doPrompt & doMid){
    hZMeas->SetTitle("nonprompt mid-rapidity");
  }
  if(!doPrompt & !doMid){
    hZMeas->SetTitle("nonprompt fow-rapidity");
  }  
  
  //  hZMeas->GetYaxis()->SetRangeUser(0.,.6);
  //  if(doMid) hZMeas->GetXaxis()->SetRangeUser(0.4,1.);
  //else hZMeas->GetXaxis()->SetRangeUser(0.2,1.);
  hZMeas->GetXaxis()->SetTitle("z");
  hZMeas->GetYaxis()->SetTitleOffset(2.);
  hZMeas->GetYaxis()->SetTitle("dN/dz (per 0.2)");

  hZMeas->SetMaximum(hZMeas->GetMaximum()*1.2);
  
  hZMeas->SetLineColor(kBlack);
  hZMeas->SetMarkerColor(kBlack);
  hZMeas->SetMarkerStyle(28);
  hZMeas->SetMarkerSize(1.5);

  hZTrue->SetLineColor(kblue);
  hZTrue->SetMarkerColor(kblue);
  hZTrue->SetMarkerStyle(20);
  hZTrue->SetMarkerSize(1.4);
    
  hZUnf_SI1->SetLineColor(kgreen);
  hZUnf_SI2->SetLineColor(kazure);
  hZUnf_SI3->SetLineColor(kred);
  hZUnf_SI4->SetLineColor(kpink);
    
  hZUnf_SI1->SetMarkerColor(kgreen);
  hZUnf_SI2->SetMarkerColor(kazure);
  hZUnf_SI3->SetMarkerColor(kred);
  hZUnf_SI4->SetMarkerColor(kpink);

  hZUnf_SI1->SetMarkerColor(kgreen);
  hZUnf_SI2->SetMarkerColor(kazure);
  hZUnf_SI3->SetMarkerColor(kred);
  hZUnf_SI4->SetMarkerColor(kpink);

  hZUnf_SI1->SetMarkerStyle(24);
  hZUnf_SI2->SetMarkerStyle(25);
  hZUnf_SI3->SetMarkerStyle(22);
  hZUnf_SI4->SetMarkerStyle(30);
  
  hZUnf_SI1->SetMarkerSize(1.5);
  hZUnf_SI2->SetMarkerSize(1.5);
  hZUnf_SI3->SetMarkerSize(1.5);
  hZUnf_SI4->SetMarkerSize(1.5);
  
  
  hZMeas->Draw("EP");
  hZTrue->Draw("EPsame");
  hZUnf_SI1->Draw("EPsame");
  hZUnf_SI2->Draw("EPsame");
  hZUnf_SI3->Draw("EPsame");
  hZUnf_SI4->Draw("EPsame");
  
  legend->AddEntry(hZMeas, "meas","ep");
  legend->AddEntry(hZTrue, "truth","ep");
  legend->AddEntry(hZUnf_SI1, "unf SI#1","ep");
  legend->AddEntry(hZUnf_SI2, "unf SI#2","ep");
  legend->AddEntry(hZUnf_SI3, "unf SI#3","ep");
  legend->AddEntry(hZUnf_SI4, "unf SI#4","ep");
  
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

  hZUnf_SI1_ratioMeas->SetStats(0);
  //  hZUnf_SI1_ratioMeas->SetTitle("ratio to measured");
  hZUnf_SI1_ratioMeas->GetYaxis()->SetTitle("ratio to measured");
  hZUnf_SI1_ratioMeas->GetYaxis()->SetTitleOffset(1.6);
  hZUnf_SI1_ratioMeas->GetYaxis()->SetRangeUser(0.,1.4);
  hZUnf_SI1_ratioMeas->GetXaxis()->SetTitle("z");

  hZTrue_ratioMeas->SetLineColor(kblue);
  hZTrue_ratioMeas->SetMarkerColor(kblue);
  hZTrue_ratioMeas->SetMarkerStyle(27);
  hZTrue_ratioMeas->SetMarkerSize(1.);
    
  hZUnf_SI1_ratioMeas->SetLineColor(kgreen);
  hZUnf_SI2_ratioMeas->SetLineColor(kazure);
  hZUnf_SI3_ratioMeas->SetLineColor(kred);
  hZUnf_SI4_ratioMeas->SetLineColor(kpink);
  
  hZTrue_ratioMeas->SetLineStyle(7);
  hZUnf_SI3_ratioMeas->SetLineStyle(1);
  hZUnf_SI4_ratioMeas->SetLineStyle(8);

  hZTrue_ratioMeas->SetLineWidth(2);
  hZUnf_SI3_ratioMeas->SetLineWidth(2);
  hZUnf_SI4_ratioMeas->SetLineWidth(2);
    
  hZUnf_SI1_ratioMeas->SetMarkerColor(kgreen);
  hZUnf_SI2_ratioMeas->SetMarkerColor(kazure);
  hZUnf_SI3_ratioMeas->SetMarkerColor(kpink);
  hZUnf_SI4_ratioMeas->SetMarkerColor(kred);
  
  hZUnf_SI1_ratioMeas->Draw("HIST");
  hZTrue_ratioMeas->Draw("HISTsame");
  hZUnf_SI2_ratioMeas->Draw("HISTsame");
  hZUnf_SI3_ratioMeas->Draw("HISTsame");
  hZUnf_SI4_ratioMeas->Draw("HISTsame");
    
  
  mycan1->Update();

  float xCoord2;
  if(doMid) xCoord2 = 0.44;
  else xCoord2 = 0.2;
    
  TLine *line1 = new TLine(xCoord2,0,xCoord2,gPad->GetUymax());
  line1->SetLineColor(kred);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");

  mycan1->Update();
    
  if(doPrompt && !doMid){
    mycan1->SaveAs("../plots/unf_mc_prompt_fwd.pdf");
  }

  if(!doPrompt && !doMid){
    mycan1->SaveAs("../plots/unf_mc_nonprompt_fwd.pdf");
  }
  
  
}

void PlotRatios_MCUnfoldedTruth_Fwd(){

  plot(true,false);
  //plot(false,false);
  
}
