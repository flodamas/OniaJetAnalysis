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

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string filenameMeas = "";
  string outputfile = "";

  if(doPrompt && doMid){
    cout << "prompt mid" << endl;
    filename1 = "../unfOutput/step1/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename2 = "../unfOutput/step2/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename3 = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename4 = "../unfOutput/step4/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filenameMeas = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    outputfile = "unfResult_prompt_mid.root";
  }
  
  if(doPrompt && !doMid){
    cout << "prompt fwd" << endl;
    filename1 ="../unfOutput/step1/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename2 ="../unfOutput/step2/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename3 ="../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename4 ="../unfOutput/step4/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filenameMeas = "../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal.root";
    outputfile = "unfResult_prompt_fwd.root";
  }
  
  if(!doPrompt && doMid){
    cout << "prompt mid" << endl;
    filename1 = "../unfOutput/step1/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename2 = "../unfOutput/step2/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename3 = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename4 = "../unfOutput/step4/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filenameMeas = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    outputfile = "unfResult_nonprompt_mid.root";
  }

  if(!doPrompt && !doMid){
    cout << "prompt fwd" << endl;
    filename1 ="../unfOutput/step1/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename2 ="../unfOutput/step2/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename3 ="../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename4 ="../unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filenameMeas = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal.root";
    outputfile = "unfResult_nonprompt_fwd.root";
  }
  
  
  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  TFile *fileMeas = new TFile(filenameMeas.c_str());
  
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)fileMeas->Get("fh2MeasData;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",2,2);
  hZMeas->Draw();
  
  TH2D *h2UnfResp1 = (TH2D*)file1->Get("hReco_Iter1;1");
  TH2D *h2UnfResp2 = (TH2D*)file2->Get("hReco_Iter1;1");
  TH2D *h2UnfResp3 = (TH2D*)file3->Get("hReco_Iter1;1");
  TH2D *h2UnfResp4 = (TH2D*)file4->Get("hReco_Iter1;1");
  
  TH1D *hZUnf_SI1 = (TH1D*)h2UnfResp1->ProjectionX("hZUnf_SI1",2,2);
  TH1D *hZUnf_SI2 = (TH1D*)h2UnfResp2->ProjectionX("hZUnf_SI2",2,2);
  TH1D *hZUnf_SI3 = (TH1D*)h2UnfResp3->ProjectionX("hZUnf_SI3",2,2);
  TH1D *hZUnf_SI4 = (TH1D*)h2UnfResp4->ProjectionX("hZUnf_SI4",2,2);

  /*
  int zBins = 5;
  cout << "after rebin" << endl;
  for(int i = 0; i < zBins; i++){
    //cout << "bin = " << i << " bin content = " << hZUnf_SI3->GetBinContent(i+1) << " bin err = " <<  hZUnf_SI3->GetBinError(i+1) << " rel uncert = " << hZUnf_SI3->GetBinError(i+1)*100/hZUnf_SI3->GetBinContent(i+1) << endl;
    cout << "bin = " << i << " rel unc (in%) = " << hZUnf_SI3->GetBinError(i+1)*100/hZUnf_SI3->GetBinContent(i+1) << endl;
  }
  */
  
  /*
  cout << "measured" << endl;
  
  for(int i = 0; i < zBins; i++){
    cout << "bin = " << i << " bin content = " << hZMeas->GetBinContent(i+1) << " bin err = " <<  hZMeas->GetBinError(i+1) << " rel uncert = " << hZMeas->GetBinError(i+1)*100/hZMeas->GetBinContent(i+1) << endl;
  }
  */

  
  /*
  TFile *file_out = new TFile(outputfile.c_str(),"RECREATE");
  hZUnf_SI3->Write("zUnf");
  hZMeas->Write("zMeas");
  file_out->Close();
  */
  
  TH1D *hZUnf_SI1_ratioMeas;
  TH1D *hZUnf_SI2_ratioMeas;
  TH1D *hZUnf_SI3_ratioMeas;
  TH1D *hZUnf_SI4_ratioMeas;
  
  hZUnf_SI1_ratioMeas = (TH1D*)hZUnf_SI1->Clone("hZUnf_SI1_ratioMeas");
  hZUnf_SI1_ratioMeas->Divide(hZMeas);

  hZUnf_SI2_ratioMeas = (TH1D*)hZUnf_SI2->Clone("hZUnf_SI2_ratioMeas");
  hZUnf_SI2_ratioMeas->Divide(hZMeas);

  hZUnf_SI3_ratioMeas = (TH1D*)hZUnf_SI3->Clone("hZUnf_SI3_ratioMeas");
  hZUnf_SI3_ratioMeas->Divide(hZMeas);

  hZUnf_SI4_ratioMeas = (TH1D*)hZUnf_SI4->Clone("hZUnf_SI4_ratioMeas");
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

  if(doPrompt && !doMid) hZMeas->SetMaximum(hZMeas->GetMaximum()*1.5);
  else hZMeas->SetMaximum(hZMeas->GetMaximum()*1.3);
  hZMeas->GetXaxis()->SetTitle("z");
  hZMeas->GetYaxis()->SetTitleOffset(1.9);
  hZMeas->GetYaxis()->SetTitle("dN/dz (per 0.2)");
  if(doMid) hZMeas->GetYaxis()->SetTitle("dN/dz (per 0.14)");

  hZMeas->SetLineColor(kBlack);
  hZMeas->SetMarkerColor(kBlack);
  hZMeas->SetMarkerStyle(28);
  hZMeas->SetMarkerSize(1.5);
  
  hZUnf_SI1->SetLineColor(kgreen);
  hZUnf_SI2->SetLineColor(kazure);
  hZUnf_SI3->SetLineColor(kred);
  hZUnf_SI4->SetLineColor(kpink);

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
  hZUnf_SI1->Draw("EPsame");
  hZUnf_SI2->Draw("EPsame");
  hZUnf_SI3->Draw("EPsame");
  hZUnf_SI4->Draw("EPsame");
  
  legend->AddEntry(hZMeas, "meas","ep");
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
  hZUnf_SI1_ratioMeas->GetYaxis()->SetTitle("ratio to measured");
  hZUnf_SI1_ratioMeas->GetYaxis()->SetRangeUser(0.,1.4);
  //  if(doMid) hZUnf_SI1_ratioMeas->GetXaxis()->SetRangeUser(0.4,1.);
  //else hZUnf_SI1_ratioMeas->GetXaxis()->SetRangeUser(0.2,1.);
  hZUnf_SI1_ratioMeas->GetXaxis()->SetTitle("z");

  hZUnf_SI1_ratioMeas->SetLineColor(kgreen);
  hZUnf_SI2_ratioMeas->SetLineColor(kazure);
  hZUnf_SI3_ratioMeas->SetLineColor(kred);
  hZUnf_SI4_ratioMeas->SetLineColor(kpink);

  hZUnf_SI1_ratioMeas->SetMarkerColor(kgreen);
  hZUnf_SI2_ratioMeas->SetMarkerColor(kazure);
  hZUnf_SI3_ratioMeas->SetMarkerColor(kred);
  hZUnf_SI4_ratioMeas->SetMarkerColor(kpink);
  
  hZUnf_SI3_ratioMeas->SetLineStyle(1);
  hZUnf_SI3_ratioMeas->SetLineWidth(2);

  hZUnf_SI4_ratioMeas->SetLineStyle(8);
  hZUnf_SI4_ratioMeas->SetLineWidth(2);
      
  hZUnf_SI1_ratioMeas->Draw("HIST");
  hZUnf_SI2_ratioMeas->Draw("HISTsame");
  hZUnf_SI3_ratioMeas->Draw("HISTsame");
  hZUnf_SI4_ratioMeas->Draw("HISTsame");
  
  mycan1->Update();

  TLine *line1 = new TLine(xCoord,0,xCoord,gPad->GetUymax());
  line1->SetLineColor(kred);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");

  mycan1->Update();

  if(doPrompt && doMid){
    mycan1->SaveAs("../plots/unf_data_prompt_mid_49z15ptBins7zMeasBins_old_nominal_Diag.pdf");
  }

  if(doPrompt && !doMid){
    mycan1->SaveAs("../plots/unf_data_prompt_fwd_50z15ptBins_old_nominal_Diag.pdf");
  }

  if(!doPrompt && doMid){
    mycan1->SaveAs("../plots/unf_data_nonprompt_mid_49z15ptBins7zMeasBins_nominal_Diag.pdf");
  }

  if(!doPrompt && !doMid){
    mycan1->SaveAs("../plots/unf_data_nonprompt_fwd_50z15ptBins_nominal_Diag.pdf");
  }
  
  
}
void PlotRatios_DataUnfolded_afterDiag(){

  //plot(true,true);
  //plot(true,false);

  //plot(false,true);
  plot(false,false);
}
