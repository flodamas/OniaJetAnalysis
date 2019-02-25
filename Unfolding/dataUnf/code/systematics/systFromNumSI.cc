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

void compute(bool doPrompt = false, bool doMid = true){

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string outputfile = "";

  if(doPrompt && doMid){
    filename1 = "../../unfOutput/step1/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename2 = "../../unfOutput/step2/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename3 = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename4 = "../../unfOutput/step4/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    outputfile = "../../unfOutput/numberSIErrs/SIerr_Prompt_Mid_newNominal.root";
  }
  
  if(doPrompt && !doMid){
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    outputfile = "../../unfOutput/numberSIErrs/SIerr_Prompt_Fwd_newNominal.root";
  }

  if(!doPrompt && doMid){
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    outputfile = "../../unfOutput/numberSIErrs/SIerr_NonPrompt_Mid.root";
  }
  
  if(!doPrompt && !doMid){
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    outputfile = "../../unfOutput/numberSIErrs/SIerr_NonPrompt_Fwd.root";
  }
  
  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  /*
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)file1->Get("fh2MeasData;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",2,2);
  */
  
  TH2D *h2UnfResp1 = (TH2D*)file1->Get("hReco_Iter1;1");
  TH2D *h2UnfResp2 = (TH2D*)file2->Get("hReco_Iter1;1");
  TH2D *h2UnfResp3 = (TH2D*)file3->Get("hReco_Iter1;1");
  TH2D *h2UnfResp4 = (TH2D*)file4->Get("hReco_Iter1;1");

  TH1D *hZUnf_SI1 = (TH1D*)h2UnfResp1->ProjectionX("hZUnf_SI1",2,2);
  TH1D *hZUnf_SI2 = (TH1D*)h2UnfResp2->ProjectionX("hZUnf_SI2",2,2);
  TH1D *hZUnf_SI3 = (TH1D*)h2UnfResp3->ProjectionX("hZUnf_SI3",2,2);
  TH1D *hZUnf_SI4 = (TH1D*)h2UnfResp4->ProjectionX("hZUnf_SI4",2,2);

  TH1D *hZUnf_SI2_SI3;
  TH1D *hZUnf_SI4_SI3;
  
  hZUnf_SI2_SI3 = (TH1D*)hZUnf_SI2->Clone("hZUnf_SI2_SI3");
  hZUnf_SI2_SI3->Divide(hZUnf_SI3);

  hZUnf_SI4_SI3 = (TH1D*)hZUnf_SI4->Clone("hZUnf_SI4_SI3");
  hZUnf_SI4_SI3->Divide(hZUnf_SI3);

  TH1D *hZUnf_SI1_scale = (TH1D*)hZUnf_SI1->Clone("hZUnf_SI1_scale");
  //  hZUnf_SI1_scale->Reset();
  
  TH1D *hZUnf_SI2_scale = (TH1D*)hZUnf_SI2->Clone("hZUnf_SI2_scale");
  //hZUnf_SI2_scale->Reset();
  
  TH1D *hZUnf_SI3_scale = (TH1D*)hZUnf_SI3->Clone("hZUnf_SI3_scale");
  //hZUnf_SI3_scale->Reset();
  
  TH1D *hZUnf_SI4_scale = (TH1D*)hZUnf_SI4->Clone("hZUnf_SI4_scale");
  //hZUnf_SI4_scale->Reset();  

  int binMin = 2;
  int binMax = 5;

  if(doMid) binMin = 4;
  if(doMid) binMax = 7;
    
  hZUnf_SI1_scale->Scale(1./hZUnf_SI1_scale->Integral(binMin,binMax));
  hZUnf_SI2_scale->Scale(1./hZUnf_SI2_scale->Integral(binMin,binMax));
  hZUnf_SI3_scale->Scale(1./hZUnf_SI3_scale->Integral(binMin,binMax));
  hZUnf_SI4_scale->Scale(1./hZUnf_SI4_scale->Integral(binMin,binMax));

  TH1D *hZUnf_SI2_SI3_scale = (TH1D*)hZUnf_SI2_scale->Clone();
  hZUnf_SI2_SI3_scale->Divide(hZUnf_SI3_scale);

  TH1D *hZUnf_SI4_SI3_scale = (TH1D*)hZUnf_SI4_scale->Clone();
  hZUnf_SI4_SI3_scale->Divide(hZUnf_SI3_scale);
    
  float binContUnf = 0;
  float ratio1 = 0;
  float ratio2 = 0;
  float binErrUnf = 0;

  float binContUnf_scale = 0;
  float ratio1_scale = 0;
  float ratio2_scale = 0;
  float binErrUnf_scale = 0;
   
  TH1D *hZUnf_SI3_newUncert = (TH1D*)hZUnf_SI3->Clone("hZUnf_SI3_newUncert");
  hZUnf_SI3_newUncert->Reset();

  TH1D *hZUnf_SI3_newUncert_scale = (TH1D*)hZUnf_SI3_scale->Clone("hZUnf_SI3_newUncert_scale");
  hZUnf_SI3_newUncert_scale->Reset();
    
  int nBins = 5;
  if(doMid) nBins = 7;
  
  for(int i = 0; i < nBins; i++){

    binContUnf = hZUnf_SI3->GetBinContent(i+1);
    ratio1 = hZUnf_SI2_SI3->GetBinContent(i+1);
    ratio2 = hZUnf_SI4_SI3->GetBinContent(i+1);

    //cout << "bin i = " << i << " binContUnf = " << binContUnf << " ; 1 - ratio1 = " << 1-ratio1 << " ; 1 - ratio2 = " << 1-ratio2 << endl;
    
    if(abs(1-ratio1) > abs(1-ratio2)) {
      binErrUnf = abs(1-ratio1)*binContUnf;
      cout << "rel err (in %) = " << (1-ratio1)*100 << endl;
    }
    else {
      binErrUnf = abs(1-ratio2)*binContUnf;
      cout << "rel err (in %) = " << (1-ratio2)*100 << endl;
    }

    hZUnf_SI3_newUncert->SetBinContent(i+1,binContUnf);
    hZUnf_SI3_newUncert->SetBinError(i+1,binErrUnf);

    binContUnf_scale = hZUnf_SI3_scale->GetBinContent(i+1);
    ratio1_scale = hZUnf_SI2_SI3_scale->GetBinContent(i+1);
    ratio2_scale = hZUnf_SI4_SI3_scale->GetBinContent(i+1);

    //    cout << "bin i = " << i << " binContUnf = " << binContUnf_scale << " ; 1 - ratio1 = " << 1-ratio1_scale << " ; 1 - ratio2 = " << 1-ratio2_scale << endl;

    cout << "i = " << i << " , bin content = " << binContUnf << endl;
    
    if(abs(1-ratio1_scale) > abs(1-ratio2_scale)) {
      binErrUnf_scale = abs(1-ratio1_scale)*binContUnf;
      cout << "scaled rel err (in %) = " << (1-ratio1_scale)*100 << endl;
    }
    else {
      binErrUnf_scale = abs(1-ratio2_scale)*binContUnf;
      cout << "scaled rel err (in %) = " << (1-ratio2_scale)*100 << endl;
    }

    hZUnf_SI3_newUncert_scale->SetBinContent(i+1,binContUnf);
    cout << "binContUnf = " << binContUnf << endl;
    hZUnf_SI3_newUncert_scale->SetBinError(i+1,binErrUnf_scale);
        
  }
  

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
  
  hZUnf_SI2_scale->SetStats(0);
  if(doPrompt & doMid){
    hZUnf_SI2_scale->SetTitle("prompt mid-rapidity");
  }
  if(doPrompt & !doMid){
    hZUnf_SI2_scale->SetTitle("prompt fwd-rapidity");
  }
  if(!doPrompt & doMid){
    hZUnf_SI2_scale->SetTitle("nonprompt mid-rapidity");
  }
  if(!doPrompt & !doMid){
    hZUnf_SI2_scale->SetTitle("nonprompt fwd-rapidity");
  }  

  if(doPrompt && !doMid) hZUnf_SI2_scale->SetMaximum(hZUnf_SI2_scale->GetMaximum()*1.5);
  else hZUnf_SI2_scale->SetMaximum(hZUnf_SI2_scale->GetMaximum()*1.3);

  hZUnf_SI2_scale->GetXaxis()->SetTitle("z");
  hZUnf_SI2_scale->GetYaxis()->SetTitleOffset(1.9);
  hZUnf_SI2_scale->GetYaxis()->SetTitle("dN/dz (per 0.2)");
  if(doMid) hZUnf_SI2_scale->GetYaxis()->SetTitle("dN/dz (per 0.14)");

  hZUnf_SI2_scale->SetMarkerColor(kazure);
  hZUnf_SI3_scale->SetMarkerColor(kred);
  hZUnf_SI4_scale->SetMarkerColor(kpink);
  

  hZUnf_SI2_scale->SetMarkerStyle(25);
  hZUnf_SI3_scale->SetMarkerStyle(22);
  hZUnf_SI4_scale->SetMarkerStyle(30);
  

  hZUnf_SI2_scale->SetMarkerSize(1.5);
  hZUnf_SI3_scale->SetMarkerSize(1.5);
  hZUnf_SI4_scale->SetMarkerSize(1.5);
  

  hZUnf_SI2_scale->Draw("EP");
  hZUnf_SI3_scale->Draw("EPsame");
  hZUnf_SI4_scale->Draw("EPsame");
  
  legend->AddEntry(hZUnf_SI2_scale, "unf SI#2","ep");
  legend->AddEntry(hZUnf_SI3_scale, "unf SI#3","ep");
  legend->AddEntry(hZUnf_SI4_scale, "unf SI#4","ep");
  
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
  
  hZUnf_SI2_SI3_scale->SetStats(0);
  hZUnf_SI2_SI3_scale->GetYaxis()->SetTitle("ratio to SI3");
  hZUnf_SI2_SI3_scale->GetYaxis()->SetRangeUser(0.,1.4);

  hZUnf_SI2_SI3_scale->GetXaxis()->SetTitle("z");

  hZUnf_SI2_SI3_scale->SetLineColor(kazure);
  hZUnf_SI4_SI3_scale->SetLineColor(kpink);

  hZUnf_SI2_SI3_scale->SetMarkerColor(kazure);
  hZUnf_SI4_SI3_scale->SetMarkerColor(kpink);
  
  hZUnf_SI4_SI3_scale->SetLineStyle(8);
  hZUnf_SI4_SI3_scale->SetLineWidth(2);
      
  hZUnf_SI2_SI3_scale->Draw("HIST");
  hZUnf_SI4_SI3_scale->Draw("HISTsame");
  
  mycan1->Update();

  TLine *line1 = new TLine(xCoord,0,xCoord,gPad->GetUymax());
  line1->SetLineColor(kred);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");

  mycan1->Update();
  
  if(doPrompt && doMid){
    mycan1->SaveAs("../../plots/unf_data_prompt_mid_SI_syst_scaled.png");
  }

  if(doPrompt && !doMid){
    mycan1->SaveAs("../../plots/unf_data_prompt_fwd_SI_syst_scaled.png");
  }

  if(!doPrompt && doMid){
    mycan1->SaveAs("../plots/unf_data_nonprompt_mid_SI_syst_scaled.png");
  }

  if(!doPrompt && !doMid){
    mycan1->SaveAs("../plots/unf_data_nonprompt_fwd_SI_syst_scaled.png");
  }
  

  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  
  hZUnf_SI2->Write("zUnfSI2");
  hZUnf_SI3->Write("zUnfSI3");
  hZUnf_SI4->Write("zUnfSI4");
  hZUnf_SI2_SI3->Write("zUnfSI2ovSI3");
  hZUnf_SI4_SI3->Write("zUnfSI4ovSI3");
  hZUnf_SI3_newUncert->Write("zUnfSI3_newUncert");
  hZUnf_SI3_newUncert_scale->Write("zUnfSI3_newUncert_scaledErr");
  
  outfile->Close();
  
  
}

void systFromNumSI(){

  //compute(true,true);
  compute(true,false);

  //compute(false,true);
  //compute(false,false);
    
  
}
