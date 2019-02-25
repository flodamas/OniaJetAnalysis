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

  string filename_newNominal = "";
  string filename_jesPr = "";
  string filename_jesNonPr = "";
  string outputfile = "";
  
  if(doPrompt && doMid){
    cout << "prompt mid" << endl;
    filename_jesPr = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename_newNominal = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename_jesNonPr = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nonPrJES_nominal_Diag.root";
    outputfile = "../unfOutput/JERJESErrs/addJESerr_Prompt_Mid_newNominal.root";
  }
  
  if(doPrompt && !doMid){
    cout << "prompt fwd" << endl;
    filename_jesPr ="../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename_newNominal = "../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename_jesNonPr = "../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nonPrJES_nominal_Diag.root";
    outputfile = "../unfOutput/JERJESErrs/addJESerr_Prompt_Fwd_newNominal.root";
  }
  
  TFile *file_jesPr = new TFile(filename_jesPr.c_str());
  TFile *file_newNominal = new TFile(filename_newNominal.c_str());
  TFile *file_jesNonPr = new TFile(filename_jesNonPr.c_str());
  
  // after errors got back
  
  TH2D *h2Unf_jesPr = (TH2D*)file_jesPr->Get("hReco_Iter1;1");
  TH2D *h2Unf_newNominal = (TH2D*)file_newNominal->Get("hReco_Iter1;1");
  TH2D *h2Unf_jesNonPr = (TH2D*)file_jesNonPr->Get("hReco_Iter1;1");

  TH1D *hZUnf_jesPr = (TH1D*)h2Unf_jesPr->ProjectionX("hZUnf_jesPr",2,2);
  TH1D *hZUnf_jesNonPr = (TH1D*)h2Unf_jesNonPr->ProjectionX("hZUnf_jesNonPr",2,2);
  TH1D *hZUnf_newNominal = (TH1D*)h2Unf_newNominal->ProjectionX("hZUnf_newNominal",2,2);

  TH1D *hZUnf_newNominal_noScale = (TH1D*)h2Unf_newNominal->ProjectionX("hZUnf_newNominal_noScale",2,2);
  TH1D *hZUnf_SI3_wJESadd = (TH1D*)hZUnf_newNominal_noScale->Clone();
  hZUnf_SI3_wJESadd->Reset();
  
  
  TH1D *hZUnf_jesPr_nominal = (TH1D*)hZUnf_jesPr->Clone("hZUnf_jesPr_nominal");
  hZUnf_jesPr_nominal->Reset();

  TH1D *hZUnf_jesNonPr_nominal = (TH1D*)hZUnf_jesNonPr->Clone("hZUnf_jesNonPr_nominal");
  hZUnf_jesNonPr_nominal->Reset();
    
  int binMin = 2;
  int binMax = 5;

  if(doMid) binMin = 4;
  if(doMid) binMax = 7;

  /*
  hZUnf_jesPr->Scale(1/hZUnf_jesPr->Integral(binMin,binMax));
  hZUnf_jesNonPr->Scale(1/hZUnf_jesNonPr->Integral(binMin,binMax));
  hZUnf_newNominal->Scale(1/hZUnf_newNominal->Integral(binMin,binMax));
  */
  
  int zBins = 5;
  if(doMid) zBins = 7;

  float int_jesPr = hZUnf_jesPr->Integral(binMin,binMax);
  float int_newNominal = hZUnf_newNominal->Integral(binMin,binMax);
  float int_jesNonPr = hZUnf_jesNonPr->Integral(binMin,binMax);

  cout << "int_jesPr = " << int_jesPr << " , int_newNominal = " << int_newNominal << " , int_jesNonPr = " << int_jesNonPr << endl;
  cout << " err 1 = " <<  (int_newNominal-int_jesPr)/int_newNominal  << " err 2 = " << (int_newNominal-int_jesNonPr)/int_newNominal << endl; 

  /*
  for(int i = 0; i < zBins; i++){

    // prompt JES
    
    float binCont_jesPr = hZUnf_jesPr->GetBinContent(i+1);
    float binCont_newNominal = hZUnf_newNominal->GetBinContent(i+1);
    float binCont_jesNonPr = hZUnf_jesNonPr->GetBinContent(i+1);

    // ratio of prompt JES and new nominal
    float ratio_jesPr = 0;
    if(binCont_newNominal>0) ratio_jesPr = (binCont_newNominal-binCont_jesPr)/binCont_newNominal;

    float ratio_jesNonPr = 0;
    if(binCont_newNominal>0) ratio_jesNonPr = (binCont_newNominal-binCont_jesNonPr)/binCont_newNominal;

    
    cout << "bin = " << i << endl;
    cout << "New nominal = " << binCont_newNominal << " , with jes from prompt mc = " << binCont_jesPr  << " , with jes matched to nonprompt mc = " << binCont_jesNonPr << endl;
    
    hZUnf_jesPr_nominal->SetBinContent(i+1,ratio_jesPr);
    hZUnf_jesNonPr_nominal->SetBinContent(i+1,ratio_jesNonPr);

    float binCont = binCont_newNominal;
    float binErrJES = 0;
    if(abs(ratio_jesPr)>abs(ratio_jesNonPr)) {
      binErrJES = abs(ratio_jesPr);
      cout << "rel error in % = " << binErrJES*100 << endl;
    }
    else {
      binErrJES = abs(ratio_jesNonPr);
      cout << "rel error in % = " << binErrJES*100 << endl;
    }

    // put no scaled bin content
    
    float noScaleBinCont = hZUnf_newNominal_noScale->GetBinContent(i+1);
    float binErrJESVal = binErrJES*noScaleBinCont;
    
    hZUnf_SI3_wJESadd->SetBinContent(i+1,noScaleBinCont);
    hZUnf_SI3_wJESadd->SetBinError(i+1,binErrJESVal);
    
  }
  
  TLegend *legend = new TLegend(0.5,0.75,.7,0.93,"","brNDC");
  legend->SetHeader("");

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.035);
    
  TCanvas * can1 = new TCanvas("can1","can1",900,600);
  can1->Divide(2,1);

  can1->cd(1);
  
  hZUnf_jesPr->SetMarkerColor(kRed);
  hZUnf_jesPr->SetMarkerSize(1.1);
  hZUnf_jesPr->SetMarkerStyle(20);
  hZUnf_jesPr->GetXaxis()->SetTitle("z");
  hZUnf_jesPr->GetYaxis()->SetTitle("1/N dN/dz");
  hZUnf_jesPr->GetYaxis()->SetTitleOffset(2.3);
  hZUnf_jesPr->SetMaximum(hZUnf_jesPr->GetMaximum()*1.5);
  
  hZUnf_newNominal->SetMarkerColor(kGreen+2);
  hZUnf_newNominal->SetMarkerSize(1.1);
  hZUnf_newNominal->SetMarkerStyle(21);

  hZUnf_jesNonPr->SetMarkerColor(kBlue);
  hZUnf_jesNonPr->SetMarkerSize(1.5);
  hZUnf_jesNonPr->SetMarkerStyle(33);
  
  legend->AddEntry(hZUnf_jesPr, "prompt JES","ep");
  legend->AddEntry(hZUnf_newNominal, "new nominal","ep");
  legend->AddEntry(hZUnf_jesNonPr, "nonprompt JES","ep");
  
  hZUnf_jesPr->Draw("EP");
  hZUnf_newNominal->Draw("EPsame");
  hZUnf_jesNonPr->Draw("EPsame");
  
  legend->Draw("same");

  float xCoord2;
  if(doMid) xCoord2 = 0.44;
  else xCoord2 = 0.2;

  TLine *line1 = new TLine(xCoord2,hZUnf_jesPr->GetMinimum(),xCoord2,hZUnf_jesPr->GetMaximum());
  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->Draw("same");
  
  can1->cd(2);

  hZUnf_jesPr_nominal->GetYaxis()->SetTitleOffset(2.);
  hZUnf_jesPr_nominal->GetXaxis()->SetTitle("z");
  hZUnf_jesPr_nominal->GetYaxis()->SetTitle("(new nominal - changed)/new nominal");
  hZUnf_jesPr_nominal->GetYaxis()->SetRangeUser(-0.15,.15);
  hZUnf_jesPr_nominal->SetLineColor(kRed);
  hZUnf_jesPr_nominal->SetLineWidth(2);
  hZUnf_jesPr_nominal->Draw("HIST");

  hZUnf_jesNonPr_nominal->SetLineColor(kBlue);
  hZUnf_jesNonPr_nominal->SetLineWidth(2);
  hZUnf_jesNonPr_nominal->SetLineStyle(9);
  hZUnf_jesNonPr_nominal->Draw("HISTsame");
    
  TLine *line2 = new TLine(xCoord2,hZUnf_jesPr_nominal->GetMinimum(),xCoord2,hZUnf_jesPr_nominal->GetMaximum());
  line2->SetLineColor(kRed);
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->Draw("same");
  */  

  /*
  if(doPrompt && doMid){
    can1->SaveAs("../plots/unf_data_prompt_mid_newNominal_noScale.pdf");
  }

  if(doPrompt && !doMid){
    can1->SaveAs("../plots/unf_data_prompt_fwd_newNominal_noScale.pdf");
  }
  */

  /*
  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  hZUnf_newNominal_noScale->Write("zUnfNominal");
  hZUnf_SI3_wJESadd->Write("zUnf_wAddJES");
  outfile->Close();
  */
  
}

void PlotRatios_DataUnfolded_JEScheck_noScale(){

  plot(true,true);
  plot(true,false);


}
