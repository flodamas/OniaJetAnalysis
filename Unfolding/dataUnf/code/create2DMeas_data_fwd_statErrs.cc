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

void create2DMeas_data_fwd_systErrs(bool doPrompt = true){

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  
  filename1 = "../data_results/unfoldingInput_lowJtPt_rap1624_statErr.root";
  filename2 = "../data_results/unfoldingInput_midJtPt_rap1624_statErr.root";
  filename3 = "../data_results/unfoldingInput_highJtPt_rap1624_statErr.root";

  TFile *file1 = new TFile(filename1.c_str(),"UPDATE");
  TFile *file2 = new TFile(filename2.c_str(),"UPDATE");
  TFile *file3 = new TFile(filename3.c_str(),"UPDATE");

  TH1F *h_lowBin;
  TH1F *h_centBin;
  TH1F *h_highBin;

  if(doPrompt){
    h_lowBin = (TH1F*)file1->Get("prHist_lowJtPt_rap1624_statErr;1");
    h_centBin = (TH1F*)file2->Get("prHist_midJtPt_rap1624_statErr;1");
    h_highBin = (TH1F*)file3->Get("prHist_highJtPt_rap1624_statErr;1");
  }
  else{
    h_lowBin = (TH1F*)file1->Get("nprHist_lowJtPt_rap1624_statErr;1");
    h_centBin = (TH1F*)file2->Get("nprHist_midJtPt_rap1624_statErr;1");
    h_highBin = (TH1F*)file3->Get("nprHist_highJtPt_rap1624_statErr;1");
  }

  TH2D *h_Meas = new TH2D("h_Meas","h_Meas",5,.0,1.0,3,15.,45.);
  h_Meas->Sumw2();
  
  int nbins = h_lowBin->GetNbinsX();

  float lowBinVals[nbins];
  float centBinVals[nbins];
  float highBinVals[nbins];

  float lowBinErrs[nbins];
  float centBinErrs[nbins];
  float highBinErrs[nbins];
    
  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_lowBin : %i ", i) << h_lowBin->GetBinContent(i) << endl;
    lowBinVals[i-1] = h_lowBin->GetBinContent(i);
    lowBinErrs[i-1] = h_lowBin->GetBinError(i);
  }

  cout << "***" << endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_centBin : %i ", i) << h_centBin->GetBinContent(i) << endl;
    centBinVals[i-1] = h_centBin->GetBinContent(i);
    centBinErrs[i-1] = h_centBin->GetBinError(i);
  }

  cout << "***"<< endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_highBin : %i ", i) << h_highBin->GetBinContent(i) << endl;
    highBinVals[i-1] = h_highBin->GetBinContent(i);
    highBinErrs[i-1] = h_highBin->GetBinError(i);
  }

  cout << "***"<< endl;

  float content = 0.0;
  float errVal = 0.0;
  
  for(int icount=0; icount < 3; icount++){
    for(int jcount=0; jcount < 5; jcount++){

      if(icount == 0) {
	content = lowBinVals[jcount];
	errVal = lowBinErrs[jcount];
      }
      if(icount == 1) {
	content = centBinVals[jcount];
	errVal = centBinErrs[jcount];
      }
      if(icount == 2) {
	content = highBinVals[jcount];
	errVal = highBinErrs[jcount];
      }
      
      h_Meas->SetBinContent(jcount+1,icount+1,content);
      h_Meas->SetBinError(jcount+1,icount+1,errVal);
    }
  }

  TCanvas * can0 = new TCanvas ("can0","can0",900,450);
  can0->Divide(3,1);

  can0->cd(1);
  h_lowBin->SetTitle("15 < jet pt < 25");
  h_lowBin->Draw();

  can0->cd(2);
  h_centBin->SetTitle("25 < jet pt < 35");
  h_centBin->Draw();

  can0->cd(3);
  h_highBin->SetTitle("35 < jet pt < 45");
  h_highBin->Draw();

  if(doPrompt) can0->SaveAs("../plots/data_z_distr_promptFwd_statErr.pdf");
  else can0->SaveAs("../plots/data_z_distr_nonpromptFwd_statErr.pdf");

  TCanvas * can1 = new TCanvas ("can1","can1",600,600);

  can1->SetRightMargin(0.2);

  if(doPrompt) h_Meas->SetTitle("prompt fow-rapidity");
  else h_Meas->SetTitle("nonprompt fow-rapidity");
  
  h_Meas->SetStats(0);
  h_Meas->GetXaxis()->SetTitle("z");
  h_Meas->GetYaxis()->SetTitle("jet p_{T}");
  h_Meas->GetYaxis()->SetTitleOffset(1.2);
  h_Meas->Draw("TEXTcolz");
  if(doPrompt) can1->SaveAs("../plots/data_meas_promptFwd_statErr.pdf");
  else can1->SaveAs("../plots/data_meas_nonpromptFwd_statErr.pdf");

  string outputfile = "";
  if(doPrompt) outputfile = "../data_results/meas_data_prompt_fwd_statErrs.root";
  else outputfile = "../data_results/meas_data_nonprompt_fwd_statErrs.root";
  
  TFile *file_data_meas = new TFile(outputfile.c_str(),"RECREATE");

  h_Meas->Write();
  h_lowBin->Write("zMeasLowBin");
  h_centBin->Write("zMeasCentBin");
  h_highBin->Write("zMeasHighBin");
  
}
