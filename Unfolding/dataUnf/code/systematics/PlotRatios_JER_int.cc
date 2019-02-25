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
    cout << "prompt mid" << endl;
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_down.root";
  }
  
  if(doPrompt && !doMid){
    cout << "prompt fwd" << endl;
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_down.root";
  }
  
  if(!doPrompt && doMid){
    cout << "nonprompt mid" << endl;
    filename_nominal = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    filename_up = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_up.root";
    filename_down = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_down.root";
  }
  
  if(!doPrompt && !doMid){
    cout << "nonprompt fwd" << endl;
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

  float int_nominal = hZUnf_SI3_nominal->Integral(binMin,binMax);
  float int_up = hZUnf_SI3_up->Integral(binMin,binMax);
  float int_down = hZUnf_SI3_down->Integral(binMin,binMax);
  
  cout << "int_nominal = " << int_nominal << " , int_up = " << int_up << " , int_down = " << int_down << endl;
  cout << "up = " << (int_nominal-int_up)/int_nominal << " , down = " << (int_nominal-int_down)/int_nominal << endl;
}

void PlotRatios_JER_int(){

  plot(true,true);
  plot(true,false);
  
  plot(false,true);
  plot(false,false);
  
}
