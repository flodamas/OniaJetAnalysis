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
  //string outputfile = "";

  if(doPrompt && doMid){
    cout << "**********" << endl;
    cout << "prompt mid" << endl;
    filename1 = "../../unfOutput/step1/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename2 = "../../unfOutput/step2/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename3 = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    filename4 = "../../unfOutput/step4/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    //outputfile = "../../unfOutput/numberSIErrs/SIerr_Prompt_Mid_newNominal.root";
  }
  
  if(doPrompt && !doMid){
    cout << "prompt fwd"<< endl;
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    //outputfile = "../../unfOutput/numberSIErrs/SIerr_Prompt_Fwd_newNominal.root";
  }

  if(!doPrompt && doMid){
    cout << "nonprompt mid"<< endl;
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    //outputfile = "../../unfOutput/numberSIErrs/SIerr_NonPrompt_Mid.root";
  }
  
  if(!doPrompt && !doMid){
    cout << "nonprompt fwd"<< endl;
    filename1 ="../../unfOutput/step1/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename2 ="../../unfOutput/step2/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename3 ="../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    filename4 ="../../unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    //outputfile = "../../unfOutput/numberSIErrs/SIerr_NonPrompt_Fwd.root";
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

  int binMin = 2;
  int binMax = 5;

  if(doMid) binMin = 4;
  if(doMid) binMax = 7;
    
  // integral :

  float nominal_int = hZUnf_SI3->Integral(binMin,binMax);
  float si2_int = hZUnf_SI2->Integral(binMin,binMax);
  float si4_int = hZUnf_SI4->Integral(binMin,binMax);

  cout << "nominal = " << nominal_int << " , SI2 = " << si2_int << " , SI_4 = " << si4_int << endl;  
  cout << "rel err 1 = " << (nominal_int-si2_int)/nominal_int << " , rel err 2 = "<< (nominal_int-si4_int)/nominal_int << endl;
}

void systFromNumSI_int(){

  compute(true,true);
  compute(true,false);
  compute(false,true);
  compute(false,false);
    
  
}
