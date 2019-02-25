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

  string filename = "";
  string outputfile = "";

  if(doPrompt && doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal_Diag.root";
    outputfile = "../../unfOutput/JERJESErrs/JERJESerr_Prompt_Mid_newNominal.root";
  }
  
  if(doPrompt && !doMid){
    filename ="../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal_Diag.root";
    outputfile = "../../unfOutput/JERJESErrs/JERJESerr_Prompt_Fwd_newNominal.root";
  }

  if(!doPrompt && doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal_Diag.root";
    outputfile = "../../unfOutput/JERJESErrs/JERJESerr_NonPrompt_Mid.root";
  }

  if(!doPrompt && !doMid){
    filename ="../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal_Diag.root";
    outputfile = "../../unfOutput/JERJESErrs/JERJESerr_NonPrompt_Fwd.root";
  }
  
  
  /*
JES 
    mid-rapidity 
    0.0<z<0.2: 0.000 0.000 
    0.2<z<0.4: 0.066 -0.063 
    0.4<z<0.6: 0.002 -0.003 
    0.6<z<0.8: -0.030 0.029 
    0.8<z<1.0: -0.035 0.036 
    forward rapidity 
    0.0<z<0.2: 0.000 -0.384 
    0.2<z<0.4: 0.061 -0.065 
    0.4<z<0.6: 0.005 -0.002 
    0.6<z<0.8: -0.028 0.026 
    0.8<z<1.0: -0.044 0.042 

   */
  
  
  TFile *file = new TFile(filename.c_str());
  
  TH2D *h2ZUnf_SI3;
  h2ZUnf_SI3=(TH2D*)file->Get("hReco_Iter1;1");

  TH1D *hZUnf_SI3;
  hZUnf_SI3=(TH1D*)h2ZUnf_SI3->ProjectionX("hZUnf_SI3",2,2);

  hZUnf_SI3->Draw();

  TH1D *hZUnf_SI3_wJES = (TH1D*)hZUnf_SI3->Clone();
  hZUnf_SI3_wJES->Reset();
  
  TH1D *hZUnf_SI3_wJER = (TH1D*)hZUnf_SI3->Clone();
  hZUnf_SI3_wJER->Reset();
  
  float binCont = 0;
  float JESval = 0;
  float JERval = 0;
  float binErrJES = 0;
  float binErrJER = 0;

  int nBins = 5;
  if(doMid) nBins = 7;
  
  for(int i = 0; i < nBins; i++){

    binCont = hZUnf_SI3->GetBinContent(i+1);
      
    if(doMid){
      if(i+1<3) {
	JESval=0.;
	JERval=0.;
      }
      if(i+1==3){
	JESval=0.066;
	if(doPrompt) JERval=0.097;
	else JERval=0.122;
      }
      if(i+1==4){
	JESval=0.004;
	if(doPrompt) JERval=0.010;
	else JERval=0.010;
      }
      if(i+1==5){
	JESval=0.032;
	if(doPrompt) JERval=0.008;
	else JERval=0.024;
      }
      if(i+1==6){
	JESval=0.053;
	if(doPrompt) JERval=0.012;
	else JERval=0.010;
      }
      if(i+1==7){
	JESval=0.056;
	if(doPrompt) JERval=0.022;
	else JERval=0.045;
      }
      
    }
    
    else{
      if(i+1==1) {
	JESval=0.;
	JERval=0.;
      }
      if(i+1==2){
	JESval=0.079;
	if(doPrompt) JERval=0.019;
	else JERval=0.058;
      }
      if(i+1==3){
	JESval=0.010;
	if(doPrompt) JERval=0.008;
	else JERval=0.010;
      }
      if(i+1==4){
	JESval=0.037;
	if(doPrompt) JERval=0.036;
	else JERval=0.045;
      }
      if(i+1==5){
	JESval=0.068;
	if(doPrompt) JERval=0.006;
	else JERval=0.085;
      }
      
    }

    binErrJES = JESval*binCont;
    binErrJER = JERval*binCont;
    
    hZUnf_SI3_wJES->SetBinContent(i+1,binCont);
    hZUnf_SI3_wJES->SetBinError(i+1,binErrJES);

    hZUnf_SI3_wJER->SetBinContent(i+1,binCont);
    hZUnf_SI3_wJER->SetBinError(i+1,binErrJER);
    
  }
  

  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  
  hZUnf_SI3->Write("zUnfNominal");
  hZUnf_SI3_wJES->Write("zUnf_wJESErr");
  hZUnf_SI3_wJER->Write("zUnf_wJERErr");
  
  outfile->Close();
  
  
}

void systFromJERJES(){

  //compute(true,true);
  //compute(true,false);
  //compute(false,true);
  compute(false,false);
    
}
