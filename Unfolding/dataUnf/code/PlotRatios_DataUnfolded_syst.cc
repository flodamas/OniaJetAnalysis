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

  string filename = "";
  string filename3 = "";

  string outputfile = "";

  if(doPrompt && doMid){
    cout << "prompt mid" << endl;
    filename3 = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_syst_nominal_Diag.root";
    filename = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_syst_nominal.root";
    outputfile = "unfResult_prompt_mid_wQSyst.root";
  }
  
  if(doPrompt && !doMid){
    cout << "prompt fwd" << endl;
    filename3 ="../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_syst_nominal_Diag.root";
    filename = "../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_syst_nominal.root";
    outputfile = "unfResult_prompt_fwd_wQSyst.root";
  }

  if(!doPrompt && doMid){
    cout << "nonprompt mid" << endl;
    filename3 = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_syst_nominal_Diag.root";
    filename = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_syst_nominal.root";
    outputfile = "unfResult_nonprompt_mid_wQSyst.root";
  }

  if(!doPrompt && !doMid){
    cout << "nonprompt fwd" << endl;
    filename3 ="../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_syst_nominal_Diag.root";
    filename = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_syst_nominal.root";
    outputfile = "unfResult_nonprompt_fwd_wQSyst.root";
  }
  

  TFile *file3 = new TFile(filename3.c_str());
  TFile *file = new TFile(filename.c_str());
  
  
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)file->Get("fh2MeasData;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",2,2);

  //unfolded result

  TH2D *h2UnfResp3 = (TH2D*)file->Get("hReco_Iter3;1");
  TH1D *hZUnf_SI3 = (TH1D*)h2UnfResp3->ProjectionX("hZUnf_SI3",6,10);
  if(!doMid) hZUnf_SI3->Rebin(10);
  else hZUnf_SI3->Rebin(7);
  
  // after errors got back
  
  TH2D *h2UnfResp3B = (TH2D*)file3->Get("hReco_Iter1;1");
  TH2D *h2UnfResp3Inv = (TH2D*)file3->Get("hReco_Invert;1");

  TH1D *hZUnf_SI3B = (TH1D*)h2UnfResp3B->ProjectionX("hZUnf_SI3B",2,2);
  TH1D *hZUnf_SI3Inv = (TH1D*)h2UnfResp3Inv->ProjectionX("hZUnf_SI3Inv",2,2);

  int zBins = 5;
  if(doMid) zBins = 7;
  
  for(int i = 0; i < zBins; i++){

    float binContBayes = hZUnf_SI3B->GetBinContent(i+1);
    float binErrBayes = hZUnf_SI3B->GetBinError(i+1);
    float relErrBayes = 0;
    if(binContBayes>0) relErrBayes = binErrBayes/binContBayes;
    
    float binContInv = hZUnf_SI3Inv->GetBinContent(i+1);
    float binErrInv = hZUnf_SI3Inv->GetBinError(i+1);
    float relErrInv = 0;
    if(binContInv>0) relErrInv = binErrInv/binContInv;
        
    float binContUnf = hZUnf_SI3->GetBinContent(i+1);
    float binErrUnf = hZUnf_SI3->GetBinError(i+1);
    float relErrUnf = 0;
    if(binContUnf>0) relErrUnf = binErrUnf/binContUnf;
        
    float binContMeas = hZMeas->GetBinContent(i+1);
    float binErrMeas = hZMeas->GetBinError(i+1);
    float relErrMeas = 0;
    if(binContMeas>0) relErrMeas = binErrMeas/binContMeas;
    
    
    cout << "bin = " << i << endl;
    cout << "Bayes :: content = " << binContBayes << " +- " << binErrBayes << " rel unc (%) = " << relErrBayes*100 << endl;
    cout << "Inverted :: content = " << binContInv << " +- " << binErrInv << " rel unc (%) = " << relErrInv*100 <<  endl;
    cout << "Unfolded :: content = " << binContUnf << " +- " << binErrUnf << " rel unc (%) = " << relErrUnf*100 <<  endl;
    cout << "Meas :: content = " << binContMeas << " +- " << binErrMeas << " rel unc (%) = " << relErrMeas*100 << endl;
  }


  /*
  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  hZUnf_SI3B->Write("zUnf");
  hZMeas->Write("zMeas");
  outfile->Close();
  */
  
}

void PlotRatios_DataUnfolded_syst(){

  //plot(true,true);
  //plot(true,false);
  //plot(false,true);
  plot(false,false);

}
