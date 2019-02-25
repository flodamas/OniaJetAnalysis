#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

using namespace std;

void unfold(bool doPrompt = true, bool doMid = true, Int_t iterMax = 3) {

  //gSystem->Load("libRooUnfold");
  
Int_t iterMin =1;
Int_t iterDef = iterMax - 2;
  
//#ifdef __CINT__
//  gSystem->Load("libRooUnfold");
//#endif

  string testInputName = "";
  string trainInputName = "";
  string outputName = "";

  string SF_name = "";
  SF_name = "_nominal";
  
  if(doPrompt && doMid){
    testInputName = Form("dataUnfNewMidBins/unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nonPrJES%s.root",SF_name.c_str());
    trainInputName = "dataUnfNewMidBins/code/diag4DMatrixResponseInvMid.root";
    outputName = Form("dataUnfNewMidBins/unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_nonPrJES%s_Diag.root",SF_name.c_str());
  }

  if(doPrompt && !doMid){
    testInputName = Form("dataUnfNewMidBins/unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nonPrJES%s.root",SF_name.c_str());
    trainInputName = "dataUnfNewMidBins/code/diag4DMatrixResponseInv.root";
    outputName = Form("dataUnfNewMidBins/unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_nonPrJES%s_Diag.root",SF_name.c_str());
  }

  if(!doPrompt && doMid){
    testInputName = Form("dataUnfNewMidBins/unfOutput/step4/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins%s.root",SF_name.c_str());
    trainInputName = "dataUnfNewMidBins/code/diag4DMatrixResponseInvMid.root";
    outputName = Form("dataUnfNewMidBins/unfOutput/step4/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins%s_Diag.root",SF_name.c_str());
  }

  if(!doPrompt && !doMid){
    testInputName = Form("dataUnfNewMidBins/unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins%s.root",SF_name.c_str());
    trainInputName = "dataUnfNewMidBins/code/diag4DMatrixResponseInv.root";
    outputName = Form("dataUnfNewMidBins/unfOutput/step4/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins%s_Diag.root",SF_name.c_str());
  }
  

  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  //unfolded as measured
  TH2D *hMeasured = (TH2D*)f_measured->Get("hReco_Iter3;1");
  TMatrixD *covmat = (TMatrixD*)f_measured->Get("covmat3;1");
  
  TFile *f = new TFile(trainInputName.c_str());
  RooUnfoldResponse *resp = (RooUnfoldResponse*)f->Get("resp;1");
  TH2F *hSmear = (TH2F*)f->Get("fh2RespDimM;1");
  TH2F *hTrue = (TH2F*)f->Get("fh2RespDimT;1");    
    
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//kNoError;//

  
  TH2D *hReco;
  RooUnfoldBayes unfold;
  
  TH2D *hRecoInvert;
  RooUnfoldInvert unfoldInvert;
        

  unfold = RooUnfoldBayes(resp, hMeasured, 1);
  unfold.SetMeasuredCov(*covmat);
  
  hReco = (TH2D*)unfold.Hreco(errorTreatment);
  hReco->SetName(Form("hReco_Iter%d",1));
    
  unfoldInvert = RooUnfoldInvert(resp, hMeasured);
  unfoldInvert.SetMeasuredCov(*covmat);

  hRecoInvert = (TH2D*)unfoldInvert.Hreco(errorTreatment);
  hRecoInvert->SetName("hReco_Invert");
  
  
  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  resp->Write();
  
  hReco->Write();
  unfold.Write(Form("unfold%d",1));

  hRecoInvert->Write();
  unfoldInvert.Write("unfoldInvert");
  
  
  fout->Write();
  fout->Close();
  
}

void unfoldStepNonPrJESDiag(){

  unfold(true,true);
  unfold(true,false);
  
}
