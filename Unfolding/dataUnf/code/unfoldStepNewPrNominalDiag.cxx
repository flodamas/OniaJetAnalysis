#include "inputParams.h"

void unfold(bool doPrompt = true, bool doPbPb = true, Int_t iterMax = 3) {
  if (!setSystTag()) return;

  Int_t iterMin =1;
  Int_t iterDef = iterMax - 2;
  
  string testInputName = "";
  string trainInputName = "";
  string outputName = "";
  
  //string SF_name = "";
  //SF_name = "_nominal";

  testInputName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%s.root",iterMax,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  trainInputName = "/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/code/diag4DMatrixResponseInv.root";
  outputName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",iterMax,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  
  
  TFile *f_measured = new TFile(testInputName.c_str());
  f_measured->ls();
  //unfolded as measured
  TH2D *hMeasured = (TH2D*)f_measured->Get(Form("hReco_Iter%d;1",nIter));
  TMatrixD *covmat = (TMatrixD*)f_measured->Get(Form("covmat%d;1",nIter));
  
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

void unfoldStepNewPrNominalDiag(int step){
  unfold(true,true,step);
  if (step<=nSIter_pp && (centShift==0))
    unfold(true,false,step);
}
