#include "inputParams.h"

void plot(bool doPrompt = true, bool doPbPb = false){
  if (!setSystTag()) return;
  gSystem->mkdir("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/");

  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string filenameMeas = "";
  string outputfile = "";

  int iterFinal = nSIter;
  if (!doPbPb) iterFinal = nSIter_pp;

  filename1 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename2 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",2,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename3 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",iterFinal-1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filename4 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_Diag%s.root",iterFinal,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  filenameMeas = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%s.root",iterFinal,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, systTag.c_str());
  outputfile = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, iterFinal, systTag.c_str());
  
  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  TFile *fileMeas = new TFile(filenameMeas.c_str());
  
  TH2D *h2ZMeas;
  h2ZMeas=(TH2D*)fileMeas->Get("fh2MeasData;1");
  TH1D *hZMeas = h2ZMeas->ProjectionX("hZMeas",3,3);
  hZMeas->Draw();

  cout <<"getting the hists"<<endl;  
  TH2D *h2UnfResp1 = (TH2D*)file1->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp2 = (TH2D*)file2->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp3 = (TH2D*)file3->Get(Form("hReco_Iter%d;1",1));
  TH2D *h2UnfResp4 = (TH2D*)file4->Get(Form("hReco_Iter%d;1",1));

  cout <<"hists found"<<endl;  

  Int_t min_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midLowerPt+0.00001);
  Int_t max_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midUpperPt-0.00001);

  cout <<"min_PriorFolded = "<<min_PriorFolded<<"; max_PriorFolded = "<<max_PriorFolded<<endl;
  
  TH1D *hZUnf_SI1 = (TH1D*)h2UnfResp1->ProjectionX("hZUnf_SI1",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI2 = (TH1D*)h2UnfResp2->ProjectionX("hZUnf_SI2",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI3 = (TH1D*)h2UnfResp3->ProjectionX("hZUnf_SI3",min_PriorFolded,max_PriorFolded);
  TH1D *hZUnf_SI4 = (TH1D*)h2UnfResp4->ProjectionX("hZUnf_SI4",min_PriorFolded,max_PriorFolded);

  cout<<"saving"<<endl;
  TFile *file_out = new TFile(outputfile.c_str(),"RECREATE");
  hZUnf_SI1->Write(Form("zUnfSI%d",1));
  hZUnf_SI2->Write(Form("zUnfSI%d",2));
  hZUnf_SI3->Write(Form("zUnfSI%d",iterFinal-1));
  hZUnf_SI4->Write(Form("zUnfSI%d",iterFinal));
  hZMeas->Write("zMeas");
  file_out->Close();
}

void finalResults_statErrs(){
  plot(true,true);
  if (centShift ==0)
  plot(true,false);
}
