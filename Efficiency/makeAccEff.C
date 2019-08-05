#define makeAccEff_cxx
#define _USE_MATH_DEFINES
//#include "makeAccEff.h"
//#include "compAccEff.C"
#include "systAccEff.C"

//Double_t ptbins []= {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 100.0};
//Double_t ybins []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

Double_t ptbinsAna []= {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 14.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};

void oniaTree::AccEffCalc()
{ 
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;
  string centTag [] = {"cent010","cent1020","cent2030","cent3040","cent4060","cent6080","cent80100","cent100140","cent140200"};

  cout<<"[INFO] Importing the numerators and denominators of the corrections."<<endl;
  //TFile*prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PbPb.root");
  TFile*prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
  if (!prAccFile_pbpb) {
    cout<<"[ERROR] pbpb prompt acceptance file not found!"<<endl;
    if (isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc();
      //prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PbPb.root");
      prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  
  //TFile*nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PbPb.root");
  TFile*nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
  if (!nprAccFile_pbpb) {
    cout<<"[ERROR] pbpb nonprompt acceptance file not found!"<<endl;
    if (isAcc && isPbPb && !isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc();
      //nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PbPb.root");
      nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  

  TFile*prEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/prEffHists_PbPb.root");
  if (!prEffFile_pbpb) {
    cout<<"[ERROR] pbpb prompt efficiency file not found!"<<endl;
    if (!isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/prEffHists_PbPb.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  
  TFile*nprEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/nprEffHists_PbPb.root");
  if (!nprEffFile_pbpb) {
    cout<<"[ERROR] pbpb nonprompt Eff file not found!"<<endl;
    if (!isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      nprEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/nprEffHists_PbPb.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  
  ////////////////////////////////////pp/////////////////////////////////
  TFile*prAccFile_pp = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
  if (!prAccFile_pp) {
    cout<<"[ERROR] pp prompt acceptance file not found!"<<endl;
    if (isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
    AccCalc();
    prAccFile_pp = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  /*
  TFile*nprAccFile_pp = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
  if (!nprAccFile_pp) {
    cout<<"[ERROR] pp nonprompt acceptance file not found!"<<endl;
    if (isAcc && !isPbPb && !isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc();
      nprAccFile_pp = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  */
  TFile*prEffFile_pp = TFile::Open("FilesAccxEff/Eff/prEffHists_PP.root");
  if (!prEffFile_pp) {
    cout<<"[ERROR] pp prompt efficiency file not found!"<<endl;
    if (!isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile_pp = TFile::Open("FilesAccxEff/Eff/prEffHists_PP.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }
  /*
  TFile*nprEffFile_pp = TFile::Open("FilesAccxEff/Eff/nprEffHists_PP.root");
  if (!nprEffFile_pp) {
    cout<<"[ERROR] pp nonprompt Eff file not found!"<<endl;
    if (!isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      nprEffFile_pp = TFile::Open("FilesAccxEff/Eff/nprEffHists_PP.root");
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  */

  TH2F *prAccNum_pbpb = (TH2F*) prAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *prAccDen_pbpb = (TH2F*) prAccFile_pbpb->Get("hdeno_2d");
  TH2F *nprAccNum_pbpb = (TH2F*) nprAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *nprAccDen_pbpb = (TH2F*) nprAccFile_pbpb->Get("hdeno_2d");
  TH2F *prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_nominal");
  TH2F *prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get("hdeno_pty");
  TH2F *nprEffNum_pbpb = (TH2F*) nprEffFile_pbpb->Get("hnum_2d_nominal");
  TH2F *nprEffDen_pbpb = (TH2F*) nprEffFile_pbpb->Get("hdeno_2d");

  TH2F *prAccNum_pp = (TH2F*) prAccFile_pp->Get("hnum_2d_nominal");
  TH2F *prAccDen_pp = (TH2F*) prAccFile_pp->Get("hdeno_2d");
  //TH2F *nprAccNum_pp = (TH2F*) nprAccFile_pp->Get("hnum_2d_nominal");
  //TH2F *nprAccDen_pp = (TH2F*) nprAccFile_pp->Get("hdeno_2d");
  TH2F *prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_nominal");
  TH2F *prEffDen_pp = (TH2F*) prEffFile_pp->Get("hdeno_pty");
  //TH2F *nprEffNum_pp = (TH2F*) nprEffFile_pp->Get("hnum_2d_nominal");
  //TH2F *nprEffDen_pp = (TH2F*) nprEffFile_pp->Get("hdeno_2d");

  cout <<"PbPb pr num"<<endl;
  prEffNum_pbpb->Multiply(prAccNum_pbpb);
  cout <<"PbPb pr den"<<endl;
  prEffDen_pbpb->Multiply(prAccDen_pbpb);
  cout <<"PbPb npr num"<<endl;
  nprEffNum_pbpb->Multiply(nprAccNum_pbpb);
  cout <<"PbPb npr den"<<endl;
  nprEffDen_pbpb->Multiply(nprAccDen_pbpb);

  cout <<"pp pr num"<<endl;
  prEffNum_pp->Multiply(prAccNum_pp);
  cout <<"pp pr den"<<endl;
  prEffDen_pp->Multiply(prAccDen_pp);
  //cout <<"pp npr num"<<endl;
  //nprEffNum_pp->Multiply(nprAccNum_pp);
  //cout <<"pp npr den"<<endl;
  //nprEffDen_pp->Multiply(nprAccDen_pp);

  TEfficiency* prCorr_pbpb = new TEfficiency("prCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pbpb->SetPassedHistogram(*prEffNum_pbpb,"f");
  prCorr_pbpb->SetTotalHistogram(*prEffDen_pbpb,"f");
  prCorr_pbpb->SetName("hcorr_Jpsi_PbPb_pr");
  
  TEfficiency* nprCorr_pbpb = new TEfficiency("nprCorr_pbpb", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  nprCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
  nprCorr_pbpb->SetPassedHistogram(*nprEffNum_pbpb,"f");
  nprCorr_pbpb->SetTotalHistogram(*nprEffDen_pbpb,"f");
  nprCorr_pbpb->SetName("hcorr_Jpsi_PbPb_npr");
  
  TEfficiency* prCorr_pp = new TEfficiency("prCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pp->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp->SetPassedHistogram(*prEffNum_pp,"f");
  prCorr_pp->SetTotalHistogram(*prEffDen_pp,"f");
  prCorr_pp->SetName("hcorr_Jpsi_PP_pr");
  /*
  TEfficiency* nprCorr_pp = new TEfficiency("nprCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  nprCorr_pp->SetStatisticOption(TEfficiency::kBBayesian);
  nprCorr_pp->SetPassedHistogram(*nprEffNum_pp,"f");
  nprCorr_pp->SetTotalHistogram(*nprEffDen_pp,"f");
  nprCorr_pp->SetName("hcorr_Jpsi_PP_npr");
  */

  TFile* fsave = new TFile("../Fitter/Input/correction_AccEff_centMaps.root","RECREATE");
  prCorr_pbpb->Write("hcorr_Jpsi_PbPb_pr");
  //nprCorr_pbpb->Write("hcorr_Jpsi_PbPb_npr");
  prCorr_pp->Write("hcorr_Jpsi_PP_pr");
  //nprCorr_pp->Write("hcorr_Jpsi_PP_npr");
  for (int i=0; i<9; i++)
    {
      prAccNum_pbpb = (TH2F*) prAccFile_pbpb->Get("hnum_2d_nominal");
      prAccDen_pbpb = (TH2F*) prAccFile_pbpb->Get("hdeno_2d");
      prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hnum_nominal_%s",centTag[i].c_str()));
      prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hdeno_pty_%s",centTag[i].c_str()));
      prEffNum_pbpb->Multiply(prAccNum_pbpb);
      prEffDen_pbpb->Multiply(prAccDen_pbpb);
      prCorr_pbpb = new TEfficiency(Form("prCorr_pp_%s",centTag[i].c_str()), "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
      prCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
      prCorr_pbpb->SetPassedHistogram(*prEffNum_pbpb,"f");
      prCorr_pbpb->SetTotalHistogram(*prEffDen_pbpb,"f");
      prCorr_pbpb->SetName(Form("hcorr_Jpsi_PbPb_pr_%s",centTag[i].c_str()));
      prCorr_pbpb->Write(Form("hcorr_Jpsi_PbPb_pr_%s",centTag[i].c_str()));

      nprAccNum_pbpb = (TH2F*) nprAccFile_pbpb->Get("hnum_2d_nominal");
      nprAccDen_pbpb = (TH2F*) nprAccFile_pbpb->Get("hdeno_2d");
      nprEffNum_pbpb = (TH2F*) nprEffFile_pbpb->Get(Form("hnum_nominal_%s",centTag[i].c_str()));
      nprEffDen_pbpb = (TH2F*) nprEffFile_pbpb->Get(Form("hdeno_pty_%s",centTag[i].c_str()));
      nprEffNum_pbpb->Multiply(nprAccNum_pbpb);
      nprEffDen_pbpb->Multiply(nprAccDen_pbpb);
      nprCorr_pbpb = new TEfficiency(Form("nprCorr_pp_%s",centTag[i].c_str()), "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
      nprCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
      nprCorr_pbpb->SetPassedHistogram(*nprEffNum_pbpb,"f");
      nprCorr_pbpb->SetTotalHistogram(*nprEffDen_pbpb,"f");
      nprCorr_pbpb->SetName(Form("hcorr_Jpsi_PbPb_npr_%s",centTag[i].c_str()));
      nprCorr_pbpb->Write(Form("hcorr_Jpsi_PbPb_npr_%s",centTag[i].c_str()));
    }
  fsave->Close();
}

void oniaTree::EffCalc () {  

  if (isAcc) {cout<<"[ERROR] you're trying to make Efficiency with Acceptance trees."<<endl; return;}
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;

  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbinsAna, ptbinsAna); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbinsAna, ptbinsAna); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs y; y; N_{reco}", nybins, ybins); hnum_y->Sumw2();
  
  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_noweights->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal->Sumw2();

  //for centrality dependance
  TH2F* hdeno_pty_cent010 = new TH2F ("hdeno_pty_cent010", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent010->Sumw2();
  TH2F* hnum_nominal_cent010 = new TH2F ("hnum_nominal_cent010", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent010->Sumw2();

  TH2F* hdeno_pty_cent1020 = new TH2F ("hdeno_pty_cent1020", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent1020->Sumw2();
  TH2F* hnum_nominal_cent1020 = new TH2F ("hnum_nominal_cent1020", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent1020->Sumw2();

  TH2F* hdeno_pty_cent2030 = new TH2F ("hdeno_pty_cent2030", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent2030->Sumw2();
  TH2F* hnum_nominal_cent2030 = new TH2F ("hnum_nominal_cent2030", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent2030->Sumw2();

  TH2F* hdeno_pty_cent3040 = new TH2F ("hdeno_pty_cent3040", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent3040->Sumw2();
  TH2F* hnum_nominal_cent3040 = new TH2F ("hnum_nominal_cent3040", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent3040->Sumw2();

  TH2F* hdeno_pty_cent4060 = new TH2F ("hdeno_pty_cent4060", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent4060->Sumw2();
  TH2F* hnum_nominal_cent4060 = new TH2F ("hnum_nominal_cent4060", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent4060->Sumw2();

  TH2F* hdeno_pty_cent6080 = new TH2F ("hdeno_pty_cent6080", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent6080->Sumw2();
  TH2F* hnum_nominal_cent6080 = new TH2F ("hnum_nominal_cent6080", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent6080->Sumw2();

  TH2F* hdeno_pty_cent80100 = new TH2F ("hdeno_pty_cent80100", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent80100->Sumw2();
  TH2F* hnum_nominal_cent80100 = new TH2F ("hnum_nominal_cent80100", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent80100->Sumw2();

  TH2F* hdeno_pty_cent100140 = new TH2F ("hdeno_pty_cent100140", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent100140->Sumw2();
  TH2F* hnum_nominal_cent100140 = new TH2F ("hnum_nominal_cent100140", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent100140->Sumw2();

  TH2F* hdeno_pty_cent140200 = new TH2F ("hdeno_pty_cent140200", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent140200->Sumw2();
  TH2F* hnum_nominal_cent140200 = new TH2F ("hnum_nominal_cent140200", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal_cent140200->Sumw2();

  //for pp
  TH2F* hnum_muidtrg_plus1sig_syst = new TH2F ("hnum_muidtrg_plus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muidtrg_plus1sig_syst->Sumw2();
  TH2F* hnum_muidtrg_minus1sig_syst = new TH2F ("hnum_muidtrg_minus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muidtrg_minus1sig_syst->Sumw2();
  TH2F* hnum_muidtrg_plus1sig_stat = new TH2F ("hnum_muidtrg_plus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muidtrg_plus1sig_stat->Sumw2();
  TH2F* hnum_muidtrg_minus1sig_stat = new TH2F ("hnum_muidtrg_minus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muidtrg_minus1sig_stat->Sumw2();

  TH2F* hnum_glb_plus1sig_syst = new TH2F ("hnum_glb_plus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_glb_plus1sig_syst->Sumw2();
  TH2F* hnum_glb_minus1sig_syst = new TH2F ("hnum_glb_minus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_glb_minus1sig_syst->Sumw2();
  TH2F* hnum_glb_plus1sig_stat = new TH2F ("hnum_glb_plus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_glb_plus1sig_stat->Sumw2();
  TH2F* hnum_glb_minus1sig_stat = new TH2F ("hnum_glb_minus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_glb_minus1sig_stat->Sumw2();

  //for PbPb
  TH2F* hnum_muid_plus1sig_syst = new TH2F ("hnum_muid_plus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_plus1sig_syst->Sumw2();
  TH2F* hnum_muid_minus1sig_syst = new TH2F ("hnum_muid_minus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_minus1sig_syst->Sumw2();
  TH2F* hnum_muid_plus1sig_stat = new TH2F ("hnum_muid_plus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_plus1sig_stat->Sumw2();
  TH2F* hnum_muid_minus1sig_stat = new TH2F ("hnum_muid_minus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_minus1sig_stat->Sumw2();

  TH2F* hnum_trg_plus1sig_syst = new TH2F ("hnum_trg_plus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trg_plus1sig_syst->Sumw2();
  TH2F* hnum_trg_minus1sig_syst = new TH2F ("hnum_trg_minus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trg_minus1sig_syst->Sumw2();
  TH2F* hnum_trg_plus1sig_stat = new TH2F ("hnum_trg_plus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trg_plus1sig_stat->Sumw2();
  TH2F* hnum_trg_minus1sig_stat = new TH2F ("hnum_trg_minus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trg_minus1sig_stat->Sumw2();

  TH2F* hnum_trk_plus1sig_syst = new TH2F ("hnum_trk_plus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_plus1sig_syst->Sumw2();
  TH2F* hnum_trk_minus1sig_syst = new TH2F ("hnum_trk_minus1sig_syst", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_minus1sig_syst->Sumw2();
  TH2F* hnum_trk_plus1sig_stat = new TH2F ("hnum_trk_plus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_plus1sig_stat->Sumw2();
  TH2F* hnum_trk_minus1sig_stat = new TH2F ("hnum_trk_minus1sig_stat", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_minus1sig_stat->Sumw2();
  

  Long64_t nentries =fChain->GetEntries();
  nentries=10;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cent = getMCHiBinFromhiHF(hiHF);
      if (!isPbPb && isPr) {
	if (pthat >= 15 && pthat < 25)  weight = 0.0247699;
	else if (pthat >= 25 && pthat < 35) weight =  0.00311931;
	else if (pthat >= 35 && pthat < 45) weight = 0.000693027;
	else if (pthat >= 45) weight = 0.000212618;
      }

      else if (isPbPb) weight = Gen_weight*findNcoll(cent); 
      else weight = 1.0;

      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{

	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m=GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();

	  if (jpsi_pt<3 || jpsi_pt>150) continue;
	  if (fabs(jpsi_rap)>2.4) continue;

	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;

	  hdeno_pty->Fill(jpsi_rap, jpsi_pt, weight);
	  if (isPbPb){
	    if (cent<=10) hdeno_pty_cent010->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=20) hdeno_pty_cent1020->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=30) hdeno_pty_cent2030->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=40) hdeno_pty_cent3040->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=60) hdeno_pty_cent4060->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=80) hdeno_pty_cent6080->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=100) hdeno_pty_cent80100->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=140) hdeno_pty_cent100140->Fill(jpsi_rap, jpsi_pt, weight);
	    else if (cent<=200) hdeno_pty_cent140200->Fill(jpsi_rap, jpsi_pt, weight);
	  }

	  hdeno_pt->Fill(jpsi_pt, weight);
	  if (jpsi_pt>6.5) 
	    hdeno_y->Fill(jpsi_rap,weight);
	 
	  int whichRec = Gen_QQ_whichRec[iQQ];
	  if (whichRec < 0) continue;
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(whichRec);
	  //cout <<"gen pt = "<<jpsi_pt<< ", gen rap = "<<jpsi_rap<<"reco pt ="<<RecoQQ4mom->Pt()<<"reco rap ="<<RecoQQ4mom->Rapidity()<<endl;
	  //cout <<"pthat = "<<Gen_pthat<<". Gen_weight = " <<Gen_weight<<", NColl = "<<findNcoll(cent)<<"cent = "<<cent<<endl;
 
	  if (!areMuonsInAcceptance2019(whichRec)) continue;
	  if (!passQualityCuts2019(whichRec)) continue;
	  //if (!isPbPb && !isTriggerMatch(whichRec, triggerIndex_PP)) continue;
	  if (Reco_QQ_sign[whichRec]!=0) continue;
	  if (RecoQQ4mom->Pt()<3 || RecoQQ4mom->Pt()>150) continue;
	  if (fabs(RecoQQ4mom->Rapidity())>2.4) continue;
	  if (RecoQQ4mom->M()<2.6 || RecoQQ4mom->M()>3.5) continue;

	  if (isPbPb){
	    if(!isTriggerMatch(whichRec, triggerIndex_PbPb)) continue;
	    if (!(pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter2Th4)) continue;
	  }
	  if (!isPbPb){
	    if (!isTriggerMatch(whichRec, triggerIndex_PP)) continue;
	    if (!(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;
	  }

	  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[whichRec]);
	  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[whichRec]);

	  double filterIdxMuPl = 0, filterIdxMuMi = 0;
	  //check trigger filters, 26=L2Filter, 27=L3Filter
	  if (isPbPb) {
	    if (isFilterMatch(Reco_QQ_mupl_idx[whichRec],27) && !isFilterMatch(Reco_QQ_mumi_idx[whichRec],27)) {
	      filterIdxMuPl=1; 
	      filterIdxMuMi=0;
	    }
	    else if (!isFilterMatch(Reco_QQ_mupl_idx[whichRec],27) && isFilterMatch(Reco_QQ_mumi_idx[whichRec],27)) {
	      filterIdxMuPl=0; 
	      filterIdxMuMi=1;
	    }
	    else { 
	      if (!(isFilterMatch(Reco_QQ_mupl_idx[whichRec],27) && isFilterMatch(Reco_QQ_mumi_idx[whichRec],27))) cout <<"the muons do not pass the filters but I will apply random SF for now"<<endl;
	      TRandom* rn = new TRandom(); 
	      filterIdxMuPl=rn->Integer(2); 
	      filterIdxMuMi=1-filterIdxMuPl;
	    }
	  }

	  tnp_weight=1.0;
	  hnum_noweights->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), weight);
	  double muplPt = RecoQQmupl4mom->Pt();
	  double muplY = RecoQQmupl4mom->Rapidity();
	  double mumiPt = RecoQQmumi4mom->Pt();
	  double mumiY = RecoQQmumi4mom->Rapidity();

	  if (isPbPb) {
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_nominal->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    //for centrality dependance
	    if (cent<=10) hnum_nominal_cent010->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=20) hnum_nominal_cent1020->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=30) hnum_nominal_cent2030->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=40) hnum_nominal_cent3040->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=60) hnum_nominal_cent4060->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=80) hnum_nominal_cent6080->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=100) hnum_nominal_cent80100->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=140) hnum_nominal_cent100140->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	    else if (cent<=200) hnum_nominal_cent140200->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,-1)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,-1)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_plus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,-2)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,-2)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_minus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,+1)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,+1)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_plus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,+2)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,+2)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_minus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,-1)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,-1);
	    hnum_trg_plus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,-2)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,-2);
	    hnum_trg_minus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,+1)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,+1);
	    hnum_trg_plus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,+2)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,+2);
	    hnum_trg_minus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,-1)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,-1)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_plus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,-2)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,-2)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_minus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,+1)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,+1)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_plus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,+2)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,+2)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_minus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	  }
	  else {
	    auto muplGlb = tnp_weight_GlobalMuon_TightAcceptance_pp(muplPt,muplY);
	    auto mumiGlb = tnp_weight_GlobalMuon_TightAcceptance_pp(mumiPt,mumiY);
	    auto muplMuIdTrg = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(muplPt,muplY);
	    auto mumiMuIdTrg = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mumiPt,mumiY);

	    tnp_weight = std::get<0>(muplGlb) * std::get<0>(muplMuIdTrg)
	      *std::get<0>(mumiGlb) * std::get<0>(mumiMuIdTrg);
	    hnum_nominal->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)+std::get<2>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)+std::get<2>(mumiMuIdTrg));
	    hnum_muidtrg_plus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)-std::get<2>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)-std::get<2>(mumiMuIdTrg));
	    hnum_muidtrg_minus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)+std::get<1>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)+std::get<1>(mumiMuIdTrg));
	    hnum_muidtrg_plus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)-std::get<1>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)-std::get<1>(mumiMuIdTrg));
	    hnum_muidtrg_minus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = (std::get<0>(muplGlb)+std::get<2>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)+std::get<2>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_plus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = (std::get<0>(muplGlb)-std::get<2>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)-std::get<2>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_minus1sig_syst->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = (std::get<0>(muplGlb)+std::get<1>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)+std::get<1>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_plus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);

	    tnp_weight = (std::get<0>(muplGlb)-std::get<1>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)-std::get<1>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_minus1sig_stat->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight*weight);
	  }
	  
	  hnum_pt->Fill(RecoQQ4mom->Pt(), weight);
	  if (RecoQQ4mom->Pt()>6.5)
	    hnum_y->Fill(RecoQQ4mom->Rapidity(),weight);
	}
    }

  gSystem->mkdir("FilesAccxEff");  
  gSystem->mkdir("FilesAccxEff/Eff");
  TFile* fsave = new TFile (Form("FilesAccxEff/Eff/%sEffHists_%s.root", isPr?"pr":"npr",isPbPb?"PbPb":"PP"), "RECREATE");
  hdeno_pty->Write();
  hnum_noweights->Write();
  hnum_nominal->Write();
  if (isPbPb) {
    hnum_muid_plus1sig_syst->Write();
    hnum_muid_minus1sig_syst->Write();
    hnum_muid_plus1sig_stat->Write();
    hnum_muid_minus1sig_stat->Write();
    hnum_trg_plus1sig_syst->Write();
    hnum_trg_minus1sig_syst->Write();
    hnum_trg_plus1sig_stat->Write();
    hnum_trg_minus1sig_stat->Write();
    hnum_trk_plus1sig_syst->Write();
    hnum_trk_minus1sig_syst->Write();
    hnum_trk_plus1sig_stat->Write();
    hnum_trk_minus1sig_stat->Write();
    hdeno_pty_cent010->Write();
    hnum_nominal_cent010->Write();
    hdeno_pty_cent1020->Write();
    hnum_nominal_cent1020->Write();
    hdeno_pty_cent2030->Write();
    hnum_nominal_cent2030->Write();
    hdeno_pty_cent3040->Write();
    hnum_nominal_cent3040->Write();
    hdeno_pty_cent4060->Write();
    hnum_nominal_cent4060->Write();
    hdeno_pty_cent6080->Write();
    hnum_nominal_cent6080->Write();
    hdeno_pty_cent80100->Write();
    hnum_nominal_cent80100->Write();
    hdeno_pty_cent100140->Write();
    hnum_nominal_cent100140->Write();
    hdeno_pty_cent140200->Write();
    hnum_nominal_cent140200->Write();
  }
  else {
    hnum_muidtrg_plus1sig_syst->Write();
    hnum_muidtrg_minus1sig_syst->Write();
    hnum_muidtrg_plus1sig_stat->Write();
    hnum_muidtrg_minus1sig_stat->Write();
    hnum_glb_plus1sig_syst->Write();
    hnum_glb_minus1sig_syst->Write();
    hnum_glb_plus1sig_stat->Write();
    hnum_glb_minus1sig_stat->Write();
  }
  hdeno_pt->Write();
  hnum_pt->Write();
  hdeno_y->Write();
  hnum_y->Write();
  
  fsave->Close();
  
  //delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_noweights; delete hnum_nominal; delete hnum_binned; delete hnum_plus1sig; delete hnum_minus1sig; delete hnum_muid_sta; delete hnum_muid; delete hnum_muid_plus1sig; delete hnum_muid_minus1sig; delete hnum_sta; delete hnum_sta_plus1sig; delete hnum_sta_minus1sig; delete hnum_trk_plus1sig; delete hnum_trk_minus1sig;
}//end of EffStep()


void oniaTree::AccCalc () {
  if (!isAcc) {cout<<"[ERROR] you're trying to make Acceptance with Efficiency trees."<<endl; return;}
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;
  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;

  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbinsAna, ptbinsAna); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbinsAna, ptbinsAna); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs y; y; N_{reco}", nybins, ybins); hnum_y->Sumw2();

  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal->Sumw2();
  
  Long64_t nentries =fChain->GetEntries();
  //nentries = 2000000;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m=GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	    
	  if (jpsi_pt<3 || jpsi_pt>150) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  hdeno_pty->Fill(jpsi_rap, jpsi_pt);

	  hdeno_pt->Fill(jpsi_pt);
	  if (jpsi_pt>6.5)
	    hdeno_y->Fill(jpsi_rap);

	  //if (!areGenMuonsInAcceptance2019(iQQ)) continue;
	  TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	  TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	  ////////mupl in acc 
	  if (!isGlobalMuonInAccept2019(GenQQmupl4mom)) continue;
	  if (!isGlobalMuonInAccept2019(GenQQmumi4mom)) continue;
 
	  hnum_nominal->Fill(jpsi_rap, jpsi_pt);
	    
	  //if (jpsi_pt<6.5) continue;
	  hnum_pt->Fill(jpsi_pt);
	  if(jpsi_pt>=6.5)
	    hnum_y->Fill(jpsi_rap);
	}//end of genQQ loop
    }//end of entries loop

  gSystem->mkdir("FilesAccxEff");
  gSystem->mkdir("FilesAccxEff/Acc");
  TFile* fsave = new TFile (Form("FilesAccxEff/Acc/%sAccHists_%s.root", isPr?"pr":"npr", isPbPb?"PbPb":"PP"), "RECREATE");
  hdeno_pty->Write("hdeno_2d");
  hnum_nominal->Write("hnum_2d_nominal");
  hdeno_pt->Write("hdeno_pt");
  hnum_pt->Write("hnum_pt");
  hdeno_y->Write("hdeno_y");
  hnum_y->Write("hnum_y");
  fsave->Close();
  
  delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_nominal;
}//end of AccCalc function

void oniaTree::ClosureTest()
{
  //int nptbins = sizeof(ptbins)/sizeof(double)-1;
  gStyle->SetOptStat(0);
  Double_t etabins []={0, 1.6, 2.4};
  TH1F* mixCount = new TH1F ("mixCount", "pt distribution at reco level with mixAccEff", 12, 6.5, 30.5);
  TH1F* genCount = new TH1F ("genCount", "pt distribution at gen level", 12, 6.5, 30.5);
  TH1F* sepCount = new TH1F ("sepCount", "pt distribution at reco level with separate AccEff", 12, 6.5, 30.5);

  TFile* corrfile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ");

  //TEfficiency* prEff = (TEfficiency*) corrfile->Get(Form("hcorr_Jpsi_%s_pr",isPbPb?"PbPb":"PP"));
  //TEfficiency* nprEff = (TEfficiency*) corrfile->Get(Form("hcorr_Jpsi_%s_npr",isPbPb?"PbPb":"PP"));

  //get the acceptance files
  TFile*accFile = TFile::Open(Form("FilesAccxEff/Acc/%sAccHists_%s.root",isPr?"pr":"npr",isPbPb?"PbPb":"PP"));
  TH2F *accNum = (TH2F*) accFile->Get("hnum_2d_nominal");
  TH2F *accDen = (TH2F*) accFile->Get("hdeno_2d");
  accNum->Divide(accDen);

  TFile*effPrFile = TFile::Open(Form("FilesAccxEff/Eff/prEffHists_%s.root",isPbPb?"PbPb":"PP"));
  TH2F *effPrNum = (TH2F*) effPrFile->Get("hnum_2d_noweights");
  TH2F *effPrDen = (TH2F*) effPrFile->Get("hdeno_2d");
  effPrNum->Divide(effPrDen);

  TFile*effNprFile = TFile::Open(Form("FilesAccxEff/Eff/nprEffHists_%s.root",isPbPb?"PbPb":"PP"));
  TH2F *effNprNum = (TH2F*) effNprFile->Get("hnum_2d_noweights");
  TH2F *effNprDen = (TH2F*) effNprFile->Get("hdeno_2d");
  effNprNum->Divide(effNprDen);

  Long64_t nentries = fChain->GetEntries();
  //nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;

  TString bfracFunc = "0.669157-0.784588*exp(-0.0731804*x)";
  if (isPbPb)
    bfracFunc = "0.611808-0.718494*exp(-0.102699*x)";
  
  TF1  *bfrac = new TF1("bfrac",bfracFunc, 3, 100);

  double bf =1.0; double prCorr = 1.0; double nprCorr = 1.0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] processing entry "<<jentry<<"/"<<nentries<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);

	  jpsi_m = GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	  if (jpsi_pt<6.5 || jpsi_pt>100) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  //apply the acceptance cuts
	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;
	 
	  weight = 1.0;
	  weight = 1.0/accNum->GetBinContent(accNum->FindBin(jpsi_rap,jpsi_pt));
	  if(isPr && !isPbPb) Gen_weight = 1;
	  weight = weight*Gen_weight;
	  genCount->Fill(jpsi_pt,weight);
	}
      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
	{
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	  jpsi_pt = RecoQQ4mom->Pt();
	  jpsi_rap = RecoQQ4mom->Rapidity();
	  jpsi_m = RecoQQ4mom->M();

	  if (!areMuonsInAcceptance2019(iQQ)) continue;
	  if (!passQualityCuts2019(iQQ)) continue;
	  if (isPbPb && !isTriggerMatch(iQQ, triggerIndex_PbPb)) continue;
	  if (!isPbPb && !isTriggerMatch(iQQ, triggerIndex_PP)) continue;
	  if (Reco_QQ_sign[iQQ]!=0) continue;
	  if (jpsi_pt<6.5 || jpsi_pt>100) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  if (jpsi_m<2.6 || jpsi_m>3.5) continue;
	  if (isPbPb && !(pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter2Th4)) continue;
	  if (!isPbPb && !(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;
	  if (Reco_QQ_whichGen[iQQ] < 0) continue;
	  if (!areGenMuonsInAcceptance2019(Reco_QQ_whichGen[iQQ])) continue;
	  if(isPr && !isPbPb) Gen_weight=1;

	  bf = bfrac->Eval(jpsi_pt);
	  //prCorr = (prEff->GetEfficiency(prEff->FindFixBin(jpsi_rap,jpsi_pt)));
	  prCorr = (effPrNum->GetBinContent(effPrNum->FindBin(jpsi_rap,jpsi_pt)))*(accNum->GetBinContent(accNum->FindBin(jpsi_rap,jpsi_pt)));
	  //nprCorr = (nprEff->GetEfficiency(nprEff->FindFixBin(jpsi_rap,jpsi_pt)));
	  nprCorr = (effNprNum->GetBinContent(effNprNum->FindBin(jpsi_rap,jpsi_pt)))*(accNum->GetBinContent(accNum->FindBin(jpsi_rap,jpsi_pt)));
	  weight = 1.0/(bf*nprCorr + (1-bf)*prCorr);

	  weight = weight*Gen_weight;
	  mixCount->Fill(jpsi_pt, weight);
	  if (isPr) weight = 1.0/prCorr;
	  else weight = 1.0/nprCorr;
	  weight = weight*Gen_weight;
	  sepCount->Fill(jpsi_pt, weight);
	}
    }
  gSystem->mkdir("FilesAccxEff");
  gSystem->mkdir("FilesAccxEff/ClosureTest");
  TFile* testfile = new TFile (Form("FilesAccxEff/ClosureTest/ClosureTest_noTnPweights_%s_%s.root",isPbPb?"PbPb":"PP",isPr?"pr":"npr"),"RECREATE");

  genCount->Write("genYields");
  mixCount->Write("mixAccEffYields");
  sepCount->Write("sepAccEffYields");
  testfile->Close();
  
  genCount->SetLineColor(kRed);
  genCount->SetLineWidth(2);
  mixCount->SetLineColor(kBlue);
  mixCount->SetLineWidth(2);
  sepCount->SetLineColor(kGreen+2);
  sepCount->SetLineWidth(2);
  sepCount->SetLineStyle(2);
  genCount->SetTitle(Form("%s Comparison", isPr?"prompt":"nonprompt"));
  TCanvas* c = new TCanvas("c","", 1000, 1000);
  //genCount->Scale(1.0,"width");
  //mixCount->Scale(1.0,"width");
  //sepCount->Scale(1.0,"width");

  genCount->Draw();
  mixCount->Draw("same");
  sepCount->Draw("same");
  TLegend* leg = new TLegend (0.6, 0.6, 0.75, 0.75);
  leg->AddEntry(genCount, "N_{gen}","lep");
  leg->AddEntry(mixCount, "N_{reco}^{mix AccEff}","lep");
  leg->AddEntry(sepCount, "N_{reco}^{sep AccEff}","lep");
  leg->SetBorderSize(0);
  leg->Draw("same");
  c->SaveAs(Form("FilesAccxEff/ClosureTest/GenRecoComparison_noTnPweights_%s_%s.pdf",isPbPb?"PbPb":"PP",isPr?"pr":"npr"));
}

void oniaTree::Plot() {cout << "[INFO] This function is empty at the moment. It can be used to make nice plots for the analysis notes"<< endl;}
