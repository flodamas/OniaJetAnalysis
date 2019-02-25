#define makeAccEff_cxx
#define _USE_MATH_DEFINES
#include "makeAccEff.h"
#include "compAccEff.C"
#include "systAccEff.C"

Double_t ptbins2D []= {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50};
Double_t ybins2D []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


void oniaTree::AccEffCalc()
{ 
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
 
  cout<<"[INFO] Importing the numerators and denominators of the corrections."<<endl;
  TFile*prAccFile = TFile::Open("FilesAccxEff/Acc/prAccHists.root");
  if (!prAccFile) {
    cout<<"[ERROR] prompt Acc file not found!"<<endl;
    if (isAcc && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
    AccCalc();
    prAccFile = TFile::Open("FilesAccxEff/Acc/prAccHists.root");
    }
  }
  if (!prAccFile) {
    cout<<"[ERROR] Please change your settings and use AccCalc()."<<endl; return;
  }

  TFile*nprAccFile = TFile::Open("FilesAccxEff/Acc/nprAccHists.root");
  if (!nprAccFile) {
    cout<<"[ERROR] nonprompt Acc file not found!"<<endl;
    if (isAcc && !isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc();
      nprAccFile = TFile::Open("FilesAccxEff/Acc/nprAccHists.root");
    }
  }
  if (!prAccFile) {
    cout<<"[ERROR] Please change your settings and use AccCalc()."<<endl; return;
  }

  TFile*prEffFile = TFile::Open("FilesAccxEff/Eff/prEffHists.root");
  if (!prEffFile) {
    cout<<"[ERROR] prompt Eff file not found!"<<endl;
    if (!isAcc && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile = TFile::Open("FilesAccxEff/Eff/prEffHists.root");
    }
  }
  if (!prEffFile) {
    cout<<"[ERROR] Please change your settings and use EffCalc()."<<endl; return;
  }

  TFile*nprEffFile = TFile::Open("FilesAccxEff/Eff/nprEffHists.root");
  if (!nprEffFile) {
    cout<<"[ERROR] nonprompt Eff file not found!"<<endl;
    if (!isAcc && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile = TFile::Open("FilesAccxEff/Eff/nprEffHists.root");
    }
  }
  if (!prEffFile) {
    cout<<"[ERROR] Please change your settings and use EffCalc()."<<endl; return;
  }

  TH2F *prAccNum = (TH2F*) prAccFile->Get("hnum_2d_nominal");
  TH2F *prAccDen = (TH2F*) prAccFile->Get("hdeno_2d");
  TH2F *nprAccNum = (TH2F*) nprAccFile->Get("hnum_2d_nominal");
  TH2F *nprAccDen = (TH2F*) nprAccFile->Get("hdeno_2d");
  TH2F *prEffNum = (TH2F*) prEffFile->Get("hnum_2d_nominal");
  TH2F *prEffDen = (TH2F*) prEffFile->Get("hdeno_2d");
  TH2F *nprEffNum = (TH2F*) nprEffFile->Get("hnum_2d_nominal");
  TH2F *nprEffDen = (TH2F*) nprEffFile->Get("hdeno_2d");

  prAccNum->Multiply(prEffNum);
  prAccDen->Multiply(prEffDen);
  nprAccNum->Multiply(nprEffNum);
  nprAccDen->Multiply(nprEffDen);

  TEfficiency* prCorr = new TEfficiency("prCorr", "AccxEff(y,pt); y; pt; eff", ny2D, ybins2D, npt2D, ptbins2D);
  prCorr->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr->SetPassedHistogram(*prAccNum,"f");
  prCorr->SetTotalHistogram(*prAccDen,"f");
  prCorr->SetName("hcorr_Jpsi_PP_pr");

  TEfficiency* nprCorr = new TEfficiency("nprCorr", "AccxEff(y,pt); y; pt; eff", ny2D, ybins2D, npt2D, ptbins2D);
  nprCorr->SetStatisticOption(TEfficiency::kBBayesian);
  nprCorr->SetPassedHistogram(*nprAccNum,"f");
  nprCorr->SetTotalHistogram(*nprAccDen,"f");
  nprCorr->SetName("hcorr_Jpsi_PP_npr");

  TFile* fsave = new TFile("../Fitter/Input/correction_AccEff.root","RECREATE");
  prCorr->Write("hcorr_Jpsi_PP_pr");
  nprCorr->Write("hcorr_Jpsi_PP_npr");
  fsave->Close();
}

void oniaTree::EffCalc () {
  if (isAcc){
    cout<< "[ERROR] check your settings. This has to be MC and not Acc"<<endl;
    return;
  }
  Double_t ptbinsRaa []= {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};
  Double_t ybinsRaa []= {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;
  int nptRaa = sizeof(ptbinsRaa)/sizeof(double)-1;
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int nyRaa = sizeof(ybinsRaa)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptRaa, ptbinsRaa); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptRaa, ptbinsRaa); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs #eta; #eta; N_{total}", nyRaa, ybinsRaa); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs #eta; #eta; N_{reco}", nyRaa, ybinsRaa); hnum_y->Sumw2();
  
  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", ny2D, ybins2D, npt2D, ptbins2D); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_noweights->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_nominal->Sumw2();
  TH2F* hnum_binned = new TH2F ("hnum_binned", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_binned->Sumw2();
  TH2F* hnum_plus1sig = new TH2F ("hnum_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_plus1sig->Sumw2();
  TH2F* hnum_minus1sig = new TH2F ("hnum_minus1sig", "NN_{reco}vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_minus1sig->Sumw2();
  TH2F* hnum_muid_sta = new TH2F ("hnum_muid_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_muid_sta->Sumw2();
  TH2F* hnum_muid = new TH2F ("hnum_muid", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_muid->Sumw2();
  TH2F* hnum_muid_plus1sig = new TH2F ("hnum_muid_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_muid_plus1sig->Sumw2();
  TH2F* hnum_muid_minus1sig = new TH2F ("hnum_muid_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_muid_minus1sig->Sumw2();
  TH2F* hnum_sta = new TH2F ("hnum_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_sta->Sumw2();
  TH2F* hnum_sta_plus1sig = new TH2F ("hnum_sta_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_sta_plus1sig->Sumw2();
  TH2F* hnum_sta_minus1sig = new TH2F ("hnum_sta_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_sta_minus1sig->Sumw2();
  TH2F* hnum_trk_plus1sig = new TH2F ("hnum_trk_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_trk_plus1sig->Sumw2();
  TH2F* hnum_trk_minus1sig = new TH2F ("hnum_trk_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_trk_minus1sig->Sumw2();
  
  vector<TObjArray *> wHistograms = ReadFileWeight(isPr);
  
  Long64_t nentries =fChain->GetEntries();
  //nentries = 200;
  
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

	  if (jpsi_pt<3 || jpsi_pt>50) continue;
	  if (abs(jpsi_rap)>=2.4) continue;
	  if (!areGenMuonsInAcceptance2015(iQQ)) continue;

	  hdeno_pty->Fill(jpsi_rap, jpsi_pt);

	  TH1D *curve;
	  if (abs(jpsi_rap)>=0 && abs(jpsi_rap)<0.6)        curve = (TH1D*) wHistograms[0]->At(0);
	  else if (abs(jpsi_rap)>=0.6 && abs(jpsi_rap)<1.2) curve = (TH1D*) wHistograms[1]->At(0);
	  else if (abs(jpsi_rap)>=1.2 && abs(jpsi_rap)<1.8) curve = (TH1D*) wHistograms[2]->At(0);
	  else if (abs(jpsi_rap)>=1.8 && abs(jpsi_rap)<2.4) curve = (TH1D*) wHistograms[3]->At(0);
	  int ptbin = 0;
	  if (jpsi_pt>=6.5)
	  ptbin = curve->FindBin(jpsi_pt);
	  double ptweight = 1.0;
	  if (jpsi_pt>6.5)
	    ptweight = curve->GetBinContent(ptbin);

	  if (jpsi_pt>=6.5)
	    {
	      hdeno_pt->Fill(jpsi_pt, ptweight);
	      hdeno_y->Fill(abs(jpsi_rap), ptweight);
	    }
	 
	  ibest = -1;
	  if (!isMatchedGenDiMuon(iQQ)) continue;//fill the num with the matching reco 
	  if (ibest < 0) continue;
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(ibest);
	  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(ibest);
	  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(ibest);
	  if (!areMuonsInAcceptance2015(ibest)) continue;
	  if (!passQualityCuts2015(ibest)) continue;
	  if (!isTriggerMatch(ibest, triggerIndex_PP)) continue;
	  if (Reco_QQ_sign[ibest]!=0) continue;
	  if (RecoQQ4mom->Pt()<3 || RecoQQ4mom->Pt()>50) continue;
	  if (abs(RecoQQ4mom->Rapidity())>=2.4) continue;
	  if (RecoQQ4mom->M()<2.6 || RecoQQ4mom->M()>3.5) continue;
	  
	  tnp_weight=1.0;
	  hnum_noweights->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  //nominal
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_nominal->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  ptweight = ptweight*tnp_weight; //to use in the 1D histograms
	  
	  //systematics
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-10) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-10) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_binned->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_plus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_minus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_muid_sta->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_muid->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_muid_plus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_muid_minus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_sta->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
        
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_sta_plus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_sta_minus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(-1) * tnp_weight_trk_pp(-1);
	  hnum_trk_plus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(-2) * tnp_weight_trk_pp(-2);
	  hnum_trk_minus1sig->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
	  if (RecoQQ4mom->Pt()<6.5) continue;
	  hnum_pt->Fill(RecoQQ4mom->Pt(), ptweight);
	  hnum_y->Fill(abs(RecoQQ4mom->Rapidity()), ptweight);
	}
    }
  gSystem->mkdir("FilesAccxEff");  
  gSystem->mkdir("FilesAccxEff/Eff");
  TFile* fsave = new TFile (Form("FilesAccxEff/Eff/%sEffHists.root", isPr?"pr":"npr"), "RECREATE");
  hdeno_pty->Write("hdeno_2d");
  hnum_nominal->Write("hnum_2d_nominal");
  hnum_binned->Write("hnum_2d_binned");
  hnum_plus1sig->Write("hnum_2d_trg_plus1sig");
  hnum_minus1sig->Write("hnum_2d_trg_minus1sig");
  hnum_muid_sta->Write("hnum_2d_muid_sta");
  hnum_muid->Write("hnum_2d_muid");
  hnum_muid_plus1sig->Write("hnum_2d_muid_plus1sig");
  hnum_muid_minus1sig->Write("hnum_2d_muid_minus1sig");
  hnum_sta->Write("hnum_2d_sta");
  hnum_sta_plus1sig->Write("hnum_2d_sta_plus1sig");
  hnum_sta_minus1sig->Write("hnum_2d_sta_minus1sig");
  hnum_trk_plus1sig->Write("hnum_2d_trk_plus1sig");
  hnum_trk_minus1sig->Write("hnum_2d_trk_minus1sig");
  hnum_noweights->Write("hnum_2d_noweights");
  hdeno_pt->Write("hdeno_pt");
  hnum_pt->Write("hnum_pt");
  hdeno_y->Write("hdeno_y");
  hnum_y->Write("hnum_y");
  fsave->Close();
  
  delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_noweights; delete hnum_nominal; delete hnum_binned; delete hnum_plus1sig; delete hnum_minus1sig; delete hnum_muid_sta; delete hnum_muid; delete hnum_muid_plus1sig; delete hnum_muid_minus1sig; delete hnum_sta; delete hnum_sta_plus1sig; delete hnum_sta_minus1sig; delete hnum_trk_plus1sig; delete hnum_trk_minus1sig;
}//end of EffStep()


void oniaTree::AccCalc () {
  if (!isAcc) {
    cout<< "[ERROR] check your settings. This has to be MC and Acc"<<endl;
    return;
  }
  Double_t ptbinsRaa []= {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};
  Double_t ybinsRaa []= {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;
  int nptRaa = sizeof(ptbinsRaa)/sizeof(double)-1;
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int nyRaa = sizeof(ybinsRaa)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptRaa, ptbinsRaa); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptRaa, ptbinsRaa); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs #eta; #eta; N_{total}", nyRaa, ybinsRaa); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs #eta; #eta; N_{reco}", nyRaa, ybinsRaa); hnum_y->Sumw2();

  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", ny2D, ybins2D, npt2D, ptbins2D); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_noweights->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", ny2D, ybins2D, npt2D, ptbins2D); hnum_nominal->Sumw2();
  
  vector<TObjArray *> wHistograms = ReadFileWeight(isPr);
  
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
	    
	  if (jpsi_pt<3 || jpsi_pt>50) continue;
	  if (abs(jpsi_rap)>2.4) continue;
	  hdeno_pty->Fill(jpsi_rap, jpsi_pt);

	  TH1D *curve;
	  if (abs(jpsi_rap)>=0 && abs(jpsi_rap)<0.6)        curve = (TH1D*) wHistograms[0]->At(0);
	  else if (abs(jpsi_rap)>=0.6 && abs(jpsi_rap)<1.2) curve = (TH1D*) wHistograms[1]->At(0);
	  else if (abs(jpsi_rap)>=1.2 && abs(jpsi_rap)<1.8) curve = (TH1D*) wHistograms[2]->At(0);
	  else if (abs(jpsi_rap)>=1.8 && abs(jpsi_rap)<2.4) curve = (TH1D*) wHistograms[3]->At(0);
	  int ptbin = 0;
	  if (jpsi_pt>=6.5)
	    ptbin = curve->FindBin(jpsi_pt);
	  double ptweight = 1.0;
	  if (jpsi_pt>=6.5)
	    ptweight = curve->GetBinContent(ptbin);

	  if (ptweight==0 && jpsi_pt > 6.5) cout<<"ptweight = 0 for genpt = "<<jpsi_pt<<endl;

	  if (jpsi_pt>=6.5) {
	    hdeno_pt->Fill(jpsi_pt, ptweight);
	    hdeno_y->Fill(abs(jpsi_rap), ptweight);
	  }

	  if (!areGenMuonsInAcceptance2015(iQQ)) continue;

	  hnum_noweights->Fill(jpsi_rap, jpsi_pt);
	  hnum_nominal->Fill(jpsi_rap, jpsi_pt);
	    
	  //if (jpsi_pt<6.5) continue;
	  if(jpsi_pt>=6.5){
	    hnum_pt->Fill(jpsi_pt, ptweight);
	    hnum_y->Fill(abs(jpsi_rap), ptweight);
	  }

	}//end of genQQ loop
    }//end of entries loop

  gSystem->mkdir("FilesAccxEff");
  gSystem->mkdir("FilesAccxEff/Acc");
  TFile* fsave = new TFile (Form("FilesAccxEff/Acc/%sAccHists.root", isPr?"pr":"npr"), "RECREATE");
  hdeno_pty->Write("hdeno_2d");
  hnum_nominal->Write("hnum_2d_nominal");
  hnum_noweights->Write("hnum_2d_noweights");
  hdeno_pt->Write("hdeno_pt");
  hnum_pt->Write("hnum_pt");
  hdeno_y->Write("hdeno_y");
  hnum_y->Write("hnum_y");
  fsave->Close();
  
  delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_noweights; delete hnum_nominal;
}//end of AccCalc function


void oniaTree::ANEffPlots()
{
}

void oniaTree::ClosureTest()
{
  if(isAcc){
    cout<<"[ERROR] change the initial settings";
    return;
  }
  gStyle->SetOptStat(0);
  Double_t etabins []={0, 1.6, 2.4};
  TH1F* mixCount = new TH1F ("mixCount", "y distribution at reco level with mixAccEff", 2, etabins);
  TH1F* genCount = new TH1F ("genCount", "y distribution at gen level", 2, etabins);
  TH1F* sepCount = new TH1F ("sepCount", "y distribution at reco level with separate AccEff", 2, etabins);
  TFile* prf (0x0);
  TFile* nprf (0x0);
  prf = TFile::Open("Utilities/pr_correction_AccEff.root","READ");
  nprf =  TFile::Open("Utilities/npr_correction_AccEff.root","READ");

  TEfficiency* prEff = (TEfficiency*) prf->Get("hcorr_Jpsi_PP");
  TEfficiency* nprEff = (TEfficiency*) nprf->Get("hcorr_Jpsi_PP");

  Long64_t nentries = fChain->GetEntries();
  //nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;

  TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
  double bf =1.0; double prW = 1.0; double nprW = 1.0;
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
	  
	  if (jpsi_pt>6.5 && jpsi_pt<35 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
	    genCount->Fill(jpsi_rap);
	}
	  for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
	    {
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	      jpsi_pt = RecoQQ4mom->Pt();
	      jpsi_rap = RecoQQ4mom->Rapidity();
	      jpsi_m = RecoQQ4mom->M();
	      bf = bfrac->Eval(jpsi_pt);
	      if (
		  jpsi_pt > 6.5  && jpsi_pt<35 &&
		  (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		  (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		  (isTriggerMatch(iQQ, triggerIndex_PP)) &&// if it matches the trigger 
		  (isMatchedRecoDiMuon(iQQ))
		  )
		{
		  if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
		    {
		      prW = (prEff->GetEfficiency(prEff->FindFixBin(jpsi_rap,jpsi_pt)));
		      nprW = (nprEff->GetEfficiency(nprEff->FindFixBin(jpsi_rap,jpsi_pt)));
		      weight = 1.0/(bf*nprW + (1-bf)*prW);
		      
		      mixCount->Fill(jpsi_rap, weight);
		      if (isPr) weight = 1.0/prW;
		      else weight = 1.0/nprW;
		      sepCount->Fill(jpsi_rap, weight);
		    }
		}
	    }
    }
  TFile* testfile (0x0);
  if (isPr)
	testfile = new TFile ("RatioStudy/prClosureTest.root","RECREATE");
  else
    testfile = new TFile ("RatioStudy/nprClosureTest.root","RECREATE");
  genCount->Write("genYields");
  mixCount->Write("mixAccEffYields");
  sepCount->Write("sepAccEffYields");
  testfile->Close();
  
  genCount->SetLineColor(kRed);
  genCount->SetLineWidth(2);
  mixCount->SetLineColor(kBlue);
  mixCount->SetLineWidth(1);
  sepCount->SetLineColor(kGreen+2);
  sepCount->SetLineWidth(1);
  sepCount->SetLineStyle(2);
  genCount->SetTitle(Form("%s Comparison", isPr?"prompt":"nonprompt"));
  TCanvas* c = new TCanvas("c","", 1000, 1000);
  genCount->Draw();
  mixCount->Draw("same");
  sepCount->Draw("same");
  TLegend* leg = new TLegend (0.1, 0.75, 0.25, 0.9);
  leg->AddEntry(genCount, "N_{gen}","lep");
  leg->AddEntry(mixCount, "N_{reco}^{mix AccEff}","lep");
  leg->AddEntry(sepCount, "N_{reco}^{sep AccEff}","lep");
  leg->SetBorderSize(0);
  leg->Draw("same");
  c->SaveAs(Form("RatioStudy/%sGenRecoComparison.root",isPr?"pr":"npr"));
}

void oniaTree::Loop() {cout << "[INFO] This function is empty at the moment!!"<< endl;}

void oniaTree::Plot() {cout << "[INFO] This function is empty at the moment!!"<< endl;}
