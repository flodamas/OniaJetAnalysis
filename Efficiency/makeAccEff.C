#define makeAccEff_cxx
#define _USE_MATH_DEFINES
#include "makeAccEff.h"
//#include "compAccEff.C"
//#include "systAccEff.C"

Double_t ptbins []= {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50};
Double_t ybins []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


void oniaTree::AccEffCalc()
{ 
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;

  cout<<"[INFO] Importing the numerators and denominators of the corrections."<<endl;
  TFile*prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PbPb.root");
  if (!prAccFile_pbpb) {
    cout<<"[ERROR] pbpb prompt acceptance file not found!"<<endl;
    if (isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
    AccCalc();
    prAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/prAccHists_PbPb.root");
    }
    else
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PbPb.root");
  if (!nprAccFile_pbpb) {
  cout<<"[ERROR] pbpb nonprompt acceptance file not found!"<<endl;
  if (isPbPb && !isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
  AccCalc();
  nprAccFile_pbpb = TFile::Open("FilesAccxEff/Acc/nprAccHists_PbPb.root");
  }
  else 
  cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*prEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/prEffHists_PbPb.root");
  if (!prEffFile_pbpb) {
    cout<<"[ERROR] pbpb prompt efficiency file not found!"<<endl;
    if (isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/prEffHists_PbPb.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*nprEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/nprEffHists_PbPb.root");
  if (!nprEffFile_pbpb) {
    cout<<"[ERROR] pbpb nonprompt Eff file not found!"<<endl;
    if (isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      nprEffFile_pbpb = TFile::Open("FilesAccxEff/Eff/nprEffHists_PbPb.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  ////////////////////////////////////pp/////////////////////////////////
  TFile*prAccFile_pp = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
  if (!prAccFile_pp) {
    cout<<"[ERROR] pp prompt acceptance file not found!"<<endl;
    if (!isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
    AccCalc();
    prAccFile_pp = TFile::Open("FilesAccxEff/Acc/prAccHists_PP.root");
    }
    else
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*nprAccFile_pp = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
  if (!nprAccFile_pp) {
    cout<<"[ERROR] pp nonprompt acceptance file not found!"<<endl;
    if (!isPbPb && !isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc();
      nprAccFile_pp = TFile::Open("FilesAccxEff/Acc/nprAccHists_PP.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*prEffFile_pp = TFile::Open("FilesAccxEff/Eff/prEffHists_PP.root");
  if (!prEffFile_pp) {
    cout<<"[ERROR] pp prompt efficiency file not found!"<<endl;
    if (!isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      prEffFile_pp = TFile::Open("FilesAccxEff/Eff/prEffHists_PP.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TFile*nprEffFile_pp = TFile::Open("FilesAccxEff/Eff/nprEffHists_PP.root");
  if (!nprEffFile_pp) {
    cout<<"[ERROR] pp nonprompt Eff file not found!"<<endl;
    if (!isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc();
      nprEffFile_pp = TFile::Open("FilesAccxEff/Eff/nprEffHists_PP.root");
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }


  TH2F *prAccNum_pbpb = (TH2F*) prAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *prAccDen_pbpb = (TH2F*) prAccFile_pbpb->Get("hdeno_2d");
  TH2F *nprAccNum_pbpb = (TH2F*) nprAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *nprAccDen_pbpb = (TH2F*) nprAccFile_pbpb->Get("hdeno_2d");
  TH2F *prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_2d_nominal");
  TH2F *prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get("hdeno_2d");
  TH2F *nprEffNum_pbpb = (TH2F*) nprEffFile_pbpb->Get("hnum_2d_nominal");
  TH2F *nprEffDen_pbpb = (TH2F*) nprEffFile_pbpb->Get("hdeno_2d");

  TH2F *prAccNum_pp = (TH2F*) prAccFile_pp->Get("hnum_2d_nominal");
  TH2F *prAccDen_pp = (TH2F*) prAccFile_pp->Get("hdeno_2d");
  TH2F *nprAccNum_pp = (TH2F*) nprAccFile_pp->Get("hnum_2d_nominal");
  TH2F *nprAccDen_pp = (TH2F*) nprAccFile_pp->Get("hdeno_2d");
  TH2F *prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_2d_nominal");
  TH2F *prEffDen_pp = (TH2F*) prEffFile_pp->Get("hdeno_2d");
  TH2F *nprEffNum_pp = (TH2F*) nprEffFile_pp->Get("hnum_2d_nominal");
  TH2F *nprEffDen_pp = (TH2F*) nprEffFile_pp->Get("hdeno_2d");

  prAccNum_pbpb->Multiply(prEffNum_pbpb);
  prAccDen_pbpb->Multiply(prEffDen_pbpb);
  nprAccNum_pbpb->Multiply(nprEffNum_pbpb);
  nprAccDen_pbpb->Multiply(nprEffDen_pbpb);

  prAccNum_pp->Multiply(prEffNum_pp);
  prAccDen_pp->Multiply(prEffDen_pp);
  nprAccNum_pp->Multiply(nprEffNum_pp);
  nprAccDen_pp->Multiply(nprEffDen_pp);

  TEfficiency* prCorr_pbpb = new TEfficiency("prCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pbpb->SetPassedHistogram(*prAccNum_pbpb,"f");
  prCorr_pbpb->SetTotalHistogram(*prAccDen_pbpb,"f");
  prCorr_pbpb->SetName("hcorr_Jpsi_PbPb_pr");

  TEfficiency* nprCorr_pbpb = new TEfficiency("nprCorr_pbpb", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  nprCorr_pbpb->SetStatisticOption(TEfficiency::kBBayesian);
  nprCorr_pbpb->SetPassedHistogram(*nprAccNum_pbpb,"f");
  nprCorr_pbpb->SetTotalHistogram(*nprAccDen_pbpb,"f");
  nprCorr_pbpb->SetName("hcorr_Jpsi_PbPb_npr");

  TEfficiency* prCorr_pp = new TEfficiency("prCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pp->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp->SetPassedHistogram(*prAccNum_pp,"f");
  prCorr_pp->SetTotalHistogram(*prAccDen_pp,"f");
  prCorr_pp->SetName("hcorr_Jpsi_PP_pr");

  TEfficiency* nprCorr_pp = new TEfficiency("nprCorr_pp", "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  nprCorr_pp->SetStatisticOption(TEfficiency::kBBayesian);
  nprCorr_pp->SetPassedHistogram(*nprAccNum_pp,"f");
  nprCorr_pp->SetTotalHistogram(*nprAccDen_pp,"f");
  nprCorr_pp->SetName("hcorr_Jpsi_PP_npr");

  TFile* fsave = new TFile("../Fitter/Input/correction_AccEff.root","RECREATE");
  prCorr_pbpb->Write("hcorr_Jpsi_PbPb_pr");
  nprCorr_pbpb->Write("hcorr_Jpsi_PbPb_npr");
  prCorr_pp->Write("hcorr_Jpsi_PP_pr");
  nprCorr_pp->Write("hcorr_Jpsi_PP_npr");
  fsave->Close();
}

void oniaTree::EffCalc () {  
  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs y; y; N_{reco}", nybins, ybins); hnum_y->Sumw2();
  
  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_noweights->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_nominal->Sumw2();
  TH2F* hnum_binned = new TH2F ("hnum_binned", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_binned->Sumw2();
  TH2F* hnum_plus1sig = new TH2F ("hnum_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_plus1sig->Sumw2();
  TH2F* hnum_minus1sig = new TH2F ("hnum_minus1sig", "NN_{reco}vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_minus1sig->Sumw2();
  TH2F* hnum_muid_sta = new TH2F ("hnum_muid_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_sta->Sumw2();
  TH2F* hnum_muid = new TH2F ("hnum_muid", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid->Sumw2();
  TH2F* hnum_muid_plus1sig = new TH2F ("hnum_muid_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_plus1sig->Sumw2();
  TH2F* hnum_muid_minus1sig = new TH2F ("hnum_muid_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_muid_minus1sig->Sumw2();
  TH2F* hnum_sta = new TH2F ("hnum_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_sta->Sumw2();
  TH2F* hnum_sta_plus1sig = new TH2F ("hnum_sta_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_sta_plus1sig->Sumw2();
  TH2F* hnum_sta_minus1sig = new TH2F ("hnum_sta_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_sta_minus1sig->Sumw2();
  TH2F* hnum_trk_plus1sig = new TH2F ("hnum_trk_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_plus1sig->Sumw2();
  TH2F* hnum_trk_minus1sig = new TH2F ("hnum_trk_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_trk_minus1sig->Sumw2();
  

  Long64_t nentries =fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      //if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
      cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  cout<< "eventNb = "<<eventNb<<endl;
	  cout << "Gen_QQ_size = " << Gen_QQ_size << ", Gen_mu_size = " << Gen_mu_size<<endl;
	  cout << "Reco_QQ_size = " << Reco_QQ_size << ", Reco_mu_size = " << Reco_mu_size<<endl; 
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m=GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();

	  if (jpsi_pt<3 || jpsi_pt>50) continue;
	  if (fabs(jpsi_rap)>=2.4) continue;

	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;

	  hdeno_pty->Fill(jpsi_rap, jpsi_pt);

	  hdeno_pt->Fill(jpsi_pt);
	  if (jpsi_pt>6.5) 
	    hdeno_y->Fill(fabs(jpsi_rap));
	 
	  int whichRec = Gen_QQ_whichRec[iQQ];
	  if (whichRec < 0) continue;
	  cout <<"Matched, Gen_QQ_whichRec[iQQ] = "<<Gen_QQ_whichRec[iQQ]<<endl;
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(whichRec);

	  if (!areMuonsInAcceptance2019(whichRec)) continue;
	  if (!passQualityCuts2019(whichRec)) continue;
	  if (!isTriggerMatch(whichRec, triggerIndex_PP)) continue;
	  if (Reco_QQ_sign[whichRec]!=0) continue;
	  if (RecoQQ4mom->Pt()<3 || RecoQQ4mom->Pt()>50) continue;
	  if (fabs(RecoQQ4mom->Rapidity())>=2.4) continue;
	  if (RecoQQ4mom->M()<2.6 || RecoQQ4mom->M()>3.5) continue;
	  if (isPbPb && !(pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter2Th4)) continue;
	  if (!isPbPb && !(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;

	  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[whichRec]);
	  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[whichRec]);

	  tnp_weight=1.0;
	  hnum_noweights->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt());
	  
	  //nominal
	  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
	    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	  hnum_nominal->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), tnp_weight);
	  
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
	  
	  hnum_pt->Fill(RecoQQ4mom->Pt());
	  if (RecoQQ4mom->Pt()>6.5)
	    hnum_y->Fill(fabs(RecoQQ4mom->Rapidity()));
	}
    }
  gSystem->mkdir("FilesAccxEff");  
  gSystem->mkdir("FilesAccxEff/Eff");
  TFile* fsave = new TFile (Form("FilesAccxEff/Eff/%sEffHists_%s.root", isPr?"pr":"npr",isPbPb?"PbPb":"PP"), "RECREATE");
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

  int nptbins = sizeof(ptbins)/sizeof(double)-1;
  int nybins = sizeof(ybins)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins); hnum_pt->Sumw2();
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","N_{reco} vs y; y; N_{reco}", nybins, ybins); hnum_y->Sumw2();

  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", nybins, ybins, nptbins, ptbins); hnum_noweights->Sumw2();
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
	    
	  if (jpsi_pt<3 || jpsi_pt>50) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  hdeno_pty->Fill(jpsi_rap, jpsi_pt);

	  hdeno_pt->Fill(jpsi_pt);
	  if (jpsi_pt>6.5)
	    hdeno_y->Fill(fabs(jpsi_rap));

	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;

	  hnum_noweights->Fill(jpsi_rap, jpsi_pt);
	  hnum_nominal->Fill(jpsi_rap, jpsi_pt);
	    
	  //if (jpsi_pt<6.5) continue;
	  hnum_pt->Fill(jpsi_pt);
	  if(jpsi_pt>=6.5)
	    hnum_y->Fill(fabs(jpsi_rap));
	}//end of genQQ loop
    }//end of entries loop

  gSystem->mkdir("FilesAccxEff");
  gSystem->mkdir("FilesAccxEff/Acc");
  TFile* fsave = new TFile (Form("FilesAccxEff/Acc/%sAccHists_%s.root", isPr?"pr":"npr", isPbPb?"PbPb":"PP"), "RECREATE");
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

void oniaTree::ClosureTest()
{
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
	  
	  if (jpsi_pt>6.5 && jpsi_pt<35 && fabs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
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
		  (areMuonsInAcceptance2019(iQQ))&&  // 2019 Global Muon Acceptance Cuts
		  (passQualityCuts2019(iQQ)) &&  // 2019 Soft Global Muon Quality Cuts
		  (isTriggerMatch(iQQ, triggerIndex_PP)) &&// if it matches the trigger 
		  (Reco_QQ_whichGen[iQQ]!=-1)
		  )
		{
		  if (Reco_QQ_sign[iQQ]==0 && fabs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
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

void oniaTree::Plot() {cout << "[INFO] This function is empty at the moment. It can be used to make nice plots for the analysis notes"<< endl;}
