#define makeAccEff_cxx
#define _USE_MATH_DEFINES
#include "systAccEff.C"

Double_t ptbinsAna []= {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};

void oniaTree::AccEffCalc(string caseTag)
{ 
  setBins(caseTag);

  //string centTag [] = {"cent010","cent1020","cent2030","cent3040","cent4060","cent6080","cent80100","cent100140","cent140200"};
  //ncentbins = sizeof(centTag)/sizeof(centTag[0]);

  cout<<"[INFO] Importing the numerators and denominators of the corrections."<<endl;
  TFile*prAccFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
  if (!prAccFile_pbpb) {
    cout<<"[ERROR] pbpb prompt acceptance file not found!"<<endl;
    if (isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc(caseTag);
      prAccFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  TFile*prEffFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PbPb.root",caseTag.c_str()));
  if (!prEffFile_pbpb) {
    cout<<"[ERROR] pbpb prompt efficiency file not found!"<<endl;
    if (!isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc(caseTag);
      prEffFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PbPb.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  ////////////////////////////////////pp/////////////////////////////////
  TFile*prAccFile_pp = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
  if (!prAccFile_pp) {
    cout<<"[ERROR] pp prompt acceptance file not found!"<<endl;
    if (isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc(caseTag);
      prAccFile_pp = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  TFile*prEffFile_pp = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PP.root",caseTag.c_str()));
  if (!prEffFile_pp) {
    cout<<"[ERROR] pp prompt efficiency file not found!"<<endl;
    if (!isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc(caseTag);
      prEffFile_pp = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PP.root",caseTag.c_str()));
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TH2F *prAccNum_pbpb = (TH2F*) prAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *prAccDen_pbpb = (TH2F*) prAccFile_pbpb->Get("hdeno_2d");
  TH2F *prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_nominal");
  if (caseTag.find("noTnpWeights")!=std::string::npos)
    prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_noweights");
  TH2F *prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get("hdeno_pty");

  TH2F *prAccNum_pp = (TH2F*) prAccFile_pp->Get("hnum_2d_nominal");
  TH2F *prAccDen_pp = (TH2F*) prAccFile_pp->Get("hdeno_2d");
  TH2F *prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_nominal");
  if (caseTag.find("noTnpWeights")!=std::string::npos)
    prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_noweights");
  TH2F *prEffDen_pp = (TH2F*) prEffFile_pp->Get("hdeno_pty");

  TH1F *prAccNum_pbpb_1D = (TH1F*) prAccFile_pbpb->Get("hnum_pt");
  TH1F *prAccDen_pbpb_1D = (TH1F*) prAccFile_pbpb->Get("hdeno_pt");
  TH1F *prEffNum_pbpb_1D = (TH1F*) prEffFile_pbpb->Get("hnum_pt");
  TH1F *prEffDen_pbpb_1D = (TH1F*) prEffFile_pbpb->Get("hdeno_pt");

  TH1F *prAccNum_pp_1D = (TH1F*) prAccFile_pp->Get("hnum_pt");
  TH1F *prAccDen_pp_1D = (TH1F*) prAccFile_pp->Get("hdeno_pt");
  TH1F *prEffNum_pp_1D = (TH1F*) prEffFile_pp->Get("hnum_pt");
  TH1F *prEffDen_pp_1D = (TH1F*) prEffFile_pp->Get("hdeno_pt");


  TEfficiency* prCorr_pbpb_Eff = new TEfficiency("prCorr_pbpb_Eff", "Eff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pbpb_Eff->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pbpb_Eff->SetPassedHistogram(*prEffNum_pbpb,"f");
  prCorr_pbpb_Eff->SetTotalHistogram(*prEffDen_pbpb,"f");
  prCorr_pbpb_Eff->SetName("hcorr_Jpsi_PbPb_pr_Eff");

  TEfficiency* prCorr_pp_Eff = new TEfficiency("prCorr_pp_Eff", "Eff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
  prCorr_pp_Eff->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp_Eff->SetPassedHistogram(*prEffNum_pp,"f");
  prCorr_pp_Eff->SetTotalHistogram(*prEffDen_pp,"f");
  prCorr_pp_Eff->SetName("hcorr_Jpsi_PP_pr_Eff");

  TEfficiency* prCorr_pp_Acc = new TEfficiency("prCorr_pp_Acc", "Acc(y,pt); y; pt; acc", nybins, ybins, nptbins, ptbins);
  prCorr_pp_Acc->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp_Acc->SetPassedHistogram(*prAccNum_pp,"f");
  prCorr_pp_Acc->SetTotalHistogram(*prAccDen_pp,"f");
  prCorr_pp_Acc->SetName("hcorr_Jpsi_PP_pr_Acc");

  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;

  TEfficiency* prCorr_pbpb_1D_Eff = new TEfficiency("prCorr_pbpb_1D_Eff", "Eff(pt); pt; eff", nptbinsAna, ptbinsAna);
  if (caseTag.find("1DfinerBins")!=std::string::npos) { prCorr_pbpb_1D_Eff = new TEfficiency("prCorr_pbpb_1D_Eff", "Eff(pt); pt; eff", nptbins, ptbins); }
  prCorr_pbpb_1D_Eff->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pbpb_1D_Eff->SetPassedHistogram(*prEffNum_pbpb_1D,"f");
  prCorr_pbpb_1D_Eff->SetTotalHistogram(*prEffDen_pbpb_1D,"f");
  prCorr_pbpb_1D_Eff->SetName("hcorr_Jpsi_PbPb_pr_1D_Eff");

  TEfficiency* prCorr_pp_1D_Eff = new TEfficiency("prCorr_pp_1D_Eff", "Eff(pt); pt; eff", nptbinsAna, ptbinsAna);
  if (caseTag.find("1DfinerBins")!=std::string::npos) { prCorr_pp_1D_Eff = new TEfficiency("prCorr_pp_1D_Eff", "Eff(pt); pt; eff", nptbins, ptbins); }
  prCorr_pp_1D_Eff->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp_1D_Eff->SetPassedHistogram(*prEffNum_pp_1D,"f");
  prCorr_pp_1D_Eff->SetTotalHistogram(*prEffDen_pp_1D,"f");
  prCorr_pp_1D_Eff->SetName("hcorr_Jpsi_PP_pr_1D_Eff");

  TEfficiency* prCorr_pp_1D_Acc = new TEfficiency("prCorr_pp_1D_Acc", "Acc(pt); pt; eff", nptbinsAna, ptbinsAna);
  if (caseTag.find("1DfinerBins")!=std::string::npos) { prCorr_pp_1D_Acc = new TEfficiency("prCorr_pp_1D_Acc", "Acc(pt); pt; eff", nptbins, ptbins); }
  prCorr_pp_1D_Acc->SetStatisticOption(TEfficiency::kBBayesian);
  prCorr_pp_1D_Acc->SetPassedHistogram(*prAccNum_pp_1D,"f");
  prCorr_pp_1D_Acc->SetTotalHistogram(*prAccDen_pp_1D,"f");
  prCorr_pp_1D_Acc->SetName("hcorr_Jpsi_PP_pr_1D_Acc");
  

  TFile* fsave = new TFile(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str()),"RECREATE");
  prCorr_pbpb_Eff->Write("hcorr_Jpsi_PbPb_pr_Eff");
  prCorr_pp_Eff->Write("hcorr_Jpsi_PP_pr_Eff");
  prCorr_pp_Acc->Write("hcorr_Jpsi_PP_pr_Acc");

  if (caseTag.find("centBins")!=std::string::npos) {
    for (int i=0; i<ncentbins; i++) {
      prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hnum_nominal_cent%d%d",centbins[i],centbins[i+1]));
      prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hdeno_pty_cent%d%d",centbins[i],centbins[i+1]));
      prCorr_pbpb_Eff = new TEfficiency(Form("prCorr_pbpb_cent%d%d",centbins[i],centbins[i+1]), "Eff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
      prCorr_pbpb_Eff->SetStatisticOption(TEfficiency::kBBayesian);
      prCorr_pbpb_Eff->SetPassedHistogram(*prEffNum_pbpb,"f");
      prCorr_pbpb_Eff->SetTotalHistogram(*prEffDen_pbpb,"f");
      prCorr_pbpb_Eff->SetName(Form("hcorr_Jpsi_PbPb_pr_Eff_cent%d%d",centbins[i],centbins[i+1]));
      prCorr_pbpb_Eff->Write(Form("hcorr_Jpsi_PbPb_pr_Eff_cent%d%d",centbins[i],centbins[i+1]));
    }
  }
  prCorr_pbpb_1D_Eff->Write("hcorr_Jpsi_PbPb_pr_1D_Eff");
  prCorr_pp_1D_Eff->Write("hcorr_Jpsi_PP_pr_1D_Eff");
  prCorr_pp_1D_Acc->Write("hcorr_Jpsi_PP_pr_1D_Acc");
  fsave->Close();
}

void oniaTree::EffCalc (string caseTag) {  
  setBins(caseTag);

  if (isAcc) {cout<<"[ERROR] you're trying to make Efficiency with Acceptance trees."<<endl; return;}

  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;
  
  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbinsAna, ptbinsAna); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "Eff vs p_{T}; p_{T}; Eff", nptbinsAna, ptbinsAna); hnum_pt->Sumw2();
  if (caseTag.find("1DfinerBins")!=std::string::npos) {
    hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); hdeno_pt->Sumw2();
    hnum_pt = new TH1F ("hnum_pt", "Eff vs p_{T}; p_{T}; Eff", nptbins, ptbins); hnum_pt->Sumw2();
  }
  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","Eff vs y; y; Eff", nybins, ybins); hnum_y->Sumw2();


  TH1F* hdeno_cent = new TH1F ("hdeno_cent", "N_{gen} vs centrality; hiBin; N_{total}", 20, 0, 200); hdeno_cent->Sumw2();
  TH1F* hnum_cent = new TH1F ("hnum_cent", "Eff vs centrality; hiBin; Eff", 20, 0, 200); hnum_cent->Sumw2();
  
  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_noweights = new TH2F ("hnum_noweights", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_noweights->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_nominal->Sumw2();

  TH2F** hdeno_pty_cent = new TH2F*[ncentbins];
  TH2F** hnum_nominal_cent = new TH2F*[ncentbins];
  for (int i=0; i<ncentbins; i++) {
    hdeno_pty_cent[i] = new TH2F (Form("hdeno_pty_cent%d%d",centbins[i],centbins[i+1]), "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty_cent[i]->Sumw2();
    hnum_nominal_cent[i] = new TH2F (Form("hnum_nominal_cent%d%d",centbins[i],centbins[i+1]), "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_nominal_cent[i]->Sumw2();
  }
  //for pp
  TH2F* hnum_muidtrg_plus1sig_syst = new TH2F ("hnum_muidtrg_plus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muidtrg_plus1sig_syst->Sumw2();
  TH2F* hnum_muidtrg_minus1sig_syst = new TH2F ("hnum_muidtrg_minus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muidtrg_minus1sig_syst->Sumw2();
  TH2F* hnum_muidtrg_plus1sig_stat = new TH2F ("hnum_muidtrg_plus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muidtrg_plus1sig_stat->Sumw2();
  TH2F* hnum_muidtrg_minus1sig_stat = new TH2F ("hnum_muidtrg_minus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muidtrg_minus1sig_stat->Sumw2();

  TH2F* hnum_glb_plus1sig_syst = new TH2F ("hnum_glb_plus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_glb_plus1sig_syst->Sumw2();
  TH2F* hnum_glb_minus1sig_syst = new TH2F ("hnum_glb_minus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_glb_minus1sig_syst->Sumw2();
  TH2F* hnum_glb_plus1sig_stat = new TH2F ("hnum_glb_plus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_glb_plus1sig_stat->Sumw2();
  TH2F* hnum_glb_minus1sig_stat = new TH2F ("hnum_glb_minus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_glb_minus1sig_stat->Sumw2();

  //for PbPb
  TH2F* hnum_muid_plus1sig_syst = new TH2F ("hnum_muid_plus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muid_plus1sig_syst->Sumw2();
  TH2F* hnum_muid_minus1sig_syst = new TH2F ("hnum_muid_minus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muid_minus1sig_syst->Sumw2();
  TH2F* hnum_muid_plus1sig_stat = new TH2F ("hnum_muid_plus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muid_plus1sig_stat->Sumw2();
  TH2F* hnum_muid_minus1sig_stat = new TH2F ("hnum_muid_minus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_muid_minus1sig_stat->Sumw2();

  TH2F* hnum_trg_plus1sig_syst = new TH2F ("hnum_trg_plus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trg_plus1sig_syst->Sumw2();
  TH2F* hnum_trg_minus1sig_syst = new TH2F ("hnum_trg_minus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trg_minus1sig_syst->Sumw2();
  TH2F* hnum_trg_plus1sig_stat = new TH2F ("hnum_trg_plus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trg_plus1sig_stat->Sumw2();
  TH2F* hnum_trg_minus1sig_stat = new TH2F ("hnum_trg_minus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trg_minus1sig_stat->Sumw2();

  TH2F* hnum_trk_plus1sig_syst = new TH2F ("hnum_trk_plus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trk_plus1sig_syst->Sumw2();
  TH2F* hnum_trk_minus1sig_syst = new TH2F ("hnum_trk_minus1sig_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trk_minus1sig_syst->Sumw2();
  TH2F* hnum_trk_plus1sig_stat = new TH2F ("hnum_trk_plus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trk_plus1sig_stat->Sumw2();
  TH2F* hnum_trk_minus1sig_stat = new TH2F ("hnum_trk_minus1sig_stat", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_trk_minus1sig_stat->Sumw2();
  
  TH2F* hnum_tag_syst = new TH2F ("hnum_tag_syst", "Eff vs p_{T} and y; y; p_{T}; Eff", nybins, ybins, nptbins, ptbins); hnum_tag_syst->Sumw2();

  vector<TObjArray *> wHistograms;
  TH1D *curve = NULL;
  TFile* weightFile = NULL;
  TH2D* prNprWeights = NULL;
  if (caseTag.find("withPtWeights")!=std::string::npos)
    wHistograms = ReadFileWeight(isPbPb,isPr);
  else if (caseTag.find("withNonPrPtWeights")!=std::string::npos) {
    weightFile = TFile::Open(Form("weightFunctDataMC/mcWeightsPrNpr/weights_PrNpr_%s_accSample.root",isPbPb?"PbPb":"PP"));
    if (caseTag.find("useEffSample")!=std::string::npos)
      weightFile = TFile::Open(Form("weightFunctDataMC/mcWeightsPrNpr/weights_PrNpr_%s_effSample.root",isPbPb?"PbPb":"PP"));
    prNprWeights = (TH2D*) weightFile->Get("weights2D");
  }

  Long64_t nentries =fChain->GetEntries();
  //nentries=10;
  
  double normF = 1.0;
  //if (isPbPb) {
  //if (isPr) normF = 0.00276328;   
  //else normF = 0.00271722;
  //} 

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;

      if (caseTag.find("OddEvents")!=std::string::npos) {
	if (jentry%2==1) continue;
      }

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
      weight = weight*normF;
      if (caseTag.find("NoPtHatWeights")!=std::string::npos) { 
	if (isPbPb) weight = findNcoll(cent);
	else weight = 1.0;
      }
      else if (caseTag.find("NoWeights")!=std::string::npos) { 
	weight = 1.0;
      }

      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{

	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m=GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	  if (caseTag.find("absEta")!=std::string::npos) jpsi_rap = fabs(GenQQ4mom->Rapidity());
	  //if (jpsi_pt<3 || jpsi_pt>150) continue;
	  //if (fabs(jpsi_rap)>2.4) continue;

	  //if (!areGenMuonsInAcceptance2019(iQQ)) continue;

	  double ptW = 1.0;

	  if (caseTag.find("withPtWeights")!=std::string::npos) {
	    if (abs(jpsi_rap)<1.2 && jpsi_pt<6.5) weight*=1;
	    else{
	      if (abs(jpsi_rap)>=0 && abs(jpsi_rap)<0.6)        curve = (TH1D*) wHistograms[0]->At(0);
	      else if (abs(jpsi_rap)>=0.6 && abs(jpsi_rap)<1.2) curve = (TH1D*) wHistograms[1]->At(0);
	      else if (abs(jpsi_rap)>=1.2 && abs(jpsi_rap)<1.8) curve = (TH1D*) wHistograms[2]->At(0);
	      else if (abs(jpsi_rap)>=1.8 && abs(jpsi_rap)<2.4) curve = (TH1D*) wHistograms[3]->At(0);
	      ptW = curve->GetBinContent(curve->FindBin(jpsi_pt));
	    }
	  }
	  else if (caseTag.find("withNonPrPtWeights")!=std::string::npos) {
	    if (fabs(jpsi_rap)<3. && jpsi_pt<100)
	      ptW = 1.0/prNprWeights->GetBinContent(prNprWeights->FindBin(jpsi_rap, jpsi_pt));
	    if (ptW<0.0001 || ptW>1000) {
	      cout <<"ptW = "<<ptW<<" for pt = "<<jpsi_pt<<", rap = "<<jpsi_rap<<" setting it to 1"<<endl;
	      ptW =1.;
	    }
	  }
	  if (areGenMuonsInAcceptance2019(iQQ)) {
	    hdeno_pty->Fill(jpsi_rap, jpsi_pt, weight*ptW);
	    
	    if (isPbPb){
	      for (int i=0; i<ncentbins; i++) {
		if (cent<=centbins[i+1]) {
		  hdeno_pty_cent[i]->Fill(jpsi_rap, jpsi_pt, weight*ptW);
		  break;
		}
	      }
	    }
	    if (fabs(jpsi_rap)<2.4)
	      hdeno_pt->Fill(jpsi_pt, weight*ptW);
	    if (jpsi_pt>6.5) 
	      hdeno_y->Fill(jpsi_rap,weight*ptW);
	    if (isPbPb && fabs(jpsi_rap)<2.4 && jpsi_pt>6.5)
	      hdeno_cent->Fill(cent, weight*ptW);
	  }
	  
	  int whichRec = Gen_QQ_whichRec[iQQ];
	  if (whichRec < 0) continue;
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(whichRec);
	  double reco_pt = RecoQQ4mom->Pt();
	  double reco_rap = RecoQQ4mom->Rapidity();
	  if (caseTag.find("absEta")!=std::string::npos) reco_rap = fabs(RecoQQ4mom->Rapidity());

	  if (!areMuonsInAcceptance2019(whichRec)) continue;
	  if (!passQualityCuts2019(whichRec)) continue;
	  //if (!isPbPb && !isTriggerMatch(whichRec, triggerIndex_PP)) continue;
	  if (Reco_QQ_sign[whichRec]!=0) continue;
	  //if (reco_pt<3 || reco_pt>150) continue;
	  if (fabs(reco_rap)>2.4) continue;
	  if (RecoQQ4mom->M()<2.6 || RecoQQ4mom->M()>3.5) continue;

	  if (isPbPb){
	    if(!isTriggerMatch(whichRec, triggerIndex_PbPb)) continue;
	    if (!(pprimaryVertexFilter && pclusterCompatibilityFilter && phfCoincFilter2Th4)) continue;
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
	    if (!isFilterMatch(Reco_QQ_mupl_idx[whichRec],26) && !isFilterMatch(Reco_QQ_mumi_idx[whichRec],26)) cout <<"none of the muons pass the L2 filter even after the J/psi passes the trigger"<<endl;
	    if (caseTag.find("onlyL2TnpSF")!=std::string::npos) {
	      filterIdxMuPl = 26;
	      filterIdxMuMi = 26;
	    }
	    else {
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
	  }

	  tnp_weight=1.0;
	  hnum_noweights->Fill(reco_rap, reco_pt, weight*ptW);
	  double muplPt = RecoQQmupl4mom->Pt();
	  double muplY = RecoQQmupl4mom->Rapidity();
	  double mumiPt = RecoQQmumi4mom->Pt();
	  double mumiY = RecoQQmumi4mom->Rapidity();

	  if (isPbPb) {
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_nominal->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	    if (caseTag.find("noTnpWeights")!=std::string::npos) {
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(reco_pt, weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(reco_rap, weight*ptW);
	      if (fabs(jpsi_rap)<2.4 && jpsi_pt>6.5)
		hnum_cent->Fill(cent, weight*ptW);
	    }
	    else if (caseTag.find("fillWithGen")!=std::string::npos){
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(jpsi_pt, tnp_weight*weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(jpsi_rap,tnp_weight*weight*ptW);
	      if (fabs(jpsi_rap)<2.4 && jpsi_pt>6.5)
		hnum_cent->Fill(cent, tnp_weight*weight*ptW);
	    }
	    else {
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(reco_pt, tnp_weight*weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(reco_rap,tnp_weight*weight*ptW);
	      if (fabs(jpsi_rap)<2.4 && jpsi_pt>6.5)
		hnum_cent->Fill(cent, tnp_weight*weight*ptW);
	    }
	    
	    //for centrality dependance
	    if (caseTag.find("noTnpWeights")!=std::string::npos) {
	      for (int i=0; i<ncentbins; i++) {
		if (cent<=centbins[i+1]) {
		  hnum_nominal_cent[i]->Fill(reco_rap, reco_pt, weight*ptW);
		  break;
		}
	      }
	    }
	    else {
	      for (int i=0; i<ncentbins; i++) {
		if (cent<=centbins[i+1]) {
		  hnum_nominal_cent[i]->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
		  break;
		}
	      }
	    }

	    tnp_weight = tnp_weight_trk_pbpb(muplY,99)*tnp_weight_muid_pbpb(muplPt,muplY,99)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,99)
              * tnp_weight_trk_pbpb(mumiY,99)*tnp_weight_muid_pbpb(mumiPt,mumiY,99)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,99);
	    hnum_tag_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,-1)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,-1)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_plus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	      
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,-2)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,-2)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_minus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	      
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,+1)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,+1)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_plus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	      
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,+2)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,+2)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_muid_minus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	      
	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,-1)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,-1);
	    hnum_trg_plus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,-2)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,-2);
	    hnum_trg_minus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,+1)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,+1);
	    hnum_trg_plus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,0)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,+2)
	      * tnp_weight_trk_pbpb(mumiY,0)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,+2);
	    hnum_trg_minus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,-1)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,-1)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_plus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,-2)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,-2)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_minus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,+1)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,+1)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_plus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = tnp_weight_trk_pbpb(muplY,+2)*tnp_weight_muid_pbpb(muplPt,muplY,0)*tnp_weight_trg_pbpb(muplPt,muplY,filterIdxMuPl,0)
	      * tnp_weight_trk_pbpb(mumiY,+2)*tnp_weight_muid_pbpb(mumiPt,mumiY,0)*tnp_weight_trg_pbpb(mumiPt,mumiY,filterIdxMuMi,0);
	    hnum_trk_minus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	  }
	  else {
	    auto muplGlb = tnp_weight_GlobalMuon_TightAcceptance_pp(muplPt,muplY);
	    auto mumiGlb = tnp_weight_GlobalMuon_TightAcceptance_pp(mumiPt,mumiY);
	    auto muplMuIdTrg = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(muplPt,muplY);
	    auto mumiMuIdTrg = tnp_weight_HybridSoftIDTrigger_TightAcceptance_pp(mumiPt,mumiY);

	    tnp_weight = std::get<0>(muplGlb) * std::get<0>(muplMuIdTrg)
	      *std::get<0>(mumiGlb) * std::get<0>(mumiMuIdTrg);
	    hnum_nominal->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	    if (caseTag.find("noTnpWeights")!=std::string::npos) {
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(reco_pt, weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(reco_rap, weight*ptW);
	    }
	    else if (caseTag.find("fillWithGen")!=std::string::npos){
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(jpsi_pt, tnp_weight*weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(jpsi_rap,tnp_weight*weight*ptW);
	    }
	    else {
	      if (fabs(reco_rap)<2.4)
		hnum_pt->Fill(reco_pt, tnp_weight*weight*ptW);
	      if (reco_pt>6.5)
		hnum_y->Fill(reco_rap,tnp_weight*weight*ptW);
	    }
	    
	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)+std::get<2>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)+std::get<2>(mumiMuIdTrg));
	    hnum_muidtrg_plus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)-std::get<2>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)-std::get<2>(mumiMuIdTrg));
	    hnum_muidtrg_minus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)+std::get<1>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)+std::get<1>(mumiMuIdTrg));
	    hnum_muidtrg_plus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = std::get<0>(muplGlb) * (std::get<0>(muplMuIdTrg)-std::get<1>(muplMuIdTrg))
	      *std::get<0>(mumiGlb) * (std::get<0>(mumiMuIdTrg)-std::get<1>(mumiMuIdTrg));
	    hnum_muidtrg_minus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = (std::get<0>(muplGlb)+std::get<2>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)+std::get<2>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_plus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = (std::get<0>(muplGlb)-std::get<2>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)-std::get<2>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_minus1sig_syst->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = (std::get<0>(muplGlb)+std::get<1>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)+std::get<1>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_plus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);

	    tnp_weight = (std::get<0>(muplGlb)-std::get<1>(muplGlb)) * std::get<0>(muplMuIdTrg)
	      *(std::get<0>(mumiGlb)-std::get<1>(mumiGlb)) * std::get<0>(mumiMuIdTrg);
	    hnum_glb_minus1sig_stat->Fill(reco_rap, reco_pt, tnp_weight*weight*ptW);
	  }
	}
    }

  gSystem->mkdir(Form("FilesAccxEff%s",caseTag.c_str()));  
  gSystem->mkdir(Form("FilesAccxEff%s/Eff",caseTag.c_str()));
  TFile* fsave = new TFile (Form("FilesAccxEff%s/Eff/%sEffHists_%s.root", caseTag.c_str(), isPr?"pr":"npr",isPbPb?"PbPb":"PP"), "RECREATE");
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
    hnum_tag_syst->Write();

    for (int i=0; i<ncentbins; i++) {
      hdeno_pty_cent[i]->Write();
      hnum_nominal_cent[i]->Write();
    }    
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
  
  if (isPbPb) {
    hdeno_cent->Write();
    hnum_cent->Write();
  }
  
  fsave->Close();
  
  //delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_noweights; delete hnum_nominal; delete hnum_binned; delete hnum_plus1sig; delete hnum_minus1sig; delete hnum_muid_sta; delete hnum_muid; delete hnum_muid_plus1sig; delete hnum_muid_minus1sig; delete hnum_sta; delete hnum_sta_plus1sig; delete hnum_sta_minus1sig; delete hnum_trk_plus1sig; delete hnum_trk_minus1sig;
}//end of EffStep()


void oniaTree::AccCalc (string caseTag) {

  setBins(caseTag);

  if (!isAcc) {cout<<"[ERROR] you're trying to make Acceptance with Efficiency trees."<<endl; return;}
  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;

  TH1F* hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbinsAna, ptbinsAna); hdeno_pt->Sumw2();
  TH1F* hnum_pt = new TH1F ("hnum_pt", "Acc vs p_{T}; p_{T}; Acc", nptbinsAna, ptbinsAna); hnum_pt->Sumw2();
  if (caseTag.find("1DfinerBins")!=std::string::npos) {
    hdeno_pt = new TH1F ("hdeno_pt", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins); hdeno_pt->Sumw2();
    hnum_pt = new TH1F ("hnum_pt", "Acc vs p_{T}; p_{T}; Acc", nptbins, ptbins); hnum_pt->Sumw2();
  }

  TH1F* hdeno_y = new TH1F ("hdeno_y","N_{gen} vs y; y; N_{total}", nybins, ybins); hdeno_y->Sumw2();
  TH1F* hnum_y = new TH1F ("hnum_y","Acc vs y; y; Acc", nybins, ybins); hnum_y->Sumw2();

  TH2F* hdeno_pty = new TH2F ("hdeno_pty", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", nybins, ybins, nptbins, ptbins); hdeno_pty->Sumw2();
  TH2F* hnum_nominal = new TH2F ("hnum_nominal", "Acc vs p_{T} and y; y; p_{T}; Acc", nybins, ybins, nptbins, ptbins); hnum_nominal->Sumw2();
  
  vector<TObjArray *> wHistograms;
  TH1D *curve = NULL;
  TFile* weightFile = NULL;
  TH2D* prNprWeights = NULL;
  if (caseTag.find("withPtWeights")!=std::string::npos)
    wHistograms = ReadFileWeight(isPbPb,isPr);
  else if (caseTag.find("withNonPrPtWeights")!=std::string::npos) {
    weightFile = TFile::Open(Form("weightFunctDataMC/mcWeightsPrNpr/weights_PrNpr_%s_accSample.root",isPbPb?"PbPb":"PP"));
    if (caseTag.find("useEffSample")!=std::string::npos)
      weightFile = TFile::Open(Form("weightFunctDataMC/mcWeightsPrNpr/weights_PrNpr_%s_effSample.root",isPbPb?"PbPb":"PP"));
    prNprWeights = (TH2D*) weightFile->Get("weights2D");
  }

  Long64_t nentries = fChain->GetEntries();
  //nentries = 2000000;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;

      if (caseTag.find("OddEvents")!=std::string::npos) {
	if (jentry%2==1) continue;
      }

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m=GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	  if (caseTag.find("absEta")!=std::string::npos) jpsi_rap = fabs(GenQQ4mom->Rapidity());

	  if (jpsi_pt<3 || jpsi_pt>150) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  weight = 1;
	  if (caseTag.find("withPtWeights")!=std::string::npos) {
	    if (abs(jpsi_rap)<1.2 && jpsi_pt<6.5) weight = 1;
	    else {
	      if (abs(jpsi_rap)>=0 && abs(jpsi_rap)<0.6)        curve = (TH1D*) wHistograms[0]->At(0);
	      else if (abs(jpsi_rap)>=0.6 && abs(jpsi_rap)<1.2) curve = (TH1D*) wHistograms[1]->At(0);
	      else if (abs(jpsi_rap)>=1.2 && abs(jpsi_rap)<1.8) curve = (TH1D*) wHistograms[2]->At(0);
	      else if (abs(jpsi_rap)>=1.8 && abs(jpsi_rap)<2.4) curve = (TH1D*) wHistograms[3]->At(0);
	      weight = curve->GetBinContent(curve->FindBin(jpsi_pt));
	    }
	  }
	  else if (caseTag.find("withNonPrPtWeights")!=std::string::npos) {
	    weight = 1.0/prNprWeights->GetBinContent(prNprWeights->FindBin(jpsi_rap, jpsi_pt));
	  }

	  hdeno_pty->Fill(jpsi_rap, jpsi_pt, weight);
	  
	  hdeno_pt->Fill(jpsi_pt,weight);
	  if (jpsi_pt>6.5)
	    hdeno_y->Fill(jpsi_rap,weight);
	  
	  //if (!areGenMuonsInAcceptance2019(iQQ)) continue;
	  TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	  TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	  ////////muons in acc 
	  if (!isGlobalMuonInAccept2019(GenQQmupl4mom)) continue;
	  if (!isGlobalMuonInAccept2019(GenQQmumi4mom)) continue;
	  
	  hnum_nominal->Fill(jpsi_rap, jpsi_pt,weight);
	  
	  //if (jpsi_pt<6.5) continue;
	  hnum_pt->Fill(jpsi_pt,weight);
	  if(jpsi_pt>=6.5)
	    hnum_y->Fill(jpsi_rap,weight);
	}//end of genQQ loop
    }//end of entries loop
  
  gSystem->mkdir(Form("FilesAccxEff%s",caseTag.c_str()));
  gSystem->mkdir(Form("FilesAccxEff%s/Acc",caseTag.c_str()));
  TFile* fsave = new TFile (Form("FilesAccxEff%s/Acc/%sAccHists_%s.root", caseTag.c_str(), isPr?"pr":"npr", isPbPb?"PbPb":"PP"), "RECREATE");
  hdeno_pty->Write("hdeno_2d");
  hnum_nominal->Write("hnum_2d_nominal");
  hdeno_pt->Write("hdeno_pt");
  hnum_pt->Write("hnum_pt");
  hdeno_y->Write("hdeno_y");
  hnum_y->Write("hnum_y");
  fsave->Close();
  
  delete fsave; delete hdeno_pt; delete hnum_pt; delete hdeno_y; delete hnum_y; delete hdeno_pty; delete hnum_nominal;
}//end of AccCalc function

void oniaTree::ClosureTest(string caseTag, bool onlyPlot) {
  setBins(caseTag);

  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;
  gStyle->SetOptStat(0);

  TH1D* ptTotal = new TH1D("ptTotal", "pt distribution for total events", nptbinsAna, ptbinsAna);
  TH1D* ptPass = new TH1D("ptPass", "pt distribution for passing events", nptbinsAna, ptbinsAna);

  TFile* corrfile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str()),"READ");

  TEfficiency* eff = (TEfficiency*) corrfile->Get(Form("hcorr_Jpsi_%s_pr_Eff",isPbPb?"PbPb":"PP"));
  TEfficiency* acc = (TEfficiency*) corrfile->Get("hcorr_Jpsi_PP_pr_Acc");

  TList* lcorr = corrfile->GetListOfKeys();
  TIter nextCorr(lcorr);

  TObjArray* corrCent = new TObjArray();
  corrCent->SetOwner(kTRUE);

  TObjString* fname(0x0);
  while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
    {
      TEfficiency* h = static_cast<TEfficiency*>(corrfile->FindObjectAny(fname->GetString().Data()));

      TString sName(h->GetName());
      if ( sName.Contains("cent") ){
	corrCent->Add(h->Clone());
      }
      cout <<"[INFO] adding "<<sName<<" to the corrction array"<<endl;
    }

  int cent=0;
  double normF = 1.0;
  if (!onlyPlot) {
  Long64_t nentries = fChain->GetEntries();
  //nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] processing entry "<<jentry<<"/"<<nentries<<endl;

      if (caseTag.find("OddEvents")!=std::string::npos) {
	if (jentry%2==0) continue;
      }

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      if(isAcc) Gen_weight = 1;
      else if (!isPbPb && isPr) {
	Gen_weight = 1;
	if (pthat >= 15 && pthat < 25)  Gen_weight = 0.0247699;
	else if (pthat >= 25 && pthat < 35) Gen_weight =  0.00311931;
	else if (pthat >= 35 && pthat < 45) Gen_weight = 0.000693027;
	else if (pthat >= 45) Gen_weight = 0.000212618;
      }
      else if (isPbPb) {
	cent = getMCHiBinFromhiHF(hiHF);
	Gen_weight = Gen_weight*findNcoll(cent);
      }
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  
	  jpsi_m = GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	  if (caseTag.find("absEta")!=std::string::npos) jpsi_rap = fabs(GenQQ4mom->Rapidity());

	  if (jpsi_pt<6.5 || jpsi_pt>100) continue;
	  if (fabs(jpsi_rap)>2.4) continue; 
	  if (isAcc) ptTotal->Fill(jpsi_pt,Gen_weight);
	  //apply the acceptance cuts

	  if (isAcc) {
	    TLorentzVector *GenQQmupl4mom = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	    TLorentzVector *GenQQmumi4mom = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	    if (!isGlobalMuonInAccept2019(GenQQmupl4mom)) continue;
	    if (!isGlobalMuonInAccept2019(GenQQmumi4mom)) continue;
	  }

	  else {
	    if (!areGenMuonsInAcceptance2019(iQQ)) continue;
	  }
	  weight = 1.0;
	  weight = 1.0/acc->GetEfficiency(acc->FindFixBin(jpsi_rap,jpsi_pt));
	  if (isAcc) ptPass->Fill(jpsi_pt,Gen_weight*weight);
	  else ptTotal->Fill(jpsi_pt,Gen_weight);
	}

      if (isAcc) continue;
      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
	{
	  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	  jpsi_pt = RecoQQ4mom->Pt();
	  jpsi_rap = RecoQQ4mom->Rapidity();
	  if (caseTag.find("absEta")!=std::string::npos) jpsi_rap = fabs(RecoQQ4mom->Rapidity());
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
	  weight = 1.0/eff->GetEfficiency(eff->FindFixBin(jpsi_rap,jpsi_pt));
	  if (isPbPb) {
	    if (caseTag.find("centBins")!=std::string::npos) {
	      string centTag="";
	      for (int i=0; i<ncentbins; i++) {
		if (cent<=centbins[i+1]) {
		  centTag = Form("cent%d%d",centbins[i],centbins[i+1]);
		  break;
		}
	      }
	      TEfficiency* effTemp = static_cast<TEfficiency*>(corrCent->FindObject(Form("hcorr_Jpsi_PbPb_pr_Eff_%s",centTag.c_str())));
	      weight = 1.0/effTemp->GetEfficiency(effTemp->FindFixBin(jpsi_rap,jpsi_pt));
	    }
	  }
	  if (weight<0.00001 || weight>1000000)
	    weight=1; 
	  ptPass->Fill(jpsi_pt,Gen_weight*weight);
	}
    }
  }

  TFile* saveFile = NULL;
  if (onlyPlot) {
    saveFile = TFile::Open(Form("FilesAccxEff%s/ClosureTest/ClosureTest_StatisticallyIdependentSamples_%s_%s_%s.root",caseTag.c_str(),isPbPb?"PbPb":"PP",isPr?"pr":"npr",isAcc?"acc":"eff"));
    ptTotal = (TH1D*) saveFile->Get("ptTotal");
    ptPass = (TH1D*) saveFile->Get("ptPass");
  }
  else {
    gSystem->mkdir(Form("FilesAccxEff%s/ClosureTest",caseTag.c_str()));
    saveFile = new TFile (Form("FilesAccxEff%s/ClosureTest/ClosureTest_StatisticallyIdependentSamples_%s_%s_%s.root",caseTag.c_str(),isPbPb?"PbPb":"PP",isPr?"pr":"npr",isAcc?"acc":"eff"),"RECREATE");
    ptTotal->Write();
    ptPass->Write();
    saveFile->Close();
  }
  
  ptTotal->Scale(1.,"width");
  ptPass->Scale(1.,"width");
  
  ptTotal->SetLineColor(kRed);
  ptTotal->SetMarkerColor(kRed);
  ptTotal->SetMarkerStyle(kFullCircle);
  
  ptPass->SetLineColor(kBlue);
  ptPass->SetMarkerColor(kBlue);
  ptPass->SetMarkerStyle(kOpenCircle);

  TCanvas* c = new TCanvas ("c","",900,1000);
  c->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.01);
  pad1->Draw();
  
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetBottomMargin(0.2);
  pad2->SetTopMargin(0.01);
  pad2->Draw();

  TLegend* leg = new TLegend(0.59, 0.6, 0.89, 0.8);
  leg->SetBorderSize(0);
  
  TLine * ly1 = new TLine(6.5, 1, 50, 1);
  ly1->SetLineColor(kRed);
  ly1->SetLineStyle(2);
  ly1->SetLineWidth(2);
 
  TH1D* pull=NULL;

  leg->AddEntry(ptTotal,"total events","lp");
  leg->AddEntry(ptPass,"passing with corrections","lp");
  
  if (ptPass->GetMaximum()>ptTotal->GetMaximum()) ptTotal->GetYaxis()->SetRangeUser(0,ptPass->GetMaximum()*1.05);

  pad1->cd();
  ptTotal->Draw();
  ptPass->Draw("same");
  leg->Draw("same");
  //pad1->SetLogy();
  pad2->cd();
  pull = makePull(ptPass,ptTotal);
  pull->GetYaxis()->SetTitle("pass/total");
  pull->GetYaxis()->SetRangeUser(0.89,1.11);
  pull->Draw();
  ly1->Draw("same");
  pull->Draw("same");
  gSystem->mkdir("Utilities/CorrBinOptimization");
  c->SaveAs(Form("Utilities/CorrBinOptimization/ClosureTest%s_StatisticallyIdependentSamples_%s_%s_%s.pdf",caseTag.c_str(),isPbPb?"PbPb":"PP",isPr?"pr":"npr",isAcc?"acc":"eff"));
}


void oniaTree::ClosureTestPtWeights()
{
  int nptbinsAna = sizeof(ptbinsAna)/sizeof(double)-1;
  gStyle->SetOptStat(0);
  Double_t etabins []={0, 1.6, 2.4};
  TH1F* hist = new TH1F ("hist", "pt distribution at reco level without pt weights", nptbinsAna, ptbinsAna);
  TH1F* hist_ptw = new TH1F ("hist_ptw", "pt distribution at reco level with pt weights", nptbinsAna, ptbinsAna);
  TH1F* genHist = new TH1F ("genHist", "pt distribution at gen level", nptbinsAna, ptbinsAna);
  TH1F* genHist_ptw = new TH1F ("genHist_ptw", "pt distribution at gen level", nptbinsAna, ptbinsAna);

  TFile* corrfile = TFile::Open("../Fitter/Input/correction_AccEff_centMaps.root","READ");
  TFile* corrfile_ptw = TFile::Open("../Fitter/Input/correction_AccEff_centMaps_withptWeights.root","READ");

  TEfficiency* eff = (TEfficiency*) corrfile->Get(Form("hcorr_Jpsi_%s_pr",isPbPb?"PbPb":"PP"));
  TEfficiency* eff_ptw = (TEfficiency*) corrfile_ptw->Get(Form("hcorr_Jpsi_%s_pr",isPbPb?"PbPb":"PP"));

  //get the acceptance files
  
  TFile*accFile = TFile::Open(Form("FilesAccxEff/Acc/%sAccHists_PP.root",isPr?"pr":"npr"));
  TH2F *accNum = (TH2F*) accFile->Get("hnum_2d_nominal");
  TH2F *accDen = (TH2F*) accFile->Get("hdeno_2d");
  accNum->Divide(accDen);

  TFile*accFile_ptw = TFile::Open(Form("FilesAccxEff_withptWeights/Acc/%sAccHists_PP.root",isPr?"pr":"npr"));
  TH2F *accNum_ptw = (TH2F*) accFile_ptw->Get("hnum_2d_nominal");
  TH2F *accDen_ptw = (TH2F*) accFile_ptw->Get("hdeno_2d");
  accNum_ptw->Divide(accDen_ptw);

  Long64_t nentries = fChain->GetEntries();
  //nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;

  double normF = 1.0;
  if (isPbPb) {
    if (isPr) normF = 0.00276328;   
    else normF = 0.00271722;
  }   
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] processing entry "<<jentry<<"/"<<nentries<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (!isPbPb && isPr) {
	if (pthat >= 15 && pthat < 25)  Gen_weight = 0.0247699;
	else if (pthat >= 25 && pthat < 35) Gen_weight =  0.00311931;
	else if (pthat >= 35 && pthat < 45) Gen_weight = 0.000693027;
	else if (pthat >= 45) Gen_weight = 0.000212618;
      }
      else if (isPbPb) Gen_weight = Gen_weight*findNcoll(cent)*normF;
      else Gen_weight=1;
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  //cout <<"Gen_weight = "<<Gen_weight<<endl;

	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  jpsi_m = GenQQ4mom->M();
	  jpsi_pt = GenQQ4mom->Pt();
	  jpsi_rap = GenQQ4mom->Rapidity();
	  if (jpsi_pt<6.5 || jpsi_pt>50) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  //apply the acceptance cuts
	  if (!areGenMuonsInAcceptance2019(iQQ)) continue;
	  if (accNum->GetBinContent(accNum->FindBin(jpsi_rap,jpsi_pt))<0.0001) weight = Gen_weight;
	  else 
	    weight = Gen_weight*1.0/accNum->GetBinContent(accNum->FindBin(jpsi_rap,jpsi_pt)); 
	  //if (weight>10000) cout <<"acc weight without ptw = "<<weight<<endl; 
	  genHist->Fill(jpsi_pt,weight);
	  if (accNum_ptw->GetBinContent(accNum_ptw->FindBin(jpsi_rap,jpsi_pt))<0.0001) weight = Gen_weight;
	  else 
	    weight = Gen_weight*1.0/accNum_ptw->GetBinContent(accNum_ptw->FindBin(jpsi_rap,jpsi_pt)); 
	  //if (weight>10000) cout <<"acc weight  with ptw = "<<weight<<endl; 
	  genHist_ptw->Fill(jpsi_pt,weight);
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
	  if (jpsi_pt<6.5 || jpsi_pt>50) continue;
	  if (fabs(jpsi_rap)>2.4) continue;
	  if (jpsi_m<2.6 || jpsi_m>3.5) continue;
	  if (isPbPb && !(pprimaryVertexFilter && pBeamScrapingFilter && phfCoincFilter2Th4)) continue;
	  if (!isPbPb && !(pPAprimaryVertexFilter && pBeamScrapingFilter)) continue;
	  if (Reco_QQ_whichGen[iQQ] < 0) continue;
	  if (!areGenMuonsInAcceptance2019(Reco_QQ_whichGen[iQQ])) continue;

	  weight = Gen_weight*1.0/eff->GetEfficiency(eff->FindFixBin(jpsi_rap,jpsi_pt));
	  hist->Fill(jpsi_pt,weight);
	  weight = Gen_weight*1.0/eff_ptw->GetEfficiency(eff_ptw->FindFixBin(jpsi_rap,jpsi_pt));
	  hist_ptw->Fill(jpsi_pt,weight);
	}
    }
  gSystem->mkdir("FilesAccxEff");
  gSystem->mkdir("FilesAccxEff/ClosureTest");
  TFile* testfile = new TFile (Form("FilesAccxEff/ClosureTest/ClosureTest_ptWeights_%s_%s.root",isPbPb?"PbPb":"PP",isPr?"pr":"npr"),"RECREATE");
  genHist->Write();
  genHist_ptw->Write();
  hist->Write();
  hist_ptw->Write();

  testfile->Close();
  

 
  genHist->SetLineColor(kRed);
  genHist_ptw->SetLineColor(kGreen+2);
  hist->SetLineColor(kBlue);
  hist_ptw->SetLineColor(kMagenta+2);
  genHist->SetLineWidth(2);
  genHist_ptw->SetLineWidth(2);
  hist->SetLineWidth(2);
  hist_ptw->SetLineWidth(2);

  genHist->SetTitle(Form("%s Comparison", isPr?"prompt":"nonprompt"));

  TCanvas* c = new TCanvas("c","", 1000, 1000);

  genHist->Draw();
  genHist_ptw->Draw("same");
  hist->Draw("same");
  hist_ptw->Draw("same");
  TLegend* leg = new TLegend (0.6, 0.6, 0.75, 0.75);
  leg->AddEntry(genHist, "N_{gen}","lep");
  leg->AddEntry(genHist_ptw, "N_{gen} with pt weights in Acc","lep");
  leg->AddEntry(hist, "N_{reco} ","lep");
  leg->AddEntry(hist_ptw, "N_{reco} with pt weights in AccxEff","lep");
  leg->SetBorderSize(0);
  leg->Draw("same");
  c->SaveAs(Form("FilesAccxEff/ClosureTest/GenRecoComparison_ptWeights_%s_%s.pdf",isPbPb?"PbPb":"PP",isPr?"pr":"npr"));
}

void oniaTree::Plot(string caseTag) {
  setBins(caseTag);
  gStyle->SetOptStat(false);

  cout<<"[INFO] Importing the numerators and denominators of the corrections."<<endl;
  TFile*prAccFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
  if (!prAccFile_pbpb) {
    cout<<"[ERROR] pbpb prompt acceptance file not found!"<<endl;
    if (isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc(caseTag);
      prAccFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  TFile*prEffFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PbPb.root",caseTag.c_str()));
  if (!prEffFile_pbpb) {
    cout<<"[ERROR] pbpb prompt efficiency file not found!"<<endl;
    if (!isAcc && isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc(caseTag);
      prEffFile_pbpb = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PbPb.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  ////////////////////////////////////pp/////////////////////////////////
  TFile*prAccFile_pp = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
  if (!prAccFile_pp) {
    cout<<"[ERROR] pp prompt acceptance file not found!"<<endl;
    if (isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      AccCalc(caseTag);
      prAccFile_pp = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_PP.root",caseTag.c_str()));
    }
    else {
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
    }
  }
  TFile*prEffFile_pp = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PP.root",caseTag.c_str()));
  if (!prEffFile_pp) {
    cout<<"[ERROR] pp prompt efficiency file not found!"<<endl;
    if (!isAcc && !isPbPb && isPr) {cout<<"[INFO] since the settings are good I will make it"<<endl;
      EffCalc(caseTag);
      prEffFile_pp = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_PP.root",caseTag.c_str()));
    }
    else 
      cout<<"[ERROR] Please change your settings and retry."<<endl; return;
  }

  TH2F *prAccNum_pbpb = (TH2F*) prAccFile_pbpb->Get("hnum_2d_nominal");
  TH2F *prAccDen_pbpb = (TH2F*) prAccFile_pbpb->Get("hdeno_2d");
  TH2F *prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_nominal");
  if (caseTag.find("noTnpWeights")!=std::string::npos)
    prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get("hnum_noweights");
  TH2F *prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get("hdeno_pty");

  TH2F *prAccNum_pp = (TH2F*) prAccFile_pp->Get("hnum_2d_nominal");
  TH2F *prAccDen_pp = (TH2F*) prAccFile_pp->Get("hdeno_2d");
  TH2F *prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_nominal");
  if (caseTag.find("noTnpWeights")!=std::string::npos)
    prEffNum_pp = (TH2F*) prEffFile_pp->Get("hnum_noweights");
  TH2F *prEffDen_pp = (TH2F*) prEffFile_pp->Get("hdeno_pty");

  TCanvas* c = new TCanvas("c","",1000,900);
  prAccNum_pbpb->Divide(prAccDen_pbpb);
  prEffNum_pbpb->Divide(prEffDen_pbpb);
  prAccNum_pp->Divide(prAccDen_pp);
  prEffNum_pp->Divide(prEffDen_pp);
  
  prAccNum_pbpb->SetTitle("Acceptance");
  prEffNum_pbpb->SetTitle("PbPb Efficiency");
  prAccNum_pp->SetTitle("Acceptance");
  prEffNum_pp->SetTitle("pp Efficiency");

  prAccNum_pbpb->GetZaxis()->SetTitle("");
  prEffNum_pbpb->GetZaxis()->SetTitle("");
  prAccNum_pp->GetZaxis()->SetTitle("");
  prEffNum_pp->GetZaxis()->SetTitle("");

  prAccNum_pbpb->GetYaxis()->SetRangeUser(6.5,100);
  prEffNum_pbpb->GetYaxis()->SetRangeUser(6.5,100);
  prAccNum_pp->GetYaxis()->SetRangeUser(6.5,100);
  prEffNum_pp->GetYaxis()->SetRangeUser(6.5,100);

  prAccNum_pbpb->GetZaxis()->SetRangeUser(0,1);
  prEffNum_pbpb->GetZaxis()->SetRangeUser(0,1);
  prAccNum_pp->GetZaxis()->SetRangeUser(0,1);
  prEffNum_pp->GetZaxis()->SetRangeUser(0,1);

  c->cd();
  c->SetLogy();
  prAccNum_pbpb->Draw("colz");
  c->SaveAs("Utilities/pbpbAcc2D.pdf");
  prEffNum_pbpb->Draw("colz");
  c->SaveAs("Utilities/pbpbEff2D.pdf");
  prAccNum_pp->Draw("colz");
  c->SaveAs("Utilities/ppAcc2D.pdf");
  prEffNum_pp->Draw("colz");
  c->SaveAs("Utilities/ppEff2D.pdf");

  if (caseTag.find("centBins")!=std::string::npos) {
    for (int i=0; i<ncentbins; i++) {
      prEffNum_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hnum_nominal_cent%d%d",centbins[i],centbins[i+1]));
      prEffDen_pbpb = (TH2F*) prEffFile_pbpb->Get(Form("hdeno_pty_cent%d%d",centbins[i],centbins[i+1]));
      prEffNum_pbpb->Divide(prEffDen_pbpb);
      /*
      int nBinX=prEffNum_pbpb->GetNbinsX();
      int nBinY=prEffNum_pbpb->GetNbinsY();
      for (int iX=0; iX<nBinX; iX++) {
	for (int iY=0; iY<nBinY; iY++) {
	  int iBin = prEffNum_pbpb->GetBin(iX,iY);
	  if (prEffNum_pbpb->GetBinContent(iBin)< 0.0001) prEffNum_pbpb->SetBinContent(iBin, 0.0001);
	  if (prEffNum_pbpb->GetBinError(iBin)< 0.001) prEffNum_pbpb->SetBinError(iBin, 0.001);
	}
      }
      */
      prEffNum_pbpb->SetTitle(Form("PbPb Efficiency hiBin %d-%d",centbins[i],centbins[i+1]));
      prEffNum_pbpb->GetZaxis()->SetTitle("");
      prEffNum_pbpb->GetYaxis()->SetRangeUser(6.5,100);
      prEffNum_pbpb->GetZaxis()->SetRangeUser(0,1);
      prEffNum_pbpb->Draw("colz");
      c->SaveAs(Form("Utilities/pbpbEff2D_cent%d%d.pdf",centbins[i],centbins[i+1]));
    }
  }
}
      
vector<TObjArray*> oniaTree::ReadFileWeight(bool ispbpb, bool isprompt) {
  string wfilePbPb_prompt[] = {"weights_JPsi_PbPb_006_prompt.root","weights_JPsi_PbPb_0612_prompt.root","weights_JPsi_PbPb_1218_prompt.root","weights_JPsi_PbPb_1824_prompt.root"};
  string wfilePP_prompt[] = {"weights_JPsi_PP_006_prompt.root","weights_JPsi_PP_0612_prompt.root","weights_JPsi_PP_1218_prompt.root","weights_JPsi_PP_1824_prompt.root"};
  string wfilePbPb_nonprompt[] = {"weights_JPsi_PbPb_006_nonprompt.root","weights_JPsi_PbPb_0612_nonprompt.root","weights_JPsi_PbPb_1218_nonprompt.root","weights_JPsi_PbPb_1824_nonprompt.root"};
  string wfilePP_nonprompt[] = {"weights_JPsi_PP_006_nonprompt.root","weights_JPsi_PP_0612_nonprompt.root","weights_JPsi_PP_1218_nonprompt.root","weights_JPsi_PP_1824_nonprompt.root"};
	
  const int nidxf = sizeof(wfilePbPb_prompt)/sizeof(string);
	
  TFile *fweight[nidxf];
  vector<TObjArray*> objarr;
   
  cout << "ReadFileWeight: " << ispbpb << " " << isprompt << endl;
  for (int idxf=0; idxf<nidxf; idxf++) {
    if (ispbpb && isprompt)
      fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePbPb_prompt[idxf].c_str()));
    else if (ispbpb && !isprompt)
      fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePbPb_nonprompt[idxf].c_str()));
    else if (!ispbpb && isprompt)
      fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePP_prompt[idxf].c_str()));
    else if (!ispbpb && !isprompt)
      fweight[idxf] = new TFile(Form("weightFunctDataMC/%s",wfilePP_nonprompt[idxf].c_str()));
    else
      cout <<"Cannot load files " << endl;
     
    cout << "weighting file : " << fweight[idxf]->GetName() << endl;
    TObjArray *objtmp = (TObjArray*)fweight[idxf]->Get("DataOverMC");
    TObjArray *obj = (TObjArray*)objtmp->Clone(Form("%s_copy",objtmp->GetName()));
    objarr.push_back(obj);
  }

  for (int idxf=0; idxf<nidxf; idxf++) {
    fweight[idxf]->Close();
  }
  return objarr;
}  

