#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"
// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

#include <vector> 
//#ifdef __MAKECINT__ 
//#pragma link C++ class vector<vector<float> >+; 
//#endif
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TArrow.h"
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"
#include <TString.h>
#include <fstream>

using namespace std;

Double_t ptbins2D []= {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50};
Double_t ybins2D []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

void GetMCResults(bool isPbPb=false,bool isPr=true, bool underflowOff=true, bool plotBpt=false, bool plotxi=false, bool cuts18012 =false);
Bool_t isGlobalMuonInAccept2019 (TLorentzVector* Muon);

void GetMCResults_all () {
  //GetMCResults(false, true,  false, false, false);
  GetMCResults(false, false,  false, false, false, true);
  GetMCResults(false, false,  false, false, false, false);
}

void GetMCResults(bool isPbPb, bool isPr, bool underflowOff, bool plotBpt, bool plotXi, bool cuts18012) {
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;

  double sefer=0;

  TFile *treeFile = NULL;
  if (isPr)
    treeFile = TFile::Open("/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_ext_merged.root");
  else 
    treeFile = TFile::Open("/data_CMS/cms/mnguyen/jPsiJet18/mc/pp/nonprompt/genR4fix/merged_HiForestAOD.root");

  TTree* oniaTree = (TTree*) treeFile->Get("hionia/myTree");
  TTree* jetTree = NULL;
  if (cuts18012)
    jetTree = (TTree*) treeFile->Get("ak4PFJetAnalyzer/t");
  else 
    jetTree = (TTree*) treeFile->Get("ak3PFJetAnalyzer/t");

  TTree* bTree(0x0);
  bTree = (TTree*) treeFile->Get("bHadronAna/hi");
  if (!bTree) cout<<"b tree not found"<<endl;
  oniaTree->AddFriend(jetTree);
  if (plotBpt)  oniaTree->AddFriend(bTree);

  TH1F* zMidHist = new TH1F ("zMidHist","", 7, 0.02, 1); zMidHist->Sumw2();
  TH1F* zHist = new TH1F ("zHist","", 6, 0.064, 1); zHist->Sumw2();
  TH1F* zFwdHist = new TH1F ("zFwdHist","", 5, 0, 1); zFwdHist->Sumw2();
  TH1F* nMidHist = new TH1F ("nMidHist","", 10, 0, 1); nMidHist->Sumw2();
  TH1F* nFwdHist = new TH1F ("nFwdHist","", 10, 0, 1); nFwdHist->Sumw2();

  cout <<"[INFO] Importing AccFiles to correct"<<endl;
  TFile *corrFile = TFile::Open(Form("../Efficiency/FilesAccxEff_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights/Acc/%sAccHists_%s.root",isPr?"pr":"npr","PP"));
  TH2F* corrNum = (TH2F*) corrFile->Get("hnum_2d_nominal");
  TH2F* corrDeno =(TH2F*) corrFile->Get("hdeno_2d");

  TEfficiency* accCorr = new TEfficiency("accCorr", "Acc(y,pt); y; pt; eff", ny2D, ybins2D, npt2D, ptbins2D);
  accCorr->SetStatisticOption(TEfficiency::kBBayesian);
  accCorr->SetPassedHistogram(*corrNum,"f");
  accCorr->SetTotalHistogram(*corrDeno,"f");

  double accWeight;
  double zed;
  double drmin;

  TTree* fChain;

  Int_t           Gen_QQ_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  TClonesArray    *Gen_mu_4mom;
  Int_t           Gen_QQ_mupl_idx[99];
  Int_t           Gen_QQ_mumi_idx[99];
  Float_t         Gen_weight;

  Int_t           ngen;
  Float_t         genpt[99];   //[ngen]
  Float_t         geneta[99];   //[ngen]
  Float_t         geny[99];   //[ngen]
  Float_t         genphi[99];   //[ngen]
  Float_t         genm[99];   //[ngen]



  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_mupl_4mom;   //!
  TBranch        *b_Gen_QQ_mumi_4mom;   //!
  TBranch        *b_Gen_mu_4mom;
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  TBranch        *b_Gen_weight;
  TBranch        *b_ngen;   //!
  TBranch        *b_genpt;   //!
  TBranch        *b_geneta;   //!
  TBranch        *b_geny;   //!
  TBranch        *b_genphi;   //!
  TBranch        *b_genm;   //!

  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;
  Gen_mu_4mom = 0;


  if (!oniaTree) { cout<<"[ERROR] no tree found"<<endl; return;}

  fChain = oniaTree;
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
  if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
  if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
  if (fChain->GetBranch("Gen_weight")) fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
  
  if (fChain->GetBranch("ngen")) fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
  if (fChain->GetBranch("genpt")) fChain->SetBranchAddress("genpt", genpt, &b_genpt);
  if (fChain->GetBranch("geneta")) fChain->SetBranchAddress("geneta", geneta, &b_geneta);
  if (fChain->GetBranch("geny")) fChain->SetBranchAddress("geny", geny, &b_geny);
  if (fChain->GetBranch("genphi")) fChain->SetBranchAddress("genphi", genphi, &b_genphi);
  if (fChain->GetBranch("genm")) fChain->SetBranchAddress("genm", genm, &b_genm);



  Long64_t nentries =fChain->GetEntries();
  //nentries = 1000000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<Form(" for (%s, %s, %s(%s))",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotXi?"xi":"z",plotBpt?"B":"Jpsi")<<endl;
      nb = fChain->GetEntry(jentry);   
      nbytes += fChain->GetEntry(jentry);
      fChain->GetEntry(jentry);
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++) {
	  zed = -1;
	  drmin = 0.5;
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iQQ]);
	  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iQQ]);

	  if (GenQQ4mom->M()<2.6 || GenQQ4mom->M()>3.5) continue;
	  if (cuts18012) {if (GenQQ4mom->Pt() < 6.5 || GenQQ4mom->Pt() > 35) continue;}
	  else {if (GenQQ4mom->Pt() < 6.5 || GenQQ4mom->Pt() > 100) continue;}
	  if (abs(GenQQ4mom->Rapidity())>2.4) continue;
	  if (!isGlobalMuonInAccept2019(GenQQmupl) || !isGlobalMuonInAccept2019(GenQQmumi)) continue;

	  accWeight = 1.0/(accCorr->GetEfficiency(accCorr->FindFixBin(GenQQ4mom->Rapidity(), GenQQ4mom->Pt())));
	  if (isPr) Gen_weight=1;
	  accWeight = accWeight * Gen_weight;

	  float jPsiPt = GenQQ4mom->Pt();

	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    nFwdHist->Fill(sefer, accWeight);
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt() > 6.5)
	    nMidHist->Fill(sefer, accWeight);

	  for (int iJet=0; iJet<ngen; iJet++) {
	    if (cuts18012) {
	      if (genpt[iJet]<25 || genpt[iJet]>35) continue;
              if (abs(geny[iJet])>2.4) continue;
	    }
	    else {
	      if (genpt[iJet]<30 || genpt[iJet]>40) continue;
	      if (abs(geny[iJet])>2.) continue;
	    }
	    TLorentzVector v_jet;
	    v_jet.SetPtEtaPhiM(genpt[iJet], geneta[iJet], genphi[iJet], genm[iJet]);
	    if (GenQQ4mom->DeltaR (v_jet)<=drmin) {
	      drmin = GenQQ4mom->DeltaR (v_jet);
	      zed = GenQQ4mom->Pt()*1.0/genpt[iJet];
	    }	      
	  }// end of gen jet loop
	  if (zed == -1) continue;
	  if (zed > 1 && zed <= 1.000001) zed = 0.9999999;
	  if (zed > 1) continue;
	  if (plotXi) zed = log(1.0/zed);
	  
	  if (underflowOff  && zed<0.2) continue;
	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    zFwdHist->Fill(zed, accWeight);
	  if (underflowOff  && zed<0.44) continue;
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt()>6.5)
	    zMidHist->Fill(zed, accWeight);
	  zHist->Fill(zed, accWeight);
	} //end of genQQ loop 
    }//end of events loop
  gSystem->mkdir("Output/MCResults");
  TFile *fsave = new TFile(Form("Output/MCResults/mcResult_%s_%s_%s%s%s%s.root","PP",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":"", plotXi?"_xi":"",cuts18012?"_cuts18012":""),"RECREATE");
  zMidHist->Write("zDist_mid");
  nMidHist->Write("Ntot_mid");
  zFwdHist->Write("zDist_fwd");
  nFwdHist->Write("Ntot_fwd");
  zHist->Write("zDist");
  fsave->Close();
}

Bool_t isGlobalMuonInAccept2019 (TLorentzVector* Muon) {
  return (fabs(Muon->Eta()) < 2.4 &&
	  ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	   (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
	   (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
}



void drawMCPlot() {
  gStyle->SetOptStat(0);
  TFile* prFile = TFile::Open("Output/MCResults/mcResult_PP_prompt_all.root");
  TH1F* prHist = (TH1F*) prFile->Get("zDist");
  TFile* nprFile = TFile::Open("Output/MCResults/mcResult_PP_nonprompt_all.root");
  TH1F* nprHist = (TH1F*) nprFile->Get("zDist");

  prHist->SetBinContent(1,0); prHist->SetBinError(1,0);
  nprHist->SetBinContent(1,0); nprHist->SetBinError(1,0);

  prHist->Scale(1.0/prHist->Integral("width"));
  nprHist->Scale(1.0/nprHist->Integral("width"));

  prHist->SetMarkerColor(46);
  prHist->SetMarkerStyle(kFullCircle);
  prHist->SetMarkerSize(0.5);
  prHist->SetLineColor(46);
  prHist->SetLineWidth(2);

  nprHist->SetMarkerColor(30);
  nprHist->SetMarkerStyle(kFullCircle);
  nprHist->SetMarkerSize(0.5);
  nprHist->SetLineColor(30);
  nprHist->SetLineWidth(2);

  //double ybins [] = {0, 1.6, 2.4};
  TH1F *axisHist = new TH1F("axisHist","", 5, 0, 1);
  axisHist->GetYaxis()->SetRangeUser(0, 7.5);
  axisHist->GetYaxis()->SetTitle("1/N dN/dz");
  axisHist->GetXaxis()->SetTitle("z");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetYaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  axisHist->GetXaxis()->SetRangeUser(0.2, 1);
  TLegend* leg = new TLegend(0.62,0.62,0.88,0.80);

  leg->AddEntry(prHist, "prompt", "lp");
  leg->AddEntry(nprHist, "nonprompt", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *  text2 = new TLatex(0.160 ,0.62, "|#eta_{Jet}| < 2");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.03);
  text2->SetLineWidth(2);

  TLatex *  text21 = new TLatex(0.160 ,0.67, "30 < p_{T,Jet} < 40 GeV");
  text21->SetNDC();
  text21->SetTextFont(42);
  text21->SetTextSize(0.03);
  text21->SetLineWidth(2);

  TLatex *  text22 = new TLatex(0.160 ,0.72, "p_{T,J/#psi} > 6.5 GeV");
  text22->SetNDC();
  text22->SetTextFont(42);
  text22->SetTextSize(0.03);
  text22->SetLineWidth(2);

  TLatex *  text = new TLatex(0.18 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.06);
  text->SetLineWidth(8);

  TLatex *  text1 = new TLatex(0.16 ,0.77,"Preliminary");
  text1->SetNDC();
  text1->SetTextFont(52);
  text1->SetTextSize(0.04);
  text1->SetLineWidth(2);

  //TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
  //text4->SetNDC();
  //text4->SetTextFont(42);
  //text4->SetTextSize(0.04);
  //text4->SetLineWidth(2);


  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  axisHist->Draw();
  prHist->Draw("same e1");
  prHist->Draw("same hist");
  nprHist->Draw("same e1");
  nprHist->Draw("same hist");
  leg->Draw("same");
  text2->Draw("same");
  text21->Draw("same");
  text22->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  gSystem->mkdir("Output/MCResults");
  c->SaveAs("Output/MCResults/mcTruthPlot.pdf");
  c->SaveAs("Output/MCResults/mcTruthPlot.png");

}


void compMeasured(bool isPrompt=true, bool isGen=false, bool MattStep=false) {
  gStyle->SetOptStat(0);

  TFile* file19007 = TFile::Open(Form("TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_jetR4_AccEff_JEC_updatedCorrForComp.root",isPrompt?"PR":"NOPR"));
  TFile* file18012 = TFile::Open(Form("~/DimuonCADIs/HIN-16-004/Fitter/TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_AccEff_JEC_updatedCorrForComp.root",isPrompt?"PR":"NOPR"));

  TTree* tree19 = (TTree*) file19007->Get("treeForUnfolding");
  TTree* tree18 = (TTree*) file18012->Get("treeForUnfolding");

  TH1D* hist19 = new TH1D("hist19","", 5, 0.3, 1);
  TH1D* hist18 = new TH1D("hist18","", 5, 0.3, 1);
  
  if (isGen) {
    tree19->Draw("gen_z>>hist19","corr_AccEff_comp*corr_ptw*(jp_gen_pt>6.5 && jp_gen_pt<35 && fabs(jp_gen_rap)<1.6 && jt_ref_pt>25 && jt_ref_pt<35 && gen_z<1 && gen_z>0.44)");
    tree18->Draw("gen_z>>hist18","corr_AccEff_comp*corr_ptw*(jp_gen_pt>6.5 && jp_gen_pt<35 && fabs(jp_gen_rap)<1.6 && jt_ref_pt>25 && jt_ref_pt<35 && gen_z<1 && gen_z>0.44)");
  }
  else if (MattStep){
    tree19->Draw("z>>hist19","corr_AccEff_comp*corr_ptw*(jp_pt>6.5 && jp_pt<35 && fabs(jp_rap)<1.6 && jt_ref_pt>25 && jt_ref_pt<35 && z<1 &&z>0.44)");
    tree18->Draw("z>>hist18","corr_AccEff_comp*corr_ptw*(jp_pt>6.5 && jp_pt<35 && fabs(jp_rap)<1.6 && jt_ref_pt>25 && jt_ref_pt<35 && z<1 &&z>0.44)");
  }
  else {
    tree19->Draw("z>>hist19","corr_AccEff_comp*corr_ptw*(jp_pt>6.5 && jp_pt<35 && fabs(jp_rap)<1.6 && jt_pt>25 && jt_pt<35 && z<1 &&z>0.44)");
    tree18->Draw("z>>hist18","corr_AccEff_comp*corr_ptw*(jp_pt>6.5 && jp_pt<35 && fabs(jp_rap)<1.6 && jt_pt>25 && jt_pt<35 && z<1 &&z>0.44)");
  }
  
  hist19->Scale(1./hist19->Integral("width"));
  hist18->Scale(1./hist18->Integral("width"));

  hist19->SetLineColor(8);
  hist19->SetMarkerColor(8);
  hist19->SetMarkerStyle(kFullCircle);

  hist18->SetLineColor(9);
  hist18->SetMarkerColor(9);
  hist18->SetMarkerStyle(kFullSquare);

  hist19->GetYaxis()->SetTitle("1/N dN/dz");
  hist19->GetXaxis()->SetTitle("z");
  hist19->GetXaxis()->SetRangeUser(0,1);
  hist19->GetYaxis()->SetRangeUser(0,1.2*hist19->GetMaximum());
  if (hist18->GetMaximum()>hist19->GetMaximum()) hist19->GetYaxis()->SetRangeUser(0,1.5*hist18->GetMaximum());
  else hist19->GetYaxis()->SetRangeUser(0,1.5*hist19->GetMaximum());

  TLegend* leg = new TLegend(0.2,0.2,0.5,0.3);
  if (isPrompt) leg = new TLegend(0.2,0.4,0.5,0.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist19,"2017 pp, R=0.4","lp");
  leg->AddEntry(hist18,"2015 pp, R=0.4","lp");
  TCanvas* c = new TCanvas("c","",900,900);
  hist19->Draw("e1");
  hist18->Draw("e1 same");
  leg->Draw("same");
  c->SaveAs(Form("Output/%sDistComp%s_18012vs19007_ppMC%s_midRap_updatedCorrForComp_noUnderflow.pdf",isGen?"Tr":"Measure",MattStep?"recoZWithRefCuts":"",isPrompt?"PR":"NOPR"));
  c->SaveAs(Form("Output/%sDistComp%s_18012vs19007_ppMC%s_midRap_updatedCorrForComp_noUnderflow.png",isGen?"Tr":"Measure",MattStep?"recoZWithRefCuts":"",isPrompt?"PR":"NOPR"));
}

void compTr(bool isPrompt = false){
  gStyle->SetOptStat(0);

  TFile *file19007 = new TFile(Form("Output/MCResults/mcResult_PP_%s_all_cuts18012.root",isPrompt?"prompt":"nonprompt"));
  TFile *file18012 = new TFile(Form("~/DimuonCADIs/HIN-16-004/Fitter/Output/MCResults/mcResult_%s_all.root",isPrompt?"prompt":"nonprompt"));

  TH1F* hist19 = (TH1F*) file19007->Get("zDist_mid");
  TH1F* hist18 = (TH1F*) file18012->Get("zDist_mid");

  int nBin = hist19->GetNbinsX()+1;
  for (int i=0; i<=nBin; i++) {
    cout <<"[INFO] bin "<<i<<", bin center = "<<hist19->GetBinCenter(i)<<", lower edge = "<<hist19->GetBinLowEdge(i)<<", bin content = "<<hist19->GetBinContent(i)<<endl;
    if (hist19->GetBinCenter(i)<0.44 || hist19->GetBinCenter(i)>1) {
      hist19->SetBinContent(i,0);
      hist19->SetBinError(i,0);
    }
  }

  nBin = hist18->GetNbinsX()+1;
  for (int i=0; i<=nBin; i++) {
    if (hist18->GetBinCenter(i)<0.44 || hist18->GetBinCenter(i)>1) {
      hist18->SetBinContent(i,0);
      hist18->SetBinError(i,0);
    }
  }

  hist19->Scale(1./hist19->Integral("width"));
  hist18->Scale(1./hist18->Integral("width"));

  hist19->SetLineColor(8);
  hist19->SetMarkerColor(8);
  hist19->SetMarkerStyle(kFullCircle);

  hist18->SetLineColor(9);
  hist18->SetMarkerColor(9);
  hist18->SetMarkerStyle(kFullSquare);

  hist19->GetYaxis()->SetTitle("1/N dN/dz");
  hist19->GetXaxis()->SetTitle("z");
  hist19->GetXaxis()->SetRangeUser(0,1);
  hist19->GetYaxis()->SetRangeUser(0,1.2*hist19->GetMaximum());
  if (hist18->GetMaximum()>hist19->GetMaximum()) hist19->GetYaxis()->SetRangeUser(0,1.2*hist18->GetMaximum());

  hist19->GetXaxis()->SetRangeUser(023,0.99);

  TLegend* leg = new TLegend(0.2,0.2,0.5,0.3);
  if (isPrompt) leg = new TLegend(0.2,0.4,0.5,0.5);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hist19,"2017 pp, R=0.4","lp");
  leg->AddEntry(hist18,"2015 pp, R=0.4","lp");
  TCanvas* c = new TCanvas("c","",900,900);
  hist19->Draw("e1");
  hist18->Draw("e1 same");
  leg->Draw("same");
  c->SaveAs(Form("Output/GenDistComp_18012vs19007_ppMC%s_midRap_noUnderflow.pdf",isPrompt?"PR":"NOPR"));
  c->SaveAs(Form("Output/GenDistComp_18012vs19007_ppMC%s_midRap_noUnderflow.png",isPrompt?"PR":"NOPR"));

}

void updatetreeCorr(bool is19, bool isPr) {

  ////////
  TFile *accFile19 = TFile::Open(Form("../Efficiency/FilesAccxEff_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights_noTnpWeights/Acc/%sAccHists_%s.root",isPr?"pr":"npr","PP"));
  TH2F* accNum19 = (TH2F*) accFile19->Get("hnum_2d_nominal");
  TH2F* accDeno19 =(TH2F*) accFile19->Get("hdeno_2d");
  accNum19->Divide(accDeno19);

  TFile *effFile19 = TFile::Open(Form("../Efficiency/FilesAccxEff_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights_noTnpWeights/Eff/%sEffHists_%s.root",isPr?"pr":"npr","PP"));
  TH2F* effNum19 = (TH2F*) effFile19->Get("hnum_noweights");
  TH2F* effDeno19 =(TH2F*) effFile19->Get("hdeno_pty");
  effNum19->Divide(effDeno19);
  cout <<"got the acceptance and efficiency files for 19-007"<<endl;
  //MyEfficiency/FilesAccxEff/Acc/prAccHists.root
  TFile *accFile18 = TFile::Open(Form("~/DimuonCADIs/HIN-16-004/MyEfficiency/FilesAccxEff/Acc/%sAccHists.root",isPr?"pr":"npr"));
  TH2F* accNum18 = (TH2F*) accFile18->Get("hnum_2d_nominal");
  TH2F* accDeno18 =(TH2F*) accFile18->Get("hdeno_2d");
  accNum18->Divide(accDeno18);
  cout <<"got the acceptance files for 18-012"<<endl;

  TFile *effFile18 = TFile::Open(Form("~/DimuonCADIs/HIN-16-004/MyEfficiency/FilesAccxEff/Eff/%sEffHists.root",isPr?"pr":"npr"));
  TH2F* effNum18 = (TH2F*) effFile18->Get("hnum_2d_noweights");
  TH2F* effDeno18 =(TH2F*) effFile18->Get("hdeno_2d");
  effNum18->Divide(effDeno18);

  cout <<"got the acceptance and efficiency files for 18-012"<<endl;
  ////////
  TFile* file19007 = TFile::Open(Form("TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_jetR4_AccEff_JEC.root",isPr?"PR":"NOPR"));
  TFile* file18012 = TFile::Open(Form("~/DimuonCADIs/HIN-16-004/Fitter/TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_AccEff_JEC.root",isPr?"PR":"NOPR"));

  TTree* told = (TTree*) file19007->Get("treeForUnfolding");
  if (!is19) told = (TTree*) file18012->Get("treeForUnfolding");

  string fileName = Form("TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_jetR4_AccEff_JEC_updatedCorrForComp.root",isPr?"PR":"NOPR");
  if (!is19) fileName = Form("~/DimuonCADIs/HIN-16-004/Fitter/TreesForUnfolding/tree_MCJPSI%s_PP_NoBkg_AccEff_JEC_updatedCorrForComp.root",isPr?"PR":"NOPR");
  TFile* fsave = new TFile(fileName.c_str(),"RECREATE");

  TTree *tnew = told->CloneTree(0);
  tnew->SetAutoSave(0);
  tnew->SetAutoFlush(0);

  float jp_pt; float jp_rap; float corr_AccEff;
  float corr_AccEff_comp=1.0;


  told->SetBranchAddress("jp_pt", &jp_pt);
  told->SetBranchAddress("jp_rap", &jp_rap);
  told->SetBranchAddress("corr_AccEff", &corr_AccEff);
  tnew->Branch("corr_AccEff_comp", &corr_AccEff_comp, "corr_AccEff_comp/F");

  int nentries = told->GetEntries();

  for (int i=0; i<nentries; i++) {
    if (i%10000 == 0 ) cout <<"[INFO] entry = "<< i <<"/"<<nentries<<endl;
    told->GetEntry(i);
    if (corr_AccEff>0.51) {
      if (is19)
	corr_AccEff_comp = 1./(accNum19->GetBinContent(accNum19->FindBin(jp_rap,jp_pt)) * effNum19->GetBinContent(effNum19->FindBin(jp_rap,jp_pt)));
      else 
	corr_AccEff_comp = 1./(accNum18->GetBinContent(accNum18->FindBin(jp_rap,jp_pt)) * effNum18->GetBinContent(effNum18->FindBin(jp_rap,jp_pt)));
      if (!corr_AccEff_comp) corr_AccEff_comp=1;
    }
    else corr_AccEff_comp=1;
    tnew->Fill();
  }
  
  tnew->Write("treeForUnfolding");
  fsave->Close();
}


void updatetreeCorrData(bool is19) {
  ////////
  //TFile* file19007 = TFile::Open("/data_CMS/cms/diab/JpsiJet/TreesForUnfolding/tree_DATA_PP_jetR4_AccEff_JEC.root");
  TFile* file19007 = TFile::Open("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_MCJPSINOPR_PP_NoBkg_jetR3_AccEff_JEC.root");
  TFile* file18012 = TFile::Open("/data_CMS/cms/diab/JpsiJet/TreesForUnfolding_2015/tree_DATA_PP_AccEff_JEC.root");

  TTree* told = (TTree*) file19007->Get("treeForUnfolding");
  if (!is19) told = (TTree*) file18012->Get("treeForUnfolding");

  //string fileName = "/data_CMS/cms/diab/JpsiJet/TreesForUnfolding/tree_DATA_PP_jetR4_AccEff_JEC_updatedCorr.root";
  string fileName = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_MCJPSINOPR_PP_NoBkg_jetR3_AccEff_JEC_updatedCorr.root";
  if (!is19) fileName = "/data_CMS/cms/diab/JpsiJet/TreesForUnfolding_2015/tree_DATA_PP_AccEff_JEC_updatedCorr.root";
  TFile* fsave = new TFile(fileName.c_str(),"RECREATE");

  TTree *tnew = told->CloneTree(0);
  tnew->SetAutoSave(0);
  tnew->SetAutoFlush(0);

  float corr_AccEff;
  float corr_AccEff_comp=1.0;

  told->SetBranchAddress("corr_AccEff", &corr_AccEff);
  tnew->Branch("corr_AccEff_comp", &corr_AccEff_comp, "corr_AccEff_comp/F");

  int nentries = told->GetEntries();

  for (int i=0; i<nentries; i++) {
    if (i%10000 == 0 ) cout <<"[INFO] entry = "<< i <<"/"<<nentries<<endl;
    told->GetEntry(i);
    if (corr_AccEff>0.51) {
      corr_AccEff_comp=corr_AccEff;
    }
    else corr_AccEff_comp=1;
    tnew->Fill();
  }
  
  tnew->Write("treeForUnfolding");
  fsave->Close();
}
