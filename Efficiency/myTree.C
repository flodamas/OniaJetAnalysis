#define myTree_cxx
#define _USE_MATH_DEFINES
#include "myTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>

#include "tnp_weight.h"

using namespace std;
using namespace  RooFit;
//using namespace RooPlot;

TLorentzVector* matchReco;
TLorentzVector* matchGen;

Float_t jpsi_m;
Float_t jpsi_pt;
Float_t jpsi_eta;
Float_t jpsi_phi;
Float_t jpsi_rap;
Float_t dr;
Float_t dphi;
Float_t dphimin;
Float_t deta;
Float_t drmin; 
Float_t z=100;
int triggerIndex_PP =0;
int k;
Float_t weight;
Float_t tnp_weight;
Double_t ptbins []={3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35, 40, 45, 50};
int nptbins = ((sizeof(ptbins)/sizeof(double))-1);
Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
int netabins = ((sizeof(etabins)/sizeof(double))-1);

Double_t ybins []={0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
int nybins = ((sizeof(ybins)/sizeof(double))-1);


void myTree::EffCalc()
{
  if (isMc)
    {
      TH1F* ptg = new TH1F ("ptg", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins);
      TH1F* ptr = new TH1F ("ptr", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins);
      TH1F* rapr = new TH1F ("rapr","N_{reco} vs #eta; #eta; N_{reco}", netabins, etabins);
      TH1F* rapg = new TH1F ("rapg","N_{gen} vs #eta; #eta; N_{total}", netabins, etabins);
      TH2F* ptrapg = new TH2F ("ptrapg", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", netabins, etabins, nptbins, ptbins);
      TH2F* ptrapr = new TH2F ("ptrapr", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);

      TH2F* num_binned = new TH2F ("num_", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_plus1sig = new TH2F ("num_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_minus1sig = new TH2F ("num_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_sta = new TH2F ("num_muid_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid = new TH2F ("num_muid", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_plus1sig = new TH2F ("num_muid_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_minus1sig = new TH2F ("num_muid_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta = new TH2F ("num_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta_plus1sig = new TH2F ("num_sta_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta_minus1sig = new TH2F ("num_sta_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_trk_plus1sig = new TH2F ("num_trk_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_trk_minus1sig = new TH2F ("num_trk_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);

      ptg->Sumw2(); ptr->Sumw2(); rapr->Sumw2(); rapg->Sumw2(); ptrapg->Sumw2(); ptrapr->Sumw2(); num_binned->Sumw2(); num_plus1sig->Sumw2(); num_minus1sig->Sumw2(); num_muid_sta->Sumw2(); num_muid->Sumw2(); num_muid_plus1sig->Sumw2(); num_muid_minus1sig->Sumw2(); num_sta->Sumw2(); num_sta_plus1sig->Sumw2(); num_sta_minus1sig->Sumw2(); num_trk_plus1sig->Sumw2(); num_trk_minus1sig->Sumw2();

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
	      if (jpsi_pt>3 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  ptg->Fill(jpsi_pt);
		  rapg->Fill(jpsi_rap);
		  ptrapg->Fill(jpsi_rap, jpsi_pt);
		}
	    }

	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m=RecoQQ4mom->M();
		  if (
		      jpsi_pt > 3  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) && // if it matches the trigger 
		      (isMatchedRecoDiMuon(iQQ))
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  /////nominal
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  ptr->Fill(jpsi_pt,tnp_weight);
			  rapr->Fill(jpsi_rap,tnp_weight);
			  ptrapr->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			  /////systematics
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-10) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-10) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_binned->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_sta->Fill(jpsi_rap, jpsi_pt, tnp_weight);
 
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(-1) * tnp_weight_trk_pp(-1);
			  num_trk_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(-2) * tnp_weight_trk_pp(-2);
			  num_trk_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight);
			}
		    }
		}
	    }
	}
      TEfficiency* gptef = new TEfficiency("gptef", "AccxEff(pt)", nptbins, ptbins);
      gptef->SetStatisticOption(TEfficiency::kBBayesian);
      gptef->SetPassedHistogram(*ptr,"f");
      gptef->SetTotalHistogram(*ptg,"f");

      TEfficiency* grapef = new TEfficiency("grapef", "AccxEff(y)", netabins, etabins);
      grapef->SetStatisticOption(TEfficiency::kBBayesian);
      grapef->SetPassedHistogram(*rapr,"f");
      grapef->SetTotalHistogram(*rapg,"f");

      /////nominal 2D
      TEfficiency* eff_nom = new TEfficiency("eff_nom", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_nom->SetStatisticOption(TEfficiency::kBBayesian);
      eff_nom->SetPassedHistogram(*ptrapr,"f");
      eff_nom->SetTotalHistogram(*ptrapg,"f");
      eff_nom->SetName("hcorr_Jpsi_PP");

      ////systematics
      TEfficiency* eff_binned = new TEfficiency("eff_binned", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_binned->SetStatisticOption(TEfficiency::kBBayesian);
      eff_binned->SetPassedHistogram(*num_binned,"f");
      eff_binned->SetTotalHistogram(*ptrapg,"f");
      eff_binned->SetName("hcorr_binned");

      TEfficiency* eff_plus1sig = new TEfficiency("eff_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_plus1sig->SetPassedHistogram(*num_plus1sig,"f");
      eff_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_plus1sig->SetName("hcorr_plus1sig");

      TEfficiency* eff_minus1sig = new TEfficiency("eff_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_minus1sig->SetPassedHistogram(*num_minus1sig,"f");
      eff_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_minus1sig->SetName("hcorr_minus1sig");

      TEfficiency* eff_muid_sta = new TEfficiency("eff_muid_sta", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_sta->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_sta->SetPassedHistogram(*num_muid_sta,"f");
      eff_muid_sta->SetTotalHistogram(*ptrapg,"f");
      eff_muid_sta->SetName("hcorr_muid_sta");

      TEfficiency* eff_muid = new TEfficiency("eff_muid", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid->SetPassedHistogram(*num_muid,"f");
      eff_muid->SetTotalHistogram(*ptrapg,"f");
      eff_muid->SetName("hcorr_muid");

      TEfficiency* eff_muid_plus1sig = new TEfficiency("eff_muid_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_plus1sig->SetPassedHistogram(*num_muid_plus1sig,"f");
      eff_muid_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_muid_plus1sig->SetName("hcorr_muid_plus1sig");

      TEfficiency* eff_muid_minus1sig = new TEfficiency("eff_muid_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_minus1sig->SetPassedHistogram(*num_muid_minus1sig,"f");
      eff_muid_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_muid_minus1sig->SetName("hcorr_muid_minus1sig");

      TEfficiency* eff_sta = new TEfficiency("eff_sta", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta->SetPassedHistogram(*num_sta,"f");
      eff_sta->SetTotalHistogram(*ptrapg,"f");
      eff_sta->SetName("hcorr_sta");

      TEfficiency* eff_sta_plus1sig = new TEfficiency("eff_sta_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta_plus1sig->SetPassedHistogram(*num_sta_plus1sig,"f");
      eff_sta_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_sta_plus1sig->SetName("hcorr_sta_plus1sig");

      TEfficiency* eff_sta_minus1sig = new TEfficiency("eff_sta_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta_minus1sig->SetPassedHistogram(*num_sta_minus1sig,"f");
      eff_sta_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_sta_minus1sig->SetName("hcorr_sta_minus1sig");

      TEfficiency* eff_trk_plus1sig = new TEfficiency("eff_trk_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_trk_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_plus1sig->SetPassedHistogram(*num_trk_plus1sig,"f");
      eff_trk_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_trk_plus1sig->SetName("hcorr_trk_plus1sig");

      TEfficiency* eff_trk_minus1sig = new TEfficiency("eff_trk_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_trk_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_minus1sig->SetPassedHistogram(*num_trk_minus1sig,"f");
      eff_trk_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_trk_minus1sig->SetName("hcorr_trk_minus1sig");


      TFile* ef (0x0);
      if (isPr)
	ef = new TFile ("Utilities/pr_correction_AccEff.root", "RECREATE");
      else
	ef = new TFile ("Utilities/npr_correction_AccEff.root", "RECREATE");


      gptef->Write("effVsPt");
      grapef->Write("effVsRap");
      ptrapg->Write("hcorr_his_deno");
      ptrapr->Write("hcorr_his_num");
      eff_nom->Write("hcorr_Jpsi_PP");

      eff_binned->Write("hcorr_binned"); 
      eff_plus1sig->Write("hcorr_plus1sig"); 
      eff_minus1sig->Write("hcorr_minus1sig"); 
      eff_muid_sta->Write("hcorr_muid_sta"); 
      eff_muid->Write("hcorr_muid"); 
      eff_muid_plus1sig->Write("hcorr_muid_plus1sig"); 
      eff_muid_minus1sig->Write("hcorr_muid_minus1sig"); 
      eff_sta->Write("hcorr_sta"); 
      eff_sta_plus1sig->Write("hcorr_sta_plus1sig"); 
      eff_sta_minus1sig->Write("hcorr_sta_minus1sig"); 
      eff_trk_plus1sig->Write("hcorr_trk_plus1sig"); 
      eff_trk_minus1sig->Write("hcorr_trk_minus1sig");

      ef->Close();
    }
  else 
    cout<< "[ERROR] this is data and not MC"<<endl;
}

void myTree::EffSyst(int jtPtRange) {
  int jtPtmin;
  int jtPtmax;
  if (jtPtRange==0) {jtPtmin=25; jtPtmax=35;}
  if (jtPtRange==-1) {jtPtmin=15; jtPtmax=25;}
  if (jtPtRange==1) {jtPtmin=35; jtPtmax=45;}

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root","MCJPSINOPR"),"READ"); //always apply the efficiency on nonprompt
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);

  cout<<"[INFO] Loading the corrections"<<endl;
  TFile* corrFile = TFile::Open("Utilities/pr_correction_AccEff.root","READ"); //always use the prompt

  string corrName [] = 
    {
      "Jpsi_PP", //nominal
      "binned",
      "plus1sig",
      "minus1sig",
      "muid_sta",
      "muid",
      "muid_plus1sig",
      "muid_minus1sig",
      "sta",
      "sta_plus1sig",
      "sta_minus1sig",
      "trk_plus1sig",
      "trk_minus1sig"
    };

  string systName [] =
    {                                                                                                                                                                              
      "binned",
      "trgSyst",
      "muidSyst",
      "staSyst",
      "trkSyst"
  };
  TObjArray *corrHis = new TObjArray(15);
  TObjArray *countHis016 = new TObjArray(15);
  TObjArray *countHis1624 = new TObjArray(15);
  TEfficiency *corrTemp = NULL; TH1F* countTemp = NULL; TH1F* countTemp2 = NULL;

  vector<double> v016;
  vector<double> v1624;

  for(int i=0; i<13; i++) {
    corrTemp = (TEfficiency*) corrFile->Get(Form("hcorr_%s",corrName[i].c_str()));
    corrHis->Add(corrTemp);
    countTemp = new TH1F (Form("his_016_%s",corrName[i].c_str()), Form("mid,  tnp %s",corrName[i].c_str()), 5, 0, 1);
    countHis016->Add(countTemp);
    countTemp = new TH1F (Form("his_1624_%s",corrName[i].c_str()), Form("fwd,  tnp %s",corrName[i].c_str()), 5, 0, 1);
    countHis1624->Add(countTemp);
  }
  //filling the histograms
  cout<<"[INFO] Filling the histograms"<<endl;
  int nentries = trNom->GetEntries();
  for (int jentry=0; jentry<nentries; jentry++) {
    if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<endl;
    trNom->GetEntry(jentry);
    
    if (jt_pt > jtPtmin && jt_pt < jtPtmax && abs(jt_rap) < 2.4 && jp_mass > 2.6 && jp_mass < 3.5) {
	if (abs(jp_rap) < 1.6 && jp_pt > 6.5 && jp_pt < 35) {
	  for (int i=0; i<13; i++) {
	    corrTemp = (TEfficiency*) corrHis->At(i);
	    countTemp = (TH1F*) countHis016->At(i);
	    countTemp->Fill(z,1.0/corrTemp->GetEfficiency(corrTemp->FindFixBin(jp_rap,jp_pt)));
	  }
	}
	else if (abs(jp_rap) > 1.6 && abs(jp_rap) < 2.4 && jp_pt > 3 && jp_pt < 35) {
	  for (int i=0; i<13; i++) {
	    corrTemp = (TEfficiency*) corrHis->At(i);
	    countTemp = (TH1F*) countHis1624->At(i);
	    countTemp->Fill(z,1.0/corrTemp->GetEfficiency(corrTemp->FindFixBin(jp_rap,jp_pt)));
	  }
	}
    }
  }
 
    cout<<"[INFO] Getting the ratios and filling the csv files"<<endl;

    ////// filling the different syst files
    for (int i = 0; i<5 ; i++)
      { 
	ofstream file016(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_tnp%s.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt",systName[i].c_str()));
	ofstream file1624(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_tnp%s.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt",systName[i].c_str()));
	file016 << "tnp " << systName[i] << endl;
	file1624 << "tnp " << systName[i] << endl;
	v016.clear();
	v1624.clear();
	for (int j=1; j<5; j++){
	  v016.clear();
	  v1624.clear();	    
	  countTemp = (TH1F*) countHis016->At(0);
	  countTemp2 = (TH1F*) countHis1624->At(0);
	  v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));

	  if (i==0) { /////binned
	    countTemp = (TH1F*) countHis016->At(1);
	    v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	    countTemp2 = (TH1F*) countHis1624->At(1);
	    v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  }
	  if (i==1) { /////trgSyst
	    countTemp = (TH1F*) countHis016->At(2);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(2);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));

	    countTemp = (TH1F*) countHis016->At(3);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(3);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  }
	  if (i==2) { ////muidSyst
	    countTemp = (TH1F*) countHis016->At(6);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(6);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));

	    countTemp = (TH1F*) countHis016->At(7);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(7);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  }
	  if (i==3) { //////staSyst
	    countTemp = (TH1F*) countHis016->At(9);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(9);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));

	    countTemp = (TH1F*) countHis016->At(10);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(10);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  }
	  if (i==4) { ///////trkSyst
	    countTemp = (TH1F*) countHis016->At(11);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(11);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));

	    countTemp = (TH1F*) countHis016->At(12);
            v016.push_back((double)(countTemp->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
            countTemp2 = (TH1F*) countHis1624->At(12);
            v1624.push_back((double)(countTemp2->GetBinContent(countTemp->FindBin(j*0.2+0.1))));
	  }
	  file016 << "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< maxdiff(v016, 1) << endl;
	  file1624 << "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< maxdiff(v1624, 1) << endl;
	}
	file016.close();
	file1624.close();
      }
}

void myTree::EffStatToy(int nToys) {
  // Randomise TEfficiency 100 times. Output will be a TObjArray with 100 TH2
  TFile* effFile = TFile::Open("Utilities/pr_correction_AccEff.root","READ"); 
  TEfficiency *eff = (TEfficiency*) effFile->Get("hcorr_Jpsi_PP");

  TObjArray *arrEffs = new TObjArray(); // Array to store efficiencies
  arrEffs->SetOwner(true);

  TRandom* rnd = new TRandom3(); // For randomisation

  for (int i=0 ; i < nToys ; i++)
    {
      if (i%10==0) cout <<"[INFO] toy "<<i<<"/"<<nToys<<endl;
      TH2* histoTot = (TH2*) eff->GetTotalHistogram()->Clone(Form("hTotToy_%d",i)); // Get corresponding histo and number of bins
      histoTot->Sumw2();
      TH2* histoPass = (TH2*) histoTot->Clone(Form("hPassToy_%d",i)); // Get corresponding histo and number of bins
      histoPass->Sumw2();
      int nBinsX = histoTot->GetNbinsX();
      int nBinsY = histoTot->GetNbinsY();
  
      for (int j = 1 ; j <= nBinsX ; j++)
	{
	  for (int k = 1 ; k <= nBinsY ; k++)
	    {
	      int bin = histoTot->GetBin(j,k);      
	      double effVal = eff->GetEfficiency(bin);
	      int ntot = histoTot->GetBinContent(bin);
	      int newPass = rnd->Binomial(ntot,effVal);
	      histoPass->SetBinContent(histoPass->GetBin(j,k),newPass);
	    }
	}
      TEfficiency* newHeff = (TEfficiency*) eff->Clone(Form("effToy_%d",i)); // Histo to store new 2D eff
      newHeff->SetPassedHistogram(*histoPass, "f");
      newHeff->SetTotalHistogram(*histoTot, "f");
      arrEffs->Add(newHeff);
    }
  effFile->Close();
  TFile* fsave = new TFile(Form("toyMC/accEff%dToys.root", nToys),"RECREATE");
  arrEffs->Write("accEffArray", TObject::kSingleKey);
  fsave->Close();
}

void myTree::EffStat(int jtPtRange){
  int jtPtmin;
  int jtPtmax;
  if (jtPtRange==0) {jtPtmin=25; jtPtmax=35;}
  if (jtPtRange==-1) {jtPtmin=15; jtPtmax=25;}
  if (jtPtRange==1) {jtPtmin=35; jtPtmax=45;}

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root","MCJPSINOPR"),"READ"); //always apply it on non prompt
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);

  TEfficiency* tempEff = NULL;

  TObjArray* hisArr016 = new TObjArray(); //hisArr016->SetOwner(kTRUE);
  TObjArray* hisArr1624 = new TObjArray(); //hisArr1624->SetOwner(kTRUE);

  TFile *corFile = TFile::Open("toyMC/accEff100Toys.root");
  TObjArray* corArr = (TObjArray*) corFile->Get("accEffArray");
  TFile *nomFile = TFile::Open("Utilities/pr_correction_AccEff.root","READ");
  TEfficiency* nomEff = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP");
  TH1F* hist016 = NULL;
  TH1F* hist1624 = NULL;

  int nentries = trNom->GetEntries();
  for (int i=0; i<=100; i++)
    {
      hist016 = new TH1F (Form("hist016_%d",i), "", 5, 0, 1);
      hist1624 = new TH1F (Form("hist1624_%d",i), "", 5, 0, 1);
      if (i==0) cout<<"[INFO] Applying the nominal AccxEff"<<endl;
      else cout<<"[INFO] Applying toy "<<i<<"/100"<<endl;
      if (i==0)
	tempEff=nomEff;
      else 
	tempEff=(TEfficiency*) corArr->At(i-1);

      for (int jentry=0; jentry<nentries; jentry++) {
	trNom->GetEntry(jentry);
	if (jt_pt > jtPtmin && jt_pt < jtPtmax && abs(jt_rap) < 2.4 && jp_mass > 2.6 && jp_mass < 3.5) {
	  if (abs(jp_rap) < 1.6 && jp_pt > 6.5 && jp_pt < 35){
	    hist016->Fill(z,1.0/tempEff->GetEfficiency(tempEff->FindFixBin(jp_rap,jp_pt)));
	  }
	  else if (abs(jp_rap) > 1.6 && abs(jp_rap) < 2.4 && jp_pt > 3 && jp_pt < 35) {
	    hist1624->Fill(z,1.0/tempEff->GetEfficiency(tempEff->FindFixBin(jp_rap,jp_pt)));
	  }
	}// end of jt pt cut
      }// end of the tree entries
      hisArr016->Add(hist016);
      hisArr1624->Add(hist1624);
    } // end of the variations

  TFile* fsave = new TFile("toyMC/AccEffStatArray.root","RECREATE"); 
  hisArr016->Write("arr016",TObject::kSingleKey);
  hisArr1624->Write("arr1624",TObject::kSingleKey);
  fsave->Close();

  cout<<"[INFO] making csv files"<<endl;
  ofstream file016(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_AccEffStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  ofstream file1624(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_AccEffStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));

  file016<<"AccxEff stat"<<endl;
  file1624<<"AccxEff stat"<<endl;
  vector<double> v016;
  vector<double> v1624;
  v016.clear();
  v1624.clear();
  TH1F* temp = NULL;
  for (int j=1; j<5; j++)
    {
      v016.clear();
      v1624.clear();
      for (int i=0; i<=100; i++)
	{
	  temp = (TH1F*) hisArr016->At(i);
	  v016.push_back(temp->GetBinContent(temp->FindBin(j*0.2+0.1)));
	  temp = (TH1F*) hisArr1624->At(i);
	  v1624.push_back(temp->GetBinContent(temp->FindBin(j*0.2+0.1)));
	}
      file016<< "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v016,true) << endl;
      file1624<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v1624,true) << endl;
    }
  file016.close();
  file1624.close();
  corFile->Close();
  nomFile->Close();
  trFile->Close();

}

void myTree::TnpToy(int min, int max) {
  if (isMc) {
    TFile* denoFile = TFile::Open("Utilities/pr_correction_AccEff.root","READ");
    TH2F *deno = (TH2F*) denoFile->Get("hcorr_his_deno");

    TObjArray *trgEff = new TObjArray(); trgEff->SetOwner(kTRUE);
    TObjArray *muidEff = new TObjArray(); muidEff->SetOwner(kTRUE);
    TObjArray *staEff = new TObjArray(); staEff->SetOwner(kTRUE);

    TEfficiency* trgTemp = NULL;
    TEfficiency* muidTemp = NULL;
    TEfficiency* staTemp = NULL;

    TH2F* trgNum = NULL;
    TH2F* muidNum = NULL;
    TH2F* staNum = NULL;
  
    for (int a=min; a<=max; a++)
      {
	trgNum = new TH2F (Form("trgNum%d",a), "", netabins, etabins, nptbins, ptbins); trgNum->Sumw2();
	muidNum = new TH2F (Form("muidNum%d",a), "", netabins, etabins, nptbins, ptbins); muidNum->Sumw2();
	staNum = new TH2F (Form("staNum%d",a), "", netabins, etabins, nptbins, ptbins);staNum->Sumw2();

	Long64_t nentries =fChain->GetEntries();
	//nentries = 2000000;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) 
	  {
	    if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<" in toy "<<a<<" (" <<min<<"->"<<max<<")"<<endl;
	    Long64_t ientry = LoadTree(jentry);
	    if (ientry < 0) break;
	    nb = fChain->GetEntry(jentry);   nbytes += nb;
	    if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	      {
		for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		  {
		    TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		    TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		    TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		    jpsi_pt = RecoQQ4mom->Pt();
		    jpsi_rap = RecoQQ4mom->Rapidity();
		    jpsi_m=RecoQQ4mom->M();
		    if (
			jpsi_pt > 3  &&
			(areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
			(passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
			(isTriggerMatch(iQQ, triggerIndex_PP)) && // if it matches the trigger 
			(isMatchedRecoDiMuon(iQQ))
			)
		      {
			if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			  {
			    //trg 
			    tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
			                 tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			    trgNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			    //muid
			    tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			                 tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
			                 tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			    muidNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			    //sta
			    tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			                 tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
			                 tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			    staNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);
			  }
		      }
		  }
	      }
	  } //end of events   
	trgTemp = new TEfficiency(Form("trgTemp%d",a), "trgToy AccxEff(y,pt); y; pt; AccxEff", netabins, etabins, nptbins, ptbins);
	trgTemp->SetStatisticOption(TEfficiency::kBBayesian);
	trgTemp->SetPassedHistogram(*trgNum,"f");
	trgTemp->SetTotalHistogram(*deno,"f");
	trgTemp->SetName(Form("trgToy%d",a));
	trgEff->Add(trgTemp);

	muidTemp = new TEfficiency(Form("muidTemp%d",a), "muidToy AccxEff(y,pt); y; pt; AccxEff", netabins, etabins, nptbins, ptbins);
	muidTemp->SetStatisticOption(TEfficiency::kBBayesian);
	muidTemp->SetPassedHistogram(*muidNum,"f");
	muidTemp->SetTotalHistogram(*deno,"f");
	muidTemp->SetName(Form("muidToy%d",a));	
	muidEff->Add(muidTemp);

	staTemp = new TEfficiency(Form("staTemp%d",a), "staToy AccxEff(y,pt); y; pt; AccxEff", netabins, etabins, nptbins, ptbins);
	staTemp->SetStatisticOption(TEfficiency::kBBayesian);
	staTemp->SetPassedHistogram(*staNum,"f");
	staTemp->SetTotalHistogram(*deno,"f");
	staTemp->SetName(Form("staToy%d",a));
	staEff->Add(staTemp);
      }//endl of 100 toy loop
    gSystem->mkdir("toyMC");
    TFile *fsave = new TFile (Form("toyMC/%stoys%d%d.root", (isPr)?"pr":"npr", min, max),"RECREATE");
    trgEff->Write(Form("trg%dtoys", max-min+1), TObject::kSingleKey);
    muidEff->Write(Form("muid%dtoys", max-min+1), TObject::kSingleKey);
    staEff->Write(Form("sta%dtoys", max-min+1), TObject::kSingleKey);
    fsave->Close();
  }
  else cout<< "this is data and not MC"<<endl;
}

void myTree::TnpStat(int jtPtRange) {
  int jtPtmin;
  int jtPtmax;
  if (jtPtRange==0) {jtPtmin=25; jtPtmax=35;}
  if (jtPtRange==-1) {jtPtmin=15; jtPtmax=25;}
  if (jtPtRange==1) {jtPtmin=35; jtPtmax=45;}

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root","MCJPSINOPR"),"READ"); //always apply it on non prompt
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);

  TH1F* temp1 = NULL;
  TH1F* temp2 = NULL;

  TObjArray* trgArr016 = new TObjArray(); trgArr016->SetOwner(kTRUE);
  TObjArray* muidArr016 = new TObjArray(); muidArr016->SetOwner(kTRUE);
  TObjArray* staArr016 = new TObjArray(); staArr016->SetOwner(kTRUE);
  
  TObjArray* trgArr1624 = new TObjArray(); trgArr1624->SetOwner(kTRUE);
  TObjArray* muidArr1624 = new TObjArray(); muidArr1624->SetOwner(kTRUE);
  TObjArray* staArr1624 = new TObjArray(); staArr1624->SetOwner(kTRUE);

  cout<<"[INFO] Creating arrays of histograms."<<endl;
  for (int i =0; i<100; i++)
    {
      temp1 = new TH1F (Form("trg016_%d",i), Form("trg %d",i), 5, 0, 1);
      trgArr016->Add(temp1);
      temp1 = new TH1F (Form("trg1624_%d",i), Form("trg %d",i), 5, 0, 1);
      trgArr1624->Add(temp1);
      temp1 = new TH1F (Form("muid016_%d",i), Form("muid %d",i), 5, 0, 1);
      muidArr016->Add(temp1);
      temp1 = new TH1F (Form("muid1624_%d",i), Form("muid %d",i), 5, 0, 1);
      muidArr1624->Add(temp1);
      temp1 = new TH1F (Form("sta016_%d",i), Form("sta %d",i), 5, 0, 1);
      staArr016->Add(temp1);
      temp1 = new TH1F (Form("sta1624_%d",i), Form("sta %d",i), 5, 0, 1);
      staArr1624->Add(temp1);
    }

  cout<<"[INFO] Loading the nominal corrections"<<endl;
  TFile* fNom = TFile::Open("Utilities/pr_correction_AccEff.root","READ"); 
  TEfficiency* corNom = (TEfficiency*) fNom->Get("hcorr_Jpsi_PP");
  TEfficiency* corTemp3 = NULL;
  TEfficiency* corTemp2 = NULL;
  TEfficiency* corTemp1 = NULL;
  TH1F* nom016 = new TH1F ("nom016", "", 5, 0, 1);
  TH1F* nom1624 = new TH1F ("nom1624", "", 5, 0, 1);

  TFile* toyFile = NULL;
  TObjArray* arrTemp1 = NULL;
  TObjArray* arrTemp2 = NULL;
  TObjArray* arrTemp3 = NULL;

  cout<<"[INFO] starting to process the tree."<<endl;

  int nentries = trNom->GetEntries();
  for (int i=0; i<10; i++)
    {
      toyFile = TFile::Open(Form("toyMC/prtoys%d%d.root", i*10, i*10+9));
      arrTemp1 = (TObjArray*) toyFile->Get("trg10toys");
      arrTemp2 = (TObjArray*) toyFile->Get("muid10toys");
      arrTemp3 = (TObjArray*) toyFile->Get("sta10toys");
      for (int j=0; j<10; j++)
	{
	  corTemp1 = (TEfficiency*) arrTemp1->At(j);
	  corTemp2 = (TEfficiency*) arrTemp2->At(j);
	  corTemp3 = (TEfficiency*) arrTemp3->At(j);
	  
	  for (int jentry=0; jentry<nentries; jentry++) {
	    if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<" in toy "<<i*10+j<<endl;
	    trNom->GetEntry(jentry);
	    if (jt_pt > jtPtmin && jt_pt < jtPtmax && abs(jt_rap) < 2.4 && jp_mass > 2.6 && jp_mass < 3.5) {
	      if (abs(jp_rap) < 1.6 && jp_pt > 6.5 && jp_pt < 35){
		if (i==0 && j==0)
		  nom016->Fill(z,1.0/corNom->GetEfficiency(corNom->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) trgArr016->At(i*10+j);
		temp1->Fill(z,1.0/corTemp1->GetEfficiency(corTemp1->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) muidArr016->At(i*10+j);
		temp1->Fill(z,1.0/corTemp2->GetEfficiency(corTemp2->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) staArr016->At(i*10+j);
		temp1->Fill(z,1.0/corTemp3->GetEfficiency(corTemp3->FindFixBin(jp_rap,jp_pt)));
	      }
	      else if (abs(jp_rap) > 1.6 && abs(jp_rap) < 2.4 && jp_pt > 3 && jp_pt < 35) {
		if (i==0 && j==0)
		  nom1624->Fill(z,1.0/corNom->GetEfficiency(corNom->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) trgArr1624->At(i*10+j);
		temp1->Fill(z,1.0/corTemp1->GetEfficiency(corTemp1->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) muidArr1624->At(i*10+j);
		temp1->Fill(z,1.0/corTemp2->GetEfficiency(corTemp2->FindFixBin(jp_rap,jp_pt)));
		temp1 = (TH1F*) staArr1624->At(i*10+j);
		temp1->Fill(z,1.0/corTemp3->GetEfficiency(corTemp3->FindFixBin(jp_rap,jp_pt)));
	      }
	    }// end of jt pt cut
	  }// end of tree entries loop
	}// end of tobjarrays inside every file
      toyFile->Close();
    }// end of files
  gStyle->SetOptStat(0);
  cout<<"[INFO] plotting the mid rapidity results."<<endl;  
  TCanvas *c1 = new TCanvas("c1", "trg plots", 1000, 1000);
  TCanvas *c2 = new  TCanvas("c2", "muid plots", 1000, 1000);
  TCanvas *c3 = new  TCanvas("c3", "sta plots", 1000, 1000);
  TLegend *leg = new TLegend(0.2, 0.6, 0.4, 0.9);

  nom016->SetLineColor(kRed);
  nom016->SetLineWidth(2);
  c1->cd();
  nom016->Draw();
  c2->cd();
  nom016->Draw();
  c3->cd();
  nom016->Draw();
  leg->AddEntry(nom016, "nominal", "lep");
  for (int i=0; i<10; i++)
    {
      c1->cd();
      temp1 = (TH1F*) trgArr016->At(i*10);
      temp1->SetLineColor(i+1); if (i==1||i==9) temp1->SetLineColor(i+12);
      temp1->SetLineWidth(1);
      temp1->Draw("same");
      leg->AddEntry(temp1, Form("var %d", i*10), "lep");

      c2->cd();
      temp1 = (TH1F*) muidArr016->At(i*10);
      temp1->SetLineColor(i+1); if (i==1||i==9) temp1->SetLineColor(i+10);
      temp1->SetLineWidth(1);
      temp1->Draw("same");

      c3->cd();
      temp1 = (TH1F*) staArr016->At(i*10);
      temp1->SetLineColor(i+1); if (i==1||i==9) temp1->SetLineColor(i+10);
      temp1->SetLineWidth(1);
      temp1->Draw("same");
    }
  c1->cd();
  nom016->Draw("same");
  leg->Draw("same");
  c2->cd();
  nom016->Draw("same");
  leg->Draw("same");
  c3->cd();
  nom016->Draw("same");
  leg->Draw("same");

  c1->SaveAs("toyMC/trgHist.root");
  c2->SaveAs("toyMC/muidHist.root");
  c3->SaveAs("toyMC/staHist.root");
  cout<<"[INFO] Getting the yields and filling the csv files"<<endl;
  ////// filling the different syst files
  
  ofstream file016_trg(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_tnptrgStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  ofstream file1624_trg(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_tnptrgStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt"));

  ofstream file016_muid(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_tnpmuidStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  ofstream file1624_muid(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_tnpmuidStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt"));

  ofstream file016_sta(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_tnpstaStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  ofstream file1624_sta(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_tnpstaStat.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt"));


  file016_trg << "tnp trg stat" << endl;
  file1624_trg << "tnp trg stat" << endl;

  file016_muid << "tnp muid stat" << endl;
  file1624_muid << "tnp muid stat" << endl;

  file016_sta << "tnp sta stat" << endl;
  file1624_sta << "tnp sta stat" << endl;
  vector<double> v016_trg;
  vector<double> v1624_trg;
  vector<double> v016_muid;
  vector<double> v1624_muid;
  vector<double> v016_sta;
  vector<double> v1624_sta;
  v016_trg.clear();
  v1624_trg.clear();
  v016_muid.clear();
  v1624_muid.clear();
  v016_sta.clear();
  v1624_sta.clear();

  for (int j=1; j<5; j++){
    v016_trg.clear();
    v1624_trg.clear();
    v016_muid.clear();
    v1624_muid.clear();
    v016_sta.clear();
    v1624_sta.clear();
    v016_trg.push_back(nom016->GetBinContent(nom016->FindBin(j*0.2+0.1)));
    v1624_trg.push_back(nom1624->GetBinContent(nom1624->FindBin(j*0.2+0.1)));
    v016_muid.push_back(nom016->GetBinContent(nom016->FindBin(j*0.2+0.1)));
    v1624_muid.push_back(nom1624->GetBinContent(nom1624->FindBin(j*0.2+0.1)));
    v016_sta.push_back(nom016->GetBinContent(nom016->FindBin(j*0.2+0.1)));
    v1624_sta.push_back(nom1624->GetBinContent(nom1624->FindBin(j*0.2+0.1)));

    for (int i = 0; i<100; i++)
      {
	temp1 = (TH1F*) trgArr016->At(i);
	v016_trg.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
	temp1 = (TH1F*) trgArr1624->At(i);
	v1624_trg.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
	temp1 = (TH1F*) muidArr016->At(i);
	v016_muid.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
	temp1 = (TH1F*) muidArr1624->At(i);
	v1624_muid.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
	temp1 = (TH1F*) staArr016->At(i);
	v016_sta.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
	temp1 = (TH1F*) staArr1624->At(i);
	v1624_sta.push_back(temp1->GetBinContent(temp1->FindBin(j*0.2+0.1)));
      }
    file016_trg<< "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v016_trg,true) << endl;
    file1624_trg<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v1624_trg, true) << endl;
    file016_muid<< "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v016_muid, true) << endl;
    file1624_muid<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v1624_muid, true) << endl;
    file016_sta<< "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v016_sta, true) << endl;
    file1624_sta<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(v1624_sta, true) << endl;
  }// end of z-range loop
  file016_trg.close();
  file1624_trg.close();
  file016_muid.close();
  file1624_muid.close();
  file016_sta.close();
  file1624_sta.close();
} //end of fct TnpStat

void myTree::EffMisMod(int jtPtRange) {
  int jtPtmin;
  int jtPtmax;
  if (jtPtRange==0) {jtPtmin=25; jtPtmax=35;}
  if (jtPtRange==-1) {jtPtmin=15; jtPtmax=25;}
  if (jtPtRange==1) {jtPtmin=35; jtPtmax=45;}

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root","MCJPSINOPR"),"READ"); //always apply it on non prompt
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);

  TFile *nomFile = TFile::Open("Utilities/pr_correction_AccEff.root","READ");
  TEfficiency* nomEff = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP");
  TFile *nprFile = TFile::Open("Utilities/npr_correction_AccEff.root","READ");
  TEfficiency* nprEff = (TEfficiency*) nprFile->Get("hcorr_Jpsi_PP");

  TH1F* hist016 = NULL;
  TH1F* hist1624 = NULL;
  TH1F* nom016 = NULL;
  TH1F* nom1624 = NULL;

  hist016 = new TH1F (Form("hist016"), "", 5, 0, 1);
  hist1624 = new TH1F (Form("hist1624"), "", 5, 0, 1);
  nom016 = new TH1F (Form("nom016"), "", 5, 0, 1);
  nom1624 =new TH1F (Form("nom1624"), "", 5, 0, 1);

  cout<<"[INFO] Processing the entries"<<endl;
  int nentries = trNom->GetEntries();
  for (int jentry=0; jentry<nentries; jentry++) {
    trNom->GetEntry(jentry);
    if (jt_pt > jtPtmin && jt_pt < jtPtmax && abs(jt_rap) < 2.4 && jp_mass > 2.6 && jp_mass < 3.5) {
      if (abs(jp_rap) < 1.6 && jp_pt > 6.5 && jp_pt < 35){
	nom016->Fill(z,1.0/nomEff->GetEfficiency(nomEff->FindFixBin(jp_rap,jp_pt)));
	hist016->Fill(z,1.0/nprEff->GetEfficiency(nprEff->FindFixBin(jp_rap,jp_pt)));
      }
      else if (abs(jp_rap) > 1.6 && abs(jp_rap) < 2.4 && jp_pt > 3 && jp_pt < 35) {
	nom1624->Fill(z,1.0/nomEff->GetEfficiency(nomEff->FindFixBin(jp_rap,jp_pt)));
	hist1624->Fill(z,1.0/nprEff->GetEfficiency(nprEff->FindFixBin(jp_rap,jp_pt)));
      }
    }// end of jt pt cut
  }// end of the tree entries
  cout<<"[INFO] filling csv files"<<endl;
  ofstream file016(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_AccEffMisMod.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  ofstream file1624(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_AccEffMisMod.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt"));
  file016<<"AccxEff mismodeling"<<endl;
  file1624<<"AccxEff mismodeling"<<endl;
  double val, var;
  for (int j=1; j<5; j++)
    {
      val=nom016->GetBinContent(nom016->FindBin(j*0.2+0.1));
      var=hist016->GetBinContent(hist016->FindBin(j*0.2+0.1));
      file016<< "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< 1.*(val-var)/val << endl;

      val=nom1624->GetBinContent(nom1624->FindBin(j*0.2+0.1));
      var=hist1624->GetBinContent(hist1624->FindBin(j*0.2+0.1));
      file1624<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< 1.*(val-var)/val << endl;
    }

  file016.close();
  file1624.close();
  nprFile->Close();
  nomFile->Close();
  trFile->Close();
} // end of the function


void myTree::FullEffSyst(int jtPtRange) {
  string systName [] = {"tnpbinned", "tnptrgSyst", "tnpmuidSyst", "tnpstaSyst", "tnptrkSyst", "tnptrgStat", "tnpmuidStat", "tnpstaStat", "AccEffStat", "AccEffMisMod"};

  ofstream file016(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_fullAccEff.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt"));
  ofstream file1624(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_fullAccEff.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt",isPr?"prompt":"nonprompt"));
  file016<<"AccxEff"<<endl;
  file1624<<"AccxEff"<<endl;
  double val1 [] = {0,0,0,0};
  double val2 [] = {0,0,0,0};
  for (int i=0; i<(sizeof(systName)/sizeof(systName[0])); i++)
    {
      vector<double> v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_016_NJpsi_%s_PP_%s.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt", systName[i].c_str()));
      vector<double> v2 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_1624_NJpsi_%s_PP_%s.csv", (jtPtRange==-1)? "lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", isPr?"prompt":"nonprompt", systName[i].c_str()));
      for (int j=0; j<4; j++)
	{
	  val1[j]=sqrt(pow(val1[j],2)+pow(v1[j],2));
	  val2[j]=sqrt(pow(val2[j],2)+pow(v2[j],2));
	}
    }
  for (int j=1; j<5; j++)
    {
      file016 << "0, 1.6, 6.5, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< val1[j-1] << endl;
      file1624 << "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<<  val2[j-1]<< endl;
    }
  file016.close();
  file1624.close();
}


double myTree::rms(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double s=0.;
  double s2=0.;
  for (unsigned int i=0; i<v.size(); i++) {
    if (v[i]==-999) continue;
    s2+=(v[i]-v[0])*(v[i]-v[0]);
    //cout<<"[INFO] v["<<i<<"] = "<< v[i] << " s = "<<s<<" s2 = "<<s2<<endl;
  }
  double ans = sqrt(s2*1.0/(v.size()-1));
  //cout<<"[INFO] s2 = "<<s2<<" ans = " <<ans/v[0] <<endl;
  if (isrelative) ans = ans/v[0];
  return ans;
}

double myTree::maxdiff(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double maxdiff=0;
  for (unsigned int i=1; i<v.size(); i++) {
    if (v[i]==-999) continue;
    maxdiff=max(maxdiff,fabs(v[i]-v[0]));
  }
  double ans = maxdiff;
  if (isrelative) ans = ans/v[0];
  return ans;
}


vector<double> myTree::readSyst(const char* systfile) {
  vector<double> ans;
    ifstream file(systfile);
    if (!(file.good())) return ans;

    string systname; getline(file,systname);

    string line;
    double zmin=0, zmax=0, rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;

    while (file.good()) {
      getline(file,line);
      if (line.size()==0) break;
      TString tline(line.c_str());
      TString t; Int_t from = 0, cnt=0;
      while (tline.Tokenize(t, from , ",")) {
	t.Strip(TString::kBoth,' ');
	value = atof(t.Data());
	if (cnt==0) rapmin = atof(t.Data());
	else if (cnt==1) rapmax = value;
	else if (cnt==2) ptmin = value;
	else if (cnt==3) ptmax = value;
	else if (cnt==4) zmin = value;
	else if (cnt==5) zmax = value;
	else if (cnt==6) centmin = value;
	else if (cnt==7) centmax = value;
	else if (cnt>8) {
	  cout << "Warning, too many fields, I'll take the last one." << endl;
	  continue;
	}
	cnt++;
      }
      ans.push_back(value);
    }
    file.close();
    return ans;
}


void myTree::ANEffPlots()
{
  if (isMc)
    {
      TH1F* Accpthist = new TH1F ("Accpthist", ";p_{T}^{#mu#mu} (GeV/c); Acceptance",nptbins,ptbins); Accpthist->Sumw2();
      TH1F* Effpthist = new TH1F("Effpthist", ";p_{T}^{#mu#mu} (GeV/c); Efficiency",nptbins,ptbins);Effpthist->Sumw2();
      TH1F* Accyhist = new TH1F("Accyhist", ";|y|^{#mu#mu}; Acceptance",nybins,ybins); Accyhist->Sumw2();
      TH1F* Effyhist = new TH1F("Effyhist", ";|y|^{#mu#mu}; Efficiency",nybins,ybins); Effyhist->Sumw2();
      TH1F* AccEffpthist = new TH1F("AccEffpthist", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff",nptbins,ptbins); AccEffpthist->Sumw2();
      TH1F* AccEffyhist = new TH1F("AccEffyhist", ";|y|^{#mu#mu}; Acc x Eff",nybins,ybins); AccEffyhist->Sumw2();

      TCanvas* c = new TCanvas("c","",1000,800);
      TH1F* Acc_pt_016_num = new TH1F ("Acc_pt_016_num", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins); Acc_pt_016_num->Sumw2();
      TH1F* Acc_pt_016_deno = new TH1F ("Acc_pt_016_deno", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);Acc_pt_016_deno->Sumw2();
      TH1F* Acc_pt_1624_num = new TH1F ("Acc_pt_1624_num", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);Acc_pt_1624_num->Sumw2();
      TH1F* Acc_pt_1624_deno = new TH1F ("Acc_pt_1624_deno", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);Acc_pt_1624_deno->Sumw2();
      TH1F* Acc_y_num = new TH1F ("Acc_y_num", ";|y|^{#mu#mu}; Acceptance", nybins, ybins);Acc_y_num->Sumw2();
      TH1F* Acc_y_deno = new TH1F ("Acc_y_deno", ";|y|^{#mu#mu}; Acceptance", nybins, ybins);Acc_y_deno->Sumw2();
      TH2F* Acc_ypt_num = new TH2F ("Acc_ypt_num", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins);Acc_ypt_num->Sumw2();
      TH2F* Acc_ypt_deno = new TH2F ("Acc_ypt_deno", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins);Acc_ypt_deno->Sumw2();
      
      TH1F* Eff_pt_016_num = new TH1F ("Eff_pt_016_num", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins); Eff_pt_016_num->Sumw2();
      TH1F* Eff_pt_016_deno = new TH1F ("Eff_pt_016_deno", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins); Eff_pt_016_deno->Sumw2();
      TH1F* Eff_pt_1624_num = new TH1F ("Eff_pt_1624_num", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);Eff_pt_1624_num->Sumw2();
      TH1F* Eff_pt_1624_deno = new TH1F ("Eff_pt_1624_deno", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);Eff_pt_1624_deno->Sumw2();
      TH1F* Eff_y_num = new TH1F ("Eff_y_num", ";|y|^{#mu#mu}; Efficiency", nybins, ybins);Eff_y_num->Sumw2();
      TH1F* Eff_y_deno = new TH1F ("Eff_y_deno", ";|y|^{#mu#mu}; Efficiency", nybins, ybins);Eff_y_deno->Sumw2();
      TH2F* Eff_ypt_num = new TH2F ("Eff_ypt_num", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins); Eff_ypt_num->Sumw2();
      TH2F* Eff_ypt_deno = new TH2F ("Eff_ypt_deno", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins); Eff_ypt_deno->Sumw2();
      
      TEfficiency* Acc_pt_016 = new TEfficiency ("Acc_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins); 
      TEfficiency* Acc_pt_1624 = new TEfficiency ("Acc_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TEfficiency* Acc_y = new TEfficiency ("Acc_y", ";|y|^{#mu#mu}; Acceptance", nybins, ybins);
      TEfficiency* Acc_ypt = new TEfficiency ("Acc_ypt", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins);
      
      TEfficiency* Eff_pt_016 = new TEfficiency ("Eff_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TEfficiency* Eff_pt_1624 = new TEfficiency ("Eff_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TEfficiency* Eff_y = new TEfficiency ("Eff_y", ";|y|^{#mu#mu}; Efficiency", nybins, ybins);
      TEfficiency* Eff_ypt = new TEfficiency ("Eff_ypt", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins);
      
      TEfficiency* AccEff_pt_016 = new TEfficiency ("AccEff_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff", nptbins, ptbins);
      TEfficiency* AccEff_pt_1624 = new TEfficiency ("AccEff_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff", nptbins, ptbins);
      TEfficiency* AccEff_y = new TEfficiency ("AccEff_y", ";|y|^{#mu#mu}; Acc x Eff", nybins, ybins);
      TEfficiency* AccEff_ypt = new TEfficiency ("AccEff_ypt", ";y^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", netabins, etabins, nptbins, ptbins);
      
      Long64_t nentries = fChain->GetEntries();
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
	      if ( abs(GenQQ4mom->Rapidity()) < 2.4 && abs (GenQQ4mom->Rapidity()) > 1.6 && GenQQ4mom->Pt() > 3)
		{
		  Acc_pt_1624_deno->Fill(GenQQ4mom->Pt());
		  if(areGenMuonsInAcceptance2015(iQQ))
		    {
		      Acc_pt_1624_num->Fill(GenQQ4mom->Pt());
		      Eff_pt_1624_deno->Fill(GenQQ4mom->Pt());
		    }
		}
	      if (abs(GenQQ4mom->Rapidity()) < 1.6 && GenQQ4mom->Pt() > 6.5)
		{
		  Acc_pt_016_deno->Fill(GenQQ4mom->Pt());
		  if(areGenMuonsInAcceptance2015(iQQ))
		    {
		      Acc_pt_016_num->Fill(GenQQ4mom->Pt());
		      Eff_pt_016_deno->Fill(GenQQ4mom->Pt());
		    }
		}
	    if (abs(GenQQ4mom->Rapidity()) < 2.4 && GenQQ4mom->Pt() > 6.5)
	      {
		Acc_y_deno->Fill(abs(GenQQ4mom->Rapidity()));
		if(areGenMuonsInAcceptance2015(iQQ))
		  {
		    Acc_y_num->Fill(abs(GenQQ4mom->Rapidity()));
                    Eff_y_deno->Fill(abs(GenQQ4mom->Rapidity()));
		  }
	      }
	    if (abs(GenQQ4mom->Rapidity()) < 2.4 && GenQQ4mom->Pt() > 3)
	      {
	       	Acc_ypt_deno->Fill(GenQQ4mom->Rapidity(), GenQQ4mom->Pt());
		if (areGenMuonsInAcceptance2015(iQQ))
		  {
		    Acc_ypt_num->Fill(GenQQ4mom->Rapidity(), GenQQ4mom->Pt());
		    Eff_ypt_deno->Fill(GenQQ4mom->Rapidity(), GenQQ4mom->Pt());
		  }
	      }
	  }
	  for(int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
	    {
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	      TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
	      TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
	      if(
		 (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts                                                                                                                
		 (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts                                                                                                                  
		 (isTriggerMatch(iQQ, triggerIndex_PP)) && // if it matches the trigger                                                                                                                
		 (isMatchedRecoDiMuon(iQQ))
		 )
		{
		    tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
		    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);

		  if (abs(RecoQQ4mom->Rapidity()) < 2.4 && abs (RecoQQ4mom->Rapidity()) > 1.6 && RecoQQ4mom->Pt() > 3)
		    Eff_pt_1624_num->Fill(RecoQQ4mom->Pt(), tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 1.6 && RecoQQ4mom->Pt() > 6.5)
		    Eff_pt_016_num->Fill(RecoQQ4mom->Pt(), tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 2.4 && RecoQQ4mom->Pt() > 6.5)
		    Eff_y_num->Fill(abs(RecoQQ4mom->Rapidity()),tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 2.4 && RecoQQ4mom->Pt() > 3)
		    Eff_ypt_num->Fill(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(),tnp_weight);
		}
	    }
	}

      Acc_pt_1624->SetStatisticOption(TEfficiency::kBBayesian);
      Acc_pt_1624->SetPassedHistogram(*Acc_pt_1624_num,"f");
      Acc_pt_1624->SetTotalHistogram(*Acc_pt_1624_deno,"f");
      
      Acc_pt_016->SetStatisticOption(TEfficiency::kBBayesian);
      Acc_pt_016->SetPassedHistogram(*Acc_pt_016_num,"f");
      Acc_pt_016->SetTotalHistogram(*Acc_pt_016_deno,"f");

      Acc_y->SetStatisticOption(TEfficiency::kBBayesian);
      Acc_y->SetPassedHistogram(*Acc_y_num,"f");
      Acc_y->SetTotalHistogram(*Acc_y_deno,"f");
      
      Acc_ypt->SetStatisticOption(TEfficiency::kBBayesian);
      Acc_ypt->SetPassedHistogram(*Acc_ypt_num,"f");
      Acc_ypt->SetTotalHistogram(*Acc_ypt_deno,"f");

      Eff_pt_1624->SetStatisticOption(TEfficiency::kBBayesian);
      Eff_pt_1624->SetPassedHistogram(*Eff_pt_1624_num,"f");
      Eff_pt_1624->SetTotalHistogram(*Eff_pt_1624_deno,"f");

      Eff_pt_016->SetStatisticOption(TEfficiency::kBBayesian);
      Eff_pt_016->SetPassedHistogram(*Eff_pt_016_num,"f");
      Eff_pt_016->SetTotalHistogram(*Eff_pt_016_deno,"f");
      
      Eff_y->SetStatisticOption(TEfficiency::kBBayesian);
      Eff_y->SetPassedHistogram(*Eff_y_num,"f");
      Eff_y->SetTotalHistogram(*Eff_y_deno,"f");

      Eff_ypt->SetStatisticOption(TEfficiency::kBBayesian);
      Eff_ypt->SetPassedHistogram(*Eff_ypt_num,"f");
      Eff_ypt->SetTotalHistogram(*Eff_ypt_deno,"f");

      AccEff_pt_1624->SetStatisticOption(TEfficiency::kBBayesian);
      AccEff_pt_1624->SetPassedHistogram(*Eff_pt_1624_num,"f");
      AccEff_pt_1624->SetTotalHistogram(*Acc_pt_1624_deno,"f");

      AccEff_pt_016->SetStatisticOption(TEfficiency::kBBayesian);
      AccEff_pt_016->SetPassedHistogram(*Eff_pt_016_num,"f");
      AccEff_pt_016->SetTotalHistogram(*Acc_pt_016_deno,"f");

      AccEff_y->SetStatisticOption(TEfficiency::kBBayesian);
      AccEff_y->SetPassedHistogram(*Eff_y_num,"f");
      AccEff_y->SetTotalHistogram(*Acc_y_deno,"f");

      AccEff_ypt->SetStatisticOption(TEfficiency::kBBayesian);
      AccEff_ypt->SetPassedHistogram(*Eff_ypt_num,"f");
      AccEff_ypt->SetTotalHistogram(*Acc_ypt_deno,"f");


      TFile* Corr = new TFile(Form("ANPlots/%sAccEffCorr.root", isPr?"pr":"npr"),"RECREATE");
      Corr->cd();
      Acc_pt_1624->Write("Acc_pt_1624");
      Acc_pt_016->Write("Acc_pt_016");
      Acc_y->Write("Acc_y");
      Acc_ypt->Write("Acc_ypt");
      Eff_pt_1624->Write("Eff_pt_1624");
      Eff_pt_016->Write("Eff_pt_016");
      Eff_y->Write("Eff_y");
      Eff_ypt->Write("Eff_ypt");
      AccEff_pt_1624->Write("AccEff_pt_1624");
      AccEff_pt_016->Write("AccEff_pt_016");
      AccEff_y->Write("AccEff_y");
      AccEff_ypt->Write("AccEff_ypt");
      Corr->Close();
      
      c->cd();
      gStyle->SetOptStat(0);
      //Acc_pt_1624->SetMinimum(0);
      //Acc_pt_1624->SetMaximum(1.2);
      //Acc_pt_016->SetMinimum(0);
      //Acc_pt_016->SetMaximum(1.2);
      //Acc_y->SetMinimum(0);
      //Acc_y->SetMaximum(1.2);
      //Eff_pt_1624->SetMinimum(0);
      //Eff_pt_1624->SetMaximum(1.2);
      //Eff_pt_016->SetMinimum(0);
      //Eff_pt_016->SetMaximum(1.2);
      //Eff_y->SetMinimum(0);
      //Eff_y->SetMaximum(1.2);

      Accpthist->SetMinimum(0);
      Accpthist->SetMaximum(1.2);
      Effpthist->SetMinimum(0);
      Effpthist->SetMaximum(1.2);
      Accyhist->SetMinimum(0);
      Accyhist->SetMaximum(1.2);
      Effyhist->SetMinimum(0);
      Effyhist->SetMaximum(1.2);
      AccEffpthist->SetMinimum(0);
      AccEffpthist->SetMaximum(1.2);
      AccEffyhist->SetMinimum(0);
      AccEffyhist->SetMaximum(1.2);

      TLine* pt1 = new TLine(3,1,50,1);
      pt1->SetLineColor(kRed);
      pt1->SetLineStyle(2);
      TLine* y1 = new TLine(0,1,2.4,1);
      y1->SetLineColor(kRed);
      y1->SetLineStyle(2);
      
      Acc_pt_1624->SetMarkerStyle(27);
      Acc_pt_1624->SetMarkerColor(40);
      Acc_pt_1624->SetMarkerSize(1);
      Acc_pt_1624->SetLineColor(39);
      Acc_pt_016->SetMarkerStyle(30);
      Acc_pt_016->SetMarkerColor(46);
      Acc_pt_016->SetMarkerSize(1);
      Acc_pt_016->SetLineColor(48);
      Acc_y->SetMarkerStyle(26);
      Acc_y->SetMarkerColor(8);
      Acc_y->SetMarkerSize(1);
      Acc_y->SetLineColor(29);

      Eff_pt_1624->SetMarkerStyle(27);
      Eff_pt_1624->SetMarkerColor(40);
      Eff_pt_1624->SetMarkerSize(1);
      Eff_pt_1624->SetLineColor(39);
      Eff_pt_016->SetMarkerStyle(30);
      Eff_pt_016->SetMarkerColor(46);
      Eff_pt_016->SetMarkerSize(1);
      Eff_pt_016->SetLineColor(48);
      Eff_y->SetMarkerStyle(26);
      Eff_y->SetMarkerColor(8);
      Eff_y->SetMarkerSize(1);
      Eff_y->SetLineColor(29);
      
      AccEff_pt_1624->SetMarkerStyle(27);
      AccEff_pt_1624->SetMarkerColor(40);
      AccEff_pt_1624->SetMarkerSize(1);
      AccEff_pt_1624->SetLineColor(39);
      AccEff_pt_016->SetMarkerStyle(30);
      AccEff_pt_016->SetMarkerColor(46);
      AccEff_pt_016->SetMarkerSize(1);
      AccEff_pt_016->SetLineColor(48);
      AccEff_y->SetMarkerStyle(26);
      AccEff_y->SetMarkerColor(8);
      AccEff_y->SetMarkerSize(1);
      AccEff_y->SetLineColor(29);

      TLegend *leg = new TLegend(0.6, 0.775, 0.8, 0.875);
      leg->AddEntry(Acc_pt_016, "|y| < 1.6", "lep");
      leg->AddEntry(Acc_pt_1624, "1.6 < |y| < 2.4", "lep");
      leg->SetBorderSize(1);
      leg->SetFillStyle(0);
      
      TPaveText *tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
      tbox->AddText(Form("%s J/#psi",isPr?"Prompt":"Nonprompt"));
      tbox->SetBorderSize(0);
      tbox->SetFillColor(0);
      tbox->SetFillStyle(0);

      Accpthist->Draw();    
      Acc_pt_1624->Draw("same");
      Acc_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccVsPt.pdf", isPr?"pr":"npr"));

      Effpthist->Draw();
      Eff_pt_1624->Draw("same");
      Eff_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sEffVsPt.pdf", isPr?"pr":"npr"));
      
      AccEffpthist->Draw();
      AccEff_pt_1624->Draw("same");
      AccEff_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccEffVsPt.pdf", isPr?"pr":"npr"));

      tbox->AddText("6.5 < p_{T} < 35 GeV/c");
      
      Accyhist->Draw();
      Acc_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccVsY.pdf", isPr?"pr":"npr"));

      Effyhist->Draw();
      Eff_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sEffVsY.pdf", isPr?"pr":"npr"));
      
      AccEffyhist->Draw();
      AccEff_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccEffVsY.pdf", isPr?"pr":"npr"));

      Acc_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sAcc2D.pdf", isPr?"pr":"npr"));
      Eff_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sEff2D.pdf", isPr?"pr":"npr"));
      AccEff_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sAccEff2D.pdf", isPr?"pr":"npr"));
    }
  else {
    cout << "[ERROR] this is data not MC" << endl;
  }
}

void myTree::ClosureTest()
{
  TH1F* rpt = new TH1F ("rpt", "pt distribution at reco level", 10, 6.5, 26.5);
  TH1F* gpt = new TH1F ("gpt", "pt distribution at gen level", 10, 6.5, 26.5);
  TFile* f (0x0);
  if (isPr)
    f= TFile::Open("prEff.root");
  else
    f= TFile::Open ("nprEff.root");
  TEfficiency* eff = (TEfficiency*) f->Get("ptrap");

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if (isMc)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m = GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  gpt->Fill(jpsi_pt);
		}
	    }
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  weight=1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			  rpt->Fill(jpsi_pt, weight);
			}
		    }
		}
	    }
	}
      TFile* testfile (0x0);
      if (isPr)
	testfile = new TFile ("prClosureTest.root","RECREATE");
      else
	testfile = new TFile ("nprClosureTest.root","RECREATE");
      rpt->Write("recopt");
      gpt->Write("genpt");
      testfile->Close();
    }
  else
    cout<< "this is data and not MC"<<endl;
}

void myTree::Loop() {cout << "[INFO] This function is empty at the moment!!"<< endl;}

void myTree::Plot() {cout << "[INFO] This function is empty at the moment!!"<< endl;}

Bool_t myTree::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit) 
  {
    Bool_t cond = true;
    cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) ); 
    cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    return cond;
  };

Bool_t myTree::isGlobalMuonInAccept2015 (TLorentzVector* Muon) 
  {
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
  };

Bool_t myTree::areMuonsInAcceptance2015 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
  };

Bool_t myTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
  {
    TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenQQ);
    TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenQQ);
    return (isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi));
  };  
  
Bool_t myTree::passQualityCuts2015 (Int_t iRecoQQ) 
  {
    Bool_t cond = true;
    cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
    cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));    
    // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);
    
    cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
    cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));
    // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  };


Double_t myTree::deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

Bool_t myTree::isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
{
  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);
  
  bool isMatched(false);
  int iGenMuon(0);
  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
  {
    TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
    TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) 
      {
	isMatched = true;
	matchGen = (TLorentzVector*) Gen_QQ_4mom->At(iGenMuon);
      }
    iGenMuon++;
  }
  return isMatched;
};



Bool_t myTree::isMatchedGenDiMuon(int iGenDiMuon, double maxDeltaR)
{
  TLorentzVector* GenMuonpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuonmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuon = (TLorentzVector*) Gen_QQ_4mom->At(iGenDiMuon);

  bool isMatched(false);
  int iRecoMuon(0);
  while ( !isMatched && (iRecoMuon < Reco_QQ_size) )
  {
    TLorentzVector *RecoMuonpl = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRecoMuon);
    TLorentzVector *RecoMuonmi = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRecoMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR) ) 
      {
	isMatched = true;
	matchReco = (TLorentzVector*) Reco_QQ_4mom->At(iRecoMuon);
      }
    iRecoMuon++;
  }
  
  return isMatched;
};

void myTree::Unfolding() {cout << "[INFO] This function is empty at the moment!!"<< endl;}

void myTree::JetPtRange()
{ 
  gSystem->mkdir("JtPtRange");
  gSystem->cd("JtPtRange");
  int imatch=0;
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "", 1000,800);
  TH1F* genptdist = new TH1F ("genptdist", ";p_{T}(jet) (GeV/c)", 100, 10, 50);
  TH1F* recoptdist = new TH1F ("recoptdist", ";p_{T}(jet) (GeV/c)", 100, 10, 50);
  TH1F* ptres = new TH1F ("ptres", ";p_{T}^{gen}-p_{T}^{reco}(jet) (GeV/c)", 40, -20, 20);
  Long64_t nentries = 1000000; //fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;  

      for (Long64_t ijet=0; ijet<nref; ijet++)
	{
	  
	  if (abs(jty[ijet])<1.6 /*&& abs(jty[ijet])<2.1*/ && jtpt[ijet]>25 && jtpt[ijet]<40)
	    {
	      recoptdist->Fill(jtpt[imatch]);
	      genptdist->Fill(refpt[imatch]);
	      ptres->Fill(refpt[imatch]-jtpt[imatch]); 
	    }
	}
    }

  TString filename;
  if (isPr)
    filename ="prJtPt2540_mid.dat";
  else
    filename ="nprJtPt2540_mid.dat";

  ofstream file_out(filename);

  for (int i=0; i<genptdist->GetSize(); i++)
    {
      file_out<< "for min pt(jet) = "<<genptdist->GetXaxis()->GetBinLowEdge(i)<< " Ngen/Nreco = "<<genptdist->Integral(i, genptdist->GetSize()-1)*1.0/recoptdist->Integral()<<endl;
      cout<< "for min pt(jet) = "<<genptdist->GetXaxis()->GetBinLowEdge(i)<< " Ngen/Nreco = "<<genptdist->Integral(i, genptdist->GetSize()-1)*1.0/recoptdist->Integral()<<endl;
    }
  file_out.close();
  c->cd();
  ptres->Draw();
  if (isPr)
    c->SaveAs("prJtPtRes2540_mid.pdf");
  else
    c->SaveAs("nprJtPtRes2540_mid.pdf");
  recoptdist->SetLineColor(4);
  recoptdist->Draw();
  genptdist->SetLineColor(2);
  genptdist->Draw("same");
  TLegend* l = new TLegend (0.65,0.7,0.85,0.8);
  l->AddEntry(genptdist, "gen", "lp");
  l->AddEntry(recoptdist, "reco", "lp");
  l->SetBorderSize(0);
  l->Draw("same");

  TPaveText* tbox0 = new TPaveText(0.6,0.5,0.8,0.7, "BRNDC");
  if (isPr)
    tbox0->AddText("prompt MC");
  else
    tbox0->AddText("nonprompt MC");
  tbox0->AddText("25 < p_{T}^{reco}(jet) < 40 GeV/c");
  tbox0->AddText("|y(jet)| < 1.6");
  tbox0->SetBorderSize(0);
  tbox0->SetFillColor(0);
  tbox0->Draw("same");
  if (isPr)
    c->SaveAs("prJtPtDist2540_mid.pdf");
  else
    c->SaveAs("nprJtPtDist2540_mid.pdf");
  gSystem->cd("..");
}
