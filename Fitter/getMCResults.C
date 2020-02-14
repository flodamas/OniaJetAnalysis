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

void GetMCResults (bool isPbPb=false,bool isPr=true, bool underflowOff=true, bool plotBpt=false, bool plotxi=false);
Bool_t isGlobalMuonInAccept2019 (TLorentzVector* Muon);

void GetMCResults_all ()
{
  //GetMCResults(true, true, false);
  //GetMCResults(true, false, false);
  GetMCResults(false, true,  true, false);
  //GetMCResults(false, false,  true, false);
  //GetMCResults(false,  false, false);

  //GetMCResults(false, true, true);
  //GetMCResults(false, false, true);
}

void GetMCResults (bool isPbPb, bool isPr, bool underflowOff, bool plotBpt, bool plotXi)
{
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;

  double sefer=0;

  //TFile *treeFile = TFile::Open(Form("/data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen"));
  //cout << Form("[INFO] Reading tree from file /data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen")<<endl;
  TFile *treeFile = TFile::Open(Form("/data_CMS/cms/diab/JpsiJet/MC/pp/%s/HiForestAOD%s_merged.root", isPr?"prompt/v3":"nonprompt/v4",isPr?"_ext":""));

  TTree* oniaTree = (TTree*) treeFile->Get("hionia/myTree");
  TTree* jetTree = (TTree*) treeFile->Get("ak4PFJetAnalyzer/t");
  TTree* bTree(0x0);
  bTree = (TTree*) treeFile->Get("bHadronAna/hi");
  if (!bTree) cout<<"b tree not found"<<endl;
  oniaTree->AddFriend(jetTree);
  if (plotBpt)  oniaTree->AddFriend(bTree);

  TH1F* zMidHist = new TH1F ("zMidHist","", 7, 0.02, 1); zMidHist->Sumw2();
  TH1F* zFwdHist = new TH1F ("zFwdHist","", 5, 0, 1); zFwdHist->Sumw2();
  TH1F* nMidHist = new TH1F ("nMidHist","", 10, 0, 1); nMidHist->Sumw2();
  TH1F* nFwdHist = new TH1F ("nFwdHist","", 10, 0, 1); nFwdHist->Sumw2();

  cout <<"[INFO] Importing AccFiles to correct"<<endl;
  TFile *corrFile = TFile::Open(Form("../Efficiency/FilesAccxEff/Acc/%sAccHists_%s.root",isPr?"pr":"npr",isPbPb?"PbPb":"PP"));
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
  //vector<vector<float> >    *jtbHadronPt; //
  std::vector<float>   *pt =0;
  std::vector<float>   *eta =0;
  std::vector<float>   *phi = 0;
  Int_t           mult;


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
  //TBranch        *b_jtbHadronPt;
  TBranch        *b_mult;
  TBranch        *b_pt = 0;
  TBranch        *b_eta = 0;
  TBranch        *b_phi = 0;

  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;
  Gen_mu_4mom = 0;
  //jtbHadronPt = 0;

  if (!oniaTree) { cout<<"[ERROR] no tree found"<<endl; return;}

  fChain = oniaTree;
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  //if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  //if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
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
  if (fChain->GetBranch("mult")) fChain->SetBranchAddress("mult", &mult, &b_mult);
  //if (fChain->GetBranch("pt")) fChain->SetBranchAddress("pt", &pt, &b_pt);
  //if (fChain->GetBranch("eta")) fChain->SetBranchAddress("eta", &eta, &b_eta);
  //if (fChain->GetBranch("phi")) fChain->SetBranchAddress("phi", &phi, &b_phi);
  cout<<"[INFO] all branch addresses set"<<endl;
  /*
  fChain->SetBranchStatus("*",0);

  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchStatus("Gen_QQ_size",1);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchStatus("Gen_QQ_4mom",1);
  if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1);
  if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1);
  if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchStatus("Gen_mu_4mom",1); 
  if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchStatus("Gen_QQ_mupl_idx",1);
  if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchStatus("Gen_QQ_mumi_idx",1);
  if (fChain->GetBranch("Gen_weight")) fChain->SetBranchStatus("Gen_weight",1);
  if (fChain->GetBranch("ngen")) fChain->SetBranchStatus("ngen",1);
  if (fChain->GetBranch("genpt")) fChain->SetBranchStatus("genpt",1);
  if (fChain->GetBranch("geneta")) fChain->SetBranchStatus("geneta",1);
  if (fChain->GetBranch("geny")) fChain->SetBranchStatus("geny",1);
  if (fChain->GetBranch("genphi")) fChain->SetBranchStatus("genphi",1);
  if (fChain->GetBranch("genm")) fChain->SetBranchStatus("genm",1);
  if (fChain->GetBranch("mult")) fChain->SetBranchStatus("mult",1);
  if (fChain->GetBranch("pt")) fChain->SetBranchStatus("pt",1);
  if (fChain->GetBranch("eta")) fChain->SetBranchStatus("eta",1);
  if (fChain->GetBranch("phi")) fChain->SetBranchStatus("phi",1);
  */

  //get the histograms for the Fonll weights

  TFile *fFonll = new TFile("/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/acc/fFONLL.root");
  TGraph *gFonllMid = (TGraph*)fFonll->Get("gFonllMid");
  TGraph *gFonllFwd = (TGraph*)fFonll->Get("gFonllFwd");


  TFile *fPythia = new TFile("/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/acc/fPythia.root");
  TH1F  *hPythiaMid = (TH1F *)fPythia->Get("hPythiaMid");
  TH1F  *hPythiaFwd = (TH1F *)fPythia->Get("hPythiaFwd");

  Long64_t nentries =fChain->GetEntries();
  //nentries = 1000000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<Form(" for (%s, %s, %s(%s))",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotXi?"xi":"z",plotBpt?"B":"Jpsi")<<endl;
      nb = fChain->GetEntry(jentry);   
      nbytes += fChain->GetEntry(jentry);
      fChain->GetEntry(jentry);
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  zed = -1;
	  drmin = 0.5;
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iQQ]);
	  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iQQ]);

	  if (GenQQ4mom->M()<2.6 || GenQQ4mom->M()>3.5) continue;
	  if (GenQQ4mom->Pt() < 3 || GenQQ4mom->Pt() > 35) continue;
	  if (abs(GenQQ4mom->Rapidity())>2.4) continue;
	  if (!isGlobalMuonInAccept2019(GenQQmupl) || !isGlobalMuonInAccept2019(GenQQmumi)) continue;

	  accWeight = 1.0/(accCorr->GetEfficiency(accCorr->FindFixBin(GenQQ4mom->Rapidity(), GenQQ4mom->Pt())));
	  if (isPr) Gen_weight=1;
	  accWeight = accWeight * Gen_weight;

	  float jPsiPt = GenQQ4mom->Pt();
	  float fonllXs = gFonllMid->Eval(jPsiPt);
	  if (abs(GenQQ4mom->Rapidity())>1.6) fonllXs = gFonllFwd->Eval(jPsiPt);

	  int iBin = hPythiaMid->FindBin(jPsiPt);
	  if (abs(GenQQ4mom->Rapidity())>1.6) iBin = hPythiaFwd->FindBin(jPsiPt);

	  float binCent = hPythiaMid->GetBinCenter(iBin);
	  if (abs(GenQQ4mom->Rapidity())>1.6) binCent = hPythiaFwd->GetBinCenter(iBin);

	  float pythiaXs = 0.;

	  if(iBin == 1 || jPsiPt > binCent)
	    {
	      float x1 = binCent;
	      float x2 = hPythiaMid->GetBinCenter(iBin+1);
	      float y1 = hPythiaMid->GetBinContent(iBin);
	      float y2 = hPythiaMid->GetBinContent(iBin+1);
	      if (abs(GenQQ4mom->Rapidity())>1.6) {
		x2 = hPythiaFwd->GetBinCenter(iBin+1);
		y1 = hPythiaFwd->GetBinContent(iBin);
		y2 = hPythiaFwd->GetBinContent(iBin+1);
	      }
	      float m = (y2-y1)/(x2-x1);
	      float b = y1 - m*x1;
	      pythiaXs = m*jPsiPt + b;
	    }
	  else{
	    float x1 = hPythiaMid->GetBinCenter(iBin-1);
	    if (abs(GenQQ4mom->Rapidity())>1.6) x1 = hPythiaFwd->GetBinCenter(iBin-1);
	    float x2 = binCent;
	    float y1 = hPythiaMid->GetBinContent(iBin-1);
	    float y2 = hPythiaMid->GetBinContent(iBin);
	    if (abs(GenQQ4mom->Rapidity())>1.6) {
	      y1 = hPythiaFwd->GetBinContent(iBin-1);
	      y2 = hPythiaFwd->GetBinContent(iBin);
	    }
	    float m = (y2-y1)/(x2-x1);
	    float b = y1 - m*x1;
	    pythiaXs = m*jPsiPt + b;
	  }

	  //if (pythiaXs==0)
	  //cout<<"[INFO] pythiaXS =0"<<endl;
	  //if (!isPr && pythiaXs>0)
	  //accWeight = accWeight*fonllXs/pythiaXs;

	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    nFwdHist->Fill(sefer, accWeight);
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt() > 6.5)
	    nMidHist->Fill(sefer, accWeight);

	  for (int iJet=0; iJet<ngen; iJet++) {
	    if (genpt[iJet]<25 || genpt[iJet]>35) continue;
	    if (abs(geny[iJet])>2.4) continue;
	    
	    TLorentzVector v_jet;
	    v_jet.SetPtEtaPhiM(genpt[iJet], geneta[iJet], genphi[iJet], genm[iJet]);
	    if (GenQQ4mom->DeltaR (v_jet)<=drmin) {
	      drmin = GenQQ4mom->DeltaR (v_jet);
	      if (!plotBpt)
		zed = GenQQ4mom->Pt()*1.0/genpt[iJet];
	      if (plotBpt) {
		int ibestB = -1;
		float drminb = 999;
		
		for (int ib = 0; ib<mult; ib++) {
		  float deta = eta->at(ib);
		  deta = eta->at(ib) - geneta[iJet];
		  float dphi = phi->at(ib) - genphi[iJet];
		  dphi = acos(cos(dphi));
		  float dR = sqrt(deta*deta + dphi*dphi);
		  if (dR < drminb) {
		    drminb = dR;
		    ibestB = ib;
		  }
		}// end of b loop 
		if (ibestB > 0)
		  zed = pt->at(ibestB)*1.0/genpt[iJet];
		if (ibestB > 0 && drminb>0.4) cout <<"[WARNING] Bbig dR(b-jet) = "<<drminb<<"for pt(b) = "<<pt->at(ibestB)<<endl;
	      }// end of b cond
	    }	      
	  }// end of gen jet loop
	  if (zed == -1) continue;
	  if (zed > 1 && zed <= 1.000001) zed = 0.9999999;
	  
	  if (plotXi) zed = log(1.0/zed);
	  
	  if (underflowOff  && zed<0.2) continue;
	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    zFwdHist->Fill(zed, accWeight);
	  if (underflowOff  && zed<0.44) continue;
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt()>6.5)
	    zMidHist->Fill(zed, accWeight);
	} //end of genQQ loop 
    }//end of events loop
  //if (underflowOff)
  //{
  //zMidHist->Scale(1.0/zMidHist->Integral("width"));
  //zFwdHist->Scale(1.0/zFwdHist->Integral("width"));
  //}
  //zHist->Scale(1.0/0.2);
  gSystem->mkdir("Output/MCResults");
  //gSystem->mkdir("Output/MCResults/fonllCorr");
  //TFile *fsave = new TFile(Form("Output/MCResults/fonllCorr/mcResult_%s_%s%s%s.root",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":"", plotXi?"_xi":""),"RECREATE");
  TFile *fsave = new TFile(Form("Output/MCResults/mcResult_%s_%s_%s%s%s.root",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":"", plotXi?"_xi":""),"RECREATE");
  zMidHist->Write("zDist_mid");
  nMidHist->Write("Ntot_mid");
  zFwdHist->Write("zDist_fwd");
  nFwdHist->Write("Ntot_fwd");
  fsave->Close();
}

Bool_t isGlobalMuonInAccept2019 (TLorentzVector* Muon)
{
  return (fabs(Muon->Eta()) < 2.4 &&
	  ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	   (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
	   (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
}


void drawPrBPlot() {
  gStyle->SetOptStat(0);
  TFile* prFile = TFile::Open("Output/MCResults/mcResult_prompt_all.root");
  TH1F* midJpsiHist = (TH1F*) prFile->Get("zDist_mid");
  TH1F* fwdJpsiHist = (TH1F*) prFile->Get("zDist_fwd");
  TFile* nprFile = TFile::Open("Output/MCResults/mcResult_nonprompt_all_bHadronPt.root");
  TH1F* midBHist = (TH1F*) nprFile->Get("zDist_mid");
  TH1F* fwdBHist = (TH1F*) nprFile->Get("zDist_fwd");

  midJpsiHist->Scale(1.0/midJpsiHist->Integral("width"));
  fwdJpsiHist->Scale(1.0/fwdJpsiHist->Integral("width"));
  midBHist->Scale(1.0/midBHist->Integral("width"));
  fwdBHist->Scale(1.0/fwdBHist->Integral("width"));

  midJpsiHist->SetMarkerColor(kOrange+2);
  midJpsiHist->SetMarkerStyle(kFullCircle);
  midJpsiHist->SetMarkerSize(1);
  midJpsiHist->SetLineColor(kOrange+2);
  midJpsiHist->SetLineWidth(2);

  fwdJpsiHist->SetMarkerColor(kOrange+2);
  fwdJpsiHist->SetMarkerStyle(kFullCircle);
  fwdJpsiHist->SetMarkerSize(1);
  fwdJpsiHist->SetLineColor(kOrange+2);
  fwdJpsiHist->SetLineWidth(2);

  midBHist->SetMarkerColor(kBlue);
  midBHist->SetMarkerStyle(kFullCircle);
  midBHist->SetMarkerSize(1);
  midBHist->SetLineColor(kBlue);
  midBHist->SetLineWidth(2);

  fwdBHist->SetMarkerColor(kBlue);
  fwdBHist->SetMarkerStyle(kFullCircle);
  fwdBHist->SetMarkerSize(1);
  fwdBHist->SetLineColor(kBlue);
  fwdBHist->SetLineWidth(2);

  double ybins [] = {0, 1.6, 2.4};
  TH1F *axisHist = new TH1F("axisHist","", 5, 0, 1);
  axisHist->GetYaxis()->SetRangeUser(0, 7.5);
  axisHist->GetYaxis()->SetTitle("1/N dN/dz");
  axisHist->GetXaxis()->SetTitle("z");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  TLegend* leg = new TLegend(0.62,0.62,0.88,0.80);

  leg->AddEntry(midJpsiHist, "b hadrons", "lp");
  leg->AddEntry(midBHist, "prompt J/#psi", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *  text2 = new TLatex(0.160 ,0.67, "6.5 < p_{T} < 35 GeV");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.03);
  text2->SetLineWidth(2);

  TLatex *  text21 = new TLatex(0.160 ,0.72, "|y| < 1.6");
  text21->SetNDC();
  text21->SetTextFont(42);
  text21->SetTextSize(0.03);
  text21->SetLineWidth(2);

  TLatex *  text3 = new TLatex(0.160 ,0.67, "3 < p_{T} < 35 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.03);
  text3->SetLineWidth(2);

  TLatex *  text31 = new TLatex(0.160 ,0.72, "|y| < 1.6");
  text31->SetNDC();
  text31->SetTextFont(42);
  text31->SetTextSize(0.03);
  text31->SetLineWidth(2);

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

  TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.04);
  text4->SetLineWidth(2);


  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  axisHist->Draw();
  midJpsiHist->Draw("same");
  midBHist->Draw("same");
  leg->Draw("same");
  text2->Draw("same");
  text21->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  gSystem->mkdir("Output/MCResults");
  c->SaveAs("Output/MCResults/midJpsiBComp.pdf");

  axisHist->Draw();
  fwdJpsiHist->Draw("same");
  fwdBHist->Draw("same");
  leg->Draw("same");
  text3->Draw("same");
  text31->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  gSystem->mkdir("Output/MCResults");
  c->SaveAs("Output/MCResults/fwdJpsiBComp.pdf");
}

void drawXiPlot() {
  gStyle->SetOptStat(0);
  TFile* prFile = TFile::Open("Output/MCResults/mcResult_prompt_all_xi.root");
  TH1F* midprHist = (TH1F*) prFile->Get("zDist_mid");
  TH1F* fwdprHist = (TH1F*) prFile->Get("zDist_fwd");
  TFile* nprFile = TFile::Open("Output/MCResults/mcResult_nonprompt_all_xi.root");
  TH1F* midnprHist = (TH1F*) nprFile->Get("zDist_mid");
  TH1F* fwdnprHist = (TH1F*) nprFile->Get("zDist_fwd");

  midprHist->Scale(1.0/midprHist->Integral("width"));
  fwdprHist->Scale(1.0/fwdprHist->Integral("width"));
  midnprHist->Scale(1.0/midnprHist->Integral("width"));
  fwdnprHist->Scale(1.0/fwdnprHist->Integral("width"));

  midprHist->SetMarkerColor(46);
  midprHist->SetMarkerStyle(kFullCircle);
  midprHist->SetMarkerSize(1);
  midprHist->SetLineColor(46);
  midprHist->SetLineWidth(2);

  fwdprHist->SetMarkerColor(46);
  fwdprHist->SetMarkerStyle(kFullCircle);
  fwdprHist->SetMarkerSize(1);
  fwdprHist->SetLineColor(46);
  fwdprHist->SetLineWidth(2);

  midnprHist->SetMarkerColor(30);
  midnprHist->SetMarkerStyle(kFullCircle);
  midnprHist->SetMarkerSize(1);
  midnprHist->SetLineColor(30);
  midnprHist->SetLineWidth(2);

  fwdnprHist->SetMarkerColor(30);
  fwdnprHist->SetMarkerStyle(kFullCircle);
  fwdnprHist->SetMarkerSize(1);
  fwdnprHist->SetLineColor(30);
  fwdnprHist->SetLineWidth(2);

  double ybins [] = {0, 1.6, 2.4};
  TH1F *axisHist = new TH1F("axisHist","", 5, 0, 1);
  axisHist->GetYaxis()->SetRangeUser(0, 5.0);
  axisHist->GetYaxis()->SetTitle("1/N dN/d#xi");
  axisHist->GetXaxis()->SetTitle("#xi(J/#psi)");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetYaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  TLegend* leg = new TLegend(0.62,0.62,0.88,0.80);

  leg->AddEntry(midprHist, "nonprompt", "lp");
  leg->AddEntry(midnprHist, "prompt", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *  text2 = new TLatex(0.160 ,0.67, "6.5 < p_{T} < 35 GeV");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.03);
  text2->SetLineWidth(2);

  TLatex *  text21 = new TLatex(0.160 ,0.72, "|y| < 1.6");
  text21->SetNDC();
  text21->SetTextFont(42);
  text21->SetTextSize(0.03);
  text21->SetLineWidth(2);

  TLatex *  text3 = new TLatex(0.160 ,0.67, "3 < p_{T} < 35 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.03);
  text3->SetLineWidth(2);

  TLatex *  text31 = new TLatex(0.160 ,0.72, "|y| < 1.6");
  text31->SetNDC();
  text31->SetTextFont(42);
  text31->SetTextSize(0.03);
  text31->SetLineWidth(2);

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

  TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.04);
  text4->SetLineWidth(2);


  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  axisHist->Draw();
  midprHist->Draw("same");
  midnprHist->Draw("same");
  leg->Draw("same");
  text2->Draw("same");
  text21->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  gSystem->mkdir("Output/MCResults");
  c->SaveAs("Output/MCResults/midXiPlot.pdf");

  axisHist->Draw();
  fwdprHist->Draw("same");
  fwdnprHist->Draw("same");
  leg->Draw("same");
  text3->Draw("same");
  text31->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  gSystem->mkdir("Output/MCResults");
  c->SaveAs("Output/MCResults/fwdXiPlot.pdf");
}


void GetJES (bool isPbPb, bool isPr)
{
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;

  double sefer=0;

  //TFile *treeFile = TFile::Open(Form("/data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen"));
  //cout << Form("[INFO] Reading tree from file /data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen")<<endl;
  TFile *treeFile = TFile::Open(Form("/data_CMS/cms/diab/JpsiJet/MC/pp/%s/HiForestAOD%s_merged.root", isPr?"prompt/v3":"nonprompt/v4",isPr?"_ext":""));

  TTree* oniaTree = (TTree*) treeFile->Get("hionia/myTree");
  TTree* jetTree = (TTree*) treeFile->Get("ak4PFJetAnalyzer/t");

  oniaTree->AddFriend(jetTree);

  TH1F* zMidHist = new TH1F ("zMidHist","", 7, 0.02, 1); zMidHist->Sumw2();
  TH1F* zFwdHist = new TH1F ("zFwdHist","", 5, 0, 1); zFwdHist->Sumw2();
  TH1F* nMidHist = new TH1F ("nMidHist","", 10, 0, 1); nMidHist->Sumw2();
  TH1F* nFwdHist = new TH1F ("nFwdHist","", 10, 0, 1); nFwdHist->Sumw2();

  
  cout <<"[INFO] Importing AccFiles to correct"<<endl;
  TFile *corrFile = TFile::Open("Input/correction_AccEff_centMaps_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights.root");

  TEfficiency* accCorr = (TEfficiency*) corrFile->Get(Form("hcorr_Jpsi_PP_%s_Acc",isPr?"pr":"npr"));
  
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
  //jtbHadronPt = 0;

  if (!oniaTree) { cout<<"[ERROR] no tree found"<<endl; return;}

  fChain = oniaTree;
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  //if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  //if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
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

  cout<<"[INFO] all branch addresses set"<<endl;

  /*
  fChain->SetBranchStatus("*",0);
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchStatus("Gen_QQ_size",1);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchStatus("Gen_QQ_4mom",1);
  if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1);
  if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1);
  if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchStatus("Gen_mu_4mom",1); 
  if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchStatus("Gen_QQ_mupl_idx",1);
  if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchStatus("Gen_QQ_mumi_idx",1);
  if (fChain->GetBranch("Gen_weight")) fChain->SetBranchStatus("Gen_weight",1);
  if (fChain->GetBranch("ngen")) fChain->SetBranchStatus("ngen",1);
  if (fChain->GetBranch("genpt")) fChain->SetBranchStatus("genpt",1);
  if (fChain->GetBranch("geneta")) fChain->SetBranchStatus("geneta",1);
  if (fChain->GetBranch("geny")) fChain->SetBranchStatus("geny",1);
  if (fChain->GetBranch("genphi")) fChain->SetBranchStatus("genphi",1);
  if (fChain->GetBranch("genm")) fChain->SetBranchStatus("genm",1);
  */


  Long64_t nentries =fChain->GetEntries();
  //nentries = 1000000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<Form(" for %s",isPr?"prompt":"nonprompt")<<endl;
      nb = fChain->GetEntry(jentry); 
      nbytes += fChain->GetEntry(jentry);
      fChain->GetEntry(jentry);
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  zed = -1;
	  drmin = 0.5;
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iQQ]);
	  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iQQ]);

	  if (GenQQ4mom->M()<2.6 || GenQQ4mom->M()>3.5) continue;
	  if (GenQQ4mom->Pt() < 3 || GenQQ4mom->Pt() > 100) continue;
	  if (abs(GenQQ4mom->Rapidity())>2.4) continue;
	  if (!isGlobalMuonInAccept2019(GenQQmupl) || !isGlobalMuonInAccept2019(GenQQmumi)) continue;

	  accWeight = 1.0/(accCorr->GetEfficiency(accCorr->FindFixBin(GenQQ4mom->Rapidity(), GenQQ4mom->Pt())));
	  if (isPr) Gen_weight=1;
	  accWeight = accWeight * Gen_weight;

	  float jPsiPt = GenQQ4mom->Pt();


	  if(iBin == 1 || jPsiPt > binCent)
	    {
	      float x1 = binCent;
	      float x2 = hPythiaMid->GetBinCenter(iBin+1);
	      float y1 = hPythiaMid->GetBinContent(iBin);
	      float y2 = hPythiaMid->GetBinContent(iBin+1);
	      if (abs(GenQQ4mom->Rapidity())>1.6) {
		x2 = hPythiaFwd->GetBinCenter(iBin+1);
		y1 = hPythiaFwd->GetBinContent(iBin);
		y2 = hPythiaFwd->GetBinContent(iBin+1);
	      }
	      float m = (y2-y1)/(x2-x1);
	      float b = y1 - m*x1;
	      pythiaXs = m*jPsiPt + b;
	    }
	  else{
	    float x1 = hPythiaMid->GetBinCenter(iBin-1);
	    if (abs(GenQQ4mom->Rapidity())>1.6) x1 = hPythiaFwd->GetBinCenter(iBin-1);
	    float x2 = binCent;
	    float y1 = hPythiaMid->GetBinContent(iBin-1);
	    float y2 = hPythiaMid->GetBinContent(iBin);
	    if (abs(GenQQ4mom->Rapidity())>1.6) {
	      y1 = hPythiaFwd->GetBinContent(iBin-1);
	      y2 = hPythiaFwd->GetBinContent(iBin);
	    }
	    float m = (y2-y1)/(x2-x1);
	    float b = y1 - m*x1;
	    pythiaXs = m*jPsiPt + b;
	  }

	  //if (pythiaXs==0)
	  //cout<<"[INFO] pythiaXS =0"<<endl;
	  //if (!isPr && pythiaXs>0)
	  //accWeight = accWeight*fonllXs/pythiaXs;

	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    nFwdHist->Fill(sefer, accWeight);
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt() > 6.5)
	    nMidHist->Fill(sefer, accWeight);

	  for (int iJet=0; iJet<ngen; iJet++) {
	    if (genpt[iJet]<25 || genpt[iJet]>35) continue;
	    if (abs(geny[iJet])>2.4) continue;
	    
	    TLorentzVector v_jet;
	    v_jet.SetPtEtaPhiM(genpt[iJet], geneta[iJet], genphi[iJet], genm[iJet]);
	    if (GenQQ4mom->DeltaR (v_jet)<=drmin) {
	      drmin = GenQQ4mom->DeltaR (v_jet);
	      if (!plotBpt)
		zed = GenQQ4mom->Pt()*1.0/genpt[iJet];
	      if (plotBpt) {
		int ibestB = -1;
		float drminb = 999;
		
		for (int ib = 0; ib<mult; ib++) {
		  float deta = eta->at(ib);
		  deta = eta->at(ib) - geneta[iJet];
		  float dphi = phi->at(ib) - genphi[iJet];
		  dphi = acos(cos(dphi));
		  float dR = sqrt(deta*deta + dphi*dphi);
		  if (dR < drminb) {
		    drminb = dR;
		    ibestB = ib;
		  }
		}// end of b loop 
		if (ibestB > 0)
		  zed = pt->at(ibestB)*1.0/genpt[iJet];
		if (ibestB > 0 && drminb>0.4) cout <<"[WARNING] Bbig dR(b-jet) = "<<drminb<<"for pt(b) = "<<pt->at(ibestB)<<endl;
	      }// end of b cond
	    }	      
	  }// end of gen jet loop
	  if (zed == -1) continue;
	  if (zed > 1 && zed <= 1.000001) zed = 0.9999999;
	  
	  if (plotXi) zed = log(1.0/zed);
	  
	  if (underflowOff  && zed<0.2) continue;
	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    zFwdHist->Fill(zed, accWeight);
	  if (underflowOff  && zed<0.44) continue;
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt()>6.5)
	    zMidHist->Fill(zed, accWeight);
	} //end of genQQ loop 
    }//end of events loop
  //if (underflowOff)
  //{
  //zMidHist->Scale(1.0/zMidHist->Integral("width"));
  //zFwdHist->Scale(1.0/zFwdHist->Integral("width"));
  //}
  //zHist->Scale(1.0/0.2);
  gSystem->mkdir("Output/MCResults");
  //gSystem->mkdir("Output/MCResults/fonllCorr");
  //TFile *fsave = new TFile(Form("Output/MCResults/fonllCorr/mcResult_%s_%s%s%s.root",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":"", plotXi?"_xi":""),"RECREATE");
  TFile *fsave = new TFile(Form("Output/MCResults/mcResult_%s_%s_%s%s%s.root",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":"", plotXi?"_xi":""),"RECREATE");
  zMidHist->Write("zDist_mid");
  nMidHist->Write("Ntot_mid");
  zFwdHist->Write("zDist_fwd");
  nFwdHist->Write("Ntot_fwd");
  fsave->Close();
}
