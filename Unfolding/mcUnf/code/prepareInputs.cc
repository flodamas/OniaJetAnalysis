#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TBranchElement.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "THnSparse.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;

void prepare(bool doPrompt = false, bool doTrain=true, bool doMid = true, Int_t stepNumber = 1){

  string filename = "";
  string outputfile = "";
  string filenamePrevStep = "";
  
  if(doPrompt && doTrain){
    filename = "../../TreesForUnfolding_corrAccEff/tree_MCJPSIPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_midRapidity_Train_49z15ptBins7zMeasBins.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_fwdRapidity_Train_50z15ptBins.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins.root",stepNumber-1);
  }
  if(doPrompt && !doTrain){
    filename = "../../TreesForUnfolding_corrAccEff/tree_MCJPSIPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_midRapidity_Test_49z15ptBins7zMeasBins.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_prompt_fwdRapidity_Test_50z15ptBins.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins.root",stepNumber-1);
  }
  
  if(!doPrompt && doTrain){
    filename = "../../TreesForUnfolding_corrAccEff/tree_MCJPSINOPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_midRapidity_Train_49z15ptBins7zMeasBins.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_fwdRapidity_Train_50z15ptBins.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins.root",stepNumber-1);
  }
  if(!doPrompt && !doTrain){
    filename = "../../TreesForUnfolding_corrAccEff/tree_MCJPSINOPR_PP_NoBkg_AccEff_JEC.root";
    
    if(doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_midRapidity_Test_49z15ptBins7zMeasBins.root",stepNumber);
    if(!doMid) outputfile = Form("../unfInput/step%i/unfolding_4D_nonprompt_fwdRapidity_Test_50z15ptBins.root",stepNumber);
    
    if(stepNumber > 1 && doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins.root",stepNumber-1);
    if(stepNumber > 1 && !doMid) filenamePrevStep = Form("../unfOutput/step%i/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins.root",stepNumber-1);
  }

  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");

  int n_entries = t_unf->GetEntries();

  /*
  int const nbins_z = 50;
  float min_z = 0.;
  float max_z = 1.0;
  */
  
  int nbins_z = 50;
  if(doMid) nbins_z = 49;

  float min_z = 0.;
  if(doMid) min_z = 0.02;
  
  float max_z = 1.0;
  
  int nbins_jtpt = 15;
  
  Double_t z_frac_JetPtBin[nbins_jtpt][nbins_z];
  
  //only for prompt fwd in steps > 1, until the stats are better in mc
  Double_t *z_frac_allJetPtBins = new Double_t[nbins_z];

  // if this is step 2, 3, 4 ... take the z unfolded and use it for the prior, else put z prior to flat
  
  if(stepNumber > 1) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());

    //    cout << "file OK" << endl;
    
    TH1D *h_z_unf_jetPtBin[nbins_jtpt];

    for(int ibin = 0; ibin < nbins_jtpt; ibin++){
      h_z_unf_jetPtBin[ibin] = (TH1D*)filePrevStep->Get(Form("hMUnf_%i_Iter3;1",ibin));
    }
    
    TH1D *h_z_allBins = (TH1D*)h_z_unf_jetPtBin[0]->Clone();
    for(int ibin = 1; ibin < nbins_jtpt; ibin++){
      h_z_allBins->Add(h_z_unf_jetPtBin[ibin]);
    }
    
    h_z_allBins->Scale(1/h_z_allBins->Integral());

    for(int ibin = 1; ibin < nbins_jtpt; ibin++){
      h_z_unf_jetPtBin[ibin]->Scale(1/h_z_unf_jetPtBin[ibin]->Integral());
    }
            
    for(int iz = 1; iz <= nbins_z; iz++){
      
      z_frac_allJetPtBins[iz-1] = h_z_allBins->GetBinContent(iz);

      //cout << "iz = " << iz << " allJtPtBins content = " << z_frac_allJetPtBins[iz-1] << endl;
      
      for(int ijtpt = 0; ijtpt < nbins_jtpt; ijtpt++){

	z_frac_JetPtBin[ijtpt][iz-1] = h_z_unf_jetPtBin[ijtpt]->GetBinContent(iz);
	//cout <<"ijtpt = " << ijtpt << " z_frac_JetPtBin[ijtpt][iz-1] = " << z_frac_JetPtBin[ijtpt][iz-1] << endl;

      //leave it here, until the stats in prompt mc are better, remove later
	/*
	if(doPrompt && !doMid){
	  z_frac_JetPtBin[ijtpt][iz-1] = z_frac_allJetPtBins[iz-1];
	}
	*/	
      }
    }
    
  }
  else{
    for(int i = 1; i <= nbins_z; i++){
      for(int ijtpt = 0; ijtpt < nbins_jtpt; ijtpt++){
	z_frac_JetPtBin[ijtpt][i-1] = 1.0;
      }
    }
    
  }
  
  //J/Psi variables:
  
  float jp_pt = 0.;
  float jp_eta = 0.;
  float jp_mass = 0.;
     
  t_unf->SetBranchAddress("jp_pt", &jp_pt);
  t_unf->SetBranchAddress("jp_eta", &jp_eta);
  t_unf->SetBranchAddress("jp_mass", &jp_mass);

  //jet reco variables:

  float jt_pt = 0.;
  float jt_eta = 0.;
  
  t_unf->SetBranchAddress("jt_pt", &jt_pt);
  t_unf->SetBranchAddress("jt_eta", &jt_eta);

  float corr_AccEff = 0.;
  t_unf->SetBranchAddress("corr_AccEff", &corr_AccEff);

  float corr_ptw = 0.;
  t_unf->SetBranchAddress("corr_ptw",&corr_ptw);

  float finalCorr = 0.;
      
  //z corrected

  float z = 0.;
  
  t_unf->SetBranchAddress("z", &z);

  //jet gen variables

  float jt_ref_pt = 0.;
  
  t_unf->SetBranchAddress("jt_ref_pt", &jt_ref_pt);

  //gen z:

  float gen_z = 0.;
  
  t_unf->SetBranchAddress("gen_z", &gen_z);

  // histograms for renormalizing response matrix:
   
  TH2F * hRespZJetPtCentBin = new TH2F ("hRespZJetPtCentBin","z_gen vs z_raw for 25 < refpt < 40",nbins_z,min_z,max_z,nbins_z,min_z,max_z);
  TH2F * hRespZJetPtLowBin = new TH2F ("hRespZJetPtLowBin","z_gen vs z_raw for 10 < refpt < 25",nbins_z,min_z,max_z,nbins_z,min_z,max_z);
  TH2F * hRespZJetPtHighBin = new TH2F ("hRespZJetPtHighBin","z_gen vs z_raw for 40 < refpt < 55",nbins_z,min_z,max_z,nbins_z,min_z,max_z);

  TH1D * h_z_gen = new TH1D ("h_z_gen","z_gen for 25 < refpt < 55",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtCentBin = new TH1D ("hGenZJetPtCentBin","z_gen for 25 < refpt < 40",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtLowBin = new TH1D ("hGenZJetPtLowBin","z_gen for 10 < refpt < 25",nbins_z,min_z,max_z);
  TH1D * hGenZJetPtHighBin = new TH1D ("hGenZJetPtHighBin","z_gen for 40 < refpt < 55",nbins_z,min_z,max_z);
    
  int nbins_jetpt = 3;
  float min_jetpt = 15.;
  float max_jetpt = 45.;

  TH2F * h_jetpPtGen_jetPtReco = new TH2F ("h_jetpPtGen_jetPtReco","gen vs reco corr jet pt",nbins_jetpt,min_jetpt,max_jetpt,nbins_jetpt,min_jetpt,max_jetpt);

  hRespZJetPtCentBin->Sumw2();
  hRespZJetPtLowBin->Sumw2();
  hRespZJetPtHighBin->Sumw2();
  h_jetpPtGen_jetPtReco->Sumw2();

  // 4D histogram creation
  
  int fDim = 4.;
  Double_t fValue[fDim];
  
  Int_t* bins_sparce = new Int_t[fDim];
  Double_t *xmin_sparce = new Double_t[fDim];
  Double_t *xmax_sparce = new Double_t[fDim];

  bins_sparce[0] = 15;
  xmin_sparce[0] = 15.0;
  xmax_sparce[0] = 45.0;

  if(doMid){
    bins_sparce[1] = 49;
    xmin_sparce[1] = 0.02;
    xmax_sparce[1] = 1.0;
  }
  else{
    bins_sparce[1] = 50;
    xmin_sparce[1] = 0.;
    xmax_sparce[1] = 1.0;
  }

  bins_sparce[2] = 3;
  xmin_sparce[2] = 15.0;
  xmax_sparce[2] = 45.0;

  if(doMid){
    bins_sparce[3] = 7;
    xmin_sparce[3] = 0.02;
    xmax_sparce[3] = 1.0;
  }
  else{
    bins_sparce[3] = 5;
    xmin_sparce[3] = 0.;
    xmax_sparce[3] = 1.0;
  }
    
  //initial not normalized 4D  
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  //4D with z gen flat
  THnSparseF * fSparse_newZNorm = new THnSparseF("hs_newZNorm", "hs_newZNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newZNorm->Sumw2();
  fSparse_newZNorm->CalculateErrors();

  //4D with z gen flat and gen pt flat
  THnSparseF * fSparse_newJetPtNorm = new THnSparseF("hs_newJetPtNorm", "hs_newJetPtNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newJetPtNorm->Sumw2();
  fSparse_newJetPtNorm->CalculateErrors();

  float jp_pt_cut = 0.;
  
  for(int nEv = 0; nEv < n_entries; nEv++){

    if(doTrain && nEv%2 == 1) continue;
    if(!doTrain && nEv%2 == 0) continue;
     
    t_unf->GetEntry(nEv);
    
    // check event content

    if(doMid) jp_pt_cut = 6.5;
    if(!doMid) jp_pt_cut = 3.0;
    
    if(jp_pt < jp_pt_cut || jp_pt > 35.0) continue;
    
    if(doMid && TMath::Abs(jp_eta) > 1.6) continue;
    
    if(!doMid){
      if(TMath::Abs(jp_eta) < 1.6 || TMath::Abs(jp_eta) > 2.4) continue;
    }
            
    if(jp_mass < 2.6 || jp_mass > 3.5) continue;

    if(jt_pt < 15.0 || jt_pt > 45.0) continue;
    if(jt_ref_pt < 15.0 || jt_ref_pt > 45.0) continue;
    if(TMath::Abs(jt_eta) > 2.4) continue;
    

    double gen_z2 = gen_z;
    if(gen_z == 1.) gen_z2 = 0.99;

    double z2 = z;
    if(z == 1.) z2 = 0.99;

    if(gen_z2 >=1. ) continue;
    if(z2 >= 1. ) continue;
        
    if(jt_ref_pt >= 25.0 && jt_ref_pt < 35.){
      hRespZJetPtCentBin->Fill(gen_z2,z2);
      hGenZJetPtCentBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= 15.0 && jt_ref_pt < 25.0){
      hRespZJetPtLowBin->Fill(gen_z2,z2);
      hGenZJetPtLowBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= 35.0 && jt_ref_pt < 45.0){
      hRespZJetPtHighBin->Fill(gen_z2,z2);
      hGenZJetPtHighBin->Fill(gen_z2);
    }
        	  
    h_jetpPtGen_jetPtReco->Fill(jt_ref_pt,jt_pt);
    h_z_gen->Fill(gen_z2);
    
    //response matrix fill
    fValue[0] = jt_ref_pt;
    fValue[1] = gen_z2;
    fValue[2] = jt_pt;
    fValue[3] = z2;

    finalCorr = corr_AccEff*corr_ptw;
    fSparse->Fill(fValue,finalCorr);
  }

  
  // check 4D content

  cout << "number of bins :" << fSparse->GetNbins() <<endl;


  ///
  int nBinZ_gen = 50;
  if(doMid) nBinZ_gen = 49;
  int nBinZ_reco = 5;
  if(doMid) nBinZ_reco = 7;
  
  int nBinJet_gen = 15;
  int nBinJet_reco = 3;

  float jetPt_gen_binWidth = 30./nBinJet_gen;
  float z_gen_binWidth = 0.; // = 1./nBinZ_gen;
  if(doMid) z_gen_binWidth = 0.98/nBinZ_gen;
  else z_gen_binWidth = 1./nBinZ_gen;
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){

      float izMid = iz*z_gen_binWidth+z_gen_binWidth/2;
      if(doMid) izMid = iz*z_gen_binWidth+z_gen_binWidth/2+0.02;
      float ijetMid = 15.+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;

      double valIntegral = 0;

      //cout << "iz = " << iz << " izMid = " << izMid << endl;
      
      double zUnfFrac = 1.;
      zUnfFrac = z_frac_JetPtBin[ijet][iz];
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  float jjetMid = 10.*jjet+20;
	  float jzMid = 0.2*jz+0.1;
	  if(doMid) jzMid = 0.14*jz+0.09;

	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};

	  int bin = fSparse->GetBin(x);
	  valIntegral+=fSparse->GetBinContent(bin);
	}
      }
      
      double scaleFactor = 0;
      //if(valIntegral>0) scaleFactor = 1./valIntegral;
      if(valIntegral>0) scaleFactor = zUnfFrac/valIntegral;
      
      
      int countBin = 0;
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){
	    float jjetMid = 10.*jjet+20;
	    float jzMid = 0.2*jz+0.1;
	    if(doMid) jzMid = 0.14*jz+0.09;

	    const double x[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin = fSparse->GetBin(x);

	    //this does not scale errors properly
	    //fSparse_newPrior->Fill(x, fSparse->GetBinContent(bin)*scaleFactor);

	    //do it bin by bin
	    
	    const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin2 = fSparse_newZNorm->GetBin(x2);

	    double binCont = fSparse_newZNorm->GetBinContent(bin2);
	    double binErr  = fSparse_newZNorm->GetBinError(bin2);

	    double newBinCont = binCont+fSparse->GetBinContent(bin)*scaleFactor;
	    double newBinErr = binErr+fSparse->GetBinError(bin)*scaleFactor;

	    fSparse_newZNorm->SetBinContent(bin2, newBinCont);
	    fSparse_newZNorm->SetBinError(bin2, newBinErr);
	    
	  }
      }
    }

  }


  
  /*
  TH1D *gTruth1DX = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(0,"E"));
  TH1D *gTruth1DY = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(1,"E"));
  
  TCanvas * can1 = new TCanvas("can1","can1",1200,600);
  can1->Divide(2,1);
  can1->cd(1);
  gTruth1DX->Draw("EP");
  can1->cd(2);
  gTruth1DY->Draw("EP");
  can1->SaveAs("matrixTest.png");
  */
  

  /// jet pt normalization -> put it to match jet pt truth

  for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      
    //float ijetMid = 10.*ijet+20;
    //float ijetMid = ijet+15.5;
    float ijetMid = 15.+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;
    double valIntegralOrig = 0;
    double valIntegralZNorm = 0;

    for(int iz = 0; iz < nBinZ_gen; iz++){
      float izMid = iz*z_gen_binWidth+z_gen_binWidth/2;
      if(doMid) izMid = iz*z_gen_binWidth+z_gen_binWidth/2+0.02;
      
      //float izMid = 0.02*iz+0.01;
           
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){

	  float jjetMid = 10.*jjet+20;
	  float jzMid = 0.2*jz+0.1;
	  if(doMid) jzMid = 0.14*jz+0.09;

	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin = fSparse->GetBin(x);
	  valIntegralOrig+=fSparse->GetBinContent(bin);

	  const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin2 = fSparse_newZNorm->GetBin(x2);
	  valIntegralZNorm+=fSparse_newZNorm->GetBinContent(bin2);
	  	  
	}

      }

    }

    double scaleFactor = 0;
    if(valIntegralZNorm>0) scaleFactor = valIntegralOrig/valIntegralZNorm;
    
      for(int iz = 0; iz < nBinZ_gen; iz++){
	//float izMid = 0.02*iz+0.01;

	float izMid = iz*z_gen_binWidth+z_gen_binWidth/2;
	if(doMid) izMid = iz*z_gen_binWidth+z_gen_binWidth/2+0.02;
	
	for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){

	    float jjetMid = 10.*jjet+20;
	    float jzMid = 0.2*jz+0.1;
	    if(doMid) jzMid = 0.14*jz+0.09;
	        
	    const double x3[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin3 = fSparse_newZNorm->GetBin(x3);

	    //does not scale errors properly
	    //fSparse_newJetPtNorm->Fill(x3, fSparse_newZNorm->GetBinContent(bin3)*scaleFactor);

	    const double x4[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin4 = fSparse_newJetPtNorm->GetBin(x4);

	    double binCont = fSparse_newJetPtNorm->GetBinContent(bin4);
	    double binErr  = fSparse_newJetPtNorm->GetBinError(bin4);

	    double newBinCont = binCont+fSparse_newZNorm->GetBinContent(bin3)*scaleFactor;
	    double newBinErr = binErr+fSparse_newZNorm->GetBinError(bin3)*scaleFactor;

	    fSparse_newJetPtNorm->SetBinContent(bin4, newBinCont);
	    fSparse_newJetPtNorm->SetBinError(bin4, newBinErr);
	    	    
	  }
	}
	
      }
        
  }


  /*
  TH2D *gTruth2D = dynamic_cast<TH2D*>(fSparse_newJetPtNorm->Projection(0,1));
  gTruth2D->Draw("colz");

  
  TH1D *gTruth1DX_jtPt = dynamic_cast<TH1D*>(fSparse_newJetPtNorm->Projection(0,"E"));
  TH1D *gTruth1DY_jtPt = dynamic_cast<TH1D*>(fSparse_newJetPtNorm->Projection(1,"E"));

  TCanvas * can3 = new TCanvas("can3","can3",1200,600);
  can3->Divide(2,1);
  can3->cd(1);
  gTruth1DX_jtPt->Draw("EP");
  can3->cd(2);
  gTruth1DY_jtPt->Draw("EP");
  can3->SaveAs("matrixTest_jtPt.png");
  */
  
  ///
    
  TFile *file_forProf_varBin = new TFile(outputfile.c_str(),"RECREATE");

  h_jetpPtGen_jetPtReco->Write();
  
  fSparse->Write();

  fSparse_newZNorm->Write();
  fSparse_newJetPtNorm->Write();
  h_z_gen->Write();

  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
  
  hRespZJetPtCentBin->Write();
  hRespZJetPtLowBin->Write();
  hRespZJetPtHighBin->Write();
  
}

void prepareInputs(Int_t step = 1){

  prepare(true,true,true,step);
  prepare(true,false,true,step);

  prepare(true,true,false,step);
  prepare(true,false,false,step);

  //

  prepare(false,true,true,step);
  prepare(false,false,true,step);

  prepare(false,true,false,step);
  prepare(false,false,false,step);


}

