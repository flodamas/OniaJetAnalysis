#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void prepare(bool doPrompt = false, bool doPbPb = true, bool doTrain=true, Int_t stepNumber = 1){
  printInput();
  if (!setCaseTag()) return;
  
  string filename = "";
  string outputfile = "";
  string filenamePrevStep = "";

  gSystem->mkdir(Form("%s/mcUnf/unfInput",unfPath.c_str()));
  gSystem->mkdir(Form("%s/mcUnf/unfInput/step%i",unfPath.c_str(),stepNumber));
  gSystem->mkdir(Form("%s/mcUnf/unfOutput",unfPath.c_str()));
  gSystem->mkdir(Form("%s/mcUnf/unfOutput/step%i",unfPath.c_str(),stepNumber));

  filename = Form("~/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_%s_%s_NoBkg%s_AccEff_JEC.root",doPrompt?"MCJPSIPR":"MCJPSINOPR",doPbPb?"PbPb":"PP",mc2015?"":Form("_jetR%d",(int)(jetR*10)));
  outputfile = Form("%s/mcUnf/unfInput/step%i/unfolding_4D_%s_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), stepNumber, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doTrain?"Train":"Test", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  
  if(stepNumber > 1) filenamePrevStep = Form("%s/mcUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), stepNumber-1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());

  cout <<"filename = "<<filename<<endl;
  cout <<"outputfile = "<<outputfile<<endl;
  cout <<"filenamePrevStep = "<<filenamePrevStep<<endl;

  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");

  int n_entries = t_unf->GetEntries();
    
  Double_t z_frac_JetPtBin_dataMC[nBinJet_reco][nBinZ_reco];
  TH2D* mcDataWeight = NULL;
  if (dataDist) {
    if (stepNumber==1) {
      TFile *fileData = TFile::Open(Form("%s/dataUnf/data_results/meas_%s_data_%s_statErrs.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"prompt")); //use prompt data always 
      TH2D *hMeasuredData = (TH2D*) fileData->Get("h_Meas;1");
      TH2D *hMeasuredMC = (TH2D*) hMeasuredData->Clone("hMeasuredMC");
      mcDataWeight = (TH2D*) hMeasuredData->Clone("mcDataWeight");
      t_unf->Draw("jt_pt:z>>hMeasuredMC",Form("corr_AccEff*corr_ptw*((centr>=%d && centr<%d) && (jp_pt>%f && jp_pt<%f) && (TMath::Abs(jp_eta)<%f) && (jp_mass>2.6 && jp_mass<3.5) && (jt_pt>%f && jt_pt<%f) && (jt_ref_pt>%f && jt_ref_pt<%f) && (TMath::Abs(jt_eta)<2.4) && gen_z <=1. && z<=1)",centShift==0?10:centShift==-1?5:15, centShift==0?190:centShift==-1?185:195, min_jp_pt, max_jp_pt, max_jp_eta, min_jetpt, max_jetpt, min_jetpt, max_jetpt));
      mcDataWeight->Divide(hMeasuredMC);
    }//end of (stepNumber==1)
    else {
      TFile *fileW = TFile::Open(Form("%s/mcUnf/unfInput/step1/unfolding_4D_%s_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root",unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doTrain?"Train":"Test", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str()));
      mcDataWeight = (TH2D*) fileW->Get("mcDataWeight");
    }
  }//end of if (dataDist)
  
  Double_t z_frac_JetPtBin[nBinJet_gen][nBinZ_gen];
  
  //only for prompt fwd in steps > 1, until the stats are better in mc
  Double_t *z_frac_allJetPtBins = new Double_t[nBinZ_gen];
  
  // if this is step 2, 3, 4 ... take the z unfolded and use it for the prior, else put z prior to flat    
  if(stepNumber > 1) {
    TFile *filePrevStep = new TFile(filenamePrevStep.c_str());
    // cout << "file OK" << endl;
    TH1D *h_z_unf_jetPtBin[nBinJet_gen];

    for(int ibin = 0; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin] = (TH1D*)filePrevStep->Get(Form("hMUnf_%i_Iter%d;1",ibin,nIter));
    }
    
    TH1D *h_z_allBins = (TH1D*)h_z_unf_jetPtBin[0]->Clone();
    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_allBins->Add(h_z_unf_jetPtBin[ibin]);
    }
    
    h_z_allBins->Scale(1/h_z_allBins->Integral());

    for(int ibin = 1; ibin < nBinJet_gen; ibin++){
      h_z_unf_jetPtBin[ibin]->Scale(1/h_z_unf_jetPtBin[ibin]->Integral());
    }
            
    for(int iz = 1; iz <= nBinZ_gen; iz++){
      z_frac_allJetPtBins[iz-1] = h_z_allBins->GetBinContent(iz);
      //cout << "iz = " << iz << " allJtPtBins content = " << z_frac_allJetPtBins[iz-1] << endl;

      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][iz-1] = h_z_unf_jetPtBin[ijtpt]->GetBinContent(iz);
	//cout <<"ijtpt = " << ijtpt << " z_frac_JetPtBin[ijtpt][iz-1] = " << z_frac_JetPtBin[ijtpt][iz-1] << endl;
      }
    }
  }
  else{
    for(int i = 1; i <= nBinZ_gen; i++){
      for(int ijtpt = 0; ijtpt < nBinJet_gen; ijtpt++){
	z_frac_JetPtBin[ijtpt][i-1] = 1.0;
      }
    }    
  }
  
  //J/Psi variables:
  float jp_pt = 0.; t_unf->SetBranchAddress("jp_pt", &jp_pt);
  float jp_eta = 0.; t_unf->SetBranchAddress("jp_eta", &jp_eta);
  float jp_mass = 0.; t_unf->SetBranchAddress("jp_mass", &jp_mass);

  //jet reco variables:
  float jt_pt = 0.; t_unf->SetBranchAddress("jt_pt", &jt_pt);
  float jt_eta = 0.; t_unf->SetBranchAddress("jt_eta", &jt_eta);

  //centrality
  int centr = 0; t_unf->SetBranchAddress("centr", &centr);
  
  //correction variables:
  float corr_AccEff = 0.; t_unf->SetBranchAddress("corr_AccEff", &corr_AccEff);
  float corr_ptw = 0.; t_unf->SetBranchAddress("corr_ptw",&corr_ptw);
  float finalCorr = 0.;
      
  //z corrected
  float z = 0.; t_unf->SetBranchAddress("z", &z);

  //jet gen variables
  float jt_ref_pt = 0.; t_unf->SetBranchAddress("jt_ref_pt", &jt_ref_pt);

  //gen z:
  float gen_z = 0.; t_unf->SetBranchAddress("gen_z", &gen_z);

  // histograms for renormalizing response matrix:
  TH2F * hRespZJetPtCentBin = new TH2F ("hRespZJetPtCentBin","z_gen vs z_raw for 30 < refpt < 40",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtLowBin = new TH2F ("hRespZJetPtLowBin","z_gen vs z_raw for 20 < refpt < 30",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);
  TH2F * hRespZJetPtHighBin = new TH2F ("hRespZJetPtHighBin","z_gen vs z_raw for 40 < refpt < 50",nBinZ_gen,min_z,max_z,nBinZ_gen,min_z,max_z);

  TH1D * h_z_gen = new TH1D ("h_z_gen","z_gen for 20 < refpt < 50",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtCentBin = new TH1D ("hGenZJetPtCentBin","z_gen for 30 < refpt < 40",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtLowBin = new TH1D ("hGenZJetPtLowBin","z_gen for 20 < refpt < 30",nBinZ_gen,min_z,max_z);
  TH1D * hGenZJetPtHighBin = new TH1D ("hGenZJetPtHighBin","z_gen for 40 < refpt < 50",nBinZ_gen,min_z,max_z);

  TH2F * h_jetpPtGen_jetPtReco = new TH2F ("h_jetpPtGen_jetPtReco","gen vs reco corr jet pt",nBinJet_reco,min_jetpt,max_jetpt,nBinJet_reco,min_jetpt,max_jetpt);

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

  bins_sparce[0] = nBinJet_gen;
  xmin_sparce[0] = min_jetpt;
  xmax_sparce[0] = max_jetpt;

  bins_sparce[1] = nBinZ_gen;
  xmin_sparce[1] = min_z;
  xmax_sparce[1] = max_z;
  
  bins_sparce[2] = nBinJet_reco;
  xmin_sparce[2] = min_jetpt;
  xmax_sparce[2] = max_jetpt;

  bins_sparce[3] = nBinZ_reco;
  xmin_sparce[3] = min_z;
  xmax_sparce[3] = max_z;
  
  //initial not normalized 4D  
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  THnSparseF * fSparseOrig = new THnSparseF("hs_orig", "hs_orig", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparseOrig->Sumw2();
  fSparseOrig->CalculateErrors();

  //4D with z gen flat
  THnSparseF * fSparse_newZNorm = new THnSparseF("hs_newZNorm", "hs_newZNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newZNorm->Sumw2();
  fSparse_newZNorm->CalculateErrors();

  //4D with z gen flat and gen pt = truth
  THnSparseF * fSparse_newJetPtNorm = new THnSparseF("hs_newJetPtNorm", "hs_newJetPtNorm", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse_newJetPtNorm->Sumw2();
  fSparse_newJetPtNorm->CalculateErrors();
  
  for(int nEv = 0; nEv < n_entries; nEv++){
        
    if (!sameSample) {
      if (doPbPb) {
	if(doTrain && nEv%10 >= 8) continue;
	if(!doTrain && nEv%10 < 8) continue;
      }
      else {
	if(doTrain && nEv%2 == 1) continue;
	if(!doTrain && nEv%2 == 0) continue;
      }
    }
    
    t_unf->GetEntry(nEv);

    //if (doPbPb && (centr<10 || centr>=190)) continue;
    if (!doTrain) {
      if (doPbPb && centShift==-1 && (centr<5 || centr>=185)) continue;
      if (doPbPb && centShift==0 && (centr<10 || centr>=190)) continue;
      if (doPbPb && centShift==1 && (centr<15 || centr>=195)) continue;
    }
    else {
      if (doPbPb && (centr<10 || centr>=190)) continue;
    }
    // check event content
    if(jp_pt < min_jp_pt || jp_pt > max_jp_pt) continue;    
    if(TMath::Abs(jp_eta) > max_jp_eta ) continue;
    if(jp_mass < 2.6 || jp_mass > 3.5) continue;
    if(jt_pt < min_jetpt || jt_pt > max_jetpt) continue;
    if(jt_ref_pt < min_jetpt || jt_ref_pt > max_jetpt) continue;
    if(TMath::Abs(jt_eta) > max_jt_eta) continue;    

    double gen_z2 = gen_z;
    if(gen_z == 1.) gen_z2 = 0.99;

    double z2 = z;
    if(z == 1.) z2 = 0.99;

    if(gen_z2 >=1. ) continue;
    if(z2 >= 1. ) continue;

    /*    
    if(jt_ref_pt >= 30.0 && jt_ref_pt < 40.){
      hRespZJetPtCentBin->Fill(gen_z2,z2);
      hGenZJetPtCentBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= min_jetpt && jt_ref_pt < 30.0){
      hRespZJetPtLowBin->Fill(gen_z2,z2);
      hGenZJetPtLowBin->Fill(gen_z2);
    }
    
    if(jt_ref_pt >= 40.0 && jt_ref_pt < max_jetpt){
      hRespZJetPtHighBin->Fill(gen_z2,z2);
      hGenZJetPtHighBin->Fill(gen_z2);
    }
        	  
    h_jetpPtGen_jetPtReco->Fill(jt_ref_pt,jt_pt);
    h_z_gen->Fill(gen_z2);
    */

    //response matrix fill
    fValue[0] = jt_ref_pt;
    fValue[1] = gen_z2;
    fValue[2] = jt_pt;
    fValue[3] = z2;

    finalCorr = corr_AccEff*corr_ptw;
    fSparseOrig->Fill(fValue,finalCorr);
    if (dataDist) finalCorr = finalCorr*mcDataWeight->GetBinContent(mcDataWeight->FindBin(z2,jt_pt));
    fSparse->Fill(fValue,finalCorr);
  }

  
  // check 4D content
  cout << "number of bins :" << fSparse->GetNbins() <<endl;
  
  ///
  for(int iz = 0; iz < nBinZ_gen; iz++){
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){

      float izMid = min_z+iz*z_gen_binWidth+z_gen_binWidth/2;
      float ijetMid = min_jetpt+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;

      double valIntegral = 0;
      //cout << "iz = " << iz << " izMid = " << izMid << endl;
      double zUnfFrac = 1.;
      zUnfFrac = z_frac_JetPtBin[ijet][iz];
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	  float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;

	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};

	  int bin = fSparse->GetBin(x);
	  valIntegral+=fSparse->GetBinContent(bin);
	} // end of jz < nBinZ_reco
      }// end of jjet < nBinJet_reco
      
      double scaleFactor = 0;
      //if(valIntegral>0) scaleFactor = 1./valIntegral;
      if(valIntegral>0) scaleFactor = zUnfFrac/valIntegral;
      int countBin = 0;

      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){
	    float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	    float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;
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
	    
	  } //end of jz < nBinZ_reco
      } //end of jjet < nBinJet_reco
    }  //end of ijet < nBinJet_gen
  } //end of iz < nBinZ_gen  
  
  TH1D *gTruth1DX = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(0,"E"));
  TH1D *gTruth1DY = dynamic_cast<TH1D*>(fSparse_newZNorm->Projection(1,"E"));
  
  TCanvas * can1 = new TCanvas("can1","can1",1200,600);
  can1->Divide(2,1);
  can1->cd(1);
  gTruth1DX->Draw("EP");
  can1->cd(2);
  gTruth1DY->Draw("EP");
  can1->SaveAs("matrixTest.png");
  
  
  /// jet pt normalization -> put it to match jet pt truth

  for(int ijet = 0; ijet < nBinJet_gen; ijet++){      
    float ijetMid = min_jetpt+ijet*jetPt_gen_binWidth+jetPt_gen_binWidth/2;
    double valIntegralOrig = 0;
    double valIntegralZNorm = 0;

    for(int iz = 0; iz < nBinZ_gen; iz++){
      float izMid = min_z+iz*z_gen_binWidth+z_gen_binWidth/2;

      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){

	  float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
	  float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;

	  const double x[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin = fSparse->GetBin(x);
	  valIntegralOrig+=fSparse->GetBinContent(bin);

	  const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin2 = fSparse_newZNorm->GetBin(x2);
	  valIntegralZNorm+=fSparse_newZNorm->GetBinContent(bin2);	  	  
	} //end of jz < nBinZ_reco
      } //end of jjet < nBinJet_reco
    } //end of iz < nBinZ_gen

    double scaleFactor = 0;
    if(valIntegralZNorm>0) scaleFactor = valIntegralOrig/valIntegralZNorm;
    
      for(int iz = 0; iz < nBinZ_gen; iz++){
	float izMid = min_z+iz*z_gen_binWidth+z_gen_binWidth/2;
	
	for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	  for(int jz = 0; jz < nBinZ_reco; jz++){

	    float jjetMid = min_jetpt+jjet*jetPt_reco_binWidth+jetPt_reco_binWidth/2;
            float jzMid = min_z+jz*z_reco_binWidth+z_reco_binWidth/2;

	    const double x3[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin3 = fSparse_newZNorm->GetBin(x3);

	    const double x4[4] = {ijetMid,izMid,jjetMid,jzMid};
	    int bin4 = fSparse_newJetPtNorm->GetBin(x4);

	    double binCont = fSparse_newJetPtNorm->GetBinContent(bin4);
	    double binErr  = fSparse_newJetPtNorm->GetBinError(bin4);

	    double newBinCont = binCont+fSparse_newZNorm->GetBinContent(bin3)*scaleFactor;
	    double newBinErr = binErr+fSparse_newZNorm->GetBinError(bin3)*scaleFactor;

	    fSparse_newJetPtNorm->SetBinContent(bin4, newBinCont);
	    fSparse_newJetPtNorm->SetBinError(bin4, newBinErr);
	  }//end of jz < nBinZ_reco
	}//end of jjet < nBinJet_reco
      }//end of iz < nBinZ_gen
  }//end of ijet < nBinJet_gen

  
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
  
  ///
    
  TFile *file_forProf_varBin = new TFile(outputfile.c_str(),"RECREATE");

  h_jetpPtGen_jetPtReco->Write();
  
  fSparse->Write();
  fSparseOrig->Write();
  fSparse_newZNorm->Write();
  fSparse_newJetPtNorm->Write();
  h_z_gen->Write();

  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
  
  hRespZJetPtCentBin->Write();
  hRespZJetPtLowBin->Write();
  hRespZJetPtHighBin->Write();

  mcDataWeight->Write();
  
  file_forProf_varBin->Close();
}

void prepareInputs(Int_t step = 1){
  //prepare(bool doPrompt, bool doPbPb, bool doTrain, Int_t stepNumber)
  prepare(true,true,true,step);
  prepare(true,true,false,step);

  if (step<=nSIter_pp && centShift == 0){
    prepare(true,false,true,step);
    prepare(true,false,false,step);
  }
  //nonprompt
  prepare(false,true,true,step);
  prepare(false,true,false,step);
  
  if (step<=nSIter_pp && centShift == 0){
    prepare(false,false,true,step);
    prepare(false,false,false,step);
  }

}
