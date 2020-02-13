#include "inputParams.h"

TGraphAsymmErrors* systUncertaintyHist(bool doPbPb, bool doPrompt, TH1D* nominalHist);

void plotFinalResults(bool doPrompt = true)
{
  gStyle->SetOptStat(0);
    
  string filePPName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_pp, "_nominal");
  
  string filePbPbName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, "_nominal");
    
  string fileOutputName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PPvsPbPbvsCentShift_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIterPbPb%i_SIterPP%i_statError.pdf",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, nSIter_pp);

  TFile* filePP = TFile::Open(filePPName.c_str());
  TFile* filePbPb = TFile::Open(filePbPbName.c_str());
  
  TH1D* histPP = (TH1D*) filePP->Get(Form("zUnfSI%d",nSIter_pp));
  TH1D* histPbPb = (TH1D*) filePbPb->Get(Form("zUnfSI%d",nSIter));

  histPP->Scale(normPP*1./z_reco_binWidth);
  histPbPb->Scale(normPbPb*1./z_reco_binWidth);

  int nBin = histPP->GetNbinsX();
  for (int iBin=0;iBin<nBin; iBin++) {
    if (histPP->GetBinCenter(iBin)<unfStart) {
      histPP->SetBinContent(iBin,0);
      histPP->SetBinError(iBin,0);
      histPbPb->SetBinContent(iBin,0);
      histPbPb->SetBinError(iBin,0);
      continue;
    }
    //x_pp[iBin] = histPP->GetBinCenter(iBin);
    //y_pp[iBin] = histPP->GetBinContent(iBin);
    //exl_pp[iBin] = z_reco_binWidth/2.;
    //exh_pp[iBin] = z_reco_binWidth/2.;
    //eyl_pp[iBin] = histPP->GetBinError(iBin);
    //eyh_pp[iBin] = histPP->GetBinError(iBin);

    //x_pbpb[iBin] = histPbPb->GetBinCenter(iBin);
    //y_pbpb[iBin] = histPbPb->GetBinContent(iBin);
    //exl_pbpb[iBin] = z_reco_binWidth/2.;
    //exh_pbpb[iBin] = z_reco_binWidth/2.;
    //eyl_pbpb[iBin] = histPbPb->GetBinError(iBin);
    //eyh_pbpb[iBin] = histPbPb->GetBinError(iBin);
  }
  
  TGraphAsymmErrors* graphPP = new TGraphAsymmErrors(histPP);//(nBin,x_pp,y_pp,exl_pp,exh_pp,eyl_pp,eyh_pp);
  TGraphAsymmErrors* graphPbPb = new TGraphAsymmErrors(histPbPb);//(nBin,x_pbpb,y_pbpb,exl_pbpb,exh_pbpb,eyl_pbpb,eyh_pbpb);

  TGraphAsymmErrors* graphPPSyst = systUncertaintyHist(false, doPrompt, histPP);//new TGraphAsymmErrors(histPPSyst);
  TGraphAsymmErrors* graphPbPbSyst = systUncertaintyHist(true, doPrompt, histPbPb);//new TGraphAsymmErrors(histPbPbSyst);
  
  graphPP->SetLineColor(kred);
  graphPP->SetMarkerColor(kred);
  graphPP->SetMarkerStyle(kFullSquare);
  graphPP->SetMarkerSize(1.5);
  graphPP->SetTitle("");
  graphPP->GetXaxis()->SetTitle("z");
  graphPP->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz}");
  //graphPP->GetYaxis()->SetRangeUser(0, graphPP->GetMaximum()*1.5);
  
  graphPbPb->SetLineColor(kgreen);
  graphPbPb->SetMarkerColor(kgreen);
  graphPbPb->SetMarkerStyle(kFullCircle);
  graphPbPb->SetMarkerSize(1.5);
  
  graphPbPbSyst->SetLineColor(kgreen);
  graphPbPbSyst->SetMarkerColor(kgreen);
  graphPbPbSyst->SetMarkerStyle(kFullCircle);
  graphPbPbSyst->SetFillColorAlpha(kgreen, 0.35);
  
  graphPPSyst->SetLineColor(kred);
  graphPPSyst->SetMarkerColor(kred);
  graphPPSyst->SetMarkerStyle(kFullSquare);
  graphPPSyst->SetFillColorAlpha(kred, 0.35);
  
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(graphPP, "pp","lp");
  leg->AddEntry(graphPbPb, "PbPb","lp");
  //leg->AddEntry(graphPPSyst, "PbPb, 2.5% cent shift","lp");
  //leg->AddEntry(graphPbPbSyst, "PbPb, 7.5% cent shift","lp");
  TCanvas* c = new TCanvas("c","",1200,900);
  c->cd();
  graphPP->Draw("AP");
  graphPbPbSyst->Draw("5");
  graphPPSyst->Draw("5");
  //graphPbPb->SetDrawOption("AP");
  graphPbPb->Draw("P");
  graphPP->Draw("P");
  leg->Draw("same");
  c->SaveAs(fileOutputName.c_str());
  
  fileOutputName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_Raa_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIterPbPb%i_SIterPP%i_statError.pdf",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, nSIter_pp);
  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  Double_t eyl_syst[5] = {0,0,0,0,0};
  Double_t eyh_syst[5] = {0,0,0,0,0};

  for (int iBin = 0; iBin < nBinZ_reco-1; iBin++) {
    int bin = histPbPb->FindBin(iBin*z_reco_binWidth+min_z+0.22);
    x[iBin] = histPbPb->GetBinCenter(bin);//graphPPSyst->GetPointX(iBin);
    y[iBin] = histPbPb->GetBinContent(bin)/histPP->GetBinContent(histPP->FindBin(x[iBin]));//graphPPSyst->GetPointY(iBin)/graphPbPbSyst->GetPointY(iBin);
    exl[iBin] = z_reco_binWidth/3.;
    exh[iBin] = z_reco_binWidth/3.;
    eyh_syst[iBin] = y[iBin]*sqrt(pow(graphPbPbSyst->GetErrorYhigh(iBin)/histPbPb->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYhigh(iBin)/histPP->GetBinContent(bin),2));
    eyl_syst[iBin] = y[iBin]*sqrt(pow(graphPbPbSyst->GetErrorYlow(iBin)/histPbPb->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYlow(iBin)/histPP->GetBinContent(bin),2));
    eyh[iBin] = y[iBin]*sqrt(pow(histPbPb->GetBinError(bin)/histPbPb->GetBinContent(bin),2)+pow(histPP->GetBinError(bin)/histPP->GetBinContent(bin),2));
    eyl[iBin] = eyh[iBin];
  }

  TGraphAsymmErrors* graphRaaSyst = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl_syst,eyh_syst);
  TGraphAsymmErrors* graphRaa = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl,eyh);
  
  graphRaaSyst->SetLineColor(kviolet);
  graphRaaSyst->SetMarkerColor(kviolet);
  graphRaaSyst->SetMarkerStyle(kOpenCircle);
  graphRaaSyst->SetMarkerSize(1.5);
  graphRaaSyst->SetTitle("");
  graphRaaSyst->GetXaxis()->SetTitle("z");
  graphRaaSyst->GetYaxis()->SetTitle("Raa");
  graphRaaSyst->SetFillColorAlpha(kviolet, 0.35);

  graphRaa->SetLineColor(kviolet);
  graphRaa->SetMarkerColor(kviolet);
  graphRaa->SetMarkerStyle(kOpenCircle);
  graphRaa->SetMarkerSize(1.5);
  graphRaa->SetTitle("");
  graphRaa->GetXaxis()->SetTitle("z");
  graphRaa->GetYaxis()->SetTitle("Raa");
  graphRaa->SetFillColorAlpha(kviolet, 0.35);

  graphRaaSyst->Draw("A5");
  graphRaa->Draw("P");
  c->SaveAs(fileOutputName.c_str());
}

TGraphAsymmErrors* systUncertaintyHist(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  TCanvas* cSyst = new TCanvas("cSyst","",800,800);
  nominalHist->SetLineColor(col[0]);
  nominalHist->SetMarkerColor(col[0]);
  nominalHist->SetMarkerStyle(markerStyle[0]);
  cSyst->cd();
  nominalHist->Draw();
  TLegend* legSyst = new TLegend(0.45,0.45,0.75,0.75);
  if (!doPbPb) legSyst = new TLegend(0.4,0.25,0.6,0.45);
  legSyst->SetBorderSize(0);
  legSyst->SetFillStyle(0);
  legSyst->AddEntry(nominalHist,"nominal","lp");
  string systNames [] = {"_JESSyst","_JERSyst","_SFSyst","_nominal_centShiftSyst"};
  int nSyst = sizeof(systNames)/sizeof(systNames[0]);
  //int nBin = nominalHist->GetNbinsX();
  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  
  for (int i=0; i<nSyst; i++) {
    string fileUpName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    string fileDownName = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    TFile* fileUp = TFile::Open(fileUpName.c_str());
    TFile* fileDown = TFile::Open(fileDownName.c_str());
    if (!fileUp || !fileDown) {cout<<"[WARNING] systematic file not found:"<<systNames[i]<<endl; continue;}
    TH1D* histUp = (TH1D*) fileUp->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
    TH1D* histDown = (TH1D*) fileDown->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
    if (!histUp || !histDown) {cout<<"[WARNING] systematic histogram not found:"<<systNames[i]<<endl; continue;}

    histUp->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
    histDown->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));

    histUp->SetLineColor(col[i+1]);
    histUp->SetMarkerColor(col[i+1]);
    histUp->SetMarkerStyle(markerStyle[i+1]);
    histUp->SetLineStyle(7);
    histDown->SetLineColor(col[i+1]);
    histDown->SetMarkerColor(col[i+1]);
    histDown->SetMarkerStyle(markerStyle[i+1]);
    histDown->SetLineStyle(8);
    histUp->Draw("same");
    histDown->Draw("same");
    legSyst->AddEntry(histUp,Form("%s Up",systNames[i].c_str()),"lp");
    legSyst->AddEntry(histDown,Form("%s Down",systNames[i].c_str()),"lp");
    for (int iBin = 0; iBin < nBinZ_reco-1; iBin++) { //-1 for the underflow
      double errUp=0;
      double errDown=0;
      int nBin = nominalHist->FindBin(iBin*z_reco_binWidth+min_z+0.22);
      if (nominalHist->GetBinContent(nBin) == 0) continue;

      if (histUp->GetBinContent(nBin) > nominalHist->GetBinContent(nBin)) {
	if (histDown->GetBinContent(nBin) < nominalHist->GetBinContent(nBin)) {
	  errUp = histUp->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	}
	else {
	  errUp = histUp->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errDown = histDown->GetBinContent(nBin) - nominalHist->GetBinContent(nBin);
	  errUp = max(errUp,errDown);
	  errDown=0;
	}
      }
      else {
	if (histDown->GetBinContent(nBin) < nominalHist->GetBinContent(nBin)) {
	  errUp = nominalHist->GetBinContent(nBin) - histUp->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	  errDown = max(errUp,errDown);
	  errUp=0;
	}
	else {
	  errUp = nominalHist->GetBinContent(nBin) - histDown->GetBinContent(nBin);
	  errDown = nominalHist->GetBinContent(nBin) - histUp->GetBinContent(nBin);
	}
      }

      //errUp =  fabs(errUp = histUp->GetBinContent(nBin) - nominalHist->GetBinContent(nBin));
      //errDown = fabs(histDown->GetBinContent(nBin) - nominalHist->GetBinContent(nBin));
      //errUp = max(errUp,errDown);
      //errDown = errUp;
      //cout <<"nominal = "<<nominalHist->GetBinContent(nBin)<<", histUp = "<<histUp->GetBinContent(nBin)<<", histDown = "<<histDown->GetBinContent(nBin)<<", for syst = "<<systNames[i]<<", for z = "<<nominalHist->GetBinCenter(nBin)<<endl;
      //if (errUp < 0 || errDown < 0) cout <<"[WARNING] something is wrong: negative systematic error! errUp = "<<errUp<<", errDown = "<<errDown<<endl;
      errUp = errUp/nominalHist->GetBinContent(nBin);
      errDown = errDown/nominalHist->GetBinContent(nBin);
      
      //errUp = pow(errUp/nominalHist->GetBinContent(nBin),2) + pow(systHist->GetErrorYhigh(nBin)/nominalHist->GetBinContent(nBin),2);
      //errDown = pow(errDown/nominalHist->GetBinContent(nBin),2) + pow(systHist->GetErrorYlow(nBin)/nominalHist->GetBinContent(nBin),2);
      //errUp = sqrt(errUp);
      //errDown = sqrt(errDown);
      //cout <<"error Up = "<<errUp<<", error Down = "<<errDown<<", for syst = "<<systNames[i]<<", for z = "<<nominalHist->GetBinCenter(nBin)<<endl;
      //errUp = errUp*nominalHist->GetBinContent(nBin);
      //errDown = errDown*nominalHist->GetBinContent(nBin);
      //systHist->SetPointEYhigh(nBin,errUp);
      //systHist->SetPointEYlow(nBin,errDown);
      
      if (i==0) {
	x[iBin] = nominalHist->GetBinCenter(nBin);
	y[iBin] = nominalHist->GetBinContent(nBin);
	exl[iBin] = z_reco_binWidth/3.;
	exh[iBin] = z_reco_binWidth/3.;
	}
      eyl[iBin] = y[iBin]*sqrt(pow(errDown,2) + pow(eyl[iBin]/y[iBin],2));
      eyh[iBin] = y[iBin]*sqrt(pow(errUp,2) + pow(eyh[iBin]/y[iBin],2));
      //eyh[iBin] = sqrt(errUp*errUp + eyh[iBin]*eyh[iBin]);
      cout << "x = "<<x[iBin]<<", y = "<<y[iBin]<<", eyl = "<<eyl[iBin]<<", eyh = "<<eyh[iBin]<<endl;
    }
  }
  legSyst->Draw("same");
  cSyst->SaveAs(Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin_SIter%i.pdf",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
  TGraphAsymmErrors* systHist = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl,eyh);
  return systHist;
}
