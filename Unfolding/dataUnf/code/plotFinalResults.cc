#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt);

void plotFinalResults(bool doPrompt = true)
{
  gStyle->SetOptStat(0);
    
  string filePPName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_pp, "_nominal");
  
  string filePbPbName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, "_nominal");
    
  string fileOutputName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PPvsPbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_PbPbnIter%inSIter%i_PPnIter%inSIter%i_statError.pdf",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nIter, nSIter, nIter_pp, nSIter_pp);

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
  }

  histPP->SetTitle("");
  histPbPb->SetTitle("");
  TGraphAsymmErrors* graphPP = new TGraphAsymmErrors(histPP);//(nBin,x_pp,y_pp,exl_pp,exh_pp,eyl_pp,eyh_pp);
  TGraphAsymmErrors* graphPbPb = new TGraphAsymmErrors(histPbPb);//(nBin,x_pbpb,y_pbpb,exl_pbpb,exh_pbpb,eyl_pbpb,eyh_pbpb);

  TGraphAsymmErrors* graphPPSyst = systUncertaintyHistAll(false, doPrompt, histPP);//new TGraphAsymmErrors(histPPSyst);
  TGraphAsymmErrors* graphPbPbSyst = systUncertaintyHistAll(true, doPrompt, histPbPb);//new TGraphAsymmErrors(histPbPbSyst);
  
  graphPP->SetLineColor(kpink);
  graphPP->SetMarkerColor(kpink);
  graphPP->SetMarkerStyle(kFullSquare);
  graphPP->SetMarkerSize(1.5);
  graphPP->SetTitle("");
  graphPP->GetXaxis()->SetTitle("z");
  graphPP->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz}");
  //graphPP->GetYaxis()->SetRangeUser(0, graphPP->GetMaximum()*1.5);
  
  graphPbPb->SetLineColor(kazure);
  graphPbPb->SetMarkerColor(kazure);
  graphPbPb->SetMarkerStyle(kFullCircle);
  graphPbPb->SetMarkerSize(1.5);
  
  graphPbPbSyst->SetLineColor(kazure);
  graphPbPbSyst->SetMarkerColor(kazure);
  graphPbPbSyst->SetMarkerStyle(kFullCircle);
  graphPbPbSyst->SetFillColorAlpha(kazureLight, 0.75);
  
  graphPPSyst->SetLineColor(kpink);
  graphPPSyst->SetMarkerColor(kpink);
  graphPPSyst->SetMarkerStyle(kFullSquare);
  graphPPSyst->SetFillColorAlpha(kpinkLight, 0.75);
  
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(graphPP, "pp","lp");
  leg->AddEntry(graphPbPb, "PbPb","lp");
  TCanvas* c = new TCanvas("c","",1200,900);
  c->cd();
  graphPP->Draw("AP");
  graphPbPbSyst->Draw("5");
  graphPPSyst->Draw("5");
  graphPbPb->Draw("P");
  graphPP->Draw("P");
  leg->Draw("same");
  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());

  fileOutputName.replace(fileOutputName.find(".png"),4,".pdf");
  fileOutputName.replace(fileOutputName.find("PPvsPbPb"),8,"Raa");

  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  Double_t eyl_syst[5] = {0,0,0,0,0};
  Double_t eyh_syst[5] = {0,0,0,0,0};

  for (int iBin = 0; iBin < nBinZ_reco-1; iBin++) {
    int bin = histPbPb->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
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
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());
}

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  systUncertaintyHistReg(doPbPb, doPrompt, nominalHist);
  systUncertaintyHistTrStat(doPbPb, doPrompt);
  TCanvas* cSyst = new TCanvas("cSyst","",800,800);
  nominalHist->SetLineColor(col[0]);
  nominalHist->SetMarkerColor(col[0]);
  nominalHist->SetMarkerStyle(markerStyle[0]);
  cSyst->cd();
  nominalHist->SetStats(0);
  nominalHist->Draw();
  double maxSyst = nominalHist->GetMaximum();
  TLegend* legSyst = new TLegend(0.7,0.6,0.9,0.9);
  //if (!doPbPb) legSyst = new TLegend(0.4,0.25,0.6,0.45);
  legSyst->SetBorderSize(0);
  legSyst->SetFillStyle(0);
  legSyst->AddEntry(nominalHist,"nominal","lp");
  string systNames [] = {"_JESSyst","_JERSyst","_SFSyst","_nominal_centShiftSyst","_Regularisation", "_TrStat"};
  int nSyst = sizeof(systNames)/sizeof(systNames[0]);
  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  
  for (int i=0; i<nSyst; i++) {
    string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    TFile* fileUp = TFile::Open(fileUpName.c_str());
    TFile* fileDown = TFile::Open(fileDownName.c_str());
    if (!fileUp || !fileDown) {cout<<"[WARNING] systematic file not found:"<<systNames[i]<<endl; continue;}
    TH1D* histUp = (TH1D*) fileUp->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
    TH1D* histDown = (TH1D*) fileDown->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
    if (!histUp || !histDown) {cout<<"[WARNING] systematic histogram not found:"<<systNames[i]<<endl; continue;}

    histUp->SetStats(0);
    histDown->SetStats(0);

    histUp->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
    histDown->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));

    if (histUp->GetMaximum() > maxSyst) maxSyst = histUp->GetMaximum();
    if (histDown->GetMaximum() > maxSyst) maxSyst = histDown->GetMaximum();

    nominalHist->GetYaxis()->SetRangeUser(0,1.5*maxSyst);

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
      int nBin = nominalHist->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
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
      errUp = errUp/nominalHist->GetBinContent(nBin);
      errDown = errDown/nominalHist->GetBinContent(nBin);
            
      if (i==0) {
	x[iBin] = nominalHist->GetBinCenter(nBin);
	y[iBin] = nominalHist->GetBinContent(nBin);
	exl[iBin] = z_reco_binWidth/3.;
	exh[iBin] = z_reco_binWidth/3.;
	}
      eyl[iBin] = y[iBin]*sqrt(pow(errDown,2) + pow(eyl[iBin]/y[iBin],2));
      eyh[iBin] = y[iBin]*sqrt(pow(errUp,2) + pow(eyh[iBin]/y[iBin],2));
      cout << "x = "<<x[iBin]<<", y = "<<y[iBin]<<", eyl = "<<eyl[iBin]<<", eyh = "<<eyh[iBin]<<endl;
    }
  }
  legSyst->Draw("same");
  cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
  TGraphAsymmErrors* systHist = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl,eyh);
  return systHist;
}

void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", doPbPb?nIter_syst:nIter_syst_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter_syst:nSIter_syst_pp, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter_syst:nSIter_syst_pp));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  TH1D* nomHist = (TH1D*) nominalHist->Clone("nomHist");
  //systHist->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
  nomHist->Scale(doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    double binCenter = systHist->GetBinCenter(i);
    double err = fabs(nomHist->GetBinContent(nomHist->FindBin(binCenter))-systHist->GetBinContent(i));
    histUp->SetBinContent(i, nomHist->GetBinContent(i)+err);
    histDown->SetBinContent(i, nomHist->GetBinContent(i)-err);
  }
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularisation");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularisation");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}

void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/matrixOper/systUnc_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  TH1D* systHist = (TH1D*) fileSyst->Get("zUnf_trMatrixSyst");
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    histUp->SetBinContent(i, histUp->GetBinContent(i)+histUp->GetBinError(i));
    histDown->SetBinContent(i, histDown->GetBinContent(i)-histDown->GetBinError(i));
  }
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_TrStat");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_TrStat");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}
