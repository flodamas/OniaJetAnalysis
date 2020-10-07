#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

bool removePrelim = false;
bool plotMC = false;
bool plotNPR = true;

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt);
void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt);
void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist);
void plotFinalResults_oneStyle(bool doPrompt=true);
void plotFinalResults() {
  plotMC = false;
  //removePrelim = false;
  //plotFinalResults_oneStyle(true);
  removePrelim = true;
  plotFinalResults_oneStyle(true);

  plotMC = true;
  //removePrelim = false;
  //plotFinalResults_oneStyle(true);
  removePrelim = true;
  plotFinalResults_oneStyle(true);

}

void plotFinalResults_oneStyle(bool doPrompt) {
  gStyle->SetOptStat(0);
    
  string filePPName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_pp, "_nominal");
  
  string filePbPbName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, "_nominal");
    
  string filePrMCName = "../../../Fitter/Output/MCResults/mcResult_PP_prompt_all.root";
  string fileNprMCName = "../../../Fitter/Output/MCResults/mcResult_PP_nonprompt_all.root";

  string fileOutputName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PPvsPbPb%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_PbPbnIter%inSIter%i_PPnIter%inSIter%i_statError%s.pdf",unfPath.c_str(),plotMC?(plotNPR?"vsMC":"vsPRMC"):"",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nIter, nSIter, nIter_pp, nSIter_pp,removePrelim?"_noPreliminaryLabel":"");

  TFile* filePP = TFile::Open(filePPName.c_str());
  TFile* filePbPb = TFile::Open(filePbPbName.c_str());
  TFile* filePrMC = TFile::Open(filePrMCName.c_str());
  TFile* fileNprMC = TFile::Open(fileNprMCName.c_str());
  
  TH1D* histPP = (TH1D*) filePP->Get(Form("zUnfSI%d",nSIter_pp));
  TH1D* histPbPb = (TH1D*) filePbPb->Get(Form("zUnfSI%d",nSIter));
  TH1F* histPrMC = (TH1F*) filePrMC->Get("zDist");
  TH1F* histNprMC = (TH1F*) fileNprMC->Get("zDist");

  //double normpp = normPP;
  if (plotMC) normPP=1./histPP->Integral();
  histPP->Scale(normPP*1./z_reco_binWidth);
  histPbPb->Scale(normPbPb*1./z_reco_binWidth);
  histPrMC->Scale(1./histPrMC->Integral()*1./z_reco_binWidth);
  histNprMC->Scale(1./histNprMC->Integral()*1./z_reco_binWidth);

  int nBin = histPP->GetNbinsX();
  for (int iBin=0;iBin<nBin; iBin++) {
    if (histPP->GetBinCenter(iBin)<unfStart) {
      histPP->SetBinContent(iBin,0);
      histPP->SetBinError(iBin,0);
      histPbPb->SetBinContent(iBin,0);
      histPbPb->SetBinError(iBin,0);
      //continue;
    }
    if (histPrMC->GetBinCenter(iBin)<unfStart) {
      histPrMC->SetBinContent(iBin,0);
      histPrMC->SetBinError(iBin,0);
      histNprMC->SetBinContent(iBin,0);
      histNprMC->SetBinError(iBin,0);
      //continue;
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
  graphPP->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz} [nb]");

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

  histPrMC->SetLineColor(kOrange+2);
  histPrMC->SetLineWidth(2);
  histPrMC->SetMarkerColor(kOrange+2);
  histPrMC->SetMarkerStyle(kFullSquare);
  histPrMC->SetMarkerSize(0);

  histNprMC->SetLineColor(kgreen);
  histNprMC->SetLineWidth(2);
  histNprMC->SetMarkerColor(kgreen);
  histNprMC->SetMarkerStyle(kFullSquare);
  histNprMC->SetMarkerSize(0);

  TLegend* leg = new TLegend(0.7,0.5,0.93,0.65);
  //if (plotMC) leg = new TLegend(0.45,0.5,0.75,0.65);
  if (plotMC) leg = new TLegend(0.4,0.52,0.76,0.65);
  //if (plotMC) leg = new TLegend(0.42,0.77,0.78,0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if (!plotMC) {
    leg->AddEntry(graphPP, "pp","lp");
    leg->AddEntry(graphPbPb, "PbPb","lp");
  }
  else {
    if (plotNPR) {
    leg->AddEntry(graphPP, "Prompt data","lp");
    leg->AddEntry(histPrMC, "Prompt PYTHIA8","lp");
    leg->AddEntry(histNprMC, "Nonprompt PYTHIA8","lp");
    }
    else {
      leg->AddEntry(graphPP, "pp data","lp");
      leg->AddEntry(histPrMC, "PYTHIA 8","lp");
    }
  }

  TLatex *  text = new TLatex(0.2 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  
  TLatex *  text0 = new TLatex(0.38 ,0.82,"Preliminary");
  text0->SetNDC();
  text0->SetTextFont(52);
  text0->SetTextSize(0.055);
  text0->SetLineWidth(2);

  TLatex *  text1 = new TLatex(0.2 , 0.76, "Prompt J/#psi");
  text1->SetNDC();
  text1->SetTextFont(42);
  text1->SetTextSize(0.044);
  text1->SetLineWidth(2);

  TLatex *  text2 = new TLatex(0.2 , 0.72, "p_{T,J/#psi} > 6.5 GeV");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.035);
  text2->SetLineWidth(2);

  TLatex *  text3 = new TLatex(0.2 , 0.67, "30 < p_{T,Jet} < 40 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.035);
  text3->SetLineWidth(2);
  
  TLatex *  text4 = new TLatex(0.2 , 0.62, "|#eta_{Jet}| < 2");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.035);
  text4->SetLineWidth(2);

  TLatex *  text5 = new TLatex(0.2, 0.57,"Cent. 0-90\%");
  text5->SetNDC();
  text5->SetTextFont(42);
  text5->SetTextSize(0.035);
  text5->SetLineWidth(2);

  TLatex *  text6 = new TLatex(0.344, 0.91,"PbPb 1.6 nb^{-1}, pp 302 pb^{-1} (5.02 TeV)");
  if (plotMC) text6 = new TLatex(0.57, 0.91,"pp 302 pb^{-1} (5.02 TeV)");
  text6->SetNDC();
  text6->SetTextFont(42);
  text6->SetTextSize(0.037);
  text6->SetLineWidth(2);

  TCanvas* c = new TCanvas("c","",1000,1000);
  c->cd();
  gPad->SetLeftMargin(0.15);
  TH1D* axisHist = new TH1D("axisHist","",6,0.22,1.);
  axisHist->SetTitle("");
  axisHist->GetXaxis()->SetTitle("z");
  axisHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz} [nb]");
  //cout <<"x axis title size"<<axisHist->GetXaxis()->GetTitleSize()<<endl;
  axisHist->GetXaxis()->SetTitleSize(0.045);
  axisHist->GetYaxis()->SetTitleSize(0.04);
  if (plotMC) axisHist->GetYaxis()->SetTitleSize(0.045);
  axisHist->GetYaxis()->SetTitleOffset(1.6);
  axisHist->GetYaxis()->SetRangeUser(0, histPP->GetMaximum()*1.2);
  if (plotMC) {
    axisHist->GetYaxis()->SetTitle("1/N dN/dz");
    axisHist->GetYaxis()->SetRangeUser(0, histPrMC->GetMaximum()*1.2);
  }
  axisHist->Draw();
  if (plotMC) {
    histPrMC->Draw("same e1");
    histPrMC->Draw("same hist");
    if(plotNPR) {
      histNprMC->Draw("same e1");
      histNprMC->Draw("same hist");
    }
  }
  
  //graphPP->GetYaxis()->SetRangeUser(1,1.2);
  graphPP->GetXaxis()->SetRangeUser(0.22,1.0);
  graphPP->Draw("P");
  if(!plotMC)
    graphPbPbSyst->Draw("5");
  graphPPSyst->Draw("5");
  if(!plotMC)
    graphPbPb->Draw("P");
  graphPP->Draw("P");
  leg->Draw("same");
  text->Draw("same");
  if (!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  if(!plotMC)
  text5->Draw("same");
  text6->Draw("same");
  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());

  fileOutputName.replace(fileOutputName.find(".png"),4,".pdf");
  fileOutputName.replace(fileOutputName.find("PPvsPbPb"),8,"Raa");
  if (plotMC) fileOutputName.replace(fileOutputName.find("vsMC"),4,"");
  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  Double_t eyl_syst[5] = {0,0,0,0,0};
  Double_t eyh_syst[5] = {0,0,0,0,0};
  Double_t exl_syst[5] = {0,0,0,0,0};
  Double_t exh_syst[5] = {0,0,0,0,0};

  for (int iBin = 0; iBin < nBinZ_reco-1; iBin++) {
    int bin = histPbPb->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
    x[iBin] = histPbPb->GetBinCenter(bin);//graphPPSyst->GetPointX(iBin);
    y[iBin] = histPbPb->GetBinContent(bin)/histPP->GetBinContent(histPP->FindBin(x[iBin]));//graphPPSyst->GetPointY(iBin)/graphPbPbSyst->GetPointY(iBin);
    exl[iBin] = z_reco_binWidth/2.;
    exh[iBin] = z_reco_binWidth/2.;
    exl_syst[iBin] = z_reco_binWidth/3.;
    exh_syst[iBin] = z_reco_binWidth/3.;
    eyh_syst[iBin] = y[iBin]*sqrt(pow(graphPbPbSyst->GetErrorYhigh(iBin)/histPbPb->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYhigh(iBin)/histPP->GetBinContent(bin),2));
    eyl_syst[iBin] = y[iBin]*sqrt(pow(graphPbPbSyst->GetErrorYlow(iBin)/histPbPb->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYlow(iBin)/histPP->GetBinContent(bin),2));
    eyh[iBin] = y[iBin]*sqrt(pow(histPbPb->GetBinError(bin)/histPbPb->GetBinContent(bin),2)+pow(histPP->GetBinError(bin)/histPP->GetBinContent(bin),2));
    eyl[iBin] = eyh[iBin];
  }

  TGraphAsymmErrors* graphRaaSyst = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl_syst,exh_syst,eyl_syst,eyh_syst);
  TGraphAsymmErrors* graphRaa = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl,eyh);

  TH1D* taaHist = new TH1D("taaHist","",16,0.22,1.);
  taaHist->GetYaxis()->SetRangeUser(0,1.6);
  taaHist->SetBinContent(taaHist->FindBin(0.99),1);
  taaHist->SetBinError(taaHist->FindBin(0.99),sqrt(taaUnc*taaUnc+lumiPPUnc*lumiPPUnc));
  taaHist->SetFillColorAlpha(kGray, 0.7);
  taaHist->SetMarkerColor(kGray);
  taaHist->SetTitle("");
  taaHist->GetXaxis()->SetTitle("z");
  taaHist->GetYaxis()->SetTitle("R_{AA}");
  taaHist->GetXaxis()->SetTitleSize(0.045);
  taaHist->GetYaxis()->SetTitleSize(0.045);
  taaHist->Draw("e2");
  /*
  TH1D* lumiHist = new TH1D("lumiHist","",16,0.22,1.04);
  lumiHist->GetYaxis()->SetRangeUser(0,1.6);
  lumiHist->SetBinContent(lumiHist->FindBin(1.03),1);
  lumiHist->SetBinError(lumiHist->FindBin(1.03),lumiPPUnc);
  lumiHist->SetFillColorAlpha(kpinkLight, 0.75);
  lumiHist->SetMarkerColor(kpinkLight);
  lumiHist->Draw("same e2");
  */
  graphRaaSyst->SetLineColor(kviolet);
  graphRaaSyst->SetMarkerColor(kviolet);
  graphRaaSyst->SetMarkerStyle(kOpenCircle);
  graphRaaSyst->SetMarkerSize(1.5);
  graphRaaSyst->SetTitle("");
  graphRaaSyst->GetXaxis()->SetTitle("z");
  graphRaaSyst->GetYaxis()->SetTitle("R_{AA}");
  graphRaaSyst->SetFillColorAlpha(kviolet, 0.35);

  graphRaa->SetLineColor(kviolet);
  graphRaa->SetMarkerColor(kviolet);
  graphRaa->SetMarkerStyle(kFullCircle);
  graphRaa->SetMarkerSize(1.5);
  graphRaa->SetTitle("");
  graphRaa->GetXaxis()->SetTitle("z");
  graphRaa->GetYaxis()->SetTitle("R_{AA}");
  graphRaa->SetFillColorAlpha(kviolet, 0.35);
  graphRaaSyst->GetYaxis()->SetRangeUser(0,1.6);
  graphRaaSyst->GetXaxis()->SetRangeUser(0.22,1.0);
  graphRaaSyst->Draw("5");
  graphRaa->Draw("P");
  text->Draw("same");
  if(!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  text5->Draw("same");
  text6->Draw("same");
  TLine* y1line = new TLine(0.22,1,1.0,1);
  y1line->SetLineColor(kGray);
  y1line->SetLineStyle(2);
  y1line->SetLineWidth(2);
  y1line->Draw("same");
  if (!plotMC) {
  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".png"),4,".root");

  TFile* fsave = new TFile(fileOutputName.c_str(),"RECREATE");
  graphRaa->SetLineColor(kBlack);
  graphRaa->SetMarkerColor(kBlack);
  graphRaa->Write("graphRaa");
  graphRaaSyst->SetLineColor(kBlack);
  graphRaaSyst->SetMarkerColor(kBlack);
  graphRaaSyst->Write("graphRaaSyst");
  graphPP->SetLineColor(kBlack);
  graphPP->SetMarkerColor(kBlack);
  graphPP->Write("graphPP");
  graphPPSyst->SetLineColor(kBlack);
  graphPPSyst->SetMarkerColor(kBlack);
  graphPPSyst->Write("graphPPSyst");
  graphPbPb->SetLineColor(kBlack);
  graphPbPb->SetMarkerColor(kBlack);
  graphPbPb->Write("graphPbPb");
  graphPbPbSyst->SetLineColor(kBlack);
  graphPbPbSyst->SetMarkerColor(kBlack);
  graphPbPbSyst->Write("graphPbPbSyst");
  fsave->Close();
  }
}

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  systUncertaintyHistReg(doPbPb, doPrompt, nominalHist);
  systUncertaintyHistTrStat(doPbPb, doPrompt);
  systUncertaintyQuarkoniaSyst(doPbPb, doPrompt);
  systUncertaintyPriorSyst(doPbPb, doPrompt, nominalHist);

  TCanvas* cSyst = new TCanvas("cSyst","",800,800);
  nominalHist->SetLineColor(col[0]);
  nominalHist->SetMarkerColor(col[0]);
  nominalHist->SetMarkerStyle(markerStyle[0]);
  nominalHist->SetStats(0);
  nominalHist->Draw();
  nominalHist->SetTitle("");
  nominalHist->GetXaxis()->SetTitle("z");
  nominalHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz} [nb]");
  if (doPbPb) 
    nominalHist->GetYaxis()->SetTitle("#frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz} [nb]");

  nominalHist->GetYaxis()->SetTitleOffset(1.7);
  cSyst->cd();
  gPad->SetLeftMargin(0.15);

  TCanvas* cSystInd = new TCanvas("cSystInd","",800,800);
  cSystInd->cd();
  gPad->SetLeftMargin(0.15);

  TCanvas* cSystRat = new TCanvas("cSystRat","",800,800);
  cSystRat->cd();
  gPad->SetLeftMargin(0.15);

  double maxSyst = nominalHist->GetMaximum();
  TLegend* legSyst = new TLegend(0.7,0.6,0.9,0.9);
  //if (!doPbPb) legSyst = new TLegend(0.4,0.25,0.6,0.45);
  legSyst->SetBorderSize(0);
  legSyst->SetFillStyle(0);
  legSyst->AddEntry(nominalHist,"nominal","lp");

  TLegend* legRat = new TLegend(0.6,0.5,0.9,0.8);
  if (removePrelim)
    legRat = new TLegend(0.5,0.55,0.9,0.88);

  legRat->SetBorderSize(0);
  legRat->SetFillStyle(0);

  //string systNames [] = {"_JESSyst","_SFSyst","_nominal_centShiftSyst","_Regularization", "_TrStat","_QuarkoniaSyst","_nprPriorSyst"};
  string systNames [] = {"_QuarkoniaSyst","_JESSyst","_SFSyst","_Regularization","_nprPriorSyst","_TrStat","_nominal_centShiftSyst"};

  //string systLabels [] = {"Jet energy scale","Jet energy resolution","Underlying event","Regularization", "Transfer matrix stat.","J/#psi signal extraction","Prior"};
  string systLabels [] = {"J/#psi signal extraction","Jet energy scale","Jet energy resolution","Regularization","Prior","Transfer matrix stat.","Underlying event"};
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
    TH1D* histUpRat = NULL;
    TH1D* histDownRat = NULL;

    if (!histUp || !histDown) {cout<<"[WARNING] systematic histogram not found:"<<systNames[i]<<endl; continue;}

    histUp->SetStats(0);
    histDown->SetStats(0);

    int nB = histUp->GetNbinsX();
    for (int iB=0;iB<nB; iB++) {
      if (histUp->GetBinCenter(iB)<unfStart) {
	histUp->SetBinContent(iB,0);
	histUp->SetBinError(iB,0);
	histDown->SetBinContent(iB,0);
	histDown->SetBinError(iB,0);
	continue;
      }
    }
    
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

    histUpRat = (TH1D*) histUp->Clone("histUpRat"); 
    if (systDiff) {
      histUpRat->Scale(-1);
      histUpRat->Add(nominalHist);
      histUpRat->Scale(-1);
      histUpRat->GetYaxis()->SetTitle("Systematic uncertainty [nb]");
      histUpRat->GetYaxis()->SetRangeUser(-0.05e-3,0.12e-3);
      if(!doPbPb) histUpRat->GetYaxis()->SetRangeUser(-0.1e-3,0.2e-3);
    }
    else {
      histUpRat->Divide(nominalHist); 
      histUpRat->GetYaxis()->SetTitle("syst"); 
      histUpRat->GetYaxis()->SetRangeUser(0.6,1.4); 
      if (doPbPb)
	histUpRat->GetYaxis()->SetRangeUser(0.6,2.); 
    }
    histUpRat->SetTitle("");
    histUpRat->GetXaxis()->SetTitle("z");
    histUpRat->GetXaxis()->SetTitleSize(0.045);
    histUpRat->GetYaxis()->SetTitleSize(0.045); 
    histDownRat = (TH1D*) histDown->Clone("histDownRat"); 
    if (systDiff) {
      histDownRat->Scale(-1);
      histDownRat->Add(nominalHist);
      histDownRat->Scale(-1);
    }
    else 
      histDownRat->Divide(nominalHist);
    histUpRat->SetLineStyle(1);
    histUpRat->SetLineWidth(2);
    histDownRat->SetLineStyle(1);
    histDownRat->SetLineWidth(2);
    cSyst->cd();
    histUp->Draw("same");
    histDown->Draw("same");
    cSystRat->cd();
    histUpRat->Draw("same hist");
    histDownRat->Draw("same hist");
    cSystInd->cd();
    nominalHist->Draw();
    histUp->Draw("same");
    histDown->Draw("same");
    if (!plotMC)
      cSystInd->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp,systNames[i].c_str()));
    
    legSyst->AddEntry(histUp,Form("%s",systLabels[i].c_str()),"lp");
    legRat->AddEntry(histUpRat,Form("%s",systLabels[i].c_str()),"l");
    //legSyst->AddEntry(histDown,Form("%s Down",systNames[i].c_str()),"lp");
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
      //cout << "x = "<<x[iBin]<<", y = "<<y[iBin]<<", eyl = "<<eyl[iBin]<<", eyh = "<<eyh[iBin]<<endl;
    }
  }
  cSyst->cd();
  legSyst->Draw("same");
  if(!plotMC) {
    cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
  cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp));
  }
  cSystRat->cd();
  TLine* y1line = new TLine(min_z,1,1.0,1);
  if (systDiff) y1line = new TLine(min_z,0,1.0,0);
  y1line->SetLineColor(kBlack);
  y1line->SetLineStyle(1);
  y1line->SetLineWidth(2);
  y1line->Draw("same");

  TLatex *  text = new TLatex(0.2 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(61);
  text->SetTextSize(0.075);
  text->SetLineWidth(8);
  text->Draw("same");  

  TLatex *  text0 = new TLatex(0.38 ,0.82,"Preliminary");
  text0->SetNDC();
  text0->SetTextFont(52);
  text0->SetTextSize(0.055);
  text0->SetLineWidth(2);
  if (!removePrelim)
    text0->Draw("same");  

  TLatex *  text1 = new TLatex(0.2 , 0.76, "Prompt J/#psi");
  text1->SetNDC();
  text1->SetTextFont(42);
  text1->SetTextSize(0.044);
  text1->SetLineWidth(2);
  text1->Draw("same");  

  TLatex *  text2 = new TLatex(0.2 , 0.72, "p_{T,J/#psi} > 6.5 GeV");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.035);
  text2->SetLineWidth(2);
  text2->Draw("same");  

  TLatex *  text3 = new TLatex(0.2 , 0.67, "30 < p_{T,Jet} < 40 GeV");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.035);
  text3->SetLineWidth(2);
  text3->Draw("same");  

  TLatex *  text4 = new TLatex(0.2 , 0.62, "|#eta_{Jet}| < 2");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.035);
  text4->SetLineWidth(2);
  text4->Draw("same");  

  TLatex *  text5 = new TLatex(0.2, 0.57,"Cent. 0-90\%");
  text5->SetNDC();
  text5->SetTextFont(42);
  text5->SetTextSize(0.035);
  text5->SetLineWidth(2);
  if(doPbPb)
    text5->Draw("same");  

  TLatex *  text6 = new TLatex(0.57, 0.91,"pp 302 pb^{-1} (5.02 TeV)");
  if (doPbPb) text6 = new TLatex(0.533, 0.91,"PbPb 1.6 nb^{-1} (5.02 TeV)");
  text6->SetNDC();
  text6->SetTextFont(42);
  text6->SetTextSize(0.037);
  text6->SetLineWidth(2);
  text6->Draw("same");

  TGraphAsymmErrors* systHist = new TGraphAsymmErrors(nBinZ_reco-1,x,y,exl,exh,eyl,eyh);

  Double_t yR[5] = {0,0,0,0,0};
  Double_t exlR[5] = {z_reco_binWidth/2.,z_reco_binWidth/2.,z_reco_binWidth/2.,z_reco_binWidth/2.,z_reco_binWidth/2.};
  TGraphAsymmErrors* systHistNorm = new TGraphAsymmErrors(nBinZ_reco-1,x,yR,exlR,exlR,eyl,eyh);
  systHistNorm->SetLineColor(kBlack);
  systHistNorm->SetLineWidth(2);
  systHistNorm->SetMarkerColor(kBlack);
  systHistNorm->SetMarkerStyle(kOpenCircle);
  systHistNorm->SetMarkerSize(0);
  systHistNorm->SetFillColorAlpha(kWhite, 0);

  systHistNorm->Draw("5");
  legRat->AddEntry(systHistNorm,"Total","l");
  legRat->Draw("same");
  if(!plotMC)
    cSystRat->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/Systematic%sUncertainty_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s.pdf",unfPath.c_str(),systDiff?"Diff":"Rel",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,doPbPb?nSIter:nSIter_pp,removePrelim?"_noPreliminaryLabel":""));
  
  return systHist;
}

void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", doPbPb?nIter_syst:nIter_syst_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter_syst:nSIter_syst_pp, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
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
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularization");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_Regularization");
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
  cout <<"reading systematics from file "<<fileSystName<<endl;
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for trM stat syst"<<endl; return;}
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
void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_systError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for quarkonia syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    histUp->SetBinContent(i, histUp->GetBinContent(i)+histUp->GetBinError(i));
    histDown->SetBinContent(i, histDown->GetBinContent(i)-histDown->GetBinError(i));
  }

  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_QuarkoniaSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_QuarkoniaSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}

void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_nprPriorSyst");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
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
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_nprPriorSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, doPbPb?nSIter:nSIter_pp, "_nprPriorSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}
