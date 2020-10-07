#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

double normCent = 2.18414e-08;//2.28119e-08;
double normPeri = 4.2996e-08;//4.63974e-08;

double taaUncCent = 0.019;
double taaUncPeri = 0.036;
bool removePrelim = false;

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri);
void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri);
void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt, bool isCent, bool isPeri);
void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt, bool isCent, bool isPeri);
void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri);
void plotCentResults_oneStyle(bool doPrompt = true);

void plotCentResults() {
  //removePrelim = false;
  //plotCentResults_oneStyle(true);
  removePrelim = true;
  plotCentResults_oneStyle(true);
}
void plotCentResults_oneStyle(bool doPrompt) {
  gStyle->SetOptStat(0);
  unfStart=0.376;

  string filePPName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PP_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_pp, "_nominal");
  TFile* filePP = TFile::Open(filePPName.c_str());
  TH1D* histPP = (TH1D*) filePP->Get(Form("zUnfSI%d",nSIter_pp));
  histPP->Scale(normPP*1./z_reco_binWidth);
  cout <<"got the pp file"<<endl;
  string filePbPbName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter, "_nominal");
  TFile* filePbPb = TFile::Open(filePbPbName.c_str());
  TH1D* histPbPb = (TH1D*) filePbPb->Get(Form("zUnfSI%d",nSIter));
  histPbPb->Scale(normPbPb*1./z_reco_binWidth);
  cout <<"got the pbpb file"<<endl;
  string fileCentName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_cent, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_cent, "_centBin");
  TFile* fileCent = TFile::Open(fileCentName.c_str());
  cout <<"got the pbpb cent file"<<endl;
  TH1D* histCent = (TH1D*) fileCent->Get(Form("zUnfSI%d",nSIter_cent));
  histCent->Scale(normCent*1./z_reco_binWidth);
  cout <<"got the pbpb cent hist"<<endl;
  string filePeriName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PbPb_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nIter_peri, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nSIter_peri, "_periBin");
  TFile* filePeri = TFile::Open(filePeriName.c_str());
  cout <<"got the pbpb peri file"<<endl;
  TH1D* histPeri = (TH1D*) filePeri->Get(Form("zUnfSI%d",nSIter_peri));
  histPeri->Scale(normPeri*1./z_reco_binWidth);
  cout <<"got the pbpb peri hist"<<endl;
  string fileOutputName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_PPvsPbPbvsCent_%s_%dz%dptBins%dz%dptMeasBin_PbPbCentnIter%inSIter%i_PbPbPerinIter%inSIter%i_PPnIter%inSIter%i_statError%s.pdf",unfPath.c_str(),doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, nIter_cent, nSIter_cent, nIter_peri, nSIter_peri, nIter_pp, nSIter_pp,removePrelim?"_noPreliminaryLabel":"");

  cout <<"got the files and the hists"<<endl;
  
  int nBin = histPP->GetNbinsX();
  for (int iBin=0;iBin<nBin; iBin++) {
    if (histPP->GetBinCenter(iBin)<unfStart) {
      histPP->SetBinContent(iBin,0);
      histPP->SetBinError(iBin,0);
      histPbPb->SetBinContent(iBin,0);
      histPbPb->SetBinError(iBin,0);
      histCent->SetBinContent(iBin,0);
      histCent->SetBinError(iBin,0);
      histPeri->SetBinContent(iBin,0);
      histPeri->SetBinError(iBin,0);
      continue;
    }
  }

  histPP->SetTitle("");
  histPbPb->SetTitle("");
  histCent->SetTitle("");
  histPeri->SetTitle("");
  TGraphAsymmErrors* graphPP = new TGraphAsymmErrors(histPP);//(nBin,x_pp,y_pp,exl_pp,exh_pp,eyl_pp,eyh_pp);
  TGraphAsymmErrors* graphPbPb = new TGraphAsymmErrors(histPbPb);//(nBin,x_pbpb,y_pbpb,exl_pbpb,exh_pbpb,eyl_pbpb,eyh_pbpb);
  TGraphAsymmErrors* graphCent = new TGraphAsymmErrors(histCent);//(nBin,x_pbpb,y_pbpb,exl_pbpb,exh_pbpb,eyl_pbpb,eyh_pbpb);
  TGraphAsymmErrors* graphPeri = new TGraphAsymmErrors(histPeri);//(nBin,x_pbpb,y_pbpb,exl_pbpb,exh_pbpb,eyl_pbpb,eyh_pbpb);

  cout <<"got all the nominal hitograms"<<endl;
  TGraphAsymmErrors* graphPPSyst = systUncertaintyHistAll(false, doPrompt, histPP, false, false);//new TGraphAsymmErrors(histPPSyst);
  TGraphAsymmErrors* graphPbPbSyst = systUncertaintyHistAll(true, doPrompt, histPbPb, false, false);//new TGraphAsymmErrors(histPbPbSyst);
  cout <<"got all the inclusive hitograms"<<endl;
  TGraphAsymmErrors* graphCentSyst = systUncertaintyHistAll(true, doPrompt, histCent, true, false);//new TGraphAsymmErrors(histPbPbSyst);
  cout <<"got all the central hitograms"<<endl;
  TGraphAsymmErrors* graphPeriSyst = systUncertaintyHistAll(true, doPrompt, histPeri, false, true);//new TGraphAsymmErrors(histPbPbSyst);
  
  cout <<"got all the hitograms"<<endl;
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

  graphCent->SetLineColor(kred);
  graphCent->SetMarkerColor(kred);
  graphCent->SetMarkerStyle(kFullCircle);
  graphCent->SetMarkerSize(1.5);

  graphPeri->SetLineColor(kyellow);
  graphPeri->SetMarkerColor(kyellow);
  graphPeri->SetMarkerStyle(kFullCircle);
  graphPeri->SetMarkerSize(1.5);
  
  graphPPSyst->SetLineColor(kpink);
  graphPPSyst->SetMarkerColor(kpink);
  graphPPSyst->SetMarkerStyle(kFullSquare);
  graphPPSyst->SetFillColorAlpha(kpinkLight, 0.35);

  graphPbPbSyst->SetLineColor(kazure);
  graphPbPbSyst->SetMarkerColor(kazure);
  graphPbPbSyst->SetMarkerStyle(kFullCircle);
  graphPbPbSyst->SetFillColorAlpha(kazureLight, 0.35);

  graphCentSyst->SetLineColor(kred);
  graphCentSyst->SetMarkerColor(kred);
  graphCentSyst->SetMarkerStyle(kFullCircle);
  graphCentSyst->SetFillColorAlpha(kredLight, 0.35);

  graphPeriSyst->SetLineColor(kyellow);
  graphPeriSyst->SetMarkerColor(kyellow);
  graphPeriSyst->SetMarkerStyle(kFullCircle);
  graphPeriSyst->SetFillColorAlpha(kyellowLight, 0.35);
  
  TLegend* leg = new TLegend(0.7,0.5,0.93,0.65);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(graphPP, "pp","lp");
  leg->AddEntry(graphPbPb, "PbPb 0-90%","lp");
  leg->AddEntry(graphCent, "PbPb 0-20%","lp");
  leg->AddEntry(graphPeri, "PbPb 20-90%","lp");

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
  text6->SetNDC();
  text6->SetTextFont(42);
  text6->SetTextSize(0.037);
  text6->SetLineWidth(2);

  TCanvas* c = new TCanvas("c","",1000,1000);
  c->cd();
  gPad->SetLeftMargin(0.15);
  TH1D* axisHist = new TH1D("axisHist","",6,0.22,1.0);
  axisHist->SetTitle("");
  axisHist->GetXaxis()->SetTitle("z");
  axisHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}, #frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz}");
  axisHist->GetYaxis()->SetTitleOffset(1.6);
  axisHist->GetYaxis()->SetRangeUser(0, histPP->GetMaximum()*1.1);
  axisHist->Draw();
  //graphPP->GetYaxis()->SetRangeUser(1,1.2);
  graphPP->GetXaxis()->SetRangeUser(0.22,1.0);
  graphPP->Draw("P");
  graphPbPbSyst->Draw("5");
  graphPPSyst->Draw("5");
  graphPbPb->Draw("P");
  graphCentSyst->Draw("5");
  graphCent->Draw("P");
  graphPeriSyst->Draw("5");
  graphPeri->Draw("P");
  graphPP->Draw("P");
  leg->Draw("same");
  text->Draw("same");
  if (!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  //text5->Draw("same");
  text6->Draw("same");
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
  Double_t exl_syst[5] = {0,0,0,0,0};
  Double_t exh_syst[5] = {0,0,0,0,0};

  for (int iBin = 0; iBin < nBinZ_reco-2; iBin++) {
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

  TGraphAsymmErrors* graphRaaSyst = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl_syst,exh_syst,eyl_syst,eyh_syst);
  TGraphAsymmErrors* graphRaa = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl,exh,eyl,eyh);

  for (int iBin = 0; iBin < nBinZ_reco-2; iBin++) {
    int bin = histCent->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
    x[iBin] = histCent->GetBinCenter(bin);//graphPPSyst->GetPointX(iBin);
    y[iBin] = histCent->GetBinContent(bin)/histPP->GetBinContent(histPP->FindBin(x[iBin]));//graphPPSyst->GetPointY(iBin)/graphPbPbSyst->GetPointY(iBin);
    exl[iBin] = z_reco_binWidth/2.;
    exh[iBin] = z_reco_binWidth/2.;
    eyh_syst[iBin] = y[iBin]*sqrt(pow(graphCentSyst->GetErrorYhigh(iBin)/histCent->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYhigh(iBin)/histPP->GetBinContent(bin),2));
    eyl_syst[iBin] = y[iBin]*sqrt(pow(graphCentSyst->GetErrorYlow(iBin)/histCent->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYlow(iBin)/histPP->GetBinContent(bin),2));
    eyh[iBin] = y[iBin]*sqrt(pow(histCent->GetBinError(bin)/histCent->GetBinContent(bin),2)+pow(histPP->GetBinError(bin)/histPP->GetBinContent(bin),2));
    eyl[iBin] = eyh[iBin];
  }
  TGraphAsymmErrors* graphRaaCentSyst = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl_syst,exh_syst,eyl_syst,eyh_syst);
  TGraphAsymmErrors* graphRaaCent = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl,exh,eyl,eyh);
  for (int iBin = 0; iBin < nBinZ_reco-2; iBin++) {
    int bin = histPeri->FindBin(iBin*z_reco_binWidth+min_z+unfStart);
    x[iBin] = histPeri->GetBinCenter(bin);//graphPPSyst->GetPointX(iBin);
    y[iBin] = histPeri->GetBinContent(bin)/histPP->GetBinContent(histPP->FindBin(x[iBin]));//graphPPSyst->GetPointY(iBin)/graphPbPbSyst->GetPointY(iBin);
    exl[iBin] = z_reco_binWidth/2.;
    exh[iBin] = z_reco_binWidth/2.;
    eyh_syst[iBin] = y[iBin]*sqrt(pow(graphPeriSyst->GetErrorYhigh(iBin)/histPeri->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYhigh(iBin)/histPP->GetBinContent(bin),2));
    eyl_syst[iBin] = y[iBin]*sqrt(pow(graphPeriSyst->GetErrorYlow(iBin)/histPeri->GetBinContent(bin),2)+pow(graphPPSyst->GetErrorYlow(iBin)/histPP->GetBinContent(bin),2));
    eyh[iBin] = y[iBin]*sqrt(pow(histPeri->GetBinError(bin)/histPeri->GetBinContent(bin),2)+pow(histPP->GetBinError(bin)/histPP->GetBinContent(bin),2));
    eyl[iBin] = eyh[iBin];
  }
  TGraphAsymmErrors* graphRaaPeriSyst = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl_syst,exh_syst,eyl_syst,eyh_syst);
  TGraphAsymmErrors* graphRaaPeri = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl,exh,eyl,eyh);

  TH1D* taaHist = new TH1D("taaHist","",16,0.22,1.0);
  taaHist->GetYaxis()->SetRangeUser(0,1.6);
  taaHist->SetBinContent(taaHist->FindBin(0.99)-1,1);
  taaHist->SetBinError(taaHist->FindBin(0.99)-1,taaUnc);
  taaHist->SetFillColorAlpha(kazureLight, 0.75);
  taaHist->SetMarkerColor(kazureLight);
  taaHist->SetTitle("");
  taaHist->GetXaxis()->SetTitle("z");
  taaHist->GetYaxis()->SetTitle("R_{AA}");
  //taaHist->Draw("e2");

  TH1D* taaHistCent = new TH1D("taaHistCent","",16,0.22,1.0);
  taaHistCent->GetYaxis()->SetRangeUser(0,1.6);
  taaHistCent->SetBinContent(taaHistCent->FindBin(0.99)-1,1);
  taaHistCent->SetBinError(taaHistCent->FindBin(0.99)-1,sqrt(taaUncCent*taaUncCent+lumiPPUnc*lumiPPUnc));
  taaHistCent->SetFillColorAlpha(kredLight, 0.75);
  taaHistCent->SetMarkerColor(kredLight);
  taaHistCent->SetTitle("");
  taaHistCent->GetXaxis()->SetTitle("z");
  taaHistCent->GetYaxis()->SetTitle("R_{AA}");
  taaHistCent->GetXaxis()->SetTitleSize(0.045);
  taaHistCent->GetYaxis()->SetTitleSize(0.045);
  taaHistCent->Draw("e2");

  TH1D* taaHistPeri = new TH1D("taaHistPeri","",16,0.22,1.);
  taaHistPeri->GetYaxis()->SetRangeUser(0,1.6);
  taaHistPeri->SetBinContent(taaHistPeri->FindBin(0.99),1);
  taaHistPeri->SetBinError(taaHistPeri->FindBin(0.99),sqrt(taaUncPeri*taaUncPeri+lumiPPUnc*lumiPPUnc));
  taaHistPeri->SetFillColorAlpha(kyellowLight, 0.75);
  taaHistPeri->SetMarkerColor(kyellowLight);
  taaHistPeri->SetTitle("");
  taaHistPeri->GetXaxis()->SetTitle("z");
  taaHistPeri->GetYaxis()->SetTitle("R_{AA}");
  taaHistPeri->Draw("e2 same");
  /*
  TH1D* lumiHist = new TH1D("lumiHist","",16,0.2,1.04);
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
  graphRaaSyst->SetFillColorAlpha(kvioletLight, 0.15);


  graphRaa->SetLineColor(kviolet);
  graphRaa->SetMarkerColor(kviolet);
  graphRaa->SetMarkerStyle(kFullCircle);
  graphRaa->SetMarkerSize(1.5);
  graphRaa->SetTitle("");
  graphRaa->GetXaxis()->SetTitle("z");
  graphRaa->GetYaxis()->SetTitle("R_{AA}");
  graphRaa->SetFillColorAlpha(kviolet, 0.35);

  graphRaaCent->SetLineColor(kred);
  graphRaaCent->SetMarkerColor(kred);
  graphRaaCent->SetMarkerStyle(kFullSquare);
  graphRaaCent->SetMarkerSize(1.5);
  graphRaaCent->SetTitle("");
  graphRaaCent->GetXaxis()->SetTitle("z");
  graphRaaCent->GetYaxis()->SetTitle("R_{AA}");
  graphRaaCent->SetFillColorAlpha(kred, 0.35);

  graphRaaCentSyst->SetLineColor(kred);
  graphRaaCentSyst->SetMarkerColor(kred);
  graphRaaCentSyst->SetMarkerStyle(kFullSquare);
  graphRaaCentSyst->SetMarkerSize(1.5);
  graphRaaCentSyst->SetTitle("");
  graphRaaCentSyst->GetXaxis()->SetTitle("z");
  graphRaaCentSyst->GetYaxis()->SetTitle("R_{AA}");
  graphRaaCentSyst->SetFillColorAlpha(kredLight, 0.15);

  graphRaaPeri->SetLineColor(kyellow);
  graphRaaPeri->SetMarkerColor(kyellow);
  graphRaaPeri->SetMarkerStyle(kFullCircle);
  graphRaaPeri->SetMarkerSize(1.5);
  graphRaaPeri->SetTitle("");
  graphRaaPeri->GetXaxis()->SetTitle("z");
  graphRaaPeri->GetYaxis()->SetTitle("R_{AA}");
  graphRaaPeri->SetFillColorAlpha(kyellow, 0.35);

  graphRaaPeriSyst->SetLineColor(kyellow);
  graphRaaPeriSyst->SetMarkerColor(kyellow);
  graphRaaPeriSyst->SetMarkerStyle(kFullCircle);
  graphRaaPeriSyst->SetMarkerSize(1.5);
  graphRaaPeriSyst->SetTitle("");
  graphRaaPeriSyst->GetXaxis()->SetTitle("z");
  graphRaaPeriSyst->GetYaxis()->SetTitle("R_{AA}");
  graphRaaPeriSyst->SetFillColorAlpha(kyellowLight, 0.15);

  TLegend* legRaa = new TLegend(0.7,0.65,0.93,0.8);
  legRaa->SetBorderSize(0);
  legRaa->SetFillStyle(0);
  //legRaa->AddEntry(graphRaa, "0-90%","lp");
  legRaa->AddEntry(graphRaaCent, "0-20%","lp");
  legRaa->AddEntry(graphRaaPeri, "20-90%","lp");

  graphRaaSyst->GetYaxis()->SetRangeUser(0,1.6);
  graphRaaSyst->GetXaxis()->SetRangeUser(0.22,1.0);
  //graphRaaSyst->Draw("5");
  //graphRaa->Draw("P");
  graphRaaCentSyst->Draw("5");
  graphRaaCent->Draw("P");
  graphRaaPeriSyst->Draw("5");
  graphRaaPeri->Draw("P");
  text->Draw("same");
  if(!removePrelim)
    text0->Draw("same");
  text1->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text4->Draw("same");
  //text5->Draw("same");
  text6->Draw("same");
  legRaa->Draw("same");
  TLine* y1line = new TLine(0.22,1,1.0,1);
  y1line->SetLineColor(kGray);
  y1line->SetLineStyle(2);
  y1line->SetLineWidth(2);
  y1line->Draw("same");

  c->SaveAs(fileOutputName.c_str());
  fileOutputName.replace(fileOutputName.find(".pdf"),4,".png");
  c->SaveAs(fileOutputName.c_str());
}

TGraphAsymmErrors* systUncertaintyHistAll(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri) {
  gStyle->SetOptStat(0);
  cout<<"reading systematics"<<endl;
  systUncertaintyHistReg(doPbPb, doPrompt, nominalHist, isCent, isPeri);
  cout<<"reading reg systematics done"<<endl;
  systUncertaintyHistTrStat(doPbPb, doPrompt, isCent, isPeri);
  cout<<"reading trStat systematics done"<<endl;
  systUncertaintyQuarkoniaSyst(doPbPb, doPrompt, isCent, isPeri);
  cout<<"reading quarkonia systematics done"<<endl;
  systUncertaintyPriorSyst(doPbPb, doPrompt, nominalHist, isCent, isPeri);
  cout<<"reading prior systematics done"<<endl;

  TCanvas* cSyst = new TCanvas("cSyst","",800,800);
  nominalHist->SetLineColor(col[0]);
  nominalHist->SetMarkerColor(col[0]);
  nominalHist->SetMarkerStyle(markerStyle[0]);
  nominalHist->SetStats(0);
  nominalHist->Draw();
  nominalHist->SetTitle("");
  nominalHist->GetXaxis()->SetTitle("z");
  nominalHist->GetYaxis()->SetTitle("#frac{d#sigma^{pp}}{dz}");
  if (doPbPb) 
    nominalHist->GetYaxis()->SetTitle("#frac{1}{N_{evt} T_{AA}} #frac{dN^{PbPb}}{dz}");

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

  TLegend* legRat = new TLegend(0.7,0.6,0.9,0.9);
  legRat->SetBorderSize(0);
  legRat->SetFillStyle(0);

  string systNames [] = {"_JESSyst","_SFSyst","_nominal_centShiftSyst","_Regularisation", "_TrStat","_QuarkoniaSyst","_nprPriorSyst"};
  string systLabels [] = {"JES","JER","cent. Shift","Regularisation", "Tr. M. Stat","Quarkonia","Prior"};
  //if (isCent) systNames = {"_centBin_JESSyst","_centBin_SFSyst","_centBin_centShiftSyst","_centBin_Regularisation", "_centBin_TrStat","_centBin_QuarkoniaSyst"};
  //else if (isPeri) systNames = {"_periBin_JESSyst","_periBin_SFSyst","_periBin_centShiftSyst","_periBin_Regularisation", "_periBin_TrStat","_periBin_QuarkoniaSyst"};

  int nSyst = sizeof(systNames)/sizeof(systNames[0]);
  Double_t x[5] = {0,0,0,0,0};
  Double_t y[5] = {0,0,0,0,0};
  Double_t exl[5] = {0,0,0,0,0};
  Double_t eyl[5] = {0,0,0,0,0};
  Double_t exh[5] = {0,0,0,0,0};
  Double_t eyh[5] = {0,0,0,0,0};
  
  for (int i=0; i<nSyst; i++) {
    if ((isCent || isPeri) && i==2) systNames[i] = "_centShiftSyst";
    if (isCent) systNames[i] = "_centBin"+systNames[i];
    if (isPeri) systNames[i] = "_periBin"+systNames[i];
    string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, systNames[i].c_str());
    TFile* fileUp = TFile::Open(fileUpName.c_str());
    TFile* fileDown = TFile::Open(fileDownName.c_str());
    if (!fileUp || !fileDown) {cout<<"[WARNING] systematic file not found:"<<systNames[i]<<endl; continue;}
    TH1D* histUp = (TH1D*) fileUp->Get(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
    TH1D* histDown = (TH1D*) fileDown->Get(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
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
    
    histUp->Scale(isCent?(normCent*1./z_reco_binWidth):isPeri?(normPeri*1./z_reco_binWidth):doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
    histDown->Scale(isCent?(normCent*1./z_reco_binWidth):isPeri?(normPeri*1./z_reco_binWidth):doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));

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
    histUpRat = (TH1D*) histUp->Clone("histUpRat"); histUpRat->Divide(nominalHist); histUpRat->GetYaxis()->SetTitle("syst"); 
    histUpRat->GetYaxis()->SetRangeUser(0.6,1.4); 
    if (doPbPb)
      histUpRat->GetYaxis()->SetRangeUser(0.6,2.); 
    histUpRat->SetTitle("");
    histUpRat->GetXaxis()->SetTitle("z");
    histDownRat = (TH1D*) histDown->Clone("histDownRat"); histDownRat->Divide(nominalHist);
    histUpRat->SetLineStyle(1);
    histDownRat->SetLineStyle(1);
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
    //cSystInd->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp,systNames[i].c_str()));
    
    legSyst->AddEntry(histUp,Form("%s",systLabels[i].c_str()),"lp");
    legRat->AddEntry(histUpRat,Form("%s",systLabels[i].c_str()),"l");
    //legSyst->AddEntry(histDown,Form("%s Down",systNames[i].c_str()),"lp");
    for (int iBin = 0; iBin < nBinZ_reco-2; iBin++) { //-1 for the underflow
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
      //cout << "systNames = "<<systNames[i]<<", x = "<<x[iBin]<<", y = "<<y[iBin]<<", eyl = "<<eyl[iBin]<<", eyh = "<<eyh[iBin]<<endl;
    }
  }
  cSyst->cd();
  legSyst->Draw("same");
  cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s_SIter%i.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,isCent?"_centBin":isPeri?"_periBin":"",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  //cSyst->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  cSystRat->cd();
  legRat->Draw("same");
  TLine* y1line = new TLine(min_z,1,1.0,1);
  y1line->SetLineColor(kBlack);
  y1line->SetLineStyle(1);
  y1line->SetLineWidth(2);
  y1line->Draw("same");

  //cSystRat->SaveAs(Form("%s/dataUnf/unfOutput/finalResults/SystematicRelUncertainty_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco,isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  TGraphAsymmErrors* systHist = new TGraphAsymmErrors(nBinZ_reco-2,x,y,exl,exh,eyl,eyh);
  return systHist;
}

void systUncertaintyHistReg(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_syst_cent:isPeri?nIter_syst_peri:doPbPb?nIter_syst:nIter_syst_pp, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_syst_cent:isPeri?nSIter_syst_peri:doPbPb?nSIter_syst:nSIter_syst_pp, isCent?"_centBin":isPeri?"_periBin":"_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",isCent?nSIter_syst_cent:isPeri?nSIter_syst_peri:doPbPb?nSIter_syst:nSIter_syst_pp));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  TH1D* nomHist = (TH1D*) nominalHist->Clone("nomHist");
  //systHist->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
  nomHist->Scale(isCent?(z_reco_binWidth*1./normCent):isPeri?(z_reco_binWidth*1./normPeri):doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    double binCenter = systHist->GetBinCenter(i);
    double err = fabs(nomHist->GetBinContent(nomHist->FindBin(binCenter))-systHist->GetBinContent(i));
    histUp->SetBinContent(i, nomHist->GetBinContent(i)+err);
    histDown->SetBinContent(i, nomHist->GetBinContent(i)-err);
  }
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_Regularisation":isPeri?"_periBin_Regularisation":"_Regularisation");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_Regularisation":isPeri?"_periBin_Regularisation":"_Regularisation");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}

void systUncertaintyHistTrStat(bool doPbPb, bool doPrompt, bool isCent, bool isPeri) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/matrixOper/systUnc_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?"_centBin":isPeri?"_periBin":"_nominal");
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
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_TrStat":isPeri?"_periBin_TrStat":"_TrStat");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_TrStat":isPeri?"_periBin_TrStat":"_TrStat");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}
void systUncertaintyQuarkoniaSyst(bool doPbPb, bool doPrompt, bool isCent, bool isPeri) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_systError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin":isPeri?"_periBin":"_nominal");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for quarkonia syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    histUp->SetBinContent(i, histUp->GetBinContent(i)+histUp->GetBinError(i));
    histDown->SetBinContent(i, histDown->GetBinContent(i)-histDown->GetBinError(i));
  }

  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_QuarkoniaSyst":isPeri?"_periBin_QuarkoniaSyst":"_QuarkoniaSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_QuarkoniaSyst":isPeri?"_periBin_QuarkoniaSyst":"_QuarkoniaSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}

void systUncertaintyPriorSyst(bool doPbPb, bool doPrompt, TH1D* nominalHist, bool isCent, bool isPeri) {
  gStyle->SetOptStat(0);
  string fileSystName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%s_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_nprPriorSyst":isPeri?"_periBin_nprPriorSyst":"_nprPriorSyst");
  TFile* fileSyst = TFile::Open(fileSystName.c_str());
  if (!fileSyst) {cout<<"[WARNING] file "<<fileSystName<<" not found for regularization syst"<<endl; return;}
  TH1D* systHist = (TH1D*) fileSyst->Get(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  TH1D* histUp = (TH1D*) systHist->Clone("histUp");
  TH1D* histDown = (TH1D*) systHist->Clone("histDown");
  TH1D* nomHist = (TH1D*) nominalHist->Clone("nomHist");
  //systHist->Scale(doPbPb?(normPbPb*1./z_reco_binWidth):(normPP*1./z_reco_binWidth));
  //nomHist->Scale(doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  nomHist->Scale(isCent?(z_reco_binWidth*1./normCent):isPeri?(z_reco_binWidth*1./normPeri):doPbPb?(z_reco_binWidth*1./normPbPb):(z_reco_binWidth*1./normPP));
  int nBin = systHist->GetNbinsX();
  for (int i=0; i<=nBin;i++) {
    double binCenter = systHist->GetBinCenter(i);
    double err = fabs(nomHist->GetBinContent(nomHist->FindBin(binCenter))-systHist->GetBinContent(i));
    histUp->SetBinContent(i, nomHist->GetBinContent(i)+err);
    histDown->SetBinContent(i, nomHist->GetBinContent(i)-err);
  }
  
  string fileUpName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sUp_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_nprPriorSyst":isPeri?"_periBin_nprPriorSyst":"_nprPriorSyst");
  string fileDownName = Form("%s/dataUnf/unfOutput/finalResults/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin_SIter%i%sDown_statError.root",unfPath.c_str(),doPbPb?"PbPb":"PP", doPrompt?"prompt":"nonprompt", isCent?nIter_cent:isPeri?nIter_peri:nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp, isCent?"_centBin_nprPriorSyst":isPeri?"_periBin_nprPriorSyst":"_nprPriorSyst");
  TFile* fileUp = new TFile(fileUpName.c_str(),"RECREATE");
  histUp->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  nomHist->Write("nominal");
  systHist->Write("syst");
  fileUp->Close();
  TFile* fileDown = new TFile(fileDownName.c_str(),"RECREATE");
  histDown->Write(Form("zUnfSI%d",isCent?nSIter_cent:isPeri?nSIter_peri:doPbPb?nSIter:nSIter_pp));
  fileDown->Close();
}
