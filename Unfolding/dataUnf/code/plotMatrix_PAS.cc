#include "inputParams.h"

void plot(bool doPrompt = true, bool doPbPb = true){
  gStyle->SetOptStat(0);
  int iterFinal=4;//nSIter;

  string filename = "";
  string outputfile = "";
  filename = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput_cetnShiftBug/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%s.root",iterFinal,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal");
    
  TFile *file = new TFile(filename.c_str());
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get("invtrmat3");
  
  TH2D *h2_UnfData  = (TH2D*)file->Get("hReco_Iter3");
  TH2D *h2_MeasData  = (TH2D*)file->Get("fh2MeasData");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  
  //TH2D *h_covMatrix = new TH2D(*covMatrix);
  h_invMatrix->Rebin2D(1,nBinJet_gen*nBinZ_gen/(nBinJet_reco*nBinZ_reco));
  //h_invMatrix->GetZaxis()->SetRangeUser(0,1);

  TLine *line0 = new TLine(6,0,6,1200);
  line0->SetLineColor(kBlack);
  line0->SetLineStyle(1);
  line0->SetLineWidth(2);

  TLine *line1 = new TLine(12,0,12,1200);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(1);
  line1->SetLineWidth(2);

  TLine *line2 = new TLine(18,0,18,1200);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);

  TLine *line3 = new TLine(24,0,24,1200);
  line3->SetLineColor(kBlack);
  line3->SetLineStyle(1);
  line3->SetLineWidth(2);

  TLine *line4 = new TLine(0,240,30,240);
  line4->SetLineColor(kBlack);
  line4->SetLineStyle(1);
  line4->SetLineWidth(2);

  TLine *line5 = new TLine(0,480,30,480);
  line5->SetLineColor(kBlack);
  line5->SetLineStyle(1);
  line5->SetLineWidth(2);

  TLine *line6 = new TLine(0,720,30,720);
  line6->SetLineColor(kBlack);
  line6->SetLineStyle(1);
  line6->SetLineWidth(2);

  TLine *line7 = new TLine(0,960,30,960);
  line7->SetLineColor(kBlack);
  line7->SetLineStyle(1);
  line7->SetLineWidth(2);

  TCanvas *can = new TCanvas("can","can",600,600);

  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.25);

  gStyle->SetPalette(56);

  
  h_invMatrix->GetXaxis()->SetTickLength(0.02);
  h_invMatrix->GetYaxis()->SetTickLength(0.02);
  

  h_invMatrix->GetXaxis()->SetNdivisions(20);
  h_invMatrix->GetYaxis()->SetNdivisions(20);
  
  h_invMatrix->GetXaxis()->SetTitleSize(0);
  h_invMatrix->GetXaxis()->SetLabelSize(0);
  h_invMatrix->GetYaxis()->SetTitleSize(0);
  h_invMatrix->GetYaxis()->SetLabelSize(0);

  h_invMatrix->Draw("colz");

  line0->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");

  TLatex *texf = new TLatex(-1.5,120,"true z");
  texf->SetTextColor(kBlack);
  texf->SetTextAlign(22);
  texf->SetTextSize(0.03);
  texf->SetTextAngle(90);
  texf->Draw();

  TLatex *texf1 = new TLatex(3.,-25.,"measured z");
  texf1->SetTextColor(kBlack);
  texf1->SetTextAlign(22);
  texf1->SetTextSize(0.03);
  texf1->Draw();

  TLatex *texf2 = new TLatex(3,1250,"10 #minus 20");
  texf2->SetTextColor(kBlack);
  texf2->SetTextAlign(22);
  texf2->SetTextSize(0.03);
  texf2->Draw();

  TLatex *texf3 = new TLatex(9.,1250,"20 #minus 30");
  texf3->SetTextColor(kBlack);
  texf3->SetTextAlign(22);
  texf3->SetTextSize(0.03);
  texf3->Draw();

  TLatex *texf4 = new TLatex(15.,1250,"30 #minus 40");
  texf4->SetTextColor(kBlack);
  texf4->SetTextAlign(22);
  texf4->SetTextSize(0.03);
  texf4->Draw();

  TLatex *texf41 = new TLatex(21.,1250,"40 #minus 50");
  texf41->SetTextColor(kBlack);
  texf41->SetTextAlign(22);
  texf41->SetTextSize(0.03);
  texf41->Draw();

  TLatex *texf42 = new TLatex(27.,1250,"50 #minus 60");
  texf42->SetTextColor(kBlack);
  texf42->SetTextAlign(22);
  texf42->SetTextSize(0.03);
  texf42->Draw();

  TLatex *texf5 = new TLatex(15.,1300,"measured jet p_{T} [GeV]");
  texf5->SetTextColor(kBlack);
  texf5->SetTextAlign(22);
  texf5->SetTextSize(0.03);
  texf5->Draw();
    
  TLatex *texf6 = new TLatex(35.,120.,"10 #minus 20");
  texf6->SetTextColor(kBlack);
  texf6->SetTextAlign(22);
  texf6->SetTextSize(0.03);
  texf6->SetTextAngle(270);
  texf6->Draw();

  TLatex *texf7 = new TLatex(35.,360.,"20 #minus 30");
  texf7->SetTextColor(kBlack);
  texf7->SetTextAlign(22);
  texf7->SetTextSize(0.03);
  texf7->SetTextAngle(270);
  texf7->Draw();

  TLatex *texf8 = new TLatex(35.,600.,"30 #minus 40");
  texf8->SetTextColor(kBlack);
  texf8->SetTextAlign(22);
  texf8->SetTextSize(0.03);
  texf8->SetTextAngle(270);
  texf8->Draw();

  TLatex *texf9 = new TLatex(35.,840.,"40 #minus 50");
  texf9->SetTextColor(kBlack);
  texf9->SetTextAlign(22);
  texf9->SetTextSize(0.03);
  texf9->SetTextAngle(270);
  texf9->Draw();

  TLatex *texf10 = new TLatex(35.,1080.,"50 #minus 60");
  texf10->SetTextColor(kBlack);
  texf10->SetTextAlign(22);
  texf10->SetTextSize(0.03);
  texf10->SetTextAngle(270);
  texf10->Draw();

  
  TLatex *texf11 = new TLatex(40.,600,"true jet p_{T} [GeV]");
  texf11->SetTextColor(kBlack);
  texf11->SetTextAlign(22);
  texf11->SetTextSize(0.03);
  texf11->SetTextAngle(270);
  texf11->Draw();

  TLatex *texf12 = new TLatex(-0.5,-0.5,"0");
  texf12->SetTextColor(kBlack);
  texf12->SetTextAlign(22);
  texf12->SetTextSize(0.03);
  //texf10->Draw();

  TLatex *texf13 = new TLatex(-0.5,5.,"1");
  texf13->SetTextColor(kBlack);
  texf13->SetTextAlign(22);
  texf13->SetTextSize(0.03);
  //texf11->Draw();

  TLatex *texf14 = new TLatex(5,-0.5,"1");
  texf14->SetTextColor(kBlack);
  texf14->SetTextAlign(22);
  texf14->SetTextSize(0.03);
  //texf12->Draw();

  TString cmsText = "CMS";
  TLatex *latex = new TLatex(1.2,14.5,cmsText);;
  latex->SetTextFont(61);
  latex->SetTextSize(0.04);
  latex->SetTextAlign(22);

  TString simText = "Simulation";
  TLatex *sim = new TLatex(1.9,13.7,simText);
  sim->SetTextFont(52);
  sim->SetTextSize(0.03);
  sim->SetTextAlign(22);
    
  TString extraText = "Preliminary";
  TLatex *extra = new TLatex(1.9,12.9,extraText);
  extra->SetTextFont(52);
  extra->SetTextSize(0.03);
  extra->SetTextAlign(22);

  //latex->Draw();
  //sim->Draw();
  //extra->Draw();

  string rapRange = "";
  //if(doMid) 
  rapRange = "0 < |y| < 1.6";
  //else rapRange = "1.6 < |y| < 2.4";
  
  TLatex *texf15 = new TLatex(12.,-2.,rapRange.c_str());
  texf15->SetTextColor(kBlack);
  texf15->SetTextAlign(22);
  texf15->SetTextSize(0.03);
  //texf15->Draw();

  string jpsi = "";
  if(doPrompt) jpsi = "prompt J/#psi";
  else jpsi = "nonprompt J/#psi";

  TLatex *texf16 = new TLatex(12.1,-1.,jpsi.c_str());
  texf16->SetTextColor(kBlack);
  texf16->SetTextAlign(22);
  texf16->SetTextSize(0.03);
  //texf16->Draw();
    
  //if(doPrompt && doMid) can->SaveAs("trMatrix_pr_mid_PAS.pdf");
  //if(doPrompt && !doMid) can->SaveAs("trMatrix_pr_fwd_PAS.pdf");
  //if(!doPrompt && doMid) can->SaveAs("trMatrix_nonpr_mid_PAS.pdf");
  //if(!doPrompt && !doMid) can->SaveAs("trMatrix_nonpr_fwd_PAS.pdf");
    
}


void plotMatrix_PAS(){
  
  //plot(true,true);
  plot(true,false);
  //plot(false,true);
  //plot(false,false);
  
}
