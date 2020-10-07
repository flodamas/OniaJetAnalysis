#include "inputParams.h"

void plot(bool doPrompt = true, bool doPbPb = true){
  if (!setSystTag(doPbPb)) return;

  gStyle->SetOptStat(0);

  string filename = "";
  string outputfile = "";
  filename = Form("%s/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBin%s.root",unfPath.c_str(), 1,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", doPbPb?nIter:nIter_pp, nBinZ_reco, nBinJet_reco, nBinZ_reco, nBinJet_reco, systTag.c_str());

  outputfile = Form("%s/dataUnf/unfOutput/finalResults/trMatrix_%s_%s%s.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",noSmearing?"_noSmearing":"");


  cout <<"filename = "<<filename <<endl;    
  TFile *file = new TFile(filename.c_str());
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get(Form("invtrmat%d",doPbPb?nIter:nIter_pp));
  
  TH2D *h2_UnfData  = (TH2D*)file->Get(Form("hReco_Iter%d",doPbPb?nIter:nIter_pp));
  TH2D *h2_MeasData  = (TH2D*)file->Get("fh2MeasData");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  
  //TH2D *h_covMatrix = new TH2D(*covMatrix);
  //h_invMatrix->Rebin2D(1,nBinJet_gen*nBinZ_gen/(nBinJet_reco*nBinZ_reco));
  //h_invMatrix->GetZaxis()->SetRangeUser(0,1);

  TLine *line0 = new TLine(6,0,6,33.2);
  line0->SetLineColor(kBlack);
  line0->SetLineStyle(1);
  line0->SetLineWidth(2);

  TLine *line01 = new TLine(6,35.4,6,36);
  line01->SetLineColor(kBlack);
  line01->SetLineStyle(1);
  line01->SetLineWidth(2);


  TLine *line1 = new TLine(12,0,12,33.2);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(1);
  line1->SetLineWidth(2);

  TLine *line11 = new TLine(12,35.4,12,36);
  line11->SetLineColor(kBlack);
  line11->SetLineStyle(1);
  line11->SetLineWidth(2);

  TLine *line2 = new TLine(18,2.7,18,36);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(1);
  line2->SetLineWidth(2);

  TLine *line3 = new TLine(24,2.7,24,36);
  line3->SetLineColor(kBlack);
  line3->SetLineStyle(1);
  line3->SetLineWidth(2);

  TLine *line4 = new TLine(30,2.7,30,36);
  line4->SetLineColor(kBlack);
  line4->SetLineStyle(1);
  line4->SetLineWidth(2);

  TLine *line5 = new TLine(0,6,36,6);
  line5->SetLineColor(kBlack);
  line5->SetLineStyle(1);
  line5->SetLineWidth(2);

  TLine *line6 = new TLine(0,12,36,12);
  line6->SetLineColor(kBlack);
  line6->SetLineStyle(1);
  line6->SetLineWidth(2);

  TLine *line7 = new TLine(0,18,36,18);
  line7->SetLineColor(kBlack);
  line7->SetLineStyle(1);
  line7->SetLineWidth(2);

  TLine *line8 = new TLine(0,24,36,24);
  line8->SetLineColor(kBlack);
  line8->SetLineStyle(1);
  line8->SetLineWidth(2);

  TLine *line9 = new TLine(0,30,36,30);
  line9->SetLineColor(kBlack);
  line9->SetLineStyle(1);
  line9->SetLineWidth(2);

  TCanvas *can = new TCanvas("can","can",1000,900);

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
  line01->Draw("same");
  line1->Draw("same");
  line11->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");
  line9->Draw("same");

  TLatex *texf = new TLatex(-1.7,3,"True z");
  texf->SetTextFont(42);
  texf->SetTextColor(kBlack);
  texf->SetTextAlign(22);
  texf->SetTextSize(0.035);
  texf->SetTextAngle(90);
  texf->Draw();

  TLatex *texf0 = new TLatex(3.,-1.9,"Measured z");
  texf0->SetTextFont(42);
  texf0->SetTextColor(kBlack);
  texf0->SetTextAlign(22);
  texf0->SetTextSize(0.035);
  texf0->Draw();

  TLatex *texf1 = new TLatex(3.,37,"6.5 #minus 10");
  texf1->SetTextFont(42);
  texf1->SetTextColor(kBlack);
  texf1->SetTextAlign(22);
  texf1->SetTextSize(0.03);
  texf1->Draw();

  TLatex *texf2 = new TLatex(9,37,"10 #minus 20");
  texf2->SetTextFont(42);
  texf2->SetTextColor(kBlack);
  texf2->SetTextAlign(22);
  texf2->SetTextSize(0.03);
  texf2->Draw();

  TLatex *texf3 = new TLatex(15.,37,"20 #minus 30");
  texf3->SetTextFont(42);
  texf3->SetTextColor(kBlack);
  texf3->SetTextAlign(22);
  texf3->SetTextSize(0.03);
  texf3->Draw();

  TLatex *texf4 = new TLatex(21.,37,"30 #minus 40");
  texf4->SetTextFont(42);
  texf4->SetTextColor(kBlack);
  texf4->SetTextAlign(22);
  texf4->SetTextSize(0.03);
  texf4->Draw();

  TLatex *texf41 = new TLatex(27.,37,"40 #minus 50");
  texf41->SetTextFont(42);
  texf41->SetTextColor(kBlack);
  texf41->SetTextAlign(22);
  texf41->SetTextSize(0.03);
  texf41->Draw();

  TLatex *texf42 = new TLatex(33.,37,"50 #minus 60");
  texf42->SetTextFont(42);
  texf42->SetTextColor(kBlack);
  texf42->SetTextAlign(22);
  texf42->SetTextSize(0.03);
  texf42->Draw();

  TLatex *texf5 = new TLatex(18.,38.6,"Measured jet p_{T} [GeV]");
  texf5->SetTextFont(42);
  texf5->SetTextColor(kBlack);
  texf5->SetTextAlign(22);
  texf5->SetTextSize(0.04);
  texf5->Draw();
    
  TLatex *texf6 = new TLatex(43.,3.,"6.5 #minus 10");
  texf6->SetTextFont(42);
  texf6->SetTextColor(kBlack);
  texf6->SetTextAlign(22);
  texf6->SetTextSize(0.03);
  texf6->SetTextAngle(270);
  texf6->Draw();

  TLatex *texf7 = new TLatex(43.,9.,"10 #minus 20");
  texf7->SetTextFont(42);
  texf7->SetTextColor(kBlack);
  texf7->SetTextAlign(22);
  texf7->SetTextSize(0.03);
  texf7->SetTextAngle(270);
  texf7->Draw();

  TLatex *texf8 = new TLatex(43.,15.,"20 #minus 30");
  texf8->SetTextFont(42);
  texf8->SetTextColor(kBlack);
  texf8->SetTextAlign(22);
  texf8->SetTextSize(0.03);
  texf8->SetTextAngle(270);
  texf8->Draw();

  TLatex *texf9 = new TLatex(43.,21.,"30 #minus 40");
  texf9->SetTextFont(42);
  texf9->SetTextColor(kBlack);
  texf9->SetTextAlign(22);
  texf9->SetTextSize(0.03);
  texf9->SetTextAngle(270);
  texf9->Draw();

  TLatex *texf10 = new TLatex(43.,27.,"40 #minus 50");
  texf10->SetTextFont(42);
  texf10->SetTextColor(kBlack);
  texf10->SetTextAlign(22);
  texf10->SetTextSize(0.03);
  texf10->SetTextAngle(270);
  texf10->Draw();

  TLatex *texf11 = new TLatex(43.,33.,"50 #minus 60");
  texf11->SetTextFont(42);
  texf11->SetTextColor(kBlack);
  texf11->SetTextAlign(22);
  texf11->SetTextSize(0.03);
  texf11->SetTextAngle(270);
  texf11->Draw();

  TLatex *texf12 = new TLatex(45.,17,"True jet p_{T} [GeV]");
  texf12->SetTextFont(42);
  texf12->SetTextColor(kBlack);
  texf12->SetTextAlign(22);
  texf12->SetTextSize(0.04);
  texf12->SetTextAngle(270);
  texf12->Draw();

  TLatex *texf13 = new TLatex(-0.6,-0.6,"0.064");
  texf13->SetTextFont(42);
  texf13->SetTextColor(kBlack);
  texf13->SetTextAlign(22);
  texf13->SetTextSize(0.03);
  texf13->Draw();

  TLatex *texf14 = new TLatex(-0.6,6.1,"1");
  texf14->SetTextFont(42);
  texf14->SetTextColor(kBlack);
  texf14->SetTextAlign(22);
  texf14->SetTextSize(0.03);
  texf14->Draw();

  TLatex *texf15 = new TLatex(6.1,-0.6,"1");
  texf15->SetTextFont(42);
  texf15->SetTextColor(kBlack);
  texf15->SetTextAlign(22);
  texf15->SetTextSize(0.03);
  texf15->Draw();

  //18.,39
  TString cmsText = "CMS";
  //TLatex *latex = new TLatex(44,39.5,cmsText);
  TLatex *latex = new TLatex(4,34.5,cmsText);
  latex->SetTextFont(61);
  latex->SetTextSize(0.055);
  latex->SetTextAlign(22);

  TString simText = "Simulation";
  //TLatex *sim = new TLatex(44,37.5,simText);
  TLatex *sim = new TLatex(12,34.2,simText);
  sim->SetTextFont(52);
  sim->SetTextSize(0.04);
  sim->SetTextAlign(22);
    
  TString extraText = "Preliminary";
  TLatex *extra = new TLatex(3,31,extraText);
  extra->SetTextFont(52);
  extra->SetTextSize(0.024);
  extra->SetTextAlign(22);

  latex->Draw();
  sim->Draw();
  //extra->Draw();


  string coll = "";
  if(doPbPb) coll = "PbPb";
  else coll = "pp";
  //TLatex *texf16 = new TLatex(1.5,35,coll.c_str());
  //if (doPbPb) texf16 = new TLatex(2.6,35,coll.c_str());
  TLatex *texf16 = new TLatex(16,1.6,coll.c_str());
  if (doPbPb) texf16 = new TLatex(15,1.9,coll.c_str());
  texf16->SetTextColor(kBlack);
  texf16->SetTextFont(42);
  texf16->SetTextAlign(22);
  texf16->SetTextSize(0.04);
  texf16->Draw();

  string jpsi = "";
  if(doPrompt) jpsi = "prompt J/#psi";
  else jpsi = "nonprompt J/#psi";

  //TLatex *texf17 = new TLatex(5.,32.75,jpsi.c_str());
  TLatex *texf17 = new TLatex(23,1.6,jpsi.c_str());
  texf17->SetTextColor(kBlack);
  texf17->SetTextFont(42);
  texf17->SetTextAlign(22);
  texf17->SetTextSize(0.04);
  texf17->Draw();

  string rapRange = "|#eta_{jet}| < 2";
  TLatex *texf18 = new TLatex(32,1.4,rapRange.c_str());
  texf18->SetTextFont(42);
  texf18->SetTextColor(kBlack);
  texf18->SetTextAlign(22);
  texf18->SetTextSize(0.04);
  texf18->Draw();

  can->SaveAs(outputfile.c_str());
    
}


void plotMatrix_PAS(){
  
  plot(true,true);
  plot(true,false);
  //plot(false,true);
  //plot(false,false);
  
}
