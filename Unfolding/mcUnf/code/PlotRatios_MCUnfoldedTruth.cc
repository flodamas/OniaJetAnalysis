#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void plot(bool doPrompt = false, bool doPbPb = true){
  if (!setCaseTag()) return;
  gSystem->mkdir(Form("%s/mcUnf/plots/",unfPath.c_str()));
  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";

  int iterFinal = nSIter;
  if (!doPbPb) iterFinal = nSIter_pp;
  if (iterFinal<4) iterFinal = 4;

  filename1 = Form("%s/mcUnf/unfOutput/step1/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  filename2 = Form("%s/mcUnf/unfOutput/step2/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  filename3 = Form("%s/mcUnf/unfOutput/step%d/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), iterFinal-1, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  filename4 = Form("%s/mcUnf/unfOutput/step%d/UnfoldedDistributions_%s_%s_%diter_%dz%dptBins%dz%dptMeasBins%s.root", unfPath.c_str(), iterFinal, doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nIter, nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str());
  
  TFile *file1 = new TFile(filename1.c_str());
  TFile *file2 = new TFile(filename2.c_str());
  TFile *file3 = new TFile(filename3.c_str());
  TFile *file4 = new TFile(filename4.c_str());

  TH1D *hZMeas;
  hZMeas=(TH1D*)file3->Get(Form("hMSme_%d;1",midLowerId));

  TH1D *hZTrue = (TH1D*)file3->Get(Form("hMTru_%d;1",midLowerId)); TH1D *hZTrueTemp = NULL;
  TH1D *hZUnf_SI1 = (TH1D*)file1->Get(Form("hMUnf_%d_Iter%d;1",midLowerId,nIter)); TH1D *hZUnf_SI1_temp = NULL;
  TH1D *hZUnf_SI2 = (TH1D*)file2->Get(Form("hMUnf_%d_Iter%d;1",midLowerId,nIter)); TH1D *hZUnf_SI2_temp = NULL;
  TH1D *hZUnf_SI3 = (TH1D*)file3->Get(Form("hMUnf_%d_Iter%d;1",midLowerId,nIter)); TH1D *hZUnf_SI3_temp = NULL;
  TH1D *hZUnf_SI4 = (TH1D*)file4->Get(Form("hMUnf_%d_Iter%d;1",midLowerId,nIter)); TH1D *hZUnf_SI4_temp = NULL;  

  cout <<"midUpperId = "<<midUpperId<<", midLowerId = "<<midLowerId<<endl;

  for (int i = 1; i <= (midUpperId-midLowerId); i++) {
    hZTrueTemp = (TH1D*)file3->Get(Form("hMTru_%d;1",midLowerId+i));
    hZTrue->Add(hZTrueTemp);

    hZUnf_SI1_temp = (TH1D*)file1->Get(Form("hMUnf_%d_Iter%d;1",midLowerId+i,nIter));
    hZUnf_SI1->Add(hZUnf_SI1_temp);

    hZUnf_SI2_temp = (TH1D*)file2->Get(Form("hMUnf_%d_Iter%d;1",midLowerId+i,nIter));
    hZUnf_SI2->Add(hZUnf_SI2_temp);

    hZUnf_SI3_temp = (TH1D*)file3->Get(Form("hMUnf_%d_Iter%d;1",midLowerId+i,nIter));
    hZUnf_SI3->Add(hZUnf_SI3_temp);

    hZUnf_SI4_temp = (TH1D*)file4->Get(Form("hMUnf_%d_Iter%d;1",midLowerId+i,nIter));
    hZUnf_SI4->Add(hZUnf_SI4_temp);
  }

  hZTrue->Rebin(nBinZ_gen/nBinZ_reco);
  hZUnf_SI1->Rebin(nBinZ_gen/nBinZ_reco);
  hZUnf_SI2->Rebin(nBinZ_gen/nBinZ_reco);
  hZUnf_SI3->Rebin(nBinZ_gen/nBinZ_reco);
  hZUnf_SI4->Rebin(nBinZ_gen/nBinZ_reco);

  TH1D *hZTrue_ratioMeas;
  TH1D *hZUnf_SI1_ratioMeas;
  TH1D *hZUnf_SI2_ratioMeas;
  TH1D *hZUnf_SI3_ratioMeas;
  TH1D *hZUnf_SI4_ratioMeas;

  hZTrue_ratioMeas = (TH1D*)hZTrue->Clone();
  hZTrue_ratioMeas->Divide(hZMeas);  
  
  hZUnf_SI1_ratioMeas = (TH1D*)hZUnf_SI1->Clone();
  hZUnf_SI1_ratioMeas->Divide(hZMeas);

  hZUnf_SI2_ratioMeas = (TH1D*)hZUnf_SI2->Clone();
  hZUnf_SI2_ratioMeas->Divide(hZMeas);

  hZUnf_SI3_ratioMeas = (TH1D*)hZUnf_SI3->Clone();
  hZUnf_SI3_ratioMeas->Divide(hZMeas);

  hZUnf_SI4_ratioMeas = (TH1D*)hZUnf_SI4->Clone();
  hZUnf_SI4_ratioMeas->Divide(hZMeas);

  TH1D *hZMeas_ratioTruth;
  TH1D *hZUnf_SI1_ratioTruth;
  TH1D *hZUnf_SI2_ratioTruth;
  TH1D *hZUnf_SI3_ratioTruth;
  TH1D *hZUnf_SI4_ratioTruth;

  hZMeas_ratioTruth = (TH1D*)hZMeas->Clone();
  hZMeas_ratioTruth->Divide(hZTrue);  
  
  hZUnf_SI1_ratioTruth = (TH1D*)hZUnf_SI1->Clone();
  hZUnf_SI1_ratioTruth->Divide(hZTrue);

  hZUnf_SI2_ratioTruth = (TH1D*)hZUnf_SI2->Clone();
  hZUnf_SI2_ratioTruth->Divide(hZTrue);

  hZUnf_SI3_ratioTruth = (TH1D*)hZUnf_SI3->Clone();
  hZUnf_SI3_ratioTruth->Divide(hZTrue);

  hZUnf_SI4_ratioTruth = (TH1D*)hZUnf_SI4->Clone();
  hZUnf_SI4_ratioTruth->Divide(hZTrue);
  
  TCanvas * mycan1 = new TCanvas("mycan1","mycan1",900,500);
  mycan1->Divide(2,1);
  
  TLegend *legend = new TLegend(0.8,0.15,.95,0.9,"","brNDC");
  legend->SetHeader("");

  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);

  mycan1->cd(1);
  //gPad->SetLogy();

  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.2);
  
  hZMeas->SetStats(0);
  hZMeas->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"PP"));
  
  hZMeas->GetXaxis()->SetTitle("z");
  hZMeas->GetYaxis()->SetTitleOffset(1.6);
  hZMeas->GetYaxis()->SetTitle(Form("dN/dz (per %f)",z_reco_binWidth));

  hZMeas->GetYaxis()->SetRangeUser(0,hZMeas->GetMaximum()*1.5);
  
  hZMeas->SetLineColor(col[0]);
  hZMeas->SetMarkerColor(col[0]);
  hZMeas->SetMarkerStyle(markerStyle[0]);
  hZMeas->SetMarkerSize(markerSize[0]);
  hZMeas->SetLineWidth(lineWidth[0]);

  hZTrue->SetLineColor(col[1]);
  hZTrue->SetMarkerColor(col[1]);
  hZTrue->SetMarkerStyle(markerStyle[1]);
  hZTrue->SetMarkerSize(markerSize[1]);
  hZTrue->SetLineWidth(lineWidth[1]);

  hZUnf_SI1->SetLineColor(col[2]);
  hZUnf_SI2->SetLineColor(col[3]);
  hZUnf_SI3->SetLineColor(col[4]);
  hZUnf_SI4->SetLineColor(col[5]);
    
  hZUnf_SI1->SetMarkerColor(col[2]);
  hZUnf_SI2->SetMarkerColor(col[3]);
  hZUnf_SI3->SetMarkerColor(col[4]);
  hZUnf_SI4->SetMarkerColor(col[5]);

  hZUnf_SI1->SetMarkerStyle(markerStyle[2]);
  hZUnf_SI2->SetMarkerStyle(markerStyle[3]);
  hZUnf_SI3->SetMarkerStyle(markerStyle[4]);
  hZUnf_SI4->SetMarkerStyle(markerStyle[5]);
  
  hZUnf_SI1->SetMarkerSize(markerSize[2]);
  hZUnf_SI2->SetMarkerSize(markerSize[3]);
  hZUnf_SI3->SetMarkerSize(markerSize[4]);
  hZUnf_SI4->SetMarkerSize(markerSize[5]);
  
  hZUnf_SI1->SetLineWidth(lineWidth[2]);
  hZUnf_SI2->SetLineWidth(lineWidth[3]);
  hZUnf_SI3->SetLineWidth(lineWidth[4]);
  hZUnf_SI4->SetLineWidth(lineWidth[5]);
  
  hZMeas->Draw("EP");
  hZTrue->Draw("EPsame");
  hZUnf_SI1->Draw("EPsame");
  hZUnf_SI2->Draw("EPsame");
  hZUnf_SI4->Draw("EPsame");
  hZUnf_SI3->Draw("EPsame");
  
  legend->AddEntry(hZMeas, "meas","ep");
  legend->AddEntry(hZTrue, "truth","ep");
  legend->AddEntry(hZUnf_SI1, "unf SI#1","ep");
  legend->AddEntry(hZUnf_SI2, "unf SI#2","ep");
  legend->AddEntry(hZUnf_SI3, Form("unf SI#%d",iterFinal-1),"ep");
  legend->AddEntry(hZUnf_SI4, Form("unf SI#%d",iterFinal),"ep");
  
  legend->Draw("same");

  mycan1->Update();

  float xCoord = unfStart;
    
  TLine *line0 = new TLine(xCoord,0,xCoord,gPad->GetUymax());
  line0->SetLineColor(kRed);
  line0->SetLineStyle(2);
  line0->SetLineWidth(1);
  //line0->Draw("same");
  
  mycan1->cd(2);
  
  //gPad->SetRightMargin(0.2);

  hZUnf_SI1_ratioMeas->SetStats(0);
  hZUnf_SI1_ratioMeas->GetYaxis()->SetTitle("ratio to measured");
  hZUnf_SI1_ratioMeas->GetYaxis()->SetTitleOffset(1.6);
  hZUnf_SI1_ratioMeas->GetYaxis()->SetRangeUser(0.,hZUnf_SI1_ratioMeas->GetMaximum()*1.3);
  hZUnf_SI1_ratioMeas->GetXaxis()->SetTitle("z");

  hZTrue_ratioMeas->SetLineColor(col[1]);
  hZTrue_ratioMeas->SetMarkerColor(col[1]);
  hZTrue_ratioMeas->SetMarkerStyle(markerStyle[1]);
  hZTrue_ratioMeas->SetMarkerSize(markerSize[1]);
  hZTrue_ratioMeas->SetLineStyle(lineStyle[1]);
  hZTrue_ratioMeas->SetLineWidth(lineWidth[1]*2);
  
  hZUnf_SI1_ratioMeas->SetLineColor(col[2]);
  hZUnf_SI2_ratioMeas->SetLineColor(col[3]);
  hZUnf_SI3_ratioMeas->SetLineColor(col[4]);
  hZUnf_SI4_ratioMeas->SetLineColor(col[5]);
  
  hZUnf_SI1_ratioMeas->SetLineStyle(lineStyle[2]);
  hZUnf_SI2_ratioMeas->SetLineStyle(lineStyle[3]);
  hZUnf_SI3_ratioMeas->SetLineStyle(lineStyle[4]);
  hZUnf_SI4_ratioMeas->SetLineStyle(lineStyle[5]);

  hZUnf_SI1_ratioMeas->SetLineWidth(lineWidth[2]*2);
  hZUnf_SI2_ratioMeas->SetLineWidth(lineWidth[3]*2);
  hZUnf_SI3_ratioMeas->SetLineWidth(lineWidth[4]*2);
  hZUnf_SI4_ratioMeas->SetLineWidth(lineWidth[5]*2);

  hZUnf_SI1_ratioMeas->SetMarkerColor(col[2]);
  hZUnf_SI2_ratioMeas->SetMarkerColor(col[3]);
  hZUnf_SI3_ratioMeas->SetMarkerColor(col[4]);
  hZUnf_SI4_ratioMeas->SetMarkerColor(col[5]);
  
  hZUnf_SI1_ratioMeas->Draw("HIST");
  hZTrue_ratioMeas->Draw("HISTsame");
  //hZUnf_SI1_ratioMeas->Draw("HISTsame");
  hZUnf_SI2_ratioMeas->Draw("HISTsame");
  hZUnf_SI4_ratioMeas->Draw("HISTsame");
  hZUnf_SI3_ratioMeas->Draw("HISTsame");
  
  mycan1->Update();

  float xCoord2 = unfStart;
    
  TLine *line1 = new TLine(xCoord2,0,xCoord2,gPad->GetUymax());
  line1->SetLineColor(kRed);
  line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw("same");

  mycan1->Update();
  
  mycan1->SaveAs(Form("%s/mcUnf/plots/unf_mc_%s_%s_jetR%d_ratioMeasured_%dz%dptBins%dz%dptMeasBins%s_%dIter%dSIter.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter,iterFinal));
  mycan1->SaveAs(Form("%s/mcUnf/plots/unf_mc_%s_%s_jetR%d_ratioMeasured_%dz%dptBins%dz%dptMeasBins%s_%dIter%dSIter.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter,iterFinal));

  hZUnf_SI1_ratioTruth->SetStats(0);
  hZUnf_SI1_ratioTruth->GetYaxis()->SetTitle("ratio to truth");
  hZUnf_SI1_ratioTruth->GetYaxis()->SetTitleOffset(1.6);
  hZUnf_SI1_ratioTruth->GetYaxis()->SetRangeUser(0.,2.4);
  hZUnf_SI1_ratioTruth->GetXaxis()->SetTitle("z");

  hZMeas_ratioTruth->SetLineColor(col[0]);
  hZMeas_ratioTruth->SetMarkerColor(col[0]);
  hZMeas_ratioTruth->SetMarkerStyle(markerStyle[0]);
  hZMeas_ratioTruth->SetMarkerSize(markerSize[0]);
  hZMeas_ratioTruth->SetLineStyle(lineStyle[0]);
  hZMeas_ratioTruth->SetLineWidth(lineWidth[0]*2);
  
  hZUnf_SI1_ratioTruth->SetLineColor(col[2]);
  hZUnf_SI2_ratioTruth->SetLineColor(col[3]);
  hZUnf_SI3_ratioTruth->SetLineColor(col[4]);
  hZUnf_SI4_ratioTruth->SetLineColor(col[5]);
  
  hZUnf_SI1_ratioTruth->SetLineStyle(lineStyle[2]);
  hZUnf_SI2_ratioTruth->SetLineStyle(lineStyle[3]);
  hZUnf_SI3_ratioTruth->SetLineStyle(lineStyle[4]);
  hZUnf_SI4_ratioTruth->SetLineStyle(lineStyle[5]);

  hZUnf_SI1_ratioTruth->SetLineWidth(lineWidth[2]*2);
  hZUnf_SI2_ratioTruth->SetLineWidth(lineWidth[3]*2);
  hZUnf_SI3_ratioTruth->SetLineWidth(lineWidth[4]*2);
  hZUnf_SI4_ratioTruth->SetLineWidth(lineWidth[5]*2);

  hZUnf_SI1_ratioTruth->SetMarkerColor(col[2]);
  hZUnf_SI2_ratioTruth->SetMarkerColor(col[3]);
  hZUnf_SI3_ratioTruth->SetMarkerColor(col[4]);
  hZUnf_SI4_ratioTruth->SetMarkerColor(col[5]);
  
  hZUnf_SI1_ratioTruth->Draw("HIST");
  hZMeas_ratioTruth->Draw("HISTsame");
  //hZUnf_SI1_ratioTruth->Draw("HISTsame");
  hZUnf_SI2_ratioTruth->Draw("HISTsame");
  hZUnf_SI4_ratioTruth->Draw("HISTsame");
  hZUnf_SI3_ratioTruth->Draw("HISTsame");
  
  mycan1->Update();
    
  line1->Draw("same");
  TLine *line2 = new TLine(min_z,1,gPad->GetUxmax(),1);
  line2->SetLineColor(kblue);
  line2->SetLineStyle(1);
  line2->SetLineWidth(1);
  line2->Draw("same");

  mycan1->Update();
  
  mycan1->SaveAs(Form("%s/mcUnf/plots/unf_mc_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter%dSIter.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter,iterFinal));
  mycan1->SaveAs(Form("%s/mcUnf/plots/unf_mc_%s_%s_jetR%d_ratioTruth_%dz%dptBins%dz%dptMeasBins%s_%dIter%dSIter.png",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, caseTag.c_str(),nIter,iterFinal));
}

void PlotRatios_MCUnfoldedTruth(){
  //plot(bool doPrompt, bool doPbPb)
  plot(true,true);
  if (centShift == 0 && nSIter_pp>0 && !doCent && !doPeri)
    plot(true,false);

  plot(false,true);
  if (centShift == 0 && nSIter_pp>0 && !doCent && !doPeri)
    plot(false,false);
}
