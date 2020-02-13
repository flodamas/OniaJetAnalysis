#include "inputParams.h"
void plot(bool doPrompt = true, bool doPbPb = false){

  Int_t colCon[] = {851,851,853,854,855,856,857,858,859,860,821,822,823,824,825,826,827,828,829,830,791,792,793,794,795,796,797,798,799,800};

  int iterFinal = nSIter;
  if (!doPbPb) iterFinal = nSIter_pp;

  string filename1 = "";
  string filename2 = "";
  //string filenameMeas = "";
  string outputfile = "";

  outputfile = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/plots/UnfoldinConvergence_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%s_statError%s.pdf",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal",(centShift==0)?"":(centShift==1)?"_centShiftSystUp":"_centShiftSystDown");

  TFile *file1=NULL;
  TFile *file2=NULL;

  TH2D *h2UnfResp1=NULL;
  TH2D *h2UnfResp2=NULL;

  Int_t min_PriorFolded=0;
  Int_t max_PriorFolded=0;

  TH1D *hZUnf_SI1 = NULL;
  TH1D *hZUnf_SI2 = NULL;
  TH1D *hZUnf_ratio = NULL;
  TCanvas * mycan1 = new TCanvas("mycan1","mycan1",800,800);
  //mycan1->Divide(2,1);
  
  TLegend *legend = new TLegend(0.4,0.7,0.5,0.9);//,"","brNDC");
  legend->SetHeader("");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.02);

  TLegend *legend1 = new TLegend(0.6,0.7,0.7,0.9);//,"","brNDC");
  legend1->SetHeader("");
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->SetTextFont(42);
  legend1->SetTextSize(0.02);

  TLegend *legend2 = new TLegend(0.8,0.7,0.9,0.9);//,"","brNDC");
  legend2->SetHeader("");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->SetTextFont(42);
  legend2->SetTextSize(0.02);

  for(int i=1; i<iterFinal; i++){
  filename1 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%sDiag%s.root",i,doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal",(centShift==0)?"":(centShift==1)?"_centShiftSystUp":"_centShiftSystDown");
  filename2 = Form("/Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/unfOutput/step%i/UnfoldedDistributions_%s_%s_8iter_%dz%dptBins%dz%dptMeasBin%sDiag%s.root",(i+1),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco, "_nominal",(centShift==0)?"":(centShift==1)?"_centShiftSystUp":"_centShiftSystDown");
  file1 = new TFile(filename1.c_str());
  file2 = new TFile(filename2.c_str());

  h2UnfResp1 = (TH2D*)file1->Get("hReco_Iter1;1");
  h2UnfResp2 = (TH2D*)file2->Get("hReco_Iter1;1");

  min_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midLowerPt+0.00001);
  max_PriorFolded = h2UnfResp1->GetYaxis()->FindBin(midUpperPt-0.00001);
  
  hZUnf_SI1 = (TH1D*)h2UnfResp1->ProjectionX(Form("hZUnf_SI%d",i),min_PriorFolded,max_PriorFolded);
  hZUnf_SI2 = (TH1D*)h2UnfResp2->ProjectionX(Form("hZUnf_SI%d",i+1),min_PriorFolded,max_PriorFolded);
  
  hZUnf_ratio = (TH1D*)hZUnf_SI2->Clone(Form("hZUnf_ratio%d",i));
  hZUnf_ratio->Add(hZUnf_SI1,-1);
  hZUnf_ratio->Divide(hZUnf_SI1);

  hZUnf_ratio->SetStats(0);
  hZUnf_ratio->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"pp"));

  hZUnf_ratio->SetMaximum(hZUnf_ratio->GetMaximum()*1.5);

  hZUnf_ratio->GetXaxis()->SetTitle("z");
  hZUnf_ratio->GetYaxis()->SetTitleOffset(1.3);
  hZUnf_ratio->GetYaxis()->SetTitle("(SI_{i+1}-SI_{i})/SI_{i}");
  if (doPbPb)
    hZUnf_ratio->SetLineColor(colCon[i]);
  else
    hZUnf_ratio->SetLineColor(colCon[10*i-9]);
  //hZUnf_ratio->SetMarkerColor(col[i]);
  //hZUnf_ratio->SetMarkerStyle(markerStyle[i]);
  //hZUnf_ratio->SetMarkerSize(markerSize[i]);
  hZUnf_ratio->SetLineWidth(1);  
  if (i==1)
    hZUnf_ratio->Draw("hist");
  else
    hZUnf_ratio->Draw("hist same");
  if (i<11)
    legend->AddEntry(hZUnf_ratio, Form("SI_{%d}-SI_{%d}",i+1,i),"l");
  else if (i<21)
    legend1->AddEntry(hZUnf_ratio, Form("SI_{%d}-SI_{%d}",i+1,i),"l");
  else
    legend2->AddEntry(hZUnf_ratio, Form("SI_{%d}-SI_{%d}",i+1,i),"l");
}

  float xCoord = unfStart;
    
  TLine *line0 = new TLine(min_z,0,max_z,0);
  line0->SetLineColor(kred);
  line0->SetLineStyle(2);
  line0->SetLineWidth(2);
  line0->Draw("same");
  legend->Draw("same");
  legend1->Draw("same");
  legend2->Draw("same");
  mycan1->SaveAs(outputfile.c_str());  
  
}
void PlotConvergence_DataUnfolded_afterDiag(){
  plot(true,true);
  if (centShift==0)
  plot(true,false);
}
