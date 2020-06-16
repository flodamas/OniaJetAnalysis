#include "inputParams.h"

void SetPullStyle(TH1D* pullHist);

void checkUnfoldingEfficiency(bool doPrompt = true, bool doPbPb = true, bool checkTruth = true) {
  gStyle->SetOptStat(false);
  
  string filename = Form("~/JpsiInJetsPbPb/Fitter/TreesForUnfolding/tree_%s_%s_NoBkg%s_AccEff_JEC.root",doPrompt?"MCJPSIPR":"MCJPSINOPR",doPbPb?"PbPb":"PP",mc2015?"":Form("_jetR%d",(int)(jetR*10)));
    
  string outputfile = Form("%s/mcUnf/plots/unf_Efficiency%s_%s_%s_jetR%d_%dz%dptBins%dz%dptMeasBins.pdf",unfPath.c_str(), checkTruth?"Truth":"Meas",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt", (int) (jetR*10), nBinZ_gen, nBinJet_gen, nBinZ_reco, nBinJet_reco);
  
  TFile *file = new TFile(filename.c_str());
  TTree *t_unf = (TTree*)file->Get("treeForUnfolding");
  
  TH1D* z_Truth = new TH1D("z_Truth","",nBinZ_reco,min_z,max_z);
  TH1D* z_UnfTruth = new TH1D("z_UnfTruth","",nBinZ_reco,min_z,max_z);
  TH1D* z_UnfTruth_up = new TH1D("z_UnfTruth_up","",nBinZ_reco,min_z,max_z);
  TH1D* z_UnfTruth_down = new TH1D("z_UnfTruth_down","",nBinZ_reco,min_z,max_z);
  TH1D* z_UnfTruth_rmProbUnderFlow = new TH1D("z_UnfTruth_rmProbUnderFlow","",nBinZ_reco,min_z,max_z);
  if (checkTruth) {
    t_unf->Draw("gen_z>>z_Truth",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && gen_z<1 && jt_ref_pt > %f && jt_ref_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,max_jt_eta));
    t_unf->Draw("gen_z>>z_UnfTruth",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && gen_z<1 && jt_ref_pt > %f && jt_ref_pt < %f && jt_pt > %f && jt_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,min_jetpt,max_jetpt,max_jt_eta));
    t_unf->Draw("gen_z>>z_UnfTruth_up",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && gen_z<1 && jt_ref_pt > %f && jt_ref_pt < %f && jt_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,max_jetpt,max_jt_eta));  
    t_unf->Draw("gen_z>>z_UnfTruth_down",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && gen_z<1 && jt_ref_pt > %f && jt_ref_pt < %f && jt_pt > %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,min_jetpt,max_jt_eta));
    t_unf->Draw("gen_z>>z_UnfTruth_rmProbUnderFlow",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && gen_z<1 && jt_ref_pt > %f && jt_ref_pt < %f && jt_pt > %f && jt_pt < %f && fabs(jt_eta)<%f && (!(z<0.688 && jt_pt<10)) )",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,min_jetpt,max_jetpt,max_jt_eta));
  }
  else {
    t_unf->Draw("z>>z_Truth",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && z<1 && jt_pt > %f && jt_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt, max_jt_eta));
    t_unf->Draw("z>>z_UnfTruth",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && z<1 && jt_pt > %f && jt_pt < %f && jt_ref_pt > %f && jt_ref_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,min_jetpt,max_jetpt,max_jt_eta));
    t_unf->Draw("z>>z_UnfTruth_up",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && z<1 && jt_pt > %f && jt_pt < %f && jt_ref_pt < %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,max_jetpt,max_jt_eta));  
    t_unf->Draw("z>>z_UnfTruth_down",Form("corr_AccEff*corr_ptw*(jp_pt>%f && fabs(jp_eta)<%f && z<1 && jt_pt > %f && jt_pt < %f && jt_ref_pt > %f && fabs(jt_eta)<%f)",min_jp_pt, max_jp_eta, midLowerPt, midUpperPt,min_jetpt,max_jt_eta));
  }
  z_Truth->SetMarkerColor(kred);
  z_Truth->SetMarkerStyle(kFullSquare);
  z_Truth->SetMarkerSize(1.5);
  z_Truth->SetLineColor(kred);
    
  z_UnfTruth->SetMarkerColor(kviolet);
  z_UnfTruth->SetMarkerStyle(kFullCircle);
  z_UnfTruth->SetMarkerSize(1);
  z_UnfTruth->SetLineColor(kviolet);

  z_UnfTruth_up->SetMarkerColor(kyellow);
  z_UnfTruth_up->SetMarkerStyle(kOpenCircle);
  z_UnfTruth_up->SetMarkerSize(1);
  z_UnfTruth_up->SetLineColor(kyellow);
  
  z_UnfTruth_down->SetMarkerColor(kgreen);
  z_UnfTruth_down->SetMarkerStyle(kOpenCircle);
  z_UnfTruth_down->SetMarkerSize(1);
  z_UnfTruth_down->SetLineColor(kgreen);

  z_UnfTruth_rmProbUnderFlow->SetMarkerColor(kGray);
  z_UnfTruth_rmProbUnderFlow->SetMarkerStyle(kOpenCircle);
  z_UnfTruth_rmProbUnderFlow->SetMarkerSize(1);
  z_UnfTruth_rmProbUnderFlow->SetLineColor(kGray);

  TCanvas* can = new TCanvas("can","",900,1000);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.01);
  pad1->Draw();
  pad1->cd();
  z_Truth->Draw("e1");
  z_UnfTruth->Draw("same e1");
  z_UnfTruth_up->Draw("same e1");
  z_UnfTruth_down->Draw("same e1");
  //z_UnfTruth_rmProbUnderFlow->Draw("same e1");
  TLegend* leg = new TLegend(0.2,0.4,0.6,0.6);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(z_Truth, "gen truth","lp");
  //if (checkTruth) {
  leg->AddEntry(z_UnfTruth, Form("unfolding truth, eff = %.1f", z_UnfTruth->Integral()*100.0/z_Truth->Integral()),"lp");
  leg->AddEntry(z_UnfTruth_up, Form("upper cut, eff = %.1f",z_UnfTruth_up->Integral()*100.0/z_Truth->Integral()),"lp");
  leg->AddEntry(z_UnfTruth_down, Form("lower cut, eff = %.1f",z_UnfTruth_down->Integral()*100.0/z_Truth->Integral()),"lp");
  //if (checkTruth) leg->AddEntry(z_UnfTruth_rmProbUnderFlow, Form("unfolding truth rm ProbFit, eff = %.1f",z_UnfTruth_rmProbUnderFlow->Integral()*100.0/z_Truth->Integral()),"lp");
  //}
  //else {
  //leg->AddEntry(z_UnfTruth, "unfolding truth","lp");
  //leg->AddEntry(z_UnfTruth_up, "upper cut","lp");
  //leg->AddEntry(z_UnfTruth_down, "lower cut","lp");
  //}
  leg->Draw("same");
  
  can->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();
  TH1D* pullHist = (TH1D*) z_UnfTruth->Clone("pullHist");
  pullHist->Divide(z_Truth);
  SetPullStyle(pullHist);
  pullHist->Draw();
  TH1D* pullHist_up = (TH1D*) z_UnfTruth_up->Clone("pullHist_up");
  pullHist_up->Divide(z_Truth);
  SetPullStyle(pullHist_up);
  //pullHist->Draw();
  TH1D* pullHist_down = (TH1D*) z_UnfTruth_down->Clone("pullHist_down");
  pullHist_down->Divide(z_Truth);
  SetPullStyle(pullHist_down);
  //pullHist->Draw();  
  TLine* l1 = new TLine(min_z,1,max_z,1);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);
  l1->Draw("same");
  pullHist->Draw("same");
  pullHist_up->Draw("same");
  pullHist_down->Draw("same");
  can->SaveAs(outputfile.c_str());
  cout <<"turth = "<<z_Truth->GetEntries()<<", unf truth = "<<z_UnfTruth->GetEntries()<<" -> eff = "<<z_UnfTruth->GetEntries()*1.0/z_Truth->GetEntries()<<", eff(Integral) = "<<z_UnfTruth->Integral()*1.0/z_Truth->Integral()<<endl;
}

void SetPullStyle(TH1D* pullHist){
  pullHist->SetTitle("");
  pullHist->GetYaxis()->SetRangeUser(0.68, 1.32);
  pullHist->GetYaxis()->SetTitle("withCuts/Truth");
  pullHist->GetYaxis()->SetNdivisions(505);
  pullHist->GetYaxis()->CenterTitle(true);
  pullHist->GetYaxis()->SetNdivisions(505);
  pullHist->GetYaxis()->SetTitleSize(25);
  pullHist->GetYaxis()->SetTitleFont(43);
  pullHist->GetYaxis()->SetTitleOffset(1.4);
  pullHist->GetYaxis()->SetLabelFont(43);
  pullHist->GetYaxis()->SetLabelSize(20);
  
  pullHist->GetXaxis()->CenterTitle(true);
  //pullHist->GetXaxis()->SetTitle("");
  pullHist->GetXaxis()->SetNdivisions(510);
  pullHist->GetXaxis()->SetTitleSize(25);
  pullHist->GetXaxis()->SetTitleFont(43);
  pullHist->GetXaxis()->SetTitleOffset(2.5);
  pullHist->GetXaxis()->SetLabelFont(43);
  pullHist->GetXaxis()->SetLabelSize(25);
  //pullHist->SetMarkerColor(kBlack);
  pullHist->SetMarkerStyle(kFullCircle);
  pullHist->SetMarkerSize(1);
  //pullHist->SetLineColor(kBlack);
  pullHist->SetLineWidth(1);
}
