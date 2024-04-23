// macro to check the AccxEff stat error
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TClonesArray.h"
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <math.h>
#include <TPaveText.h>
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include <fstream>
#include <TLegend.h>
//#include "systAccEff.C"

using namespace std;

void customHist(TH1D* hist, int iColor, int iShape);
TH1D* makePull(TH1D* oldHist, TH1D* newHist, bool systErr = false);
int getCentFromName(TString sName);

void testAccEffBins(bool isData, bool isPbPb, int dimension, bool applyPtW, string caseTag) {
	int* centbins = new int[100];
	int ncentbins = 0;

	if (isData) applyPtW = false;
	if (!isPbPb && dimension > 2) {
		cout << "[WARNING] can't do centrality in PP! checking 2D corrections" << endl;
		dimension = 2;
	}
	gStyle->SetOptStat(false);
	Double_t ptbinsAna[] = {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};
	int nBin = sizeof(ptbinsAna) / sizeof(double) - 1;

	cout << "[INFO] Loading the tree" << endl;
	TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s%s_jetR3_AccEff_JEC.root", isData ? "DATA" : "MCJPSINOPR", isPbPb ? "PbPb" : "PP", isData ? "" : "_NoBkg"), "READ");
	TTree* trNom = (TTree*)trFile->Get("treeForUnfolding");
	float z;
	float jp_pt;
	float jp_rap;
	float jp_mass;
	float jt_pt;
	float jt_rap;
	float corr_ptw;
	int centr;

	if (trNom->GetBranch("z")) trNom->SetBranchAddress("z", &z);
	if (trNom->GetBranch("jp_pt")) trNom->SetBranchAddress("jp_pt", &jp_pt);
	if (trNom->GetBranch("jp_rap")) trNom->SetBranchAddress("jp_rap", &jp_rap);
	if (trNom->GetBranch("jp_mass")) trNom->SetBranchAddress("jp_mass", &jp_mass);
	if (trNom->GetBranch("jt_pt")) trNom->SetBranchAddress("jt_pt", &jt_pt);
	if (trNom->GetBranch("jt_rap")) trNom->SetBranchAddress("jt_rap", &jt_rap);
	if (trNom->GetBranch("corr_ptw")) trNom->SetBranchAddress("corr_ptw", &corr_ptw);
	if (trNom->GetBranch("centr")) trNom->SetBranchAddress("centr", &centr);

	TFile* nomFile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root", caseTag.c_str()), "READ");
	TEfficiency* prEff = (TEfficiency*)nomFile->Get(Form("hcorr_Jpsi_%s_pr_%sEff", isPbPb ? "PbPb" : "PP", (dimension == 1) ? "1D_" : ""));
	TEfficiency* prAcc = (TEfficiency*)nomFile->Get(Form("hcorr_Jpsi_PP_pr_%sAcc", (dimension == 1) ? "1D_" : ""));

	TList* lcorr = nomFile->GetListOfKeys();
	TIter nextCorr(lcorr);

	TObjArray* corrCent = new TObjArray();
	corrCent->SetOwner(kTRUE);

	TObjString* fname(0x0);
	while ((fname = static_cast<TObjString*>(nextCorr.Next()))) {
		TEfficiency* h = static_cast<TEfficiency*>(nomFile->FindObjectAny(fname->GetString().Data()));

		TString sName(h->GetName());
		if (sName.Contains("cent")) {
			corrCent->Add(h->Clone());
			centbins[ncentbins] = getCentFromName(sName);
			ncentbins++;
		}

		cout << "adding " << sName << " to the corrction array and adding " << centbins[ncentbins - 1] << " to the centrality array" << endl;
	}

	centbins[ncentbins] = 200;
	TH1D* hist = new TH1D("hist", ";p_{T};N_{J/#psi}", nBin, ptbinsAna);
	TH2D* zhist2D = new TH2D("zhist2D", ";z;p_{T,jet}", 6, 0.064, 1, 5, 10, 60);

	double eff = 1.0;
	double weight = 0;

	int nentries = trNom->GetEntries();
	for (int jentry = 0; jentry < nentries; jentry++) {
		trNom->GetEntry(jentry);
		if (jentry % 100000 == 0) cout << "[INFO] entry = " << jentry << "/" << nentries << endl;
		if (isData) corr_ptw = 1.;
		if (fabs(jp_rap) > 2.4) continue;
		if (jp_pt < 6.5) continue;
		if (jp_mass < 2.6 || jp_mass > 3.5) continue;

		if (caseTag.find("absEta") != std::string::npos) jp_rap = fabs(jp_rap);

		if (dimension == 1)
			eff = (prEff->GetEfficiency(prEff->FindFixBin(jp_pt)) * prAcc->GetEfficiency(prAcc->FindFixBin(jp_pt)));
		else if (dimension == 2)
			eff = (prEff->GetEfficiency(prEff->FindFixBin(jp_rap, jp_pt)) * prAcc->GetEfficiency(prAcc->FindFixBin(jp_rap, jp_pt)));
		else {
			string centTag = "";
			for (int i = 0; i < ncentbins; i++) {
				if (centr <= centbins[i + 1]) {
					centTag = Form("cent%d%d", centbins[i], centbins[i + 1]);
					break;
				}
			}
			TEfficiency* prEffTemp = static_cast<TEfficiency*>(corrCent->FindObject(Form("hcorr_Jpsi_PbPb_pr_Eff_%s", centTag.c_str())));
			eff = (prEffTemp->GetEfficiency(prEffTemp->FindFixBin(jp_rap, jp_pt)) * prAcc->GetEfficiency(prAcc->FindFixBin(jp_rap, jp_pt)));
		}

		if (eff < 0.0001 || eff > 1.)
			weight = 1.;
		else
			weight = (1. / eff);

		if (applyPtW) weight = weight * corr_ptw;
		hist->Fill(jp_pt, weight);
		if (z < 0 || z > 1) continue;
		zhist2D->Fill(z, jt_pt, weight);
	}

	TFile* fsave = new TFile(Form("FilesAccxEff%s/%sPtSpectrum%sWith%sPrCorr_%s.root", caseTag.c_str(), isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", (dimension == 1) ? "" : Form("%dD", dimension), isPbPb ? "PbPb" : "PP"), "RECREATE");

	hist->Write();
	zhist2D->Write();
	int nx = zhist2D->GetNbinsY();
	for (int i = 1; i <= nx; i++) {
		TH1D* zhist1D = (TH1D*)zhist2D->ProjectionX(Form("zhist1D_%d", i), i, i);
		zhist1D->Write(Form("zhist1D_%d%d", i * 10, (i + 1) * 10));
	}
	fsave->Close();
	trFile->Close();
	delete trFile;
}

void drawComparison(bool isData, bool isPbPb, bool applyPtW, string case1, int dimension1, string label1, string case2, int dimension2, string label2) {
	gStyle->SetOptStat(false);

	if (isData) applyPtW = false;

	TFile* file1 = TFile::Open(Form("FilesAccxEff%s/%sPtSpectrum%sWith%sPrCorr_%s.root", case1.c_str(), isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", (dimension1 == 1) ? "" : Form("%dD", dimension1), isPbPb ? "PbPb" : "PP"));
	TFile* file2 = TFile::Open(Form("FilesAccxEff%s/%sPtSpectrum%sWith%sPrCorr_%s.root", case2.c_str(), isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", (dimension2 == 1) ? "" : Form("%dD", dimension2), isPbPb ? "PbPb" : "PP"));

	TH1D* hist1 = (TH1D*)file1->Get("hist");
	TH1D* hist2 = (TH1D*)file2->Get("hist");

	customHist(hist1, kRed, kFullCircle);
	customHist(hist2, kBlue, kOpenCircle);

	TCanvas* c = new TCanvas("c", "", 900, 1000);
	c->cd();
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->Draw();

	c->cd();
	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetBottomMargin(0.2);
	pad2->SetTopMargin(0.01);
	pad2->Draw();

	TLegend* leg = new TLegend(0.59, 0.6, 0.89, 0.8);
	leg->SetBorderSize(0);

	TLine* ly1 = new TLine(6.5, 1, 50, 1);
	ly1->SetLineColor(kRed);
	ly1->SetLineStyle(2);
	ly1->SetLineWidth(2);

	TH1D* pull = NULL;
	TH1D* pull2 = NULL;

	leg->AddEntry(hist1, Form("%s (x)", label1.c_str()), "lp");
	leg->AddEntry(hist2, Form("%s (ref)", label2.c_str()), "lp");

	if (hist2->GetMaximum() > hist1->GetMaximum()) hist1->GetYaxis()->SetRangeUser(0, hist2->GetMaximum() * 1.05);

	pad1->cd();
	hist1->Draw();
	hist2->Draw("same");
	leg->Draw("same");
	//pad1->SetLogy();
	pad2->cd();
	pull = makePull(hist1, hist2);
	pull->GetYaxis()->SetTitle("x/ref");
	pull->Draw();
	ly1->Draw("same");
	pull->Draw("same");
	gSystem->mkdir("Utilities/CorrBinOptimization");
	c->SaveAs(Form("Utilities/CorrBinOptimization/%sPtSpectrum%sWithPrCorr_%s_Comparison%s%s_Vs%s%s.pdf", isData ? "data" : "npr", applyPtW ? "" : "NoPtWeights", isPbPb ? "PbPb" : "PP", case1.c_str(), Form("%dD", dimension1), case2.c_str(), Form("%dD", dimension2)));
	/////////////////// draw the z distributions //////////

	hist1 = (TH1D*)file1->Get("zhist1D_3040");
	hist2 = (TH1D*)file2->Get("zhist1D_3040");

	customHist(hist1, kRed, kFullCircle);
	customHist(hist2, kBlue, kOpenCircle);

	leg = new TLegend(0.59, 0.6, 0.89, 0.8);
	leg->SetBorderSize(0);

	ly1 = new TLine(0.064, 1, 1, 1);
	ly1->SetLineColor(kRed);
	ly1->SetLineStyle(2);
	ly1->SetLineWidth(2);

	leg->AddEntry(hist1, Form("%s (x)", label1.c_str()), "lp");
	leg->AddEntry(hist2, Form("%s (ref)", label2.c_str()), "lp");

	if (hist2->GetMaximum() > hist1->GetMaximum()) hist1->GetYaxis()->SetRangeUser(0, hist2->GetMaximum() * 1.05);

	pad1->cd();
	hist1->Draw();
	hist2->Draw("same");
	leg->Draw("same");
	//pad1->SetLogy();
	pad2->cd();
	pull = makePull(hist1, hist2);
	pull->GetYaxis()->SetTitle("x/ref");
	pull->Draw();
	ly1->Draw("same");
	pull->Draw("same");
	c->SaveAs(Form("Utilities/CorrBinOptimization/%sZDist%sWithPrCorr_%s_Comparison%s%s_Vs%s%s.pdf", isData ? "data" : "npr", applyPtW ? "" : "NoPtWeights", isPbPb ? "PbPb" : "PP", case1.c_str(), Form("%dD", dimension1), case2.c_str(), Form("%dD", dimension2)));
}

void drawComparison(bool isData, bool isPbPb, bool applyPtW, string inputFileTag) {
	int col[] = {1, 46, 9, 30, 28, kMagenta + 1, kBlue + 1, kRed + 1, kCyan + 2, kGreen + 1, kOrange, kYellow + 2};
	int shape[] = {kDot, kFullCircle, kOpenCircle, kFullSquare, kOpenSquare, kFullCross, kOpenCross, kFullStar, kOpenStar, kFullTriangleUp, kOpenTriangleUp, kMultiply};
	gStyle->SetOptStat(false);

	if (isData) applyPtW = false;

	TCanvas* c = new TCanvas("c", "", 900, 1000);
	c->cd();
	TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(0.01);
	pad1->Draw();

	c->cd();
	TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
	pad2->SetBottomMargin(0.2);
	pad2->SetTopMargin(0.01);
	pad2->Draw();

	TLegend* leg = new TLegend(0.59, 0.6, 0.89, 0.8);
	leg->SetBorderSize(0);

	TLine* ly1 = new TLine(6.5, 1, 50, 1);
	ly1->SetLineColor(kRed);
	ly1->SetLineStyle(2);
	ly1->SetLineWidth(2);

	TH1D* histRef = new TH1D[100];

	TCanvas* zc = new TCanvas("zc", "", 900, 1000);
	zc->cd();
	TPad* zpad1 = new TPad("zpad1", "zpad1", 0, 0.3, 1, 1.0);
	zpad1->SetBottomMargin(0.01);
	zpad1->Draw();

	zc->cd();
	TPad* zpad2 = new TPad("zpad2", "zpad2", 0, 0.05, 1, 0.3);
	zpad2->SetBottomMargin(0.2);
	zpad2->SetTopMargin(0.01);
	zpad2->Draw();

	TLegend* zleg = new TLegend(0.59, 0.6, 0.89, 0.8);
	zleg->SetBorderSize(0);

	TLine* zly1 = new TLine(0.064, 1, 1, 1);
	zly1->SetLineColor(kRed);
	zly1->SetLineStyle(2);
	zly1->SetLineWidth(2);

	TH1D* zhistRef = new TH1D[100];

	ifstream inputFile;
	inputFile.open(Form("PlotingInputFiles/%s.txt", inputFileTag.c_str()));
	Int_t nlines = 0;

	string caseTag, dimension, legLabel, pullTitle;

	TH1D* histTemp = NULL;
	TH1D* zhistTemp = NULL;
	TH1D* hist = NULL;
	TH1D* zhist = NULL;

	while (1) {
		if (!inputFile.good()) break;
		inputFile >> caseTag >> dimension >> legLabel >> pullTitle;
		cout << "caseTag = " << caseTag << ", dimension = " << dimension << ", legLabel = " << legLabel << ", pullTitle = " << pullTitle << endl;
		if (caseTag == "end") break;

		if (nlines == 0) {
			nlines++;
			continue;
		}

		TFile* file = TFile::Open(Form("FilesAccxEff%s/%sPtSpectrum%sWith%sPrCorr_%s.root", caseTag.c_str(), isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", (dimension == "1") ? "" : Form("%sD", dimension.c_str()), isPbPb ? "PbPb" : "PP"));
		histTemp = (TH1D*)file->Get("hist");
		zhistTemp = (TH1D*)file->Get("zhist1D_3040");

		hist = (TH1D*)histTemp->Clone(Form("hist_%d", nlines));
		zhist = (TH1D*)zhistTemp->Clone(Form("zhist_%d", nlines));

		customHist(hist, col[nlines], shape[nlines]);
		customHist(zhist, col[nlines], shape[nlines]);
		leg->AddEntry(hist, Form("%s", legLabel.c_str()), "lp");

		if (nlines == 1) {
			histRef = (TH1D*)hist->Clone("histRef");
			zhistRef = (TH1D*)zhist->Clone("zhistRef");
			c->cd();
			pad1->cd();
			histRef->Draw();
			zc->cd();
			zpad1->cd();
			zhistRef->Draw();
		}

		else {
			c->cd();
			pad1->cd();
			if (hist->GetMaximum() > histRef->GetMaximum()) histRef->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.2);
			hist->Draw("same");
			pad2->cd();
			TH1D* pull = makePull(hist, histRef);
			pull->GetYaxis()->SetTitle(pullTitle.c_str());
			pull->Draw("same");

			zc->cd();
			zpad1->cd();
			if (zhist->GetMaximum() > zhistRef->GetMaximum()) zhistRef->GetYaxis()->SetRangeUser(0, zhist->GetMaximum() * 1.2);
			zhist->Draw("same");
			zpad2->cd();
			TH1D* zpull = makePull(zhist, zhistRef);
			zpull->GetYaxis()->SetTitle(pullTitle.c_str());
			zpull->Draw("same");
		}
		nlines++;
	}

	pad1->cd();
	leg->Draw("same");
	pad2->cd();
	ly1->Draw("same");
	gSystem->mkdir("Utilities/CorrBinOptimization");
	c->SaveAs(Form("Utilities/CorrBinOptimization/%sPtSpectrum%sWithPrCorr_%s_Comparison%s.pdf", isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", isPbPb ? "PbPb" : "PP", inputFileTag.c_str()));

	zpad1->cd();
	leg->Draw("same");
	zpad2->cd();
	zly1->Draw("same");
	zc->SaveAs(Form("Utilities/CorrBinOptimization/%sZDist%sWithPrCorr_%s_Comparison%s.pdf", isData ? "data" : "npr", applyPtW ? "" : isData ? "" : "NoPtWeights", isPbPb ? "PbPb" : "PP", inputFileTag.c_str()));
}

void customHist(TH1D* hist, int iColor, int iShape) {
	hist->SetLineColor(iColor);
	hist->SetFillColorAlpha(iColor - 5, 0.5);
	hist->SetMarkerColor(iColor);
	hist->SetMarkerStyle(iShape);
	hist->SetMarkerSize(1);
	//hist->GetXaxis()->SetNdivisions(505);
	hist->GetXaxis()->CenterTitle(true);
	hist->GetYaxis()->CenterTitle(true);
}

TH1D* makePull(TH1D* oldHist, TH1D* newHist, bool systErr) {
	TH1D* pullHist = (TH1D*)oldHist->Clone("pullHist");

	pullHist->Divide(newHist);

	if (systErr) {
		int nbin = pullHist->GetNbinsX();
		for (int i = 0; i <= nbin; i++) {
			pullHist->SetBinError(i, oldHist->GetBinError(i) / newHist->GetBinContent(i));
			//cout <<"bin "<<i<<" error "<<pullHist->GetBinError(i)<<"("<<100.*pullHist->GetBinError(i)/pullHist->GetBinContent(i)<<"\%)"<<endl;
		}
	}

	pullHist->SetTitle("");
	pullHist->GetYaxis()->SetRangeUser(0.78, 1.22);
	pullHist->GetYaxis()->SetTitle("old/new");
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
	//pullHist->SetMarkerStyle(kFullCircle);
	pullHist->SetMarkerSize(1);
	//pullHist->SetLineColor(kBlack);
	pullHist->SetLineWidth(2);
	return pullHist;
}

int getCentFromName(TString sName) {
	if (sName.Contains("cent0")) return 0;
	for (int i = 2; i < 200; i++)
		if (sName.Contains(Form("cent%d", i)) && !sName.Contains(Form("cent%d", i * 10))) return i;
	return 200;
}
