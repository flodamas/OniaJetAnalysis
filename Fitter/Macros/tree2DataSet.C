// -*- C++ -*-

//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program creates two RooDataSets (opposite-sign and same-sign dimuons) from an Onia Tree.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 26 19:08 CET 2016
//
//

#include "Utilities/initTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"
#include "TEfficiency.h"
#include "TObjArray.h"
#include "Utilities/JetCorrector.h"
#include "Utilities/JetUncertainty.h"

map<int, double> fCentMap; // map for centrality-Ncoll mapping
double fCentBinning[200];
int fCentBins;
int* centbins = new int[100];
int ncentbins = 0;
TObjArray* fcorrArray = NULL; // Array with the 2D correction for weighting
TH2F* corrHist = NULL;
double drmin = 0.5;
int matchGR = 0;
string findMyTree(string FileName);
string findJetTree(string FileName, double jetR);
string findSkimTree(string FileName);
string findCentTree(string FileName);
bool getTChain(TChain* fChain, TChain* jChain, TChain* sChain, TChain* cChain, vector<string> FileNames);
void iniBranch(TChain* fChain, bool isMC = false);
bool checkDS(RooDataSet* DS, string DSName);
double deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
bool isMatchedDiMuon(int iRecoDiMuon, double maxDeltaR = 0.03);
double getNColl(int centr, bool isPP);
double getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP, int cent);
bool readCorrection(const char* file);
float jecCorr(double jtPt, double rawPt, double jpsiPt);
float zjecCorr(double jtPt, double rawPt, double z);
float jeuCorr(double jtPt, double z, double jeu);
bool jetsInHCALHole(double jtEta, double jtPhi);
int getCentFromName(TString sName);

bool tree2DataSet(RooWorkspace& Workspace, vector<string> InputFileNames, string DSName, string OutputFileName, bool UpdateDS = false, bool DebugMode = false) {
	RooDataSet* dataOS = NULL;
	RooDataSet* dataSS = NULL;
	RooDataSet* dataOSNoBkg = NULL;
	bool isMC = (DSName.find("MC") != std::string::npos);

	bool isPbPb = (DSName.find("PbPb") != std::string::npos);
	bool isPrompt = (DSName.find("JPSIPR") != std::string::npos);
	//	cout << "[INFO] isPrompt = " << isPrompt << endl;

	int triggerIndex_PP = PP::HLT_HIL1DoubleMu0_v1;
	int triggerIndex_PbPb = HI::HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;
	int CentFactor = 1;

	bool applyWeight = false;
	if (isMC) applyWeight = true;
	//if (isMC && isPbPb) applyWeight = true;
	//else if (isMC && !isPbPb && DSName.find("JPSIPR")!=std::string::npos) applyWeight = true;

	bool isPureSDataset = (OutputFileName.find("_PureS") != std::string::npos);

	bool applyWeight_Corr = ((OutputFileName.find("_AccEff") != std::string::npos) || (OutputFileName.find("_lJpsiEff") != std::string::npos));

	TString corrName = "";
	TString corrFileName = "";
	if (OutputFileName.find("_AccEff") != std::string::npos) {
		corrFileName = "correction_AccEff_centMaps_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights.root"; //"correction_AccEff_centMaps_binEvtComparison_1DfinerBins.root"; //"correction_AccEff_centMaps_binEvtComparison.root";
		corrName = "AccEff";
	} else if (OutputFileName.find("_lJpsiEff") != std::string::npos) {
		corrFileName = "correction_lJpsiEff.root";
		corrName = "lJpsiEff";
	}

	bool applyJEC = (OutputFileName.find("_JEC") != std::string::npos);

	double jetR = 0.4;
	if (OutputFileName.find("_jetR3") != std::string::npos) jetR = 0.3;
	if (OutputFileName.find("_jetR5") != std::string::npos) jetR = 0.5;

	bool createDS = (gSystem->AccessPathName(OutputFileName.c_str()) || UpdateDS);

	if (!gSystem->AccessPathName(OutputFileName.c_str())) {
		cout << "[INFO] Loading RooDataSet from " << OutputFileName << endl;

		TFile* DBFile = TFile::Open(OutputFileName.c_str(), "READ");
		if (isMC && isPureSDataset) {
			if (applyWeight_Corr) {
				dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")));
				if (checkDS(dataOSNoBkg, DSName) == false) { createDS = true; }
			} else {
				dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
				if (checkDS(dataOSNoBkg, DSName) == false) { createDS = true; }
			}
		} else if (applyWeight_Corr) {
			dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")));
			if (checkDS(dataOS, DSName) == false) { createDS = true; }
		} else {
			dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
			if (checkDS(dataOS, DSName) == false) { createDS = true; }
			dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
			if (checkDS(dataSS, DSName) == false) { createDS = true; }
		}
		DBFile->Close();
		delete DBFile;
	}

	if (createDS) {
		cout << "[INFO] Creating " << (isPureSDataset ? "pure signal " : "") << "RooDataSet for " << DSName << endl;
		TreeName = findMyTree(InputFileNames[0]);
		if (TreeName == "") { return false; }
		jetTreeName = findJetTree(InputFileNames[0], jetR);
		if (jetTreeName == "") { return false; }
		skimTreeName = findSkimTree(InputFileNames[0]);
		if (skimTreeName == "") { return false; }
		centTreeName = findCentTree(InputFileNames[0]);
		if (centTreeName == "") { return false; }
		TChain* theTree = new TChain(TreeName.c_str(), "");
		TChain* jetTree = new TChain(jetTreeName.c_str(), "");
		TChain* skimTree = new TChain(skimTreeName.c_str(), "");
		TChain* centTree = new TChain(centTreeName.c_str(), "");
		if (!getTChain(theTree, jetTree, skimTree, centTree, InputFileNames)) { return false; } // Import files to TChain
		initTree(theTree);                                                                      // Initialize the Tree
		iniBranch(theTree, isMC);                                                               // Initialize the Branches

		RooRealVar* mass = new RooRealVar("invMass", "#mu#mu mass", 1.0, 6.0, "GeV/c^{2}");
		RooRealVar* zed = new RooRealVar("zed", "z_{J/#psi}", -0.1, 2);
		RooRealVar* ctau = new RooRealVar("ctau", "c_{#tau}", -100000.0, 100000.0, "mm");
		RooRealVar* ctauN = new RooRealVar("ctauN", "c_{#tau}", -100000.0, 100000.0, "");
		RooRealVar* ctauTrue = new RooRealVar("ctauTrue", "c_{#tau}", -100000.0, 100000.0, "mm");
		RooRealVar* ctauNRes = new RooRealVar("ctauNRes", "c_{#tau}", -100000.0, 100000.0, "");
		RooRealVar* ctauRes = new RooRealVar("ctauRes", "c_{#tau}", -100000.0, 100000.0, "");
		RooRealVar* ctauErr = new RooRealVar("ctauErr", "#sigma_{c#tau}", -100000.0, 100000.0, "mm");
		RooRealVar* ptQQ = new RooRealVar("pt", "#mu#mu p_{T}", -1.0, 10000.0, "GeV/c");
		RooRealVar* rapQQ = new RooRealVar("rap", "#mu#mu y", -2.5, 2.5, "");
		RooRealVar* ptJet = new RooRealVar("jetpt", "Jet p_{T}", -1.0, 10000.0, "GeV/c");
		RooRealVar* rapJet = new RooRealVar("jetrap", "Jet y", -2.5, 2.5, "");
		RooRealVar* cent = new RooRealVar("cent", "centrality", -1.0, 1000.0, "");
		RooRealVar* weight = new RooRealVar("weight", "MC weight", 0.0, 10000000.0, "");
		RooRealVar* weightCorr = new RooRealVar("weightCorr", "Data correction weight", 0.0, 10000000.0, "");
		RooArgSet* cols = NULL;

		if (applyWeight && !applyWeight_Corr) {
			if (isMC) {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weight);
				cols->add(*ctauNRes);
				cols->add(*ctauRes);
				cols->add(*ptJet);
				cols->add(*rapJet);
			} else {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weight);
				cols->add(*ctauN);
				cols->add(*ptJet);
				cols->add(*rapJet);
			}
			dataOS = new RooDataSet(Form("dOS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dOS", *cols, WeightVar(*weight), StoreAsymError(*mass));
			dataSS = new RooDataSet(Form("dSS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dSS", *cols, WeightVar(*weight), StoreAsymError(*mass));
			if (isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dOSNoBkg", *cols, WeightVar(*weight), StoreAsymError(*mass));
		} else if (applyWeight_Corr) {
			if (isMC) {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weightCorr);
				cols->add(*ctauNRes);
				cols->add(*ctauRes);
				cols->add(*ptJet);
				cols->add(*rapJet);
			} else {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weightCorr);
				cols->add(*ctauN);
				cols->add(*ptJet);
				cols->add(*rapJet);
			}
			if (!readCorrection(Form("%s/Input/%s", gSystem->ExpandPathName(gSystem->pwd()), corrFileName.Data()))) { return false; }
			dataOS = new RooDataSet(Form("dOS_%s_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")), "dOS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
			if (isMC && isPureSDataset)
				dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")), "dOSNoBkg", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
			cout << "[INFO] " << corrName << " applied!" << endl;
		} else {
			if (isMC) {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent);
				cols->add(*ctauNRes);
				cols->add(*ctauRes);
				cols->add(*ptJet);
				cols->add(*rapJet);
			} else {
				cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent);
				cols->add(*ctauN);
				cols->add(*ptJet);
				cols->add(*rapJet);
			}
			dataOS = new RooDataSet(Form("dOS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dOS", *cols, StoreAsymError(*mass));
			dataSS = new RooDataSet(Form("dSS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dSS", *cols, StoreAsymError(*mass));
			if (isMC && isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")), "dOSNoBkg", *cols, StoreAsymError(*mass));
		}
		if (applyWeight) cout << "[INFO] pt weights applied!" << endl;

		////////////////////////////////////
		Long64_t nentries = theTree->GetEntries();

		if (DebugMode) nentries /= 100;

		// creating the tree to use in the unfolding
		string fl = "";

		vector<string> jecFileName;
		string jeuFileName;

		jecFileName.push_back(Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_L2Relative_AK%dPF.txt", isMC ? "MC" : "DATA", isMC ? "MC" : "DATA", (int)(jetR * 10)));
		jeuFileName = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_Uncertainty_AK%dPF.txt", /*isMC?"MC":*/ "DATA", /*isMC?"MC":*/ "DATA", (int)(jetR * 10));
		if (!isMC)
			jecFileName.push_back(Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Input/JECDatabase/textFiles/Spring18_ppRef5TeV_V4_%s/Spring18_ppRef5TeV_V4_%s_L2L3Residual_AK%dPF.txt", isMC ? "MC" : "DATA", isMC ? "MC" : "DATA", (int)(jetR * 10)));

		JetCorrector JEC(jecFileName);
		JetUncertainty JEU(jeuFileName);
		cout << "[INFO] reading jet energy corrections from JEC file: " << jecFileName[0] << endl;
		cout << "[INFO] reading uncertainties from JEU file: " << jeuFileName << endl;
		gSystem->mkdir("TreesForUnfolding");
		string trUnfFileName = Form("TreesForUnfolding/tree_%s%s_jetR%d%s%s%s.root", DSName.c_str(), (isPureSDataset ? "_NoBkg" : ""), (int)(jetR * 10), (applyWeight_Corr ? Form("_%s", corrName.Data()) : ""), (applyJEC ? "_JEC" : ""), (applyWeight ? fl.c_str() : ""));
		TFile* trUnfFile = new TFile(trUnfFileName.c_str(), "RECREATE");
		//trUnfFile->cd();
		TTree* trUnf = new TTree("treeForUnfolding", "tree used for the unfolding");
		Int_t evtNb;
		Int_t centr;
		Float_t jp_pt;
		Float_t jp_rap;
		Float_t jp_eta;
		Float_t jp_mass;
		Float_t jp_phi;
		Float_t jp_l;
		Float_t jp_gen_pt;
		Float_t jp_gen_rap;
		Float_t jp_gen_eta;
		Float_t jp_gen_phi;
		Float_t jt_pt;
		Float_t jt_pt_noZJEC;
		Float_t jt_pt_jettyCorr;
		Float_t jt_pt_JEU_Up;
		Float_t jt_pt_JEU_Down;
		Float_t jt_rap;
		Float_t jt_eta;
		Float_t jt_phi;
		Float_t jt_CHF;
		Float_t jt_NHF;
		Float_t jt_CEF;
		Float_t jt_NEF;
		Float_t jt_MUF;
		Float_t jt_CHM;
		Float_t jt_NHM;
		Float_t jt_CEM;
		Float_t jt_NEM;
		Float_t jt_MUM;
		Float_t jt_ref_pt;
		Float_t jt_ref_rap;
		Float_t jt_ref_eta;
		Float_t jt_ref_phi;
		Float_t jt_pt_genZJEC;
		Float_t jt_gen_chargedSum;
		Float_t jt_gen_hardSum;
		Float_t jt_signal_chargedSum;
		Float_t jt_signal_hardSum;
		Float_t z;
		Float_t gen_z;
		Float_t z_jettyCorr;
		Float_t corr_AccEff;
		Float_t pt_hat;
		Float_t corr_ptw;

		cout << "[INFO] Creating the tree to use in the unfolding" << endl;
		trUnf->Branch("evtNb", &evtNb, "evtNb/I");
		trUnf->Branch("centr", &centr, "centr/I");
		trUnf->Branch("jp_pt", &jp_pt, "jp_pt/F");
		trUnf->Branch("jp_rap", &jp_rap, "jp_rap/F");
		trUnf->Branch("jp_eta", &jp_eta, "jp_eta/F");
		trUnf->Branch("jp_phi", &jp_phi, "jp_phi/F");
		trUnf->Branch("jp_mass", &jp_mass, "jp_mass/F");
		trUnf->Branch("jp_l", &jp_l, "jp_l/F");
		trUnf->Branch("jt_pt", &jt_pt, "jt_pt/F");
		trUnf->Branch("jt_pt_noZJEC", &jt_pt_noZJEC, "jt_pt_noZJEC/F");
		trUnf->Branch("jt_pt_jettyCorr", &jt_pt_jettyCorr, "jt_pt_jettyCorr/F");
		trUnf->Branch("jt_pt_JEU_Up", &jt_pt_JEU_Up, "jt_pt_JEU_Up/F");
		trUnf->Branch("jt_pt_JEU_Down", &jt_pt_JEU_Down, "jt_pt_JEU_Down/F");
		trUnf->Branch("jt_rap", &jt_rap, "jt_rap/F");
		trUnf->Branch("jt_eta", &jt_eta, "jt_eta/F");
		trUnf->Branch("jt_phi", &jt_phi, "jt_phi/F");
		trUnf->Branch("z", &z, "z/F");
		trUnf->Branch("z_jettyCorr", &z_jettyCorr, "z_jettyCorr/F");
		trUnf->Branch("corr_AccEff", &corr_AccEff, "corr_AccEff/F");
		if (isMC) {
			trUnf->Branch("corr_ptw", &corr_ptw, "corr_ptw/F");
			trUnf->Branch("pt_hat", &pt_hat, "pt_hat/F");
			trUnf->Branch("jp_gen_pt", &jp_gen_pt, "jp_gen_pt/F");
			trUnf->Branch("jp_gen_rap", &jp_gen_rap, "jp_gen_rap/F");
			trUnf->Branch("jp_gen_eta", &jp_gen_eta, "jp_gen_eta/F");
			trUnf->Branch("jp_gen_phi", &jp_gen_phi, "jp_gen_phi/F");
			trUnf->Branch("jt_ref_pt", &jt_ref_pt, "jt_ref_pt/F");
			trUnf->Branch("jt_ref_rap", &jt_ref_rap, "jt_ref_rap/F");
			trUnf->Branch("jt_ref_eta", &jt_ref_eta, "jt_ref_eta/F");
			trUnf->Branch("jt_ref_phi", &jt_ref_phi, "jt_ref_phi/F");
			trUnf->Branch("jt_pt_genZJEC", &jt_pt_genZJEC, "jt_pt_genZJEC/F");
			trUnf->Branch("gen_z", &gen_z, "gen_z/F");
		}

		cout << "[INFO] Starting to process " << nentries << " nentries" << endl;

		vector<Long64_t> debugList;
		Long64_t nCandidates = 0;

		for (Long64_t jentry = 0; jentry < nentries; jentry++) {
			if (jentry % 1000000 == 0) cout << "\r[INFO] " << jentry << "/" << nentries << flush;

			pPAprimaryVertexFilter = 1;
			pBeamScrapingFilter = 1;

			if (theTree->LoadTree(jentry) < 0) break;
			if (theTree->GetTreeNumber() != fCurrent) {
				fCurrent = theTree->GetTreeNumber();
				//cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
			}
			Reco_QQ_4mom->Clear();
			Reco_mu_4mom->Clear();

			if (isMC) {
				Gen_QQ_4mom->Clear();
				Gen_mu_4mom->Clear();
			}
			theTree->GetEntry(jentry);

			bool jetFound = false;

			if (!(isPbPb ? (pprimaryVertexFilter && pclusterCompatibilityFilter && phfCoincFilter2Th4) : (pPAprimaryVertexFilter && pBeamScrapingFilter))) continue;

			// reco'ed dimuons
			for (int iQQ = 0; iQQ < Reco_QQ_size; iQQ++) {
				int iMupl = Reco_QQ_mupl_idx[iQQ];
				int iMumi = Reco_QQ_mumi_idx[iQQ];

				if (iMupl < 0 || iMumi < 0) {
					//cout << "[WARNING] entry " << jentry << " contains a reconstructed dimuon index " << iQQ << " with negative muon indices!!! Will not consider this reco'ed candidate...\n";
					debugList.push_back(jentry);
					continue;
				}

				if (!RecoQQ::areMuonsInAcceptance2019(iQQ)) continue; // 2019 Global Muon Acceptance Cuts

				if (!RecoQQ::passQualityCuts2019(iQQ)) continue; // 2019 Soft Global Muon Quality Cuts

				//if (!RecoQQ::isTriggerMatch(iQQ, isPbPb ? triggerIndex_PbPb : triggerIndex_PP)) continue; // none of the QQ candidates are matching the trigger right now, to be investigated...

				drmin = 0.5;
				zed->setVal(100);
				z = 100;
				ptJet->setVal(1);
				rapJet->setVal(100);

				TLorentzVector* RecoQQ4mom = (TLorentzVector*)Reco_QQ_4mom->At(iQQ);

				mass->setVal(RecoQQ4mom->M());
				if (theTree->GetBranch("Reco_QQ_ctau3D")) {
					ctau->setVal(Reco_QQ_ctau3D[iQQ]);
				} else if (theTree->GetBranch("Reco_QQ_ctau")) {
					ctau->setVal(Reco_QQ_ctau[iQQ]);
				} else {
					cout << "[ERROR] No ctau information found in the Onia Tree" << endl;
				}
				if (theTree->GetBranch("Reco_QQ_ctauErr3D")) {
					ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]);
				} else if (theTree->GetBranch("Reco_QQ_ctauErr")) {
					ctauErr->setVal(Reco_QQ_ctauErr[iQQ]);
				} else {
					cout << "[ERROR] No ctauErr information found in the Onia Tree" << endl;
				}

				ctauN->setVal(ctau->getVal() / ctauErr->getVal());
				ptQQ->setVal(RecoQQ4mom->Pt());
				rapQQ->setVal(RecoQQ4mom->Rapidity());
				if (isMC) {
					cent->setVal(HI::getMCHiBinFromhiHF(hiHF) * CentFactor);
					centr = HI::getMCHiBinFromhiHF(hiHF); //hiBin;
				} else {
					cent->setVal(HI::getHiBinFromhiHF(hiHF) * CentFactor);
					centr = HI::getHiBinFromhiHF(hiHF); //hiBin;
				}
				jp_pt = RecoQQ4mom->Pt();
				jp_rap = RecoQQ4mom->Rapidity();
				jp_eta = RecoQQ4mom->Eta();
				jp_phi = RecoQQ4mom->Phi();
				jp_mass = RecoQQ4mom->M();
				jp_l = Reco_QQ_ctau[iQQ];
				pt_hat = pthat;
				corr_ptw = 1;
				double rawPt = 0;

				for (Long64_t ijet = 0; ijet < nref; ijet++) {
					if (jetFound) continue;
					if (isPbPb && jetsInHCALHole(jteta[ijet], jtphi[ijet])) continue;
					TLorentzVector v_jet;
					JEC.SetJetPT(rawpt[ijet]);
					JEC.SetJetEta(jteta[ijet]);
					JEC.SetJetPhi(jtphi[ijet]);
					JEC.SetRho(1);
					JEC.SetJetArea(jtarea[ijet]);
					jt_pt_noZJEC = JEC.GetCorrectedPT();

					JEU.SetJetPT(jt_pt_noZJEC);
					JEU.SetJetEta(jteta[ijet]);
					JEU.SetJetPhi(jtphi[ijet]);

					v_jet.SetPtEtaPhiM(jt_pt_noZJEC, jteta[ijet], jtphi[ijet], jtm[ijet]);

					if (RecoQQ4mom->DeltaR(v_jet) <= drmin) {
						jetFound = true;

						drmin = RecoQQ4mom->DeltaR(v_jet);
						if (applyJEC) {
							jt_pt = jecCorr(jt_pt_noZJEC, rawpt[ijet], jp_pt);
							zed->setVal(jp_pt / jt_pt);
							if (zed->getVal() > 1 && zed->getVal() <= 1.000001) zed->setVal(0.9999999);
							ptJet->setVal(jt_pt);
							z = jp_pt / jt_pt;
							if (z > 1 && z <= 1.000001) z = 0.9999999;
							jt_pt_JEU_Down = jeuCorr(jt_pt, z, -1 * (JEU.GetUncertainty().first)); //jt_pt * (1 - JEU.GetUncertainty().first);
							jt_pt_JEU_Up = jeuCorr(jt_pt, z, JEU.GetUncertainty().second);         //jt_pt * (1 + JEU.GetUncertainty().second);
						} else {
							zed->setVal(jp_pt / jt_pt_noZJEC);
							if (zed->getVal() > 1 && zed->getVal() <= 1.000001) zed->setVal(0.9999999);
							jt_pt = jt_pt_noZJEC;
							jt_pt_JEU_Down = jt_pt * (1 - JEU.GetUncertainty().first);
							jt_pt_JEU_Up = jt_pt * (1 + JEU.GetUncertainty().second);
							ptJet->setVal(jt_pt);
							z = jp_pt / jt_pt;
							if (z > 1 && z <= 1.000001) z = 0.9999999;
						}
						rapJet->setVal(jty[ijet]);
						jt_rap = jty[ijet];
						jt_eta = jteta[ijet];
						jt_phi = jtphi[ijet];
						rawPt = rawpt[ijet];

						// just read the jet pT branches for now, TO BE FIXED
						jt_pt = jtpt[ijet]; // supposed to be corrected already, TBC by Matt
						z = jp_pt / jt_pt;

						if (DebugMode) cout << "\n[Debug] Jpsi-in-jet candidate: Jpsi pT = " << jp_pt << ", (raw) jet pT = (" << rawPt << ") " << jt_pt << ", z = " << z << endl;

						if (rawpt[ijet] > jp_pt) {
							JEC.SetJetPT(rawpt[ijet] - jp_pt);
							jt_pt_jettyCorr = jp_pt + JEC.GetCorrectedPT();
						} else
							jt_pt_jettyCorr = jp_pt;
						z_jettyCorr = jp_pt / jt_pt_jettyCorr;

						//jt_pt_jettyCorr = jecCorr(jt_pt_jettyCorr, rawpt[ijet], RecoQQ4mom->Pt());

						if (isMC) {
							jt_ref_pt = refpt[ijet];
							jt_ref_rap = refy[ijet];
							jt_ref_eta = refeta[ijet];
							jt_ref_phi = refphi[ijet];
						}
					}
					//continue;
				} // end of jet loop

				if (isMC) {
					if (isMatchedDiMuon(iQQ)) {
						if (theTree->GetBranch("Gen_QQ_ctau3D")) {
							ctauTrue->setVal(Gen_QQ_ctau3D[Reco_QQ_whichGen[iQQ]]);
						} else if (theTree->GetBranch("Gen_QQ_ctau")) {
							ctauTrue->setVal(Gen_QQ_ctau[Reco_QQ_whichGen[iQQ]]);
						} else {
							cout << "[ERROR] No ctauTrue information found in the Onia Tree" << endl;
						}
						ctauNRes->setVal((ctau->getValV() - ctauTrue->getValV()) / (ctauErr->getValV()));
						ctauRes->setVal((ctau->getValV() - ctauTrue->getValV()));
					}
				}
				if (applyWeight && !applyWeight_Corr) {
					double w = Gen_weight;
					weight->setVal(w);
				} else if (applyWeight_Corr) {
					double wCorr = 1.0 / getCorr(RecoQQ4mom->Rapidity(), RecoQQ4mom->Pt(), RecoQQ4mom->M(), !isPbPb, centr);
					//// add pt weights for prompt MC
					if (applyWeight) {
						if (!isMC) cout << "[WARNING] this is data and you are applying MC weights" << endl;
						if (!isPrompt)
							corr_ptw = Gen_weight; //cout<<"corr_ptw = "<<corr_ptw<<endl;
						else {
							if (pthat >= 15 && pthat < 25)
								corr_ptw = 0.0247699;
							else if (pthat >= 25 && pthat < 35)
								corr_ptw = 0.00311931;
							else if (pthat >= 35 && pthat < 45)
								corr_ptw = 0.000693027;
							else if (pthat >= 45)
								corr_ptw = 0.000212618;
						}
					}
					corr_AccEff = wCorr;
					weightCorr->setVal(wCorr * corr_ptw);
				}

				if (Reco_QQ_sign[iQQ] == 0) { // Opposite-Sign dimuons
					if (isMC && isPureSDataset && isMatchedDiMuon(iQQ)) {
						if (applyWeight_Corr)
							dataOSNoBkg->add(*cols, weightCorr->getVal()); //Signal-only dimuons
						else
							dataOSNoBkg->add(*cols, (applyWeight ? weight->getVal() : 1.0)); // Signal-only dimuons
						TLorentzVector* GenQQ4mom = (TLorentzVector*)Gen_QQ_4mom->At(Reco_QQ_whichGen[iQQ]);
						jp_gen_pt = GenQQ4mom->Pt();
						jp_gen_rap = GenQQ4mom->Rapidity();
						jp_gen_eta = GenQQ4mom->Eta();
						jp_gen_phi = GenQQ4mom->Phi();
						gen_z = jp_gen_pt / jt_ref_pt;
						if (gen_z > 1 && gen_z <= 1.000001) gen_z = 0.9999999;
						jt_pt_genZJEC = zjecCorr(jt_pt_noZJEC, rawPt, gen_z);
					} else if (isMC && isPureSDataset && !isMatchedDiMuon(iQQ))
						gen_z = -1;

					//if (isPbPb && !isMC && z<0.376 && jp_pt>6.5 && jt_pt>10 && jt_pt<20 && fabs(jt_eta)<2 && corr_AccEff>1000) {cout <<"found a problematic event with w = "<<corr_AccEff<<endl; continue;}
					//if (jp_pt>6.5) cout <<"less than 6.5"<<endl;

					if (z < 100 && z > 1 && jt_pt > 15)
						cout << "[WARNING] z = " << z << " in entry " << jentry << endl;
					else if (applyWeight_Corr)
						dataOS->add(*cols, weightCorr->getVal()); //Signal and background dimuons
					else
						dataOS->add(*cols, (applyWeight ? weight->getVal() : 1.0)); //Signal and background dimuons
					evtNb = jentry;
					if (isMC && isPureSDataset && gen_z >= 0 && z < 100)
						trUnf->Fill();
					else if (!isPureSDataset) // && z < 100)
						trUnf->Fill();

				} else { // Like-Sign dimuons
					if (!isPureSDataset && !applyWeight_Corr) dataSS->add(*cols, (applyWeight ? weight->getVal() : 1.0));
				}

				nCandidates++;

			} // end of reco QQ loop
		}

		cout << "\n[INFO] End of the event loop! Found " << nCandidates << " J/psi-in-jet candidates\n";

		//cout << "\n[DEBUG] List of entries for which the reconstructed dimuon has negative muon indices:\n";
		//for (auto i : debugList) cout << i << ",";

		// Close the TChain and all its pointers
		delete Reco_QQ_4mom;
		delete Reco_mu_4mom;
		delete Gen_mu_4mom;
		theTree->Reset();
		delete theTree;
		// Save the tree for the unfolding
		cout << "[INFO] Saving the tree for the unfolding" << endl;
		//trUnf->Write("treeForUnfolding");
		trUnfFile->cd();
		trUnf->Write();
		trUnfFile->Close();
		delete trUnfFile;
		// Save all the datasets
		TFile* DBFile = TFile::Open(OutputFileName.c_str(), "RECREATE");
		DBFile->cd();
		if (isMC && isPureSDataset && applyWeight_Corr) {
			//cout<<"[INFO] saving name = "<<Form("dOS_%s_NoBkg_jetR%d_%s%s", DSName.c_str(), (int)(jetR*10), corrName.Data(),(applyJEC?"_JEC":""))<<endl;
			dataOSNoBkg->Write(Form("dOS_%s_NoBkg_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")));
		} else if (isMC && isPureSDataset && !applyWeight_Corr) {
			//cout<<"[INFO] saving name = "<<Form("dOS_%s_NoBkg_jetR%d%s", DSName.c_str(),(int)(jetR*10), (applyJEC?"_JEC":""))<<endl;
			dataOSNoBkg->Write(Form("dOS_%s_NoBkg_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
		} else if (!isPureSDataset && applyWeight_Corr) {
			dataOS->Write(Form("dOS_%s_jetR%d_%s%s", DSName.c_str(), (int)(jetR * 10), corrName.Data(), (applyJEC ? "_JEC" : "")));
		} else {
			dataOS->Write(Form("dOS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
			dataSS->Write(Form("dSS_%s_jetR%d%s", DSName.c_str(), (int)(jetR * 10), (applyJEC ? "_JEC" : "")));
		}
		DBFile->Write();
		DBFile->Close();
		delete DBFile;
	}

	// Import datasets to workspace
	if (isMC && isPureSDataset) {
		if (!dataOSNoBkg) {
			cout << "[ERROR] " << DSName << "_NoBkg was not found" << endl;
			return false;
		}
		Workspace.import(*dataOSNoBkg);
	} else if (applyWeight_Corr) {
		if (!dataOS) {
			cout << "[ERROR] " << DSName << "_" << corrName.Data() << " was not found" << endl;
			return false;
		}
		Workspace.import(*dataOS);
	} else {
		if (!dataOS || !dataSS) {
			cout << "[ERROR] " << DSName << " was not found" << endl;
			return false;
		}
		Workspace.import(*dataOS);
		Workspace.import(*dataSS);
	}

	// delete the local datasets
	delete dataSS;
	delete dataOS;
	delete dataOSNoBkg;

	// delete the correction array
	if (fcorrArray) delete fcorrArray;
	delete corrHist;
	return true;
};

string findMyTree(string FileName) {
	TFile* f = TFile::Open(FileName.c_str(), "READ");
	string name = "";
	if (f->GetListOfKeys()->Contains("hionia"))
		name = "hionia/myTree";
	else if (f->GetListOfKeys()->Contains("myTree"))
		name = "myTree";
	else {
		cout << "[ERROR] myTree was not found in: " << FileName << endl;
	}
	f->Close();
	delete f;
	return name;
};

string findJetTree(string FileName, double jetR) {
	TFile* f = TFile::Open(FileName.c_str(), "READ");
	string name = "";
	if (jetR == 0.4) {
		if (f->GetListOfKeys()->Contains("ak4PFJetAnalyzer")) name = "ak4PFJetAnalyzer/t";
		//else if (f->GetListOfKeys()->Contains("akCs4PFJetAnalyzer")) name = "akCs4PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("akFlowPuCs4PFJetAnalyzer"))
			name = "akFlowPuCs4PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("t"))
			name = "t";
		else {
			cout << "[ERROR] t was not found in: " << FileName << endl;
		}
	} else if (jetR == 0.3) {
		if (f->GetListOfKeys()->Contains("ak3PFJetAnalyzer")) name = "ak3PFJetAnalyzer/t";
		//else if (f->GetListOfKeys()->Contains("akCs3PFJetAnalyzer")) name = "akCs3PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("akFlowPuCs3PFJetAnalyzer"))
			name = "akFlowPuCs3PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("t"))
			name = "t";
		else {
			cout << "[ERROR] t was not found in: " << FileName << endl;
		}
	} else if (jetR == 0.5) {
		if (f->GetListOfKeys()->Contains("ak5PFJetAnalyzer")) name = "ak4PFJetAnalyzer/t";
		//else if (f->GetListOfKeys()->Contains("akCs5PFJetAnalyzer")) name = "akCs5PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("akFlowPuCs5PFJetAnalyzer"))
			name = "akFlowPuCs5PFJetAnalyzer/t";
		else if (f->GetListOfKeys()->Contains("t"))
			name = "t";
		else {
			cout << "[ERROR] t was not found in: " << FileName << endl;
		}
	}
	f->Close();
	delete f;
	return name;
};

string findSkimTree(string FileName) {
	TFile* f = TFile::Open(FileName.c_str(), "READ");
	string name = "";
	if (f->GetListOfKeys()->Contains("skimanalysis"))
		name = "skimanalysis/HltTree";
	else if (f->GetListOfKeys()->Contains("hltanalysis"))
		name = "hltanalysis/HltTree";
	else {
		cout << "[ERROR] HltTree was not found in: " << FileName << endl;
	}
	f->Close();
	delete f;
	return name;
};

string findCentTree(string FileName) {
	TFile* f = TFile::Open(FileName.c_str(), "READ");
	string name = "";
	if (f->GetListOfKeys()->Contains("hiEvtAnalyzer"))
		name = "hiEvtAnalyzer/HiTree";
	else if (f->GetListOfKeys()->Contains("HltTree"))
		name = "HiTree";
	else {
		cout << "[ERROR] HiTree was not found in: " << FileName << endl;
	}
	f->Close();
	delete f;
	return name;
};

bool getTChain(TChain* fChain, TChain* jChain, TChain* sChain, TChain* cChain, vector<string> FileNames) {
	cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;
	for (vector<string>::iterator FileName = FileNames.begin(); FileName != FileNames.end(); ++FileName) {
		cout << "[INFO] Adding TFile " << FileName->c_str() << endl;
		fChain->Add(Form("%s/%s", FileName->c_str(), TreeName.c_str()));
		jChain->Add(Form("%s/%s", FileName->c_str(), jetTreeName.c_str()));
		sChain->Add(Form("%s/%s", FileName->c_str(), skimTreeName.c_str()));
		cChain->Add(Form("%s/%s", FileName->c_str(), centTreeName.c_str()));
	}
	if (jChain)
		fChain->AddFriend(jChain);
	if (sChain)
		fChain->AddFriend(sChain);
	if (cChain)
		fChain->AddFriend(cChain);

	if (!fChain) {
		cout << "[ERROR] fChain was not created, some input files are missing" << endl;
		return false;
	}
	return true;
};

void iniBranch(TChain* fChain, bool isMC) {
	cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
	if (fChain->GetBranch("Reco_QQ_4mom")) { fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false); }
	if (fChain->GetBranch("Reco_mu_4mom")) { fChain->GetBranch("Reco_mu_4mom")->SetAutoDelete(false); }
	if (isMC) {
		if (fChain->GetBranch("Gen_QQ_4mom")) { fChain->GetBranch("Gen_QQ_4mom")->SetAutoDelete(false); }
		if (fChain->GetBranch("Gen_mu_4mom")) { fChain->GetBranch("Gen_mu_4mom")->SetAutoDelete(false); }
	}
	fChain->SetBranchStatus("*", 0);
	RecoQQ::iniBranches(fChain);
	if (fChain->GetBranch("runNb")) { fChain->SetBranchStatus("runNb", 1); }
	if (fChain->GetBranch("Centrality")) { fChain->SetBranchStatus("Centrality", 1); }
	if (fChain->GetBranch("Reco_QQ_size")) { fChain->SetBranchStatus("Reco_QQ_size", 1); }
	if (fChain->GetBranch("Reco_QQ_sign")) { fChain->SetBranchStatus("Reco_QQ_sign", 1); }
	if (fChain->GetBranch("Reco_QQ_4mom")) { fChain->SetBranchStatus("Reco_QQ_4mom", 1); }
	if (fChain->GetBranch("Reco_QQ_mupl_idx")) { fChain->SetBranchStatus("Reco_QQ_mupl_idx", 1); }
	if (fChain->GetBranch("Reco_QQ_mumi_idx")) { fChain->SetBranchStatus("Reco_QQ_mumi_idx", 1); }
	if (fChain->GetBranch("Reco_mu_4mom")) { fChain->SetBranchStatus("Reco_mu_4mom", 1); }
	if (fChain->GetBranch("Reco_QQ_ctau3D")) { fChain->SetBranchStatus("Reco_QQ_ctau3D", 1); }
	if (fChain->GetBranch("Reco_QQ_ctauErr3D")) { fChain->SetBranchStatus("Reco_QQ_ctauErr3D", 1); }
	if (fChain->GetBranch("Reco_QQ_ctau")) { fChain->SetBranchStatus("Reco_QQ_ctau", 1); }
	if (fChain->GetBranch("Reco_QQ_ctauErr")) { fChain->SetBranchStatus("Reco_QQ_ctauErr", 1); }
	if (fChain->GetBranch("nref")) { fChain->SetBranchStatus("nref", 1); }
	if (fChain->GetBranch("jtpt")) { fChain->SetBranchStatus("jtpt", 1); }
	if (fChain->GetBranch("rawpt")) { fChain->SetBranchStatus("rawpt", 1); }
	if (fChain->GetBranch("jteta")) { fChain->SetBranchStatus("jteta", 1); }
	if (fChain->GetBranch("jty")) { fChain->SetBranchStatus("jty", 1); }
	if (fChain->GetBranch("jtphi")) { fChain->SetBranchStatus("jtphi", 1); }
	if (fChain->GetBranch("jtm")) { fChain->SetBranchStatus("jtm", 1); }
	if (fChain->GetBranch("jtarea")) { fChain->SetBranchStatus("jtarea", 1); }
	//if (fChain->GetBranch("jtPfCHF"))               { fChain->SetBranchStatus("jtPfCHF",1);               }
	//if (fChain->GetBranch("jtPfNHF"))               { fChain->SetBranchStatus("jtPfNHF",1);               }
	//if (fChain->GetBranch("jtPfCEF"))               { fChain->SetBranchStatus("jtPfCEF",1);               }
	//if (fChain->GetBranch("jtPfNEF"))               { fChain->SetBranchStatus("jtPfNEF",1);               }
	//if (fChain->GetBranch("jtPfMUF"))               { fChain->SetBranchStatus("jtPfMUF",1);               }
	//if (fChain->GetBranch("jtPfCHM"))               { fChain->SetBranchStatus("jtPfCHM",1);               }
	//if (fChain->GetBranch("jtPfNHM"))               { fChain->SetBranchStatus("jtPfNHM",1);               }
	//if (fChain->GetBranch("jtPfCEM"))               { fChain->SetBranchStatus("jtPfCEM",1);               }
	//if (fChain->GetBranch("jtPfNEM"))               { fChain->SetBranchStatus("jtPfNEM",1);               }
	//if (fChain->GetBranch("jtPfMUM"))               { fChain->SetBranchStatus("jtPfMUM",1);               }
	if (fChain->GetBranch("pPAprimaryVertexFilter")) { fChain->SetBranchStatus("pPAprimaryVertexFilter", 1); }
	if (fChain->GetBranch("pBeamScrapingFilter")) { fChain->SetBranchStatus("pBeamScrapingFilter", 1); }
	if (fChain->GetBranch("pprimaryVertexFilter")) { fChain->SetBranchStatus("pprimaryVertexFilter", 1); }
	if (fChain->GetBranch("pclusterCompatibilityFilter")) { fChain->SetBranchStatus("pclusterCompatibilityFilter", 1); }
	if (fChain->GetBranch("phfCoincFilter2Th4")) { fChain->SetBranchStatus("phfCoincFilter2Th4", 1); }
	if (fChain->GetBranch("hiBin")) { fChain->SetBranchStatus("hiBin", 1); }
	if (fChain->GetBranch("hiHF")) { fChain->SetBranchStatus("hiHF", 1); }

	if (isMC) {
		if (fChain->GetBranch("pthat")) { fChain->SetBranchStatus("pthat", 1); }
		if (fChain->GetBranch("Gen_QQ_4mom")) { fChain->SetBranchStatus("Gen_QQ_4mom", 1); }
		if (fChain->GetBranch("Gen_QQ_size")) { fChain->SetBranchStatus("Gen_QQ_size", 1); }
		if (fChain->GetBranch("Gen_QQ_mupl_idx")) { fChain->SetBranchStatus("Gen_QQ_mupl_idx", 1); }
		if (fChain->GetBranch("Gen_QQ_mumi_idx")) { fChain->SetBranchStatus("Gen_QQ_mumi_idx", 1); }
		if (fChain->GetBranch("Gen_QQ_ctau3D")) { fChain->SetBranchStatus("Gen_QQ_ctau3D", 1); }
		if (fChain->GetBranch("Gen_QQ_ctau")) { fChain->SetBranchStatus("Gen_QQ_ctau", 1); }
		if (fChain->GetBranch("Reco_QQ_whichGen")) { fChain->SetBranchStatus("Reco_QQ_whichGen", 1); }
		if (fChain->GetBranch("Gen_weight")) { fChain->SetBranchStatus("Gen_weight", 1); }
		if (fChain->GetBranch("refpt")) { fChain->SetBranchStatus("refpt", 1); }
		if (fChain->GetBranch("refeta")) { fChain->SetBranchStatus("refeta", 1); }
		if (fChain->GetBranch("refy")) { fChain->SetBranchStatus("refy", 1); }
		if (fChain->GetBranch("refphi")) { fChain->SetBranchStatus("refphi", 1); }
		if (fChain->GetBranch("refm")) { fChain->SetBranchStatus("refm", 1); }
		//if (fChain->GetBranch("genChargedSum"))       { fChain->SetBranchStatus("genChargedSum",1);         }
		//if (fChain->GetBranch("genHardSum"))          { fChain->SetBranchStatus("genHardSum",1);            }
		//if (fChain->GetBranch("signalChargedSum"))    { fChain->SetBranchStatus("signalChargedSum",1);      }
		//if (fChain->GetBranch("signalHardSum"))       { fChain->SetBranchStatus("signalHardSum",1);         }
	}
};

bool checkDS(RooDataSet* DS, string DSName) {
	bool incCent = (DSName.find("PbPb") != std::string::npos);
	bool incCtauTrue = (DSName.find("MC") != std::string::npos);
	const RooArgSet* row = DS->get();
	if (
	  (row->find("invMass") != 0) &&
	  (row->find("pt") != 0) &&
	  (row->find("ctau") != 0) &&
	  (row->find("ctauErr") != 0) &&
	  (row->find("zed") != 0) &&
	  (row->find("jetpt") != 0) &&
	  (incCent ? row->find("cent") != 0 : true) &&
	  (incCtauTrue ? row->find("ctauTrue") != 0 : true) &&
	  (incCtauTrue ? row->find("ctauRes") != 0 : true) &&
	  (incCtauTrue ? row->find("ctauNRes") != 0 : true) &&
	  (incCtauTrue ? true : row->find("ctauN") != 0)) {
		return true;
	} else {
		cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << endl;
	}

	return false;
};

double deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon) {
	double dEta = RecoMuon->Eta() - GenMuon->Eta();
	double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
	return ((double)TMath::Sqrt((dEta * dEta) + (dPhi * dPhi)));
};

bool isMatchedDiMuon(int iRecoDiMuon, double maxDeltaR) {
	if (Reco_QQ_whichGen[iRecoDiMuon] != -1) return true;
	return false;
};

double getNColl(int cen, bool isPP) {
	// Returns the corresponding Ncoll value to the "cen" centrality bin

	if (isPP) return 1.;

	int normCent = TMath::Nint(cen / 2.);

	int lcent = 0;
	int ucent = 0;
	for (int i = 0; i < fCentBins; i++) {
		ucent = fCentBinning[i];
		if ((normCent >= lcent) && (normCent < ucent))
			return fCentMap[ucent];
		else
			lcent = ucent;
	}
	return 1.;
};

bool readCorrection(const char* file) {
	TFile* froot = new TFile(file, "READ");
	if (!froot) {
		cout << "[ERROR] File " << file << " for correction of events not found" << endl;
		return false;
	}
	TList* lcorr = froot->GetListOfKeys();
	TIter nextCorr(lcorr);

	fcorrArray = new TObjArray();
	fcorrArray->SetOwner(kTRUE);

	TObjString* fname(0x0);
	while ((fname = static_cast<TObjString*>(nextCorr.Next()))) {
		TEfficiency* h = static_cast<TEfficiency*>(froot->FindObjectAny(fname->GetString().Data()));

		TString sName(h->GetName());
		if (sName.Contains("hcorr")) {
			fcorrArray->Add(h->Clone());
			if (sName.Contains("cent")) {
				centbins[ncentbins] = getCentFromName(sName);
				ncentbins++;
			}

			cout << "[INFO] Adding " << sName << " to the correction array" << endl;
		} else
			cout << "[WARNING] Correction histo " << sName.Data() << " not according to naming convention. Not included in correction array" << endl;
	}

	if (!(fcorrArray->GetSize() > 0)) {
		cout << "[ERROR] Correction array empty: No corrections found." << endl;
		return false;
	}
	centbins[ncentbins] = 200;
	delete lcorr;
	//froot->Close(); delete froot;
	return true;
};

double getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP, int cent) {
	const char* collName = "PbPb";
	const char* massName = "Jpsi";
	if (isPP)
		collName = "PP";

	if (!fcorrArray) {
		cout << "[ERROR] No correction array exist" << endl;
		return 0;
	}

	Double_t corr = 1.;
	string centTag = "";
	if (!isPP) {
		for (int iCent = 0; iCent < ncentbins; iCent++) {
			if (cent <= centbins[iCent + 1]) {
				centTag = Form("_cent%d%d", centbins[iCent], centbins[iCent + 1]);
				break;
			}
		}
	}
	TString bfracFunc = "0.611808-0.718494*exp(-0.102699*x)";
	if (isPP)
		bfracFunc = "0.669157-0.784588*exp(-0.0731804*x)";
	//TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
	TF1* bfrac = new TF1("bfrac", bfracFunc, 3, 100);
	TEfficiency* prEffHisto = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s_pr_Eff%s", collName, centTag.c_str())));
	TEfficiency* prAccHisto = static_cast<TEfficiency*>(fcorrArray->FindObject("hcorr_Jpsi_PP_pr_Acc"));

	//TEfficiency* nprcorrHisto = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s_npr",collName)));
	double prEff = 1.0;
	double prAcc = 1.0;
	double nprEff = 1.0;
	double bf = 1.0;
	if (!prEffHisto || !prAccHisto) {
		std::cout << "[Error] acc or eff histogram not provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
		return 1.;
	}

	if (pt > 3 && pt < 100 && fabs(rapidity) < 2.4) {
		prAcc = prAccHisto->GetEfficiency(prAccHisto->FindFixBin(rapidity, pt));
		prEff = prEffHisto->GetEfficiency(prEffHisto->FindFixBin(rapidity, pt));
		corr = prEff * prAcc;
		if (!(corr < 1 && corr > 0) && pt > 6.5) {
			cout << "pt = " << pt << ", rapidity = " << rapidity << endl;
			cout << "prEff = " << prEff << endl;
			cout << "nprEff = " << nprEff << endl;
			cout << "bf = " << bf << endl;
			cout << "corr = " << corr << ", 1/corr = " << 1.0 / corr << endl;
		}
	}

	if (!corr || corr < 0.00001) corr = 1.0;
	return corr;
};

/*
double getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP, int cent)
{
  const char* collName = "PbPb";
  const char* massName = "Jpsi";
  if (isPP)
    collName = "PP";

    if (!fcorrArray)
    {
    cout << "[ERROR] No correction array exist" << endl;
    return 0;
    }

  Double_t corr = 1.;

    TEfficiency* prEffHisto = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s_pr_Eff",collName)));
    TEfficiency* prAccHisto = static_cast<TEfficiency*>(fcorrArray->FindObject("hcorr_Jpsi_PP_pr_Acc"));
    double prEff = 1.0;
    double prAcc = 1.0;
    if (!prEffHisto || !prAccHisto)
      {
	std::cout << "[Error] pr or histogram not provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
	return 1.;
      }

    if (pt > 3 && pt < 50 && fabs(rapidity) < 2.4){
      prEff = prEffHisto->GetEfficiency(prEffHisto->FindFixBin(rapidity, pt));
      prAcc = prAccHisto->GetEfficiency(prAccHisto->FindFixBin(rapidity, pt));
      corr = prEff*prAcc;
      if (!(corr<1 && corr>0) && pt > 6.5) {
	cout<< "pt = "<<pt<<", rapidity = "<<rapidity<<endl;
	cout<<"corr = "<< corr <<", 1/corr = "<<1.0/corr<<endl;
	}
    }

  if(!corr || corr<0.00001) corr=1.0;
  return corr;
};
*/

float jecCorr(double jtPt, double rawPt, double jpsiPt) {
	return ((1 - (jpsiPt / rawPt)) * jtPt + ((jpsiPt / rawPt) / (jtPt / rawPt)) * jtPt);
	//return ( (1-z)*jtPt + z*rawPt );
}

float zjecCorr(double jtPt, double rawPt, double z) {
	//return ( (1-(jpsiPt/rawPt))*jtPt + ((jpsiPt/rawPt)/(jtPt/rawPt))*jtPt );
	return ((1 - z) * jtPt + z * rawPt);
}

float jeuCorr(double jtPt, double z, double jeu) {
	return ((1 - z) * (1 + jeu) * jtPt + z * jtPt);
}

bool jetsInHCALHole(double jtEta, double jtPhi) {
	if ((jtEta > -3.0) && (jtEta < -1.392) && (jtPhi > -1.57) && (jtPhi < -0.87)) return true;
	return false;
}

int getCentFromName(TString sName) {
	if (sName.Contains("cent0")) return 0;
	for (int i = 2; i < 200; i++)
		if (sName.Contains(Form("cent%d", i)) && !sName.Contains(Form("cent%d", i * 10))) return i;
	return 200;
}
