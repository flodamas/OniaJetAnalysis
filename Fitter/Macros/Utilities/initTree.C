//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 16 22:46:58 2016 by ROOT version 5.34/32
// from TTree myTree/My TTree of dimuons
// found on file: OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262620_263511_noCUT.root
//////////////////////////////////////////////////////////

#ifndef initTree_C
#define initTree_C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

Int_t fCurrent; //!current Tree number in a TChain

// Declaration of leaf types
Float_t zVtx;
Int_t trigPrescale[26];
ULong64_t HLTriggers;
Short_t Reco_QQ_size;
Short_t Reco_QQ_type[99]; //[Reco_QQ_size]
Short_t Reco_QQ_sign[99]; //[Reco_QQ_size]
TClonesArray* Reco_QQ_4mom;
Short_t Reco_QQ_mupl_idx[99];  //[Reco_QQ_size]
Short_t Reco_QQ_mumi_idx[99];  //[Reco_QQ_size]
ULong64_t Reco_QQ_trig[99];    //[Reco_QQ_size]
Float_t Reco_QQ_ctau[99];      //[Reco_QQ_size]
Float_t Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
Float_t Reco_QQ_ctau3D[99];    //[Reco_QQ_size]
Float_t Reco_QQ_ctauErr3D[99]; //[Reco_QQ_size]
Short_t Reco_QQ_whichGen[99];  //[Reco_QQ_size]
Float_t Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
Float_t Reco_QQ_dca[99];       //[Reco_QQ_size]
Short_t Reco_mu_size;
Short_t Reco_mu_type[99];        //[Reco_mu_size]
Short_t Reco_mu_whichGen[99];    //[Reco_mu_size]
Int_t Reco_mu_SelectionType[99]; //[Reco_mu_size]
Short_t Reco_mu_charge[99];      //[Reco_mu_size]
TClonesArray* Reco_mu_4mom;
Bool_t Reco_mu_highPurity[99];       //[Reco_mu_size]
Short_t Reco_mu_candType[99];        //[Reco_mu_size]
Int_t Reco_mu_nPixValHits[99];       //[Reco_mu_size]
Int_t Reco_mu_nMuValHits[99];        //[Reco_mu_size]
Int_t Reco_mu_nTrkHits[99];          //[Reco_mu_size]
Float_t Reco_mu_normChi2_inner[99];  //[Reco_mu_size]
Float_t Reco_mu_normChi2_global[99]; //[Reco_mu_size]
Int_t Reco_mu_nPixWMea[99];          //[Reco_mu_size]
Int_t Reco_mu_nTrkWMea[99];          //[Reco_mu_size]
Float_t Reco_mu_dxy[99];             //[Reco_mu_size]
Float_t Reco_mu_dz[99];              //[Reco_mu_size]
Float_t Gen_weight;
Float_t Gen_pthat;
Short_t Gen_QQ_size;
TClonesArray* Gen_QQ_4mom;
Float_t Gen_QQ_ctau[3];     //[Gen_QQ_size]
Float_t Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
Short_t Gen_QQ_mupl_idx[3]; //[Gen_QQ_size]
Short_t Gen_QQ_mumi_idx[3]; //[Gen_QQ_size]
Short_t Gen_QQ_whichRec[3]; //[Gen_QQ_size]
Short_t Gen_mu_size;
Short_t Gen_mu_charge[9]; //[Gen_mu_size]
TClonesArray* Gen_mu_4mom;
Short_t Gen_mu_whichRec[9]; //[Gen_mu_size]

Int_t evt;
Float_t b;
Int_t nref;
Float_t rawpt[99];   //[nref]
Float_t jtpt[99];    //[nref]
Float_t jteta[99];   //[nref]
Float_t jty[99];     //[nref]
Float_t jtphi[99];   //[nref]
Float_t jtm[99];     //[nref]
Float_t jtarea[99];  //[nref]
Float_t jtPfCHF[99]; //[nref]
Float_t jtPfNHF[99]; //[nref]
Float_t jtPfCEF[99]; //[nref]
Float_t jtPfNEF[99]; //[nref]
Float_t jtPfMUF[99]; //[nref]
Int_t jtPfCHM[99];   //[nref]
Int_t jtPfNHM[99];   //[nref]
Int_t jtPfCEM[99];   //[nref]
Int_t jtPfNEM[99];   //[nref]
Int_t jtPfMUM[99];   //[nref]
Float_t jttau1[99];  //[nref]
Float_t jttau2[99];  //[nref]
Float_t jttau3[99];  //[nref]
Float_t pthat;
Float_t refpt[99];              //[nref]
Float_t refeta[99];             //[nref]
Float_t refy[99];               //[nref]
Float_t refphi[99];             //[nref]
Float_t refm[99];               //[nref]
Float_t refarea[99];            //[nref]
Float_t reftau1[99];            //[nref]
Float_t reftau2[99];            //[nref]
Float_t reftau3[99];            //[nref]
Float_t refdphijt[99];          //[nref]
Float_t refdrjt[99];            //[nref]
Float_t refparton_pt[99];       //[nref]
Int_t refparton_flavor[99];     //[nref]
Int_t refparton_flavorForB[99]; //[nref]
Float_t genChargedSum[99];      //[nref]
Float_t genHardSum[99];         //[nref]
Float_t signalChargedSum[99];   //[nref]
Float_t signalHardSum[99];      //[nref]
Int_t subid[99];                //[nref]
Int_t ngen;
Int_t genmatchindex[99]; //[ngen]
Float_t genpt[99];       //[ngen]
Float_t geneta[99];      //[ngen]
Float_t geny[99];        //[ngen]
Float_t gentau1[99];     //[ngen]
Float_t gentau2[99];     //[ngen]
Float_t gentau3[99];     //[ngen]
Float_t genphi[99];      //[ngen]
Float_t genm[99];        //[ngen]
Float_t gendphijt[99];   //[ngen]
Float_t gendrjt[99];     //[ngen]
Int_t gensubid[99];      //[ngen]

Int_t pclusterCompatibilityFilter;
Int_t pprimaryVertexFilter;   // for PbPb
Int_t pPAprimaryVertexFilter; // for pp
Int_t pBeamScrapingFilter;

// List of branches
TBranch* b_zVtx;                    //!
TBranch* b_trigPrescale;            //!
TBranch* b_HLTriggers;              //!
TBranch* b_Reco_QQ_size;            //!
TBranch* b_Reco_QQ_type;            //!
TBranch* b_Reco_QQ_sign;            //!
TBranch* b_Reco_QQ_4mom;            //!
TBranch* b_Reco_QQ_mupl_idx;        //!
TBranch* b_Reco_QQ_mumi_idx;        //!
TBranch* b_Reco_QQ_trig;            //!
TBranch* b_Reco_QQ_ctau;            //!
TBranch* b_Reco_QQ_ctauErr;         //!
TBranch* b_Reco_QQ_ctau3D;          //!
TBranch* b_Reco_QQ_ctauErr3D;       //!
TBranch* b_Reco_QQ_whichGen;        //!
TBranch* b_Reco_QQ_VtxProb;         //!
TBranch* b_Reco_QQ_dca;             //!
TBranch* b_Reco_QQ_vtx;             //!
TBranch* b_Reco_mu_size;            //!
TBranch* b_Reco_mu_type;            //!
TBranch* b_Reco_mu_whichGen;        //!
TBranch* b_Reco_mu_SelectionType;   //!
TBranch* b_Reco_mu_charge;          //!
TBranch* b_Reco_mu_4mom;            //!
TBranch* b_Reco_mu_trig;            //!
TBranch* b_Reco_mu_highPurity;      //!
TBranch* b_Reco_mu_candType;        //!
TBranch* b_Reco_mu_nPixValHits;     //!
TBranch* b_Reco_mu_nMuValHits;      //!
TBranch* b_Reco_mu_nTrkHits;        //!
TBranch* b_Reco_mu_normChi2_inner;  //!
TBranch* b_Reco_mu_normChi2_global; //!
TBranch* b_Reco_mu_nPixWMea;        //!
TBranch* b_Reco_mu_nTrkWMea;        //!
TBranch* b_Reco_mu_dxy;             //!
TBranch* b_Reco_mu_dz;              //!
TBranch* b_Gen_weight;              //!
TBranch* b_Gen_pthat;               //!
TBranch* b_Gen_QQ_size;             //!
TBranch* b_Gen_QQ_4mom;             //!
TBranch* b_Gen_QQ_ctau;             //!
TBranch* b_Gen_QQ_ctau3D;           //!
TBranch* b_Gen_QQ_mupl_idx;         //!
TBranch* b_Gen_QQ_mumi_idx;         //!
TBranch* b_Gen_QQ_whichRec;         //!
TBranch* b_Gen_mu_size;             //!
TBranch* b_Gen_mu_charge;           //!
TBranch* b_Gen_mu_4mom;             //!
TBranch* b_Gen_mu_whichRec;         //!

TBranch* b_evt;                  //!
TBranch* b_b;                    //!
TBranch* b_nref;                 //!
TBranch* b_rawpt;                //!
TBranch* b_jtpt;                 //!
TBranch* b_jteta;                //!
TBranch* b_jty;                  //!
TBranch* b_jtphi;                //!
TBranch* b_jtm;                  //!
TBranch* b_jtarea;               //!
TBranch* b_jttau1;               //!
TBranch* b_jttau2;               //!
TBranch* b_jttau3;               //!
TBranch* b_pthat;                //!
TBranch* b_refpt;                //!
TBranch* b_refeta;               //!
TBranch* b_refy;                 //!
TBranch* b_refphi;               //!
TBranch* b_refm;                 //!
TBranch* b_refarea;              //!
TBranch* b_reftau1;              //!
TBranch* b_reftau2;              //!
TBranch* b_reftau3;              //!
TBranch* b_refdphijt;            //!
TBranch* b_refdrjt;              //!
TBranch* b_refparton_pt;         //!
TBranch* b_refparton_flavor;     //!
TBranch* b_refparton_flavorForB; //!
TBranch* b_genChargedSum;        //!
TBranch* b_genHardSum;           //!
TBranch* b_signalChargedSum;     //!
TBranch* b_signalHardSum;        //!
TBranch* b_subid;                //!
TBranch* b_ngen;                 //!
TBranch* b_genmatchindex;        //!
TBranch* b_genpt;                //!
TBranch* b_geneta;               //!
TBranch* b_geny;                 //!
TBranch* b_gentau1;              //!
TBranch* b_gentau2;              //!
TBranch* b_gentau3;              //!
TBranch* b_genphi;               //!
TBranch* b_genm;                 //!
TBranch* b_gendphijt;            //!
TBranch* b_gendrjt;              //!
TBranch* b_gensubid;             //!

TBranch* b_pclusterCompatibilityFilter; //!
TBranch* b_pprimaryVertexFilter;        //!
TBranch* b_pPAprimaryVertexFilter;      //!
TBranch* b_pBeamScrapingFilter;         //!

string TreeName("hionia/myTree");
string jetTreeName("ak4PFJetAnalyzer/t");
string skimTreeName("skimanalysis/HltTree");

TTree* htr;
TTree* jtr;
TTree* str;

void initTree(TChain* tree) {
	std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;

	TChain* fChain; //!pointer to the analyzed TTree or TChain

	// Set object pointer
	Reco_QQ_4mom = 0;
	Reco_mu_4mom = 0;
	Gen_QQ_4mom = 0;
	Gen_mu_4mom = 0;
	// Set branch addresses and branch pointers

	if (!tree) return;
	fChain = tree;
	fCurrent = -1;

	fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
	fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
	fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
	fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
	fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
	fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
	fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
	fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
	fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
	if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
	fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
	fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
	fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
	fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
	if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
	fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
	fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
	fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
	fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
	fChain->SetBranchAddress("Reco_mu_candType", Reco_mu_candType, &b_Reco_mu_candType);
	fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
	fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
	fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
	fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
	fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
	fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
	fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
	fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
	fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
	if (fChain->GetBranch("Gen_weight")) fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
	if (fChain->GetBranch("Gen_pthat")) fChain->SetBranchAddress("Gen_pthat", &Gen_pthat, &b_Gen_pthat);
	if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
	if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
	if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
	if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
	if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
	if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
	if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
	if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
	if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
	if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
	if (fChain->GetBranch("Gen_mu_whichRec")) fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);

	if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
	if (fChain->GetBranch("b")) fChain->SetBranchAddress("b", &b, &b_b);
	if (fChain->GetBranch("nref")) fChain->SetBranchAddress("nref", &nref, &b_nref);
	if (fChain->GetBranch("rawpt")) fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
	if (fChain->GetBranch("jtpt")) fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
	if (fChain->GetBranch("jteta")) fChain->SetBranchAddress("jteta", jteta, &b_jteta);
	if (fChain->GetBranch("jty")) fChain->SetBranchAddress("jty", jty, &b_jty);
	if (fChain->GetBranch("jtphi")) fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
	if (fChain->GetBranch("jtm")) fChain->SetBranchAddress("jtm", jtm, &b_jtm);
	if (fChain->GetBranch("jtarea")) fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
	if (fChain->GetBranch("jttau1")) fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
	if (fChain->GetBranch("jttau2")) fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
	if (fChain->GetBranch("jttau3")) fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
	if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
	if (fChain->GetBranch("refpt")) fChain->SetBranchAddress("refpt", refpt, &b_refpt);
	if (fChain->GetBranch("refeta")) fChain->SetBranchAddress("refeta", refeta, &b_refeta);
	if (fChain->GetBranch("refy")) fChain->SetBranchAddress("refy", refy, &b_refy);
	if (fChain->GetBranch("refphi")) fChain->SetBranchAddress("refphi", refphi, &b_refphi);
	if (fChain->GetBranch("refm")) fChain->SetBranchAddress("refm", refm, &b_refm);
	if (fChain->GetBranch("refarea")) fChain->SetBranchAddress("refarea", refarea, &b_refarea);
	if (fChain->GetBranch("reftau1")) fChain->SetBranchAddress("reftau1", reftau1, &b_reftau1);
	if (fChain->GetBranch("reftau2")) fChain->SetBranchAddress("reftau2", reftau2, &b_reftau2);
	if (fChain->GetBranch("reftau3")) fChain->SetBranchAddress("reftau3", reftau3, &b_reftau3);
	if (fChain->GetBranch("refdphijt")) fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
	if (fChain->GetBranch("refdrjt")) fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
	if (fChain->GetBranch("refparton_pt")) fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
	if (fChain->GetBranch("refparton_flavor")) fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
	if (fChain->GetBranch("refparton_flavorForB")) fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
	if (fChain->GetBranch("genChargedSum")) fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
	if (fChain->GetBranch("genHardSum")) fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
	if (fChain->GetBranch("signalChargedSum")) fChain->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
	if (fChain->GetBranch("signalHardSum")) fChain->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
	if (fChain->GetBranch("subid")) fChain->SetBranchAddress("subid", subid, &b_subid);
	if (fChain->GetBranch("ngen")) fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
	if (fChain->GetBranch("genmatchindex")) fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
	if (fChain->GetBranch("genpt")) fChain->SetBranchAddress("genpt", genpt, &b_genpt);
	if (fChain->GetBranch("geneta")) fChain->SetBranchAddress("geneta", geneta, &b_geneta);
	if (fChain->GetBranch("geny")) fChain->SetBranchAddress("geny", geny, &b_geny);
	if (fChain->GetBranch("gentau1")) fChain->SetBranchAddress("gentau1", gentau1, &b_gentau1);
	if (fChain->GetBranch("gentau2")) fChain->SetBranchAddress("gentau2", gentau2, &b_gentau2);
	if (fChain->GetBranch("gentau3")) fChain->SetBranchAddress("gentau3", gentau3, &b_gentau3);
	if (fChain->GetBranch("genphi")) fChain->SetBranchAddress("genphi", genphi, &b_genphi);
	if (fChain->GetBranch("genm")) fChain->SetBranchAddress("genm", genm, &b_genm);
	if (fChain->GetBranch("gendphijt")) fChain->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
	if (fChain->GetBranch("gendrjt")) fChain->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);
	if (fChain->GetBranch("gensubid")) fChain->SetBranchAddress("gensubid", gensubid, &b_gensubid);

	if (fChain->GetBranch("pclusterCompatibilityFilter")) fChain->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter, &b_pclusterCompatibilityFilter);
	if (fChain->GetBranch("pprimaryVertexFilter")) fChain->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
	if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
	if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
}
#endif // #ifndef initTree_C
