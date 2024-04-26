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
Short_t nPV;
Int_t trigPrescale[26];
ULong64_t HLTriggers;
Short_t Reco_QQ_size;
Short_t Reco_QQ_type[99]; //[Reco_QQ_size]
Short_t Reco_QQ_sign[99]; //[Reco_QQ_size]
TClonesArray* Reco_QQ_4mom;
Short_t Reco_QQ_mupl_idx[99];   //[Reco_QQ_size]
Short_t Reco_QQ_mumi_idx[99];   //[Reco_QQ_size]
ULong64_t Reco_QQ_trig[99];     //[Reco_QQ_size]
Bool_t Reco_QQ_isCowboy[99];    //[Reco_QQ_size]
Float_t Reco_QQ_ctau[99];       //[Reco_QQ_size]
Float_t Reco_QQ_ctauErr[99];    //[Reco_QQ_size]
Float_t Reco_QQ_cosAlpha[99];   //[Reco_QQ_size]
Float_t Reco_QQ_ctau3D[99];     //[Reco_QQ_size]
Float_t Reco_QQ_ctauErr3D[99];  //[Reco_QQ_size]
Float_t Reco_QQ_cosAlpha3D[99]; //[Reco_QQ_size]
Short_t Reco_QQ_whichGen[99];   //[Reco_QQ_size]
Float_t Reco_QQ_VtxProb[99];    //[Reco_QQ_size]
Float_t Reco_QQ_dca[99];        //[Reco_QQ_size]
Float_t Reco_QQ_MassErr[99];    //[Reco_QQ_size]
TClonesArray* Reco_QQ_vtx;
Short_t Reco_mu_size;
Short_t Reco_mu_type[99];        //[Reco_mu_size]
Short_t Reco_mu_whichGen[99];    //[Reco_mu_size]
Int_t Reco_mu_SelectionType[99]; //[Reco_mu_size]
Short_t Reco_mu_charge[99];      //[Reco_mu_size]
TClonesArray* Reco_mu_4mom;
ULong64_t Reco_mu_trig[99];          //[Reco_mu_size]
Bool_t Reco_mu_InTightAcc[99];       //[Reco_mu_size]
Bool_t Reco_mu_InLooseAcc[99];       //[Reco_mu_size]
Bool_t Reco_mu_highPurity[99];       //[Reco_mu_size]
Bool_t Reco_mu_isPF[99];             //[Reco_mu_size]
Short_t Reco_mu_candType[99];        //[Reco_mu_size]
Int_t Reco_mu_nPixValHits[99];       //[Reco_mu_size]
Int_t Reco_mu_nMuValHits[99];        //[Reco_mu_size]
Int_t Reco_mu_nTrkHits[99];          //[Reco_mu_size]
Float_t Reco_mu_normChi2_inner[99];  //[Reco_mu_size]
Float_t Reco_mu_normChi2_global[99]; //[Reco_mu_size]
Int_t Reco_mu_nPixWMea[99];          //[Reco_mu_size]
Int_t Reco_mu_nTrkWMea[99];          //[Reco_mu_size]
Int_t Reco_mu_StationsMatched[99];   //[Reco_mu_size]
Float_t Reco_mu_dxy[99];             //[Reco_mu_size]
Float_t Reco_mu_dxyErr[99];          //[Reco_mu_size]
Float_t Reco_mu_dz[99];              //[Reco_mu_size]
Float_t Reco_mu_dzErr[99];           //[Reco_mu_size]
Float_t Reco_mu_ptErr_inner[99];     //[Reco_mu_size]
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
Float_t rawpt[99];            //[nref]
Float_t jtpt[99];             //[nref]
Float_t jteta[99];            //[nref]
Float_t jty[99];              //[nref]
Float_t jtphi[99];            //[nref]
Float_t jtpu[99];             //[nref]
Float_t jtm[99];              //[nref]
Float_t jtarea[99];           //[nref]
Float_t jtPfCHF[99];          //[nref]
Float_t jtPfNHF[99];          //[nref]
Float_t jtPfCEF[99];          //[nref]
Float_t jtPfNEF[99];          //[nref]
Float_t jtPfMUF[99];          //[nref]
Int_t jtPfCHM[99];            //[nref]
Int_t jtPfNHM[99];            //[nref]
Int_t jtPfCEM[99];            //[nref]
Int_t jtPfNEM[99];            //[nref]
Int_t jtPfMUM[99];            //[nref]
Float_t jttau1[99];           //[nref]
Float_t jttau2[99];           //[nref]
Float_t jttau3[99];           //[nref]
Float_t discr_jetID_cuts[99]; //[nref]
Float_t discr_jetID_bdt[99];  //[nref]
Float_t discr_fr01[99];       //[nref]
Float_t trackMax[99];         //[nref]
Float_t trackSum[99];         //[nref]
Int_t trackN[99];             //[nref]
Float_t trackHardSum[99];     //[nref]
Int_t trackHardN[99];         //[nref]
Float_t chargedMax[99];       //[nref]
Float_t chargedSum[99];       //[nref]
Int_t chargedN[99];           //[nref]
Float_t chargedHardSum[99];   //[nref]
Int_t chargedHardN[99];       //[nref]
Float_t photonMax[99];        //[nref]
Float_t photonSum[99];        //[nref]
Int_t photonN[99];            //[nref]
Float_t photonHardSum[99];    //[nref]
Int_t photonHardN[99];        //[nref]
Float_t neutralMax[99];       //[nref]
Float_t neutralSum[99];       //[nref]
Int_t neutralN[99];           //[nref]
Float_t hcalSum[99];          //[nref]
Float_t ecalSum[99];          //[nref]
Float_t eMax[99];             //[nref]
Float_t eSum[99];             //[nref]
Int_t eN[99];                 //[nref]
Float_t muMax[99];            //[nref]
Float_t muSum[99];            //[nref]
Int_t muN[99];                //[nref]
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

Int_t Onia2MuMuPAT;
Int_t ana_step;
Int_t pclusterCompatibilityFilter;
Int_t pprimaryVertexFilter;   // for PbPb
Int_t pPAprimaryVertexFilter; // for pp
Int_t pBeamScrapingFilter;
Int_t collisionEventSelectionAOD;
Int_t collisionEventSelectionAODv2;
Int_t phfCoincFilter1Th3;
Int_t phfCoincFilter2Th3;
Int_t phfCoincFilter3Th3;
Int_t phfCoincFilter4Th3;
Int_t phfCoincFilter5Th3;
Int_t phfCoincFilter1Th4;
Int_t phfCoincFilter2Th4;
Int_t phfCoincFilter3Th4;
Int_t phfCoincFilter4Th4;
Int_t phfCoincFilter5Th4;
Int_t phfCoincFilter1Th5;
Int_t phfCoincFilter4Th2;
Int_t pVertexFilterCutG;
Int_t pVertexFilterCutGloose;
Int_t pVertexFilterCutGtight;
Int_t pVertexFilterCutGplus;
Int_t pVertexFilterCutE;
Int_t pVertexFilterCutEandG;
Int_t pHBHENoiseFilterResultProducer;
Int_t HBHENoiseFilterResult;
Int_t HBHENoiseFilterResultRun1;
Int_t HBHENoiseFilterResultRun2Loose;
Int_t HBHENoiseFilterResultRun2Tight;
Int_t HBHEIsoNoiseFilterResult;
Int_t superFilterPath;

Int_t run;
//ULong64_t       evt;
//UInt_t          lumi;
Float_t vx;
Float_t vy;
Float_t vz;
Float_t Npart;
Float_t Ncoll;
Float_t Nhard;
Float_t phi0;
//Float_t         b;
Int_t ProcessID;
//Float_t         pthat;
Float_t weight;
Float_t alphaQCD;
Float_t alphaQED;
Float_t qScale;
Int_t nMEPartons;
Int_t nMEPartonsFiltered;
//pair<int,int>   *pdfID;
Int_t pdfID_first;
Int_t pdfID_second;
//pair<float,float> *pdfX;
Float_t pdfX_first;
Float_t pdfX_second;
//pair<float,float> *pdfXpdf;
Float_t pdfXpdf_first;
Float_t pdfXpdf_second;
vector<float>* ttbar_w;
vector<int>* npus;
vector<float>* tnpus;
Int_t hiBin;
Float_t hiHF;
Int_t hiNevtPlane;
Float_t hiEvtPlanes[5]; //[hiNevtPlane]

// List of branches
TBranch* b_zVtx;                    //!
TBranch* b_nPV;                     //!
TBranch* b_trigPrescale;            //!
TBranch* b_HLTriggers;              //!
TBranch* b_Reco_QQ_size;            //!
TBranch* b_Reco_QQ_type;            //!
TBranch* b_Reco_QQ_sign;            //!
TBranch* b_Reco_QQ_4mom;            //!
TBranch* b_Reco_QQ_mupl_idx;        //!
TBranch* b_Reco_QQ_mumi_idx;        //!
TBranch* b_Reco_QQ_trig;            //!
TBranch* b_Reco_QQ_isCowboy;        //!
TBranch* b_Reco_QQ_ctau;            //!
TBranch* b_Reco_QQ_ctauErr;         //!
TBranch* b_Reco_QQ_cosAlpha;        //!
TBranch* b_Reco_QQ_ctau3D;          //!
TBranch* b_Reco_QQ_ctauErr3D;       //!
TBranch* b_Reco_QQ_cosAlpha3D;      //!
TBranch* b_Reco_QQ_whichGen;        //!
TBranch* b_Reco_QQ_VtxProb;         //!
TBranch* b_Reco_QQ_dca;             //!
TBranch* b_Reco_QQ_MassErr;         //!
TBranch* b_Reco_QQ_vtx;             //!
TBranch* b_Reco_mu_size;            //!
TBranch* b_Reco_mu_type;            //!
TBranch* b_Reco_mu_whichGen;        //!
TBranch* b_Reco_mu_SelectionType;   //!
TBranch* b_Reco_mu_charge;          //!
TBranch* b_Reco_mu_4mom;            //!
TBranch* b_Reco_mu_trig;            //!
TBranch* b_Reco_mu_InTightAcc;      //!
TBranch* b_Reco_mu_InLooseAcc;      //!
TBranch* b_Reco_mu_highPurity;      //!
TBranch* b_Reco_mu_isPF;            //!
TBranch* b_Reco_mu_candType;        //!
TBranch* b_Reco_mu_nPixValHits;     //!
TBranch* b_Reco_mu_nMuValHits;      //!
TBranch* b_Reco_mu_nTrkHits;        //!
TBranch* b_Reco_mu_normChi2_inner;  //!
TBranch* b_Reco_mu_normChi2_global; //!
TBranch* b_Reco_mu_nPixWMea;        //!
TBranch* b_Reco_mu_nTrkWMea;        //!
TBranch* b_Reco_mu_StationsMatched; //!
TBranch* b_Reco_mu_dxy;             //!
TBranch* b_Reco_mu_dxyErr;          //!
TBranch* b_Reco_mu_dz;              //!
TBranch* b_Reco_mu_dzErr;           //!
TBranch* b_Reco_mu_ptErr_inner;     //!
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
TBranch* b_jtpu;                 //!
TBranch* b_jtm;                  //!
TBranch* b_jtarea;               //!
TBranch* b_jtPfCHF;              //!
TBranch* b_jtPfNHF;              //!
TBranch* b_jtPfCEF;              //!
TBranch* b_jtPfNEF;              //!
TBranch* b_jtPfMUF;              //!
TBranch* b_jtPfCHM;              //!
TBranch* b_jtPfNHM;              //!
TBranch* b_jtPfCEM;              //!
TBranch* b_jtPfNEM;              //!
TBranch* b_jtPfMUM;              //!
TBranch* b_jttau1;               //!
TBranch* b_jttau2;               //!
TBranch* b_jttau3;               //!
TBranch* b_discr_jetID_cuts;     //!
TBranch* b_discr_jetID_bdt;      //!
TBranch* b_discr_fr01;           //!
TBranch* b_trackMax;             //!
TBranch* b_trackSum;             //!
TBranch* b_trackN;               //!
TBranch* b_trackHardSum;         //!
TBranch* b_trackHardN;           //!
TBranch* b_chargedMax;           //!
TBranch* b_chargedSum;           //!
TBranch* b_chargedN;             //!
TBranch* b_chargedHardSum;       //!
TBranch* b_chargedHardN;         //!
TBranch* b_photonMax;            //!
TBranch* b_photonSum;            //!
TBranch* b_photonN;              //!
TBranch* b_photonHardSum;        //!
TBranch* b_photonHardN;          //!
TBranch* b_neutralMax;           //!
TBranch* b_neutralSum;           //!
TBranch* b_neutralN;             //!
TBranch* b_hcalSum;              //!
TBranch* b_ecalSum;              //!
TBranch* b_eMax;                 //!
TBranch* b_eSum;                 //!
TBranch* b_eN;                   //!
TBranch* b_muMax;                //!
TBranch* b_muSum;                //!
TBranch* b_muN;                  //!
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

TBranch* b_Onia2MuMuPAT;                   //!
TBranch* b_ana_step;                       //!
TBranch* b_pclusterCompatibilityFilter;    //!
TBranch* b_pprimaryVertexFilter;           //!
TBranch* b_pPAprimaryVertexFilter;         //!
TBranch* b_pBeamScrapingFilter;            //!
TBranch* b_collisionEventSelectionAOD;     //!
TBranch* b_collisionEventSelectionAODv2;   //!
TBranch* b_phfCoincFilter1Th3;             //!
TBranch* b_phfCoincFilter2Th3;             //!
TBranch* b_phfCoincFilter3Th3;             //!
TBranch* b_phfCoincFilter4Th3;             //!
TBranch* b_phfCoincFilter5Th3;             //!
TBranch* b_phfCoincFilter1Th4;             //!
TBranch* b_phfCoincFilter2Th4;             //!
TBranch* b_phfCoincFilter3Th4;             //!
TBranch* b_phfCoincFilter4Th4;             //!
TBranch* b_phfCoincFilter5Th4;             //!
TBranch* b_phfCoincFilter1Th5;             //!
TBranch* b_phfCoincFilter4Th2;             //!
TBranch* b_pVertexFilterCutG;              //!
TBranch* b_pVertexFilterCutGloose;         //!
TBranch* b_pVertexFilterCutGtight;         //!
TBranch* b_pVertexFilterCutGplus;          //!
TBranch* b_pVertexFilterCutE;              //!
TBranch* b_pVertexFilterCutEandG;          //!
TBranch* b_pHBHENoiseFilterResultProducer; //!
TBranch* b_HBHENoiseFilterResult;          //!
TBranch* b_HBHENoiseFilterResultRun1;      //!
TBranch* b_HBHENoiseFilterResultRun2Loose; //!
TBranch* b_HBHENoiseFilterResultRun2Tight; //!
TBranch* b_HBHEIsoNoiseFilterResult;       //!
TBranch* b_superFilterPath;                //!

TBranch* b_run; //!
//TBranch        *b_evt;   //!
//TBranch        *b_lumi;   //!
TBranch* b_vx;    //!
TBranch* b_vy;    //!
TBranch* b_vz;    //!
TBranch* b_Npart; //!
TBranch* b_Ncoll; //!
TBranch* b_Nhard; //!
TBranch* b_NPhi0; //!
//TBranch        *b_b;   //!
TBranch* b_ProcessID; //!
//TBranch        *b_pthat;   //!
TBranch* b_weight;             //!
TBranch* b_alphaQCD;           //!
TBranch* b_alphaQED;           //!
TBranch* b_qScale;             //!
TBranch* b_nMEPartons;         //!
TBranch* b_nMEPartonsFiltered; //!
TBranch* b_pdfID_first;        //!
TBranch* b_pdfID_second;       //!
TBranch* b_pdfX_first;         //!
TBranch* b_pdfX_second;        //!
TBranch* b_pdfXpdf_first;      //!
TBranch* b_pdfXpdf_second;     //!
TBranch* b_ttbar_w;            //!
TBranch* b_npus;               //!
TBranch* b_tnpus;              //!
TBranch* b_hiBin;              //!
TBranch* b_hiHF;               //!
TBranch* b_hiNevtPlane;        //!
TBranch* b_hiEvtPlanes;        //!

string TreeName("hionia/myTree");
string jetTreeName("ak4PFJetAnalyzer/t");
string skimTreeName("skimanalysis/HltTree");
string centTreeName("hiEvtAnalyzer/HiTree");

TTree* htr;
TTree* jtr;
TTree* str;
TTree* ctr;

void initTree(TChain* tree) {
	std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;

	TChain* fChain; //!pointer to the analyzed TTree or TChain

	// Set object pointer
	Reco_QQ_4mom = 0;
	Reco_QQ_vtx = 0;
	Reco_mu_4mom = 0;
	Gen_QQ_4mom = 0;
	Gen_mu_4mom = 0;
	// Set branch addresses and branch pointers

	if (!tree) return;
	fChain = tree;
	fCurrent = -1;

	fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
	fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
	fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
	fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
	fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
	fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
	fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
	fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
	fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
	fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
	fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
	fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
	fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
	fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
	fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
	fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
	fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
	fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
	if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
	fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
	fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
	fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
	fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
	fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
	fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
	if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
	fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
	fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
	fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
	fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
	fChain->SetBranchAddress("Reco_mu_InTightAcc", Reco_mu_InTightAcc, &b_Reco_mu_InTightAcc);
	fChain->SetBranchAddress("Reco_mu_InLooseAcc", Reco_mu_InLooseAcc, &b_Reco_mu_InLooseAcc);
	fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
	fChain->SetBranchAddress("Reco_mu_isPF", Reco_mu_isPF, &b_Reco_mu_isPF);
	fChain->SetBranchAddress("Reco_mu_candType", Reco_mu_candType, &b_Reco_mu_candType);
	fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
	fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
	fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
	fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
	fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
	fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
	fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
	fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
	fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
	fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
	fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
	fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
	fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
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
	if (fChain->GetBranch("jtpu")) fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
	if (fChain->GetBranch("jtm")) fChain->SetBranchAddress("jtm", jtm, &b_jtm);
	if (fChain->GetBranch("jtarea")) fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
	if (fChain->GetBranch("jtPfCHF")) fChain->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
	if (fChain->GetBranch("jtPfNHF")) fChain->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
	if (fChain->GetBranch("jtPfCEF")) fChain->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
	if (fChain->GetBranch("jtPfNEF")) fChain->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
	if (fChain->GetBranch("jtPfMUF")) fChain->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
	if (fChain->GetBranch("jtPfCHM")) fChain->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
	if (fChain->GetBranch("jtPfNHM")) fChain->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
	if (fChain->GetBranch("jtPfCEM")) fChain->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
	if (fChain->GetBranch("jtPfNEM")) fChain->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
	if (fChain->GetBranch("jtPfMUM")) fChain->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
	if (fChain->GetBranch("jttau1")) fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
	if (fChain->GetBranch("jttau2")) fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
	if (fChain->GetBranch("jttau3")) fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
	if (fChain->GetBranch("discr_jetID_cuts")) fChain->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts, &b_discr_jetID_cuts);
	if (fChain->GetBranch("discr_jetID_bdt")) fChain->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt, &b_discr_jetID_bdt);
	if (fChain->GetBranch("discr_fr01")) fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
	if (fChain->GetBranch("trackMax")) fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
	if (fChain->GetBranch("trackSum")) fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
	if (fChain->GetBranch("trackN")) fChain->SetBranchAddress("trackN", trackN, &b_trackN);
	if (fChain->GetBranch("trackHardSum")) fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
	if (fChain->GetBranch("trackHardN")) fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
	if (fChain->GetBranch("chargedMax")) fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
	if (fChain->GetBranch("chargedSum")) fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
	if (fChain->GetBranch("chargedN")) fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
	if (fChain->GetBranch("chargedHardSum")) fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
	if (fChain->GetBranch("chargedHardN")) fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
	if (fChain->GetBranch("photonMax")) fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
	if (fChain->GetBranch("photonSum")) fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
	if (fChain->GetBranch("photonN")) fChain->SetBranchAddress("photonN", photonN, &b_photonN);
	if (fChain->GetBranch("photonHardSum")) fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
	if (fChain->GetBranch("photonHardN")) fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
	if (fChain->GetBranch("neutralMax")) fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
	if (fChain->GetBranch("neutralSum")) fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
	if (fChain->GetBranch("neutralN")) fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
	if (fChain->GetBranch("hcalSum")) fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
	if (fChain->GetBranch("ecalSum")) fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
	if (fChain->GetBranch("eMax")) fChain->SetBranchAddress("eMax", eMax, &b_eMax);
	if (fChain->GetBranch("eSum")) fChain->SetBranchAddress("eSum", eSum, &b_eSum);
	if (fChain->GetBranch("eN")) fChain->SetBranchAddress("eN", eN, &b_eN);
	if (fChain->GetBranch("muMax")) fChain->SetBranchAddress("muMax", muMax, &b_muMax);
	if (fChain->GetBranch("muSum")) fChain->SetBranchAddress("muSum", muSum, &b_muSum);
	if (fChain->GetBranch("muN")) fChain->SetBranchAddress("muN", muN, &b_muN);
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

	if (fChain->GetBranch("Onia2MuMuPAT")) fChain->SetBranchAddress("Onia2MuMuPAT", &Onia2MuMuPAT, &b_Onia2MuMuPAT);
	if (fChain->GetBranch("ana_step")) fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
	if (fChain->GetBranch("pclusterCompatibilityFilter")) fChain->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter, &b_pclusterCompatibilityFilter);
	if (fChain->GetBranch("pprimaryVertexFilter")) fChain->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
	if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
	if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
	if (fChain->GetBranch("collisionEventSelectionAOD")) fChain->SetBranchAddress("collisionEventSelectionAOD", &collisionEventSelectionAOD, &b_collisionEventSelectionAOD);
	if (fChain->GetBranch("collisionEventSelectionAODv2")) fChain->SetBranchAddress("collisionEventSelectionAODv2", &collisionEventSelectionAODv2, &b_collisionEventSelectionAODv2);
	if (fChain->GetBranch("phfCoincFilter1Th3")) fChain->SetBranchAddress("phfCoincFilter1Th3", &phfCoincFilter1Th3, &b_phfCoincFilter1Th3);
	if (fChain->GetBranch("phfCoincFilter2Th3")) fChain->SetBranchAddress("phfCoincFilter2Th3", &phfCoincFilter2Th3, &b_phfCoincFilter2Th3);
	if (fChain->GetBranch("phfCoincFilter3Th3")) fChain->SetBranchAddress("phfCoincFilter3Th3", &phfCoincFilter3Th3, &b_phfCoincFilter3Th3);
	if (fChain->GetBranch("phfCoincFilter4Th3")) fChain->SetBranchAddress("phfCoincFilter4Th3", &phfCoincFilter4Th3, &b_phfCoincFilter4Th3);
	if (fChain->GetBranch("phfCoincFilter5Th3")) fChain->SetBranchAddress("phfCoincFilter5Th3", &phfCoincFilter5Th3, &b_phfCoincFilter5Th3);
	if (fChain->GetBranch("phfCoincFilter1Th4")) fChain->SetBranchAddress("phfCoincFilter1Th4", &phfCoincFilter1Th4, &b_phfCoincFilter1Th4);
	if (fChain->GetBranch("phfCoincFilter2Th4")) fChain->SetBranchAddress("phfCoincFilter2Th4", &phfCoincFilter2Th4, &b_phfCoincFilter2Th4);
	if (fChain->GetBranch("phfCoincFilter3Th4")) fChain->SetBranchAddress("phfCoincFilter3Th4", &phfCoincFilter3Th4, &b_phfCoincFilter3Th4);
	if (fChain->GetBranch("phfCoincFilter4Th4")) fChain->SetBranchAddress("phfCoincFilter4Th4", &phfCoincFilter4Th4, &b_phfCoincFilter4Th4);
	if (fChain->GetBranch("phfCoincFilter5Th4")) fChain->SetBranchAddress("phfCoincFilter5Th4", &phfCoincFilter5Th4, &b_phfCoincFilter5Th4);
	if (fChain->GetBranch("phfCoincFilter1Th5")) fChain->SetBranchAddress("phfCoincFilter1Th5", &phfCoincFilter1Th5, &b_phfCoincFilter1Th5);
	if (fChain->GetBranch("phfCoincFilter4Th2")) fChain->SetBranchAddress("phfCoincFilter4Th2", &phfCoincFilter4Th2, &b_phfCoincFilter4Th2);
	if (fChain->GetBranch("pVertexFilterCutG")) fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
	if (fChain->GetBranch("pVertexFilterCutGloose")) fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
	if (fChain->GetBranch("pVertexFilterCutGtight")) fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
	if (fChain->GetBranch("pVertexFilterCutGplus")) fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
	if (fChain->GetBranch("pVertexFilterCutE")) fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
	if (fChain->GetBranch("pVertexFilterCutEandG")) fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);
	if (fChain->GetBranch("pHBHENoiseFilterResultProducer")) fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
	if (fChain->GetBranch("HBHENoiseFilterResult")) fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
	if (fChain->GetBranch("HBHENoiseFilterResultRun1")) fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
	if (fChain->GetBranch("HBHENoiseFilterResultRun2Loose")) fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
	if (fChain->GetBranch("HBHENoiseFilterResultRun2Tight")) fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
	if (fChain->GetBranch("HBHEIsoNoiseFilterResult")) fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
	if (fChain->GetBranch("superFilterPath")) fChain->SetBranchAddress("superFilterPath", &superFilterPath, &b_superFilterPath);

	if (fChain->GetBranch("run")) fChain->SetBranchAddress("run", &run, &b_run);
	//if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
	//if (fChain->GetBranch("lumi")) fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
	if (fChain->GetBranch("vx")) fChain->SetBranchAddress("vx", &vx, &b_vx);
	if (fChain->GetBranch("vy")) fChain->SetBranchAddress("vy", &vy, &b_vy);
	if (fChain->GetBranch("vz")) fChain->SetBranchAddress("vz", &vz, &b_vz);
	if (fChain->GetBranch("Npart")) fChain->SetBranchAddress("Npart", &Npart, &b_Npart);
	if (fChain->GetBranch("Ncoll")) fChain->SetBranchAddress("Ncoll", &Ncoll, &b_Ncoll);
	if (fChain->GetBranch("Nhard")) fChain->SetBranchAddress("Nhard", &Nhard, &b_Nhard);
	if (fChain->GetBranch("phi0")) fChain->SetBranchAddress("phi0", &phi0, &b_NPhi0);
	//if (fChain->GetBranch("b")) fChain->SetBranchAddress("b", &b, &b_b);
	if (fChain->GetBranch("ProcessID")) fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
	//if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
	if (fChain->GetBranch("weight")) fChain->SetBranchAddress("weight", &weight, &b_weight);
	if (fChain->GetBranch("alphaQCD")) fChain->SetBranchAddress("alphaQCD", &alphaQCD, &b_alphaQCD);
	if (fChain->GetBranch("alphaQED")) fChain->SetBranchAddress("alphaQED", &alphaQED, &b_alphaQED);
	if (fChain->GetBranch("qScale")) fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
	if (fChain->GetBranch("nMEPartons")) fChain->SetBranchAddress("nMEPartons", &nMEPartons, &b_nMEPartons);
	if (fChain->GetBranch("nMEPartonsFiltered")) fChain->SetBranchAddress("nMEPartonsFiltered", &nMEPartonsFiltered, &b_nMEPartonsFiltered);
	if (fChain->GetBranch("first")) fChain->SetBranchAddress("first", &pdfID_first, &b_pdfID_first);
	if (fChain->GetBranch("second")) fChain->SetBranchAddress("second", &pdfID_second, &b_pdfID_second);
	//   if (fChain->GetBranch("first")) fChain->SetBranchAddress("first", &first, &b_pdfX_first);
	//   if (fChain->GetBranch("second")) fChain->SetBranchAddress("second", &second, &b_pdfX_second);
	//   if (fChain->GetBranch("first")) fChain->SetBranchAddress("first", &first, &b_pdfXpdf_first);
	//   if (fChain->GetBranch("second")) fChain->SetBranchAddress("second", &second, &b_pdfXpdf_second);
	if (fChain->GetBranch("ttbar_w")) fChain->SetBranchAddress("ttbar_w", &ttbar_w, &b_ttbar_w);
	if (fChain->GetBranch("npus")) fChain->SetBranchAddress("npus", &npus, &b_npus);
	if (fChain->GetBranch("tnpus")) fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
	if (fChain->GetBranch("hiBin")) fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
	if (fChain->GetBranch("hiHF")) fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
	if (fChain->GetBranch("hiNevtPlane")) fChain->SetBranchAddress("hiNevtPlane", &hiNevtPlane, &b_hiNevtPlane);
	if (fChain->GetBranch("hiEvtPlanes")) fChain->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes, &b_hiEvtPlanes);
}
#endif // #ifndef initTree_C
