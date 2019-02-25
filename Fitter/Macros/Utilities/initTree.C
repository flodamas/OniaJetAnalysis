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

Int_t           fCurrent; //!current Tree number in a TChain


//////////////////////////////////////////
// temporary solution for gen leaf type //
//////////////////////////////////////////
Int_t           Gen_QQ_size;
Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_4mom;
Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
Int_t           Gen_mu_size;
Int_t           Gen_mu_type[99];   //[Gen_mu_size]
Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
TClonesArray    *Gen_mu_4mom;
Int_t           Reco_QQ_whichGen[99];

Float_t         pthat;
Float_t         refpt[56];   //[nref]
Float_t         refeta[56];   //[nref]
Float_t         refy[56];   //[nref]
Float_t         refphi[56];   //[nref]
Float_t         refm[56];   //[nref]
Float_t         refarea[56];   //[nref]
Float_t         reftau1[56];   //[nref]
Float_t         reftau2[56];   //[nref]
Float_t         reftau3[56];   //[nref]
Float_t         refdphijt[56];   //[nref]
Float_t         refdrjt[56];   //[nref]
Float_t         refparton_pt[56];   //[nref]
Int_t           refparton_flavor[56];   //[nref]
Int_t           refparton_flavorForB[56];   //[nref]
Float_t         genChargedSum[56];   //[nref]
Float_t         genHardSum[56];   //[nref]
Float_t         signalChargedSum[56];   //[nref]
Float_t         signalHardSum[56];   //[nref]
Int_t           subid[56];   //[nref]
Int_t           ngen;
Int_t           genmatchindex[28];   //[ngen]
Float_t         genpt[28];   //[ngen]
Float_t         geneta[28];   //[ngen]
Float_t         geny[28];   //[ngen]
Float_t         gentau1[28];   //[ngen]
Float_t         gentau2[28];   //[ngen]
Float_t         gentau3[28];   //[ngen]
Float_t         genphi[28];   //[ngen]
Float_t         genm[28];   //[ngen]
Float_t         gendphijt[28];   //[ngen]
Float_t         gendrjt[28];   //[ngen]
Int_t           gensubid[28];   //[ngen]

//////////////////////////////////////////

// Declaration of leaf types
UInt_t          eventNb;
UInt_t          runNb;
UInt_t          LS;
Float_t         zVtx;
Float_t         nPV;
Int_t           Centrality;
Int_t           Npix;
Int_t           NpixelTracks;
Int_t           Ntracks;
Int_t           nTrig;
Int_t           trigPrescale[30];
ULong64_t       HLTriggers;
Int_t           Reco_QQ_size;
Int_t           Reco_QQ_type[30];
Int_t           Reco_QQ_sign[30];
TClonesArray    *Reco_QQ_4mom;
Int_t           Reco_QQ_mupl_idx[30];
Int_t           Reco_QQ_mumi_idx[30];
ULong64_t       Reco_QQ_trig[30];
Bool_t          Reco_QQ_isCowboy[30];
Float_t         Reco_QQ_ctau[30];
Float_t         Reco_QQ_ctauErr[30];
Float_t         Reco_QQ_cosAlpha[30];
Float_t         Reco_QQ_ctau3D[30];
Float_t         Reco_QQ_ctauErr3D[30];
Float_t         Reco_QQ_cosAlpha3D[30];
Float_t         Reco_QQ_VtxProb[30];
Float_t         Reco_QQ_dca[30];
Float_t         Reco_QQ_MassErr[30];
TClonesArray    *Reco_QQ_vtx;
Int_t           Reco_QQ_Ntrk[30];
Int_t           Reco_mu_size;
Int_t           Reco_mu_type[20];
Int_t           Reco_mu_SelectionType[20];
Int_t           Reco_mu_charge[20];
TClonesArray    *Reco_mu_4mom;
ULong64_t       Reco_mu_trig[20];
Bool_t          Reco_mu_highPurity[20];
Bool_t          Reco_mu_TrkMuArb[20];
Bool_t          Reco_mu_TMOneStaTight[20];
Int_t           Reco_mu_nPixValHits[20];
Int_t           Reco_mu_nMuValHits[20];
Int_t           Reco_mu_nTrkHits[20];
Float_t         Reco_mu_normChi2_inner[20];
Float_t         Reco_mu_normChi2_global[20];
Int_t           Reco_mu_nPixWMea[20];
Int_t           Reco_mu_nTrkWMea[20];
Int_t           Reco_mu_StationsMatched[20];
Float_t         Reco_mu_dxy[20];
Float_t         Reco_mu_dxyErr[20];
Float_t         Reco_mu_dz[20];
Float_t         Reco_mu_dzErr[20];
Float_t         Reco_mu_pt_inner[20];
Float_t         Reco_mu_pt_global[20];
Float_t         Reco_mu_ptErr_inner[20];
Float_t         Reco_mu_ptErr_global[20];

//Int_t           evt;
Float_t         b;
Int_t           nref;
Float_t         rawpt[10];   //[nref]
Float_t         jtpt[10];   //[nref]
Float_t         jteta[10];   //[nref]
Float_t         jty[10];   //[nref]
Float_t         jtphi[10];   //[nref]
Float_t         jtpu[10];   //[nref]
Float_t         jtm[10];   //[nref]
Float_t         jtarea[10];   //[nref]
Float_t         jtPfCHF[10];   //[nref]
Float_t         jtPfNHF[10];   //[nref]
Float_t         jtPfCEF[10];   //[nref]
Float_t         jtPfNEF[10];   //[nref]
Float_t         jtPfMUF[10];   //[nref]
Int_t           jtPfCHM[10];   //[nref]
Int_t           jtPfNHM[10];   //[nref]
Int_t           jtPfCEM[10];   //[nref]
Int_t           jtPfNEM[10];   //[nref]
Int_t           jtPfMUM[10];   //[nref]
Float_t         jttau1[10];   //[nref]
Float_t         jttau2[10];   //[nref]
Float_t         jttau3[10];   //[nref]
Float_t         discr_jetID_cuts[10];   //[nref]
Float_t         discr_jetID_bdt[10];   //[nref]
Float_t         discr_fr01[10];   //[nref]
Float_t         trackMax[10];   //[nref]
Float_t         trackSum[10];   //[nref]
Int_t           trackN[10];   //[nref]
Float_t         trackHardSum[10];   //[nref]
Int_t           trackHardN[10];   //[nref]
Float_t         chargedMax[10];   //[nref]
Float_t         chargedSum[10];   //[nref]
Int_t           chargedN[10];   //[nref]
Float_t         chargedHardSum[10];   //[nref]
Int_t           chargedHardN[10];   //[nref]
Float_t         photonMax[10];   //[nref]
Float_t         photonSum[10];   //[nref]
Int_t           photonN[10];   //[nref]
Float_t         photonHardSum[10];   //[nref]
Int_t           photonHardN[10];   //[nref]
Float_t         neutralMax[10];   //[nref]
Float_t         neutralSum[10];   //[nref]
Int_t           neutralN[10];   //[nref]
Float_t         hcalSum[10];   //[nref]
Float_t         ecalSum[10];   //[nref]
Float_t         eMax[10];   //[nref]
Float_t         eSum[10];   //[nref]
Int_t           eN[10];   //[nref]
Float_t         muMax[10];   //[nref]
Float_t         muSum[10];   //[nref]
Int_t           muN[10];   //[nref]

Int_t           Onia2MuMuPAT;
Int_t           ana_step;
Int_t           pclusterCompatibilityFilter;
Int_t           pprimaryVertexFilter; // for PbPb
Int_t           pPAprimaryVertexFilter; // for pp
Int_t           pBeamScrapingFilter;
Int_t           collisionEventSelectionAOD;
Int_t           collisionEventSelectionAODv2;
Int_t           phfCoincFilter1Th3;
Int_t           phfCoincFilter2Th3;
Int_t           phfCoincFilter3Th3;
Int_t           phfCoincFilter4Th3;
Int_t           phfCoincFilter5Th3;
Int_t           phfCoincFilter1Th4;
Int_t           phfCoincFilter2Th4;
Int_t           phfCoincFilter3Th4;
Int_t           phfCoincFilter4Th4;
Int_t           phfCoincFilter5Th4;
Int_t           phfCoincFilter1Th5;
Int_t           phfCoincFilter4Th2;
Int_t           pVertexFilterCutG;
Int_t           pVertexFilterCutGloose;
Int_t           pVertexFilterCutGtight;
Int_t           pVertexFilterCutGplus;
Int_t           pVertexFilterCutE;
Int_t           pVertexFilterCutEandG;
Int_t           pHBHENoiseFilterResultProducer;
Int_t           HBHENoiseFilterResult;
Int_t           HBHENoiseFilterResultRun1;
Int_t           HBHENoiseFilterResultRun2Loose;
Int_t           HBHENoiseFilterResultRun2Tight;
Int_t           HBHEIsoNoiseFilterResult;
Int_t           superFilterPath;

UInt_t          run;
ULong64_t       evt;
//UInt_t          lumi;
Float_t         vx;
Float_t         vy;
Float_t         vz;
Int_t           hiBin;
Float_t         hiHF;
Int_t           hiNevtPlane;
Float_t         hiEvtPlanes[5];   //[hiNevtPlane]


// List of branches
TBranch        *b_eventNb;   //!
TBranch        *b_runNb;   //!
TBranch        *b_LS;   //!
TBranch        *b_zVtx;   //!
TBranch        *b_nPV;   //!
TBranch        *b_Centrality;   //!
TBranch        *b_Npix;   //!
TBranch        *b_NpixelTracks;   //!
TBranch        *b_Ntracks;   //!
TBranch        *b_nTrig;   //!
TBranch        *b_trigPrescale;   //!
TBranch        *b_HLTriggers;   //!
TBranch        *b_Reco_QQ_size;   //!
TBranch        *b_Reco_QQ_type;   //!
TBranch        *b_Reco_QQ_sign;   //!
TBranch        *b_Reco_QQ_4mom;   //!
TBranch        *b_Reco_QQ_mupl_idx;   //!
TBranch        *b_Reco_QQ_mumi_idx;   //!
TBranch        *b_Reco_QQ_trig;   //!
TBranch        *b_Reco_QQ_isCowboy;   //!
TBranch        *b_Reco_QQ_ctau;   //!
TBranch        *b_Reco_QQ_ctauErr;   //!
TBranch        *b_Reco_QQ_cosAlpha;   //!
TBranch        *b_Reco_QQ_ctau3D;   //!
TBranch        *b_Reco_QQ_ctauErr3D;   //!
TBranch        *b_Reco_QQ_cosAlpha3D;   //!
TBranch        *b_Reco_QQ_VtxProb;   //!
TBranch        *b_Reco_QQ_dca;   //!
TBranch        *b_Reco_QQ_MassErr;   //!
TBranch        *b_Reco_QQ_vtx;   //!
TBranch        *b_Reco_QQ_Ntrk;   //!
TBranch        *b_Reco_mu_size;   //!
TBranch        *b_Reco_mu_type;   //!
TBranch        *b_Reco_mu_SelectionType;   //!
TBranch        *b_Reco_mu_charge;   //!
TBranch        *b_Reco_mu_4mom;   //!
TBranch        *b_Reco_mu_trig;   //!
TBranch        *b_Reco_mu_highPurity;   //!
TBranch        *b_Reco_mu_TrkMuArb;   //!
TBranch        *b_Reco_mu_TMOneStaTight;   //!
TBranch        *b_Reco_mu_nPixValHits;   //!
TBranch        *b_Reco_mu_nMuValHits;   //!
TBranch        *b_Reco_mu_nTrkHits;   //!
TBranch        *b_Reco_mu_normChi2_inner;   //!
TBranch        *b_Reco_mu_normChi2_global;   //!
TBranch        *b_Reco_mu_nPixWMea;   //!
TBranch        *b_Reco_mu_nTrkWMea;   //!
TBranch        *b_Reco_mu_StationsMatched;   //!
TBranch        *b_Reco_mu_dxy;   //!
TBranch        *b_Reco_mu_dxyErr;   //!
TBranch        *b_Reco_mu_dz;   //!
TBranch        *b_Reco_mu_dzErr;   //!
TBranch        *b_Reco_mu_pt_inner;   //!
TBranch        *b_Reco_mu_pt_global;   //!
TBranch        *b_Reco_mu_ptErr_inner;   //!
TBranch        *b_Reco_mu_ptErr_global;   //!

//TBranch        *b_evt;   //!
TBranch        *b_b;   //!
TBranch        *b_nref;   //!
TBranch        *b_rawpt;   //!
TBranch        *b_jtpt;   //!
TBranch        *b_jteta;   //!
TBranch        *b_jty;   //!
TBranch        *b_jtphi;   //!
TBranch        *b_jtpu;   //!
TBranch        *b_jtm;   //!
TBranch        *b_jtarea;   //!
TBranch        *b_jtPfCHF;   //!
TBranch        *b_jtPfNHF;   //!
TBranch        *b_jtPfCEF;   //!
TBranch        *b_jtPfNEF;   //!
TBranch        *b_jtPfMUF;   //!
TBranch        *b_jtPfCHM;   //!
TBranch        *b_jtPfNHM;   //!
TBranch        *b_jtPfCEM;   //!
TBranch        *b_jtPfNEM;   //!
TBranch        *b_jtPfMUM;   //!
TBranch        *b_jttau1;   //!
TBranch        *b_jttau2;   //!
TBranch        *b_jttau3;   //!
TBranch        *b_discr_jetID_cuts;   //!
TBranch        *b_discr_jetID_bdt;   //!
TBranch        *b_discr_fr01;   //!
TBranch        *b_trackMax;   //!
TBranch        *b_trackSum;   //!
TBranch        *b_trackN;   //!
TBranch        *b_trackHardSum;   //!
TBranch        *b_trackHardN;   //!
TBranch        *b_chargedMax;   //!
TBranch        *b_chargedSum;   //!
TBranch        *b_chargedN;   //!
TBranch        *b_chargedHardSum;   //!
TBranch        *b_chargedHardN;   //!
TBranch        *b_photonMax;   //!
TBranch        *b_photonSum;   //!
TBranch        *b_photonN;   //!
TBranch        *b_photonHardSum;   //!
TBranch        *b_photonHardN;   //!
TBranch        *b_neutralMax;   //!
TBranch        *b_neutralSum;   //!
TBranch        *b_neutralN;   //!
TBranch        *b_hcalSum;   //!
TBranch        *b_ecalSum;   //!
TBranch        *b_eMax;   //!
TBranch        *b_eSum;   //!
TBranch        *b_eN;   //!
TBranch        *b_muMax;   //!
TBranch        *b_muSum;   //!
TBranch        *b_muN;   //!

TBranch        *b_Onia2MuMuPAT;   //!
TBranch        *b_ana_step;   //!
TBranch        *b_pclusterCompatibilityFilter;   //!
TBranch        *b_pprimaryVertexFilter;   //!
TBranch        *b_pPAprimaryVertexFilter;   //!
TBranch        *b_pBeamScrapingFilter;   //!
TBranch        *b_collisionEventSelectionAOD;   //!
TBranch        *b_collisionEventSelectionAODv2;   //!
TBranch        *b_phfCoincFilter1Th3;   //!
TBranch        *b_phfCoincFilter2Th3;   //!
TBranch        *b_phfCoincFilter3Th3;   //!
TBranch        *b_phfCoincFilter4Th3;   //!
TBranch        *b_phfCoincFilter5Th3;   //!
TBranch        *b_phfCoincFilter1Th4;   //!
TBranch        *b_phfCoincFilter2Th4;   //!
TBranch        *b_phfCoincFilter3Th4;   //!
TBranch        *b_phfCoincFilter4Th4;   //!
TBranch        *b_phfCoincFilter5Th4;   //!
TBranch        *b_phfCoincFilter1Th5;   //!
TBranch        *b_phfCoincFilter4Th2;   //!
TBranch        *b_pVertexFilterCutG;   //!
TBranch        *b_pVertexFilterCutGloose;   //!
TBranch        *b_pVertexFilterCutGtight;   //!
TBranch        *b_pVertexFilterCutGplus;   //!
TBranch        *b_pVertexFilterCutE;   //!
TBranch        *b_pVertexFilterCutEandG;   //!
TBranch        *b_pHBHENoiseFilterResultProducer;   //!
TBranch        *b_HBHENoiseFilterResult;   //!
TBranch        *b_HBHENoiseFilterResultRun1;   //!
TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
TBranch        *b_HBHEIsoNoiseFilterResult;   //!
TBranch        *b_superFilterPath;   //!

TBranch        *b_run;   //!
TBranch        *b_evt;   //!
//TBranch        *b_lumi;   //!
TBranch        *b_vx;   //!
TBranch        *b_vy;   //!
TBranch        *b_vz;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_hiHF;   //!
TBranch        *b_hiNevtPlane;   //!
TBranch        *b_hiEvtPlanes;   //!


string TreeName("hionia/myTree");
string jetTreeName("ak4PFJetAnalyzer/t");
string skimTreeName("skimanalysis/HltTree");
string centTreeName("hiEvtAnalyzer/HiTree");

TTree *htr;
TTree *jtr;
TTree *str;
TTree *ctr;

void initTree(TChain *tree)
{
  std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed). 

  TChain   *fChain;   //!pointer to the analyzed TTree or TChain

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

  if (fChain->GetBranch("eventNb")) fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  if (fChain->GetBranch("runNb")) fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
  if (fChain->GetBranch("LS")) fChain->SetBranchAddress("LS", &LS, &b_LS);
  if (fChain->GetBranch("zVtx")) fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  if (fChain->GetBranch("nPV")) fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
  if (fChain->GetBranch("Centrality")) fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  if (fChain->GetBranch("Npix")) fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
  if (fChain->GetBranch("NpixelTracks")) fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
  if (fChain->GetBranch("Ntracks")) fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
  if (fChain->GetBranch("nTrig")) fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
  if (fChain->GetBranch("trigPrescale")) fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
  if (fChain->GetBranch("HLTriggers")) fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  if (fChain->GetBranch("Reco_QQ_type")) fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  if (fChain->GetBranch("Reco_QQ_mupl_idx")) fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
  if (fChain->GetBranch("Reco_QQ_mumi_idx")) fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
  if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  if (fChain->GetBranch("Reco_QQ_isCowboy")) fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
  if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  if (fChain->GetBranch("Reco_QQ_cosAlpha")) fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
  if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  if (fChain->GetBranch("Reco_QQ_cosAlpha3D")) fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
  if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
  if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
  if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
  if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
  if (fChain->GetBranch("Reco_mu_SelectionType")) fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
  if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
  if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
  if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
  if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
  
  //if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
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
  if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
  //if (fChain->GetBranch("lumi")) fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
  if (fChain->GetBranch("vx")) fChain->SetBranchAddress("vx", &vx, &b_vx);
  if (fChain->GetBranch("vy")) fChain->SetBranchAddress("vy", &vy, &b_vy);
  if (fChain->GetBranch("vz")) fChain->SetBranchAddress("vz", &vz, &b_vz);
  if (fChain->GetBranch("hiBin")) fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
  if (fChain->GetBranch("hiHF")) fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
  if (fChain->GetBranch("hiNevtPlane")) fChain->SetBranchAddress("hiNevtPlane", &hiNevtPlane, &b_hiNevtPlane);
  if (fChain->GetBranch("hiEvtPlanes")) fChain->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes, &b_hiEvtPlanes);

}
#endif // #ifndef initTree_C
