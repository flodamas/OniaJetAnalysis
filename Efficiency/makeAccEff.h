#ifndef makeAccEff_h
#define makeAccEFf_h

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
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include <fstream>
#include <TLegend.h>

//#include "tnp_weight.h"
#include "tnp_weight_lowptpp.h"
#include "tnp_weight_lowptPbPb.h"

using namespace std;
using namespace  RooFit;

class oniaTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
   bool          isPbPb;
   bool          isPr;
   bool          isAcc;

   Double_t *ybins = new Double_t[100];
   int nybins;
   Double_t *ptbins =  new Double_t[100]; 
   int nptbins;
   int *centbins = new int[100];
   int ncentbins;
   Double_t *zbins = new Double_t[100];
   int nzbins;
   Double_t *jtptbins = new Double_t[100];
   int njtptbins;

   int centCut=40;
   //my additional variables
   
   int cent;
   Float_t jpsi_m;
   Float_t jpsi_pt;
   Float_t jpsi_eta;
   Float_t jpsi_phi;
   Float_t jpsi_rap;
   Float_t dr;
   Float_t dphi;
   Float_t dphimin;
   Float_t deta;
   Float_t drmin;
   Float_t z=100;
   int triggerIndex_PP = 3; // HLT_HIL1DoubleMu0_v1
   int triggerIndex_PbPb = 12; // HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1

   Float_t weight;
   Float_t tnp_weight;

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
   Int_t           trigPrescale[30];   //[nTrig]
   ULong64_t       HLTriggers;
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   Int_t           Reco_QQ_mupl_idx[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_idx[99];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_cosAlpha3D[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_whichGen[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[99];   //[Reco_mu_size]
   Int_t           Reco_mu_whichGen[99];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[99];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[99];   //[Reco_mu_size]
   Bool_t          Reco_mu_TrkMuArb[99];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[99];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[99];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[99];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[99];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[99];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[99];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[99];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[99];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[99];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[99];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[99];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[99];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[99];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_inner[99];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_global[99];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[99];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_global[99];   //[Reco_mu_size]
   Float_t         Gen_weight;
   Float_t         Gen_pthat;
   Int_t           Gen_QQ_size;
   Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_4mom;
   Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
   Int_t           Gen_QQ_mupl_idx[99];   //[Gen_QQ_size]
   Int_t           Gen_QQ_mumi_idx[99];   //[Gen_QQ_size]
   Int_t           Gen_QQ_whichRec[99];   //[Gen_QQ_size]
   Int_t           Gen_mu_size;
   Int_t           Gen_mu_type[99];   //[Gen_mu_size]
   Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;
   Int_t           Gen_mu_whichRec[99];   //[Gen_mu_size]


   ////////////////////////// for Acc calculation
   TClonesArray    *Gen_QQ_mupl_4mom;
   TClonesArray    *Gen_QQ_mumi_4mom; 

   Int_t           pclusterCompatibilityFilter;
   Int_t           pprimaryVertexFilter;   // for PbPb
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

   //centrality
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         pthat;

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
   TBranch        *b_Reco_QQ_whichGen;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_Ntrk;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_whichGen;   //!
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
   TBranch        *b_Gen_weight;   //!
   TBranch        *b_Gen_pthat;   //!
   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_type;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_momId;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_idx;   //!
   TBranch        *b_Gen_QQ_mumi_idx;   //!
   TBranch        *b_Gen_QQ_whichRec;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_type;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!
   TBranch        *b_Gen_mu_whichRec;   //!
   TBranch        *b_Reco_QQ_vtx;

   ////////////////////////// for Acc calculation
   TBranch        *b_Gen_QQ_mupl_4mom;
   TBranch        *b_Gen_QQ_mumi_4mom; 


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

   TBranch        *b_hiBin;
   TBranch        *b_hiHF;
   TBranch        *b_pthat;

   oniaTree(Bool_t pbpb = true, Bool_t pr = false, Bool_t acc = false);
   virtual ~oniaTree();

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     AccEffCalc(string caseTag="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     EffCalc(string caseTag="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     AccCalc(string caseTag="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     Plot(string caseTag="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     ClosureTest(string caseTag="",bool onlyPlot=false);
   virtual void     ClosureTestPtWeights();
   virtual void     TnpSyst(string caseTag ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   //virtual void     AccEffStatToy(int nToys=100);
   virtual void     AccEffStatToy_Acc(int nToys=100,string caseLabel ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     AccEffStatToy_Eff(int nToys=100,string caseLabel ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   //virtual void     AccEffStatToy_1D(int nToys=100,string caseLabel ="");
   virtual void     AccEffStat(string caseLabel ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     TnpToy(int min=0, int max=100);
   virtual void     TnpStat(string caseLabel ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   virtual void     AccEffMisMod(string caseLabel ="_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
   //virtual void     FullAccEffSyst(string caseLabel ="");
   virtual void     AccEffSyst_all();
   virtual double   readSyst(const char* systfile, double zedmin, double zedmax, double ymin, double ymax, int centrmin, int centrmax);
   virtual double   rms(vector<double> v, bool isrelative);
   virtual double   maxdiff(vector<double> v, bool isrelative);
   virtual Bool_t   isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit);
   virtual Bool_t   isFilterMatch (Int_t iRecoMu, Int_t TriggerBit);
   virtual Bool_t   isGlobalMuonInAccept2019 (TLorentzVector* Muon);
   virtual Bool_t   areMuonsInAcceptance2019 (Int_t iRecoQQ);
   virtual Bool_t   areGenMuonsInAcceptance2019 (Int_t iGenQQ);
   virtual Bool_t   passQualityCuts2019 (Int_t iRecoQQ);
   virtual Double_t deltaR (TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
   //virtual void AccEffComp();
   virtual Bool_t   isGlobalMuonInAccept2015 (TLorentzVector* Muon);
   virtual Bool_t   areMuonsInAcceptance2015 (Int_t iRecoQQ);
   virtual Bool_t   areGenMuonsInAcceptance2015 (Int_t iGenQQ);
   virtual Bool_t   passQualityCuts2015 (Int_t iRecoQQ);
   virtual void     setBins(string caseTag ="");

   virtual Int_t    getMCHiBinFromhiHF(const Double_t hiHF);
   virtual Double_t findNcoll(int hiBin);
   virtual vector<TObjArray*> ReadFileWeight(bool ispbpb, bool isprompt);
};

#endif

#ifdef makeAccEff_cxx
oniaTree::oniaTree(Bool_t pbpb, Bool_t pr, Bool_t acc) : fChain(0)
{
  TFile* f(0x0);
  isPbPb = pbpb;
  isPr = pr;
  isAcc = acc;

  TString inputFiles [8] = {
    //Acc files
    //"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root", //pp prompt 
    "/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/acc/merged_acc.root",
    //"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root", //pp nonprompt
    "/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/acc/merged_acc.root",
    //"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root", //PbPb prompt
    "/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/acc/merged_acc.root",
    //"root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root", //PbPb  nonprompt
    "/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/acc/merged_acc.root",
    //Eff files
    "/data_CMS/cms/diab/JpsiJet/MC/pp/prompt/v3/HiForestAOD_merged.root", //pp prompt
    "/data_CMS/cms/diab/JpsiJet/MC/pp/nonprompt/v5/HiForestAOD_merged.root", //pp nonprompt 
    "/data_CMS/cms/diab/JpsiJet/MC/PbPb/prompt/v7/HiForestAOD_merged.root", //PbPb prompt
    "/data_CMS/cms/diab/JpsiJet/MC/PbPb/nonprompt/v4/HiForestAOD_merged.root" //PbPb nonprompt
  };

  int inin = 0; // inin = input index

  if (isAcc && !isPbPb && isPr) inin = 0;
  else if (isAcc && !isPbPb && !isPr) inin = 1;
  else if (isAcc && isPbPb && isPr) inin = 2;
  else if (isAcc && isPbPb && !isPr) inin = 3;
  else if (!isAcc && !isPbPb && isPr) inin = 4;
  else if (!isAcc && !isPbPb && !isPr) inin = 5;
  else if (!isAcc && isPbPb && isPr) inin = 6;
  else if (!isAcc && isPbPb && !isPr) inin = 7;


  f = TFile::Open(inputFiles[inin]);

  cout<<Form("[INFO] %s %s tree in file ",(isPbPb?"PbPb":"pp"),(isPr?"prompt MC":"nonprompt MC"))<<f->GetName()<<endl;
  
  TTree * tree (0x0); TTree * skimtree (0x0); TTree * centtree (0x0);
  tree = (TTree*) f->Get("hionia/myTree");
  skimtree = (TTree*) f->Get("skimanalysis/HltTree");
  centtree = (TTree*) f->Get("hiEvtAnalyzer/HiTree");
  if(!tree){
    cout <<"[ERROR] Cannot find the onia tree"<<endl;
    return;
  }
  if(!skimtree)
    cout <<"[WARNING] Cannot find the skimanalysis tree"<<endl;
  else
    tree->AddFriend(skimtree);
  if (!centtree)
    cout<<"[WARNING] Cannot find the hiEvtAnalyzer tree"<<endl;
  else
    tree->AddFriend(centtree);

  Init(tree);
}

oniaTree::~oniaTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t oniaTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t oniaTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void oniaTree::Init(TTree *tree)
{
   Reco_QQ_4mom = 0;
   Reco_mu_4mom = 0;
   Reco_QQ_vtx = 0;

   Gen_QQ_4mom = 0;
   Gen_mu_4mom = 0;

   Gen_QQ_mupl_4mom = 0;
   Gen_QQ_mumi_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);


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
   if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchAddress("Reco_QQ_whichGen", Reco_QQ_whichGen, &b_Reco_QQ_whichGen);
   if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
   if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   if (fChain->GetBranch("Reco_mu_whichGen")) fChain->SetBranchAddress("Reco_mu_whichGen", Reco_mu_whichGen, &b_Reco_mu_whichGen);
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
   if (fChain->GetBranch("Gen_weight")) fChain->SetBranchAddress("Gen_weight", &Gen_weight, &b_Gen_weight);
   if (fChain->GetBranch("Gen_pthat")) fChain->SetBranchAddress("Gen_pthat", &Gen_pthat, &b_Gen_pthat);
   if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   if (fChain->GetBranch("Gen_QQ_type")) fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   if (fChain->GetBranch("Gen_QQ_momId")) fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchAddress("Gen_QQ_mupl_idx", Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
   if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchAddress("Gen_QQ_mumi_idx", Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
   if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchAddress("Gen_QQ_whichRec", Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
   if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   if (fChain->GetBranch("Gen_mu_type")) fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
   if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   if (fChain->GetBranch("Gen_mu_whichRec")) fChain->SetBranchAddress("Gen_mu_whichRec", Gen_mu_whichRec, &b_Gen_mu_whichRec);

   /////////////////for Acc calc
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);


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
   
   if (fChain->GetBranch("hiBin")) fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   if (fChain->GetBranch("hiHF")) fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);


   if (fChain->GetBranch("Reco_QQ_4mom")) fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false);
   if (fChain->GetBranch("Reco_mu_4mom")) fChain->GetBranch("Reco_mu_4mom")->SetAutoDelete(false);
   if (fChain->GetBranch("Gen_QQ_4mom")) fChain->GetBranch("Gen_QQ_4mom")->SetAutoDelete(false);
   if (fChain->GetBranch("Gen_mu_4mom")) fChain->GetBranch("Gen_mu_4mom")->SetAutoDelete(false);
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->GetBranch("Gen_QQ_mupl_4mom")->SetAutoDelete(false);
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->GetBranch("Gen_QQ_mumi_4mom")->SetAutoDelete(false);

   fChain->SetBranchStatus("*",0);
   if (fChain->GetBranch("HLTriggers")) fChain->SetBranchStatus("HLTriggers",1);
   if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchStatus("Reco_QQ_trig",1);
   if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchStatus("Reco_mu_trig",1);
   if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchStatus("Reco_QQ_VtxProb",1);
   
   if (fChain->GetBranch("Reco_mu_SelectionType")) fChain->SetBranchStatus("Reco_mu_SelectionType",1);
   if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchStatus("Reco_mu_TMOneStaTight",1);
   if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchStatus("Reco_mu_nTrkWMea",1);
   if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchStatus("Reco_mu_nPixWMea",1);
   if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchStatus("Reco_mu_dxy",1);
   if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchStatus("Reco_mu_dz",1);
   if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchStatus("Reco_mu_nTrkHits",1);
   if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchStatus("Reco_mu_normChi2_global",1);
   if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchStatus("Reco_mu_normChi2_inner",1);
   if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchStatus("Reco_mu_TrkMuArb",1);
   if (fChain->GetBranch("runNb")) fChain->SetBranchStatus("runNb",1);
   if (fChain->GetBranch("Centrality")) fChain->SetBranchStatus("Centrality",1);
   if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchStatus("Reco_QQ_size",1);
   if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchStatus("Reco_QQ_sign",1);
   if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchStatus("Reco_QQ_4mom",1);
   if (fChain->GetBranch("Reco_QQ_mupl_idx")) fChain->SetBranchStatus("Reco_QQ_mupl_idx",1);
   if (fChain->GetBranch("Reco_QQ_mumi_idx")) fChain->SetBranchStatus("Reco_QQ_mumi_idx",1);
   if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchStatus("Reco_QQ_ctau3D",1);
   if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchStatus("Reco_QQ_ctauErr3D",1);
   if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchStatus("Reco_QQ_ctau",1);
   if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchStatus("Reco_QQ_ctauErr",1);
   if (fChain->GetBranch("Reco_QQ_whichGen")) fChain->SetBranchStatus("Reco_QQ_whichGen",1);
   if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchStatus("Reco_mu_4mom",1);
   if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchStatus("Reco_mu_size",1);
   if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchStatus("pPAprimaryVertexFilter",1);
   if (fChain->GetBranch("pprimaryVertexFilter")) fChain->SetBranchStatus("pprimaryVertexFilter",1);
   if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchStatus("pBeamScrapingFilter",1);
   if (fChain->GetBranch("pclusterCompatibilityFilter")) fChain->SetBranchStatus("pclusterCompatibilityFilter",1);
   if (fChain->GetBranch("phfCoincFilter2Th4")) fChain->SetBranchStatus("phfCoincFilter2Th4",1);
   if (fChain->GetBranch("hiBin")) fChain->SetBranchStatus("hiBin",1);
   if (fChain->GetBranch("Gen_weight")) fChain->SetBranchStatus("Gen_weight",1);
   if (fChain->GetBranch("Gen_pthat")) fChain->SetBranchStatus("Gen_pthat",1);
   if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchStatus("Gen_QQ_4mom",1);
   if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchStatus("Gen_QQ_size",1);
   if (fChain->GetBranch("Gen_QQ_mupl_idx")) fChain->SetBranchStatus("Gen_QQ_mupl_idx",1);
   if (fChain->GetBranch("Gen_QQ_mumi_idx")) fChain->SetBranchStatus("Gen_QQ_mumi_idx",1);
   if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchStatus("Gen_QQ_ctau3D",1);
   if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchStatus("Gen_QQ_ctau",1);
   if (fChain->GetBranch("Gen_QQ_whichRec")) fChain->SetBranchStatus("Gen_QQ_whichRec",1); 
   if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchStatus("Gen_mu_4mom",1);
   if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchStatus("Gen_mu_size",1);
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1);
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1);

   if (fChain->GetBranch("hiBin")) fChain->SetBranchStatus("hiBin",1);
   if (fChain->GetBranch("hiHF")) fChain->SetBranchStatus("hiHF",1);
   if (fChain->GetBranch("pthat")) fChain->SetBranchStatus("pthat", 1);

   Notify();
}

Bool_t oniaTree::Notify()
{
   return kTRUE;
}

void oniaTree::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t oniaTree::Cut(Long64_t entry)
{
   return 1;
}


Bool_t oniaTree::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit)
{
  Bool_t cond = true;
  cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  return cond;
};

Bool_t oniaTree::isFilterMatch (Int_t iRecoMu, Int_t TriggerBit)
{
  Bool_t cond = true;
  //cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  cond = cond && ( (Reco_mu_trig[iRecoMu]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
  return cond;
};

Bool_t oniaTree::isGlobalMuonInAccept2019 (TLorentzVector* Muon)
{
  return ( fabs(Muon->Eta()) < 2.4 &&
	   ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	    (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.47-1.89*fabs(Muon->Eta())) ||
	    (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.5)));
};


Bool_t oniaTree::areMuonsInAcceptance2019 (Int_t iRecoQQ)
{
  TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[iRecoQQ]);
  TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[iRecoQQ]);
  
  return ( isGlobalMuonInAccept2019(RecoQQmupl) && isGlobalMuonInAccept2019(RecoQQmumi) );
};

Bool_t oniaTree::areGenMuonsInAcceptance2019 (Int_t iGenQQ)
{
  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGenQQ]);
  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGenQQ]);

  return ( isGlobalMuonInAccept2019(GenQQmupl) && isGlobalMuonInAccept2019(GenQQmumi) );
};

Bool_t oniaTree::passQualityCuts2019 (Int_t iRecoQQ)
{
  int iMupl = Reco_QQ_mupl_idx[iRecoQQ];
  int iMumi = Reco_QQ_mumi_idx[iRecoQQ];
  
  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  if ( ! (Reco_mu_nTrkWMea[iMumi] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMumi] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMumi]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMumi]) < 20.) ) return false;
  
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  if ( ! (Reco_mu_nTrkWMea[iMupl] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMupl] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMupl]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMupl]) < 20.) ) return false;
  
  if ( ! (Reco_QQ_VtxProb[iRecoQQ] > 0.01) ) return false;
  
  return true;
  
};

Double_t oniaTree::deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

double oniaTree::rms(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double s=0.;
  double s2=0.;
  for (unsigned int i=0; i<v.size(); i++) {
    if (v[i]==-999) continue;
    s2+=(v[i]-v[0])*(v[i]-v[0]);
    //cout<<"[INFO] v["<<i<<"] = "<< v[i] << " s = "<<s<<" s2 = "<<s2<<endl;
  }
  double ans = sqrt(s2*1.0/(v.size()-1));
  //cout<<"[INFO] s2 = "<<s2<<" ans = " <<ans/v[0] <<endl;
  if (isrelative) ans = ans/v[0];
  return ans;
}

double oniaTree::maxdiff(vector<double> v, bool isrelative) {
  if (v.size()==0 || v[0]==0) return 0;
  double maxdiff=0;
  for (unsigned int i=1; i<v.size(); i++) {
    if (v[i]==-999) continue;
    maxdiff=max(maxdiff,fabs(v[i]-v[0]));
  }
  double ans = maxdiff;
  if (isrelative) ans = ans/v[0];
  return ans;
}

double oniaTree::readSyst(const char* systfile, double zedmin, double zedmax, double ymin, double ymax, int centrmin, int centrmax) {
  double ans = 0;
  ifstream file(systfile);
  if (!(file.good())) return ans;
  
  string systname; getline(file,systname);
  
  string line;
  double zmin=0, zmax=0, rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;
  
  while (file.good()) {
    getline(file,line);
    if (line.size()==0) break;
    TString tline(line.c_str());
    TString t; Int_t from = 0, cnt=0;
    while (tline.Tokenize(t, from , ",")) {
      t.Strip(TString::kBoth,' ');
      value = atof(t.Data());
      if (cnt==0) rapmin = atof(t.Data());
      else if (cnt==1) rapmax = value;
      else if (cnt==2) ptmin = value;
      else if (cnt==3) ptmax = value;
      else if (cnt==4) zmin = value;
      else if (cnt==5) zmax = value;
      else if (cnt==6) centmin = value;
      else if (cnt==7) centmax = value;
      else if (cnt>8) {
	cout << "Warning, too many fields, I'll take the last one." << endl;
	continue;
      }
      cnt++;
    }
    if (zmin<zedmin+0.001 && zmin>zedmin-0.001 && zmax<zedmax+0.001 && zmax>zedmax-0.001 && rapmin<ymin+0.001 && rapmin>ymin-0.001 && rapmax<ymax+0.001 && rapmax>ymax-0.001 && centmin<centrmin+0.001 && centmin>centrmin-0.001 && centmax<centrmax+0.001 && centmax>centrmax-0.001) {
      cout<<"zmin = "<<zmin<<", zmax = "<<zmax<<", syst = "<<value<<endl;
      ans = value;
    }
  }
  file.close();
  if (ans==0) cout<<"[WARNING] syst =0";
  return ans;
}

Bool_t oniaTree::isGlobalMuonInAccept2015 (TLorentzVector* Muon)
{
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
};

Bool_t oniaTree::areMuonsInAcceptance2015 (Int_t iRecoQQ)
{
  TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[iRecoQQ]);
  TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[iRecoQQ]);
  
  return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
};

Bool_t oniaTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
{
  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iGenQQ]);
  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iGenQQ]);

  return ( isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi) );
};

Bool_t oniaTree::passQualityCuts2015 (Int_t iRecoQQ)
{
  int iMupl = Reco_QQ_mupl_idx[iRecoQQ];
  int iMumi = Reco_QQ_mumi_idx[iRecoQQ];
  
  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMumi]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  // if ( ! (Reco_mu_highPurity[iMumi]) ) return false;
  if ( ! (Reco_mu_TMOneStaTight[iMumi]==1) ) return false; // = isGoodMuon
  if ( ! (Reco_mu_nTrkWMea[iMumi] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMumi] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMumi]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMumi]) < 20.) ) return false;
  
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 1))) ) return false; // require the muons to be global muons
  if ( ! (Reco_mu_SelectionType[iMupl]&((ULong64_t)pow(2, 3))) ) return false; // require the muons to be tracker muons
  // if ( ! (Reco_mu_highPurity[iMupl]) ) return false;
  if ( ! (Reco_mu_TMOneStaTight[iMupl]==1) ) return false; // = isGoodMuon
  if ( ! (Reco_mu_nTrkWMea[iMupl] > 5) ) return false;
  if ( ! (Reco_mu_nPixWMea[iMupl] > 0) ) return false;
  if ( ! (fabs(Reco_mu_dxy[iMupl]) < 0.3) ) return false;
  if ( ! (fabs(Reco_mu_dz[iMupl]) < 20.) ) return false;
  
  if ( ! (Reco_QQ_VtxProb[iRecoQQ] > 0.01) ) return false;
  
  return true;
  
};

Int_t oniaTree::getMCHiBinFromhiHF(const Double_t hiHF) {
  const Int_t nBins = 200; // table of bin edges
  const Double_t binTable[nBins+1] = {0, 12.2187, 13.0371, 13.7674, 14.5129, 15.2603, 16.0086, 16.7623, 17.5335, 18.3283, 19.1596, 19.9989, 20.8532, 21.7297, 22.6773, 23.6313, 24.6208, 25.6155, 26.6585, 27.7223, 28.8632, 30.041, 31.2865, 32.5431, 33.8655, 35.2539, 36.6912, 38.2064, 39.7876, 41.4818, 43.2416, 45.0605, 46.9652, 48.9918, 51.1, 53.2417, 55.5094, 57.9209, 60.3817, 62.9778, 65.6099, 68.4352, 71.3543, 74.4154, 77.6252, 80.8425, 84.1611, 87.7395, 91.3973, 95.1286, 99.0571, 103.185, 107.482, 111.929, 116.45, 121.178, 126.081, 130.995, 136.171, 141.612, 147.298, 153.139, 159.419, 165.633, 172.114, 178.881, 185.844, 192.845, 200.244, 207.83, 215.529, 223.489, 231.878, 240.254, 249.319, 258.303, 267.508, 277.037, 286.729, 296.845, 307.458, 317.882, 328.787, 340.074, 351.295, 362.979, 375.125, 387.197, 399.604, 412.516, 425.683, 439.001, 452.667, 466.816, 481.007, 495.679, 510.588, 526.138, 541.782, 557.641, 574.141, 591.071, 608.379, 626.068, 643.616, 661.885, 680.288, 699.449, 718.925, 738.968, 758.983, 779.459, 800.376, 821.638, 843.555, 865.771, 888.339, 911.031, 934.979, 958.56, 982.582, 1007.02, 1031.9, 1057.81, 1084.01, 1111.71, 1138.21, 1165.72, 1193.73, 1221.65, 1251.51, 1281.23, 1311.01, 1341.1, 1372.4, 1404.29, 1436.52, 1468.65, 1501.91, 1535.56, 1569.69, 1604.69, 1640.65, 1676.05, 1712.62, 1749.28, 1787.43, 1825.89, 1866.07, 1906.58, 1947.84, 1989.66, 2031.4, 2072.8, 2115.32, 2159.5, 2205.23, 2252.68, 2298.58, 2345.65, 2393.36, 2442.87, 2491.45, 2541.04, 2592.81, 2645.52, 2699.1, 2753.29, 2807.93, 2864.37, 2922.6, 2979.42, 3038.68, 3098.72, 3159.29, 3221.66, 3285.9, 3350.95, 3415.81, 3482.69, 3552.62, 3623.61, 3694.63, 3767.25, 3840.28, 3917.04, 3993.66, 4073.36, 4154.33, 4238.13, 4322.21, 4409.83, 4498.89, 4589.72, 4681.56, 4777.09, 4877.95, 4987.05, 5113.04, 5279.58, 6242.82};
  
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
      }
  }
  binPos = nBins - 1 - binPos;
  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

Double_t oniaTree::findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};

TH1D* makePull(TH1D* oldHist, TH1D* newHist, bool systErr=false){
  TH1D* pullHist = (TH1D*)oldHist->Clone("pullHist");

  pullHist->Divide(newHist);

  if (systErr) {
    int nbin = pullHist->GetNbinsX();
    for (int i=0;i<=nbin;i++) {
      pullHist->SetBinError(i, oldHist->GetBinError(i)/newHist->GetBinContent(i));
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

#endif // #ifdef makeAccEff_cxx

