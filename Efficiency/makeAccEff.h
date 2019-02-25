#ifndef makeAccEff_h
#define makeAccEFf_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TClonesArray.h"
#include <TH2.h>
#include <TStyle.h>
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
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>

#include "tnp_weight.h"

using namespace std;
using namespace  RooFit;

class oniaTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
   bool          isPr;
   bool          isAcc;

   //my additional variables

   TLorentzVector* matchReco;
   TLorentzVector* matchGen;
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
   int triggerIndex_PP =0;
   int k;
   int ibest; //reco-gen matching                                                                                                                                                                          
   Float_t weight;
   Float_t tnp_weight;

   // Declaration of leaf types
   UInt_t          eventNb;
   UInt_t          runNb;
   UInt_t          LS;
   Float_t         zVtx;
   Float_t         nPV;
   Int_t           Centrality;
   Int_t           nTrig;
   Int_t           trigPrescale[99];   //[nTrig]
   ULong64_t       HLTriggers;
   Int_t           Npix;
   Int_t           NpixelTracks;
   Int_t           Ntracks;
   Float_t         SumET_HF;
   Float_t         SumET_HFplus;
   Float_t         SumET_HFminus;
   Float_t         SumET_HFplusEta4;
   Float_t         SumET_HFminusEta4;
   Float_t         SumET_ET;
   Float_t         SumET_EE;
   Float_t         SumET_EB;
   Float_t         SumET_EEplus;
   Float_t         SumET_EEminus;
   Float_t         SumET_ZDC;
   Float_t         SumET_ZDCplus;
   Float_t         SumET_ZDCminus;
   Int_t           Reco_QQ_size;
   Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   TClonesArray    *Reco_QQ_mupl_4mom;
   TClonesArray    *Reco_QQ_mumi_4mom;
   ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mupl_trig[99];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mumi_trig[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_SelectionType[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_SelectionType[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_isGoodMuon[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_isGoodMuon[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_highPurity[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_highPurity[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TrkMuArb[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TrkMuArb[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TMOneStaTight[99];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TMOneStaTight[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixValHits[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixValHits[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nMuValHits[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nMuValHits[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkHits[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkHits[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_global[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_global[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixWMea[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixWMea[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkWMea[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkWMea[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_StationsMatched[99];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_StationsMatched[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxy[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxy[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxyErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxyErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dz[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dz[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dzErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dzErr[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_global[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_global[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_inner[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_global[99];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_global[99];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[99];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[99];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
   Bool_t          Reco_mu_isGoodMuon[99];   //[Reco_mu_size]
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


   Int_t           Gen_QQ_size;
   Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_4mom;
   Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_mupl_4mom;
   TClonesArray    *Gen_QQ_mumi_4mom;
   Int_t           Gen_mu_size;
   Int_t           Gen_mu_type[99];   //[Gen_mu_size]
   Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;


   // List of branches
   TBranch        *b_eventNb;   //!
   TBranch        *b_runNb;   //!
   TBranch        *b_LS;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_Centrality;   //!
   TBranch        *b_nTrig;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_HLTriggers;   //!
   TBranch        *b_Npix;   //!
   TBranch        *b_NpixelTracks;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_SumET_HF;   //!
   TBranch        *b_SumET_HFplus;   //!
   TBranch        *b_SumET_HFminus;   //!
   TBranch        *b_SumET_HFplusEta4;   //!
   TBranch        *b_SumET_HFminusEta4;   //!
   TBranch        *b_SumET_ET;   //!
   TBranch        *b_SumET_EE;   //!
   TBranch        *b_SumET_EB;   //!
   TBranch        *b_SumET_EEplus;   //!
   TBranch        *b_SumET_EEminus;   //!
   TBranch        *b_SumET_ZDC;   //!
   TBranch        *b_SumET_ZDCplus;   //!
   TBranch        *b_SumET_ZDCminus;   //!
   TBranch        *b_Reco_QQ_size;   //!
   TBranch        *b_Reco_QQ_type;   //!
   TBranch        *b_Reco_QQ_sign;   //!
   TBranch        *b_Reco_QQ_4mom;   //!
   TBranch        *b_Reco_QQ_mupl_4mom;   //!
   TBranch        *b_Reco_QQ_mumi_4mom;   //!
   TBranch        *b_Reco_QQ_trig;   //!
   TBranch        *b_Reco_QQ_mupl_trig;   //!
   TBranch        *b_Reco_QQ_mumi_trig;   //!
   TBranch        *b_Reco_QQ_isCowboy;   //!
   TBranch        *b_Reco_QQ_ctau;   //!
   TBranch        *b_Reco_QQ_ctauErr;   //!
   TBranch        *b_Reco_QQ_ctau3D;   //!
   TBranch        *b_Reco_QQ_ctauErr3D;   //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_QQ_dca;   //!
   TBranch        *b_Reco_QQ_MassErr;   //!
   TBranch        *b_Reco_QQ_vtx;   //!
   TBranch        *b_Reco_QQ_Ntrk;   //!
   TBranch        *b_Reco_QQ_mupl_SelectionType;   //!
   TBranch        *b_Reco_QQ_mumi_SelectionType;   //!
   TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
   TBranch        *b_Reco_QQ_mupl_highPurity;   //!
   TBranch        *b_Reco_QQ_mumi_highPurity;   //!
   TBranch        *b_Reco_QQ_mupl_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mumi_TrkMuArb;   //!
   TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
   TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_inner;   //!
   TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
   TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
   TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
   TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
   TBranch        *b_Reco_QQ_mupl_dxy;   //!
   TBranch        *b_Reco_QQ_mumi_dxy;   //!
   TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
   TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
   TBranch        *b_Reco_QQ_mupl_dz;   //!
   TBranch        *b_Reco_QQ_mumi_dz;   //!
   TBranch        *b_Reco_QQ_mupl_dzErr;   //!
   TBranch        *b_Reco_QQ_mumi_dzErr;   //!
   TBranch        *b_Reco_QQ_mupl_pt_inner;   //!
   TBranch        *b_Reco_QQ_mumi_pt_inner;   //!
   TBranch        *b_Reco_QQ_mupl_pt_global;   //!
   TBranch        *b_Reco_QQ_mumi_pt_global;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_inner;   //!
   TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
   TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
   TBranch        *b_Reco_mu_size;   //!
   TBranch        *b_Reco_mu_type;   //!
   TBranch        *b_Reco_mu_SelectionType;   //!
   TBranch        *b_Reco_mu_charge;   //!
   TBranch        *b_Reco_mu_4mom;   //!
   TBranch        *b_Reco_mu_trig;   //!
   TBranch        *b_Reco_mu_isGoodMuon;   //!
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


   TBranch        *b_Gen_QQ_size;   //!
   TBranch        *b_Gen_QQ_type;   //!
   TBranch        *b_Gen_QQ_4mom;   //!
   TBranch        *b_Gen_QQ_momId;   //!
   TBranch        *b_Gen_QQ_ctau;   //!
   TBranch        *b_Gen_QQ_ctau3D;   //!
   TBranch        *b_Gen_QQ_mupl_4mom;   //!
   TBranch        *b_Gen_QQ_mumi_4mom;   //!
   TBranch        *b_Gen_mu_size;   //!
   TBranch        *b_Gen_mu_type;   //!
   TBranch        *b_Gen_mu_charge;   //!
   TBranch        *b_Gen_mu_4mom;   //!



   oniaTree(Bool_t pr = true, Bool_t acc = false);
   virtual ~oniaTree();

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     AccEffCalc();
   virtual void     EffCalc();
   virtual void     AccCalc();
   virtual void     ANEffPlots();
   virtual void     Plot();
   virtual void     ClosureTest();
   virtual void     TnpSyst(string caseLabel ="");
   virtual void     AccEffStatToy(int nToys=100);
   virtual void     AccEffStat(string caseLabel ="");
   virtual void     TnpToy(int min=0, int max=100);
   virtual void     TnpStat(string caseLabel ="");
   virtual void     AccEffMisMod(string caseLabel ="");
   virtual void     FullAccEffSyst(string caseLabel ="");
   virtual void     AccEffSyst_all();
   virtual double   readSyst(const char* systfile, double zedmin, double zedmax, double ymin, double ymax);
   virtual double   rms(vector<double> v, bool isrelative);
   virtual double   maxdiff(vector<double> v, bool isrelative);
   virtual Bool_t   isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit);
   virtual Bool_t   isGlobalMuonInAccept2015 (TLorentzVector* Muon);
   virtual Bool_t   areMuonsInAcceptance2015 (Int_t iRecoQQ);
   virtual Bool_t   areGenMuonsInAcceptance2015 (Int_t iGenQQ);
   virtual Bool_t   passQualityCuts2015 (Int_t iRecoQQ);
   virtual Bool_t   isMatchedRecoDiMuon (int iRecoDiMuon, double maxDeltaR=0.03);
   virtual Bool_t   isMatchedGenDiMuon (int iGenDiMuon, double maxDeltaR=0.03);
   virtual Double_t deltaR (TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
   virtual void AccEffComp();
   virtual vector<TObjArray*> ReadFileWeight(bool isprompt);
};

#endif

#ifdef makeAccEff_cxx
oniaTree::oniaTree(Bool_t pr = true, Bool_t acc = false) : fChain(0)
{
  TFile* f(0x0);
  isPr=pr;
  isAcc=acc;
  f = TFile::Open(Form("/data_CMS/cms/mnguyen/jPsiJet/mc/%s/%s.root", isPr?"prompt":"nonprompt",isAcc?"acc/merged_acc":"v9_ext/merged_HiForestAOD"));
  //if (isPr) {
  //if (!isAcc){
  //f = TFile::Open("/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/v9_ext/merged_HiForestAOD.root");}
  //else {
  //f = TFile::Open("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_JpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");}
  //}
  //else {
  //if (!isAcc) {
  //f = TFile::Open("/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/v9_ext/merged_HiForestAOD.root");}
  //else {
  //f = TFile::Open("root://xrootd.unl.edu//store/group/phys_heavyions/dileptons/MC2015/pp502TeV/TTrees/OniaTree_BJpsiMM_5p02TeV_TuneCUETP8M1_nofilter_pp502Fall15-MCRUN2_71_V1-v1_GENONLY.root");}
  //}
  cout<<Form("[INFO] %s from %s tree in file ",(isPr?"prompt MC":"nonprompt MC"), isAcc?"acc":"eff")<<f->GetName()<<endl;
  
  TTree * tree (0x0);
  tree = (TTree*) f->Get("hionia/myTree");
  if(!tree){
    cout <<"error in the onia tree"<<endl;
    return;
  }
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
   Reco_QQ_mupl_4mom = 0;
   Reco_QQ_mumi_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;

   Gen_QQ_4mom = 0;
   Gen_QQ_mupl_4mom = 0;
   Gen_QQ_mumi_4mom = 0;
   Gen_mu_4mom = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("LS", &LS, &b_LS);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   //fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
   fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
   fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
   fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
   fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
   fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
   fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
   fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
   fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
   fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
   fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
   fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
   fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
   fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
   fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
   if (!isAcc){
     fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
     fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
     fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
     fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
     fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
     fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
     fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
     fChain->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
     fChain->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
     fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
     fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
     fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
     fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
     fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
     fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
     fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
     fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
     fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
     fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
     fChain->SetBranchAddress("Reco_QQ_mupl_SelectionType", Reco_QQ_mupl_SelectionType, &b_Reco_QQ_mupl_SelectionType);
     fChain->SetBranchAddress("Reco_QQ_mumi_SelectionType", Reco_QQ_mumi_SelectionType, &b_Reco_QQ_mumi_SelectionType);
     fChain->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
     fChain->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
     fChain->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
     fChain->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
     fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
     fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
     fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
     fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
     fChain->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
     fChain->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
     fChain->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
     fChain->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
     fChain->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
     fChain->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
     fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
     fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
     fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
     fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
     fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
     fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
     fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
     fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
     fChain->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
     fChain->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
     fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
     fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
     fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
     fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
     fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
     fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
     fChain->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
     fChain->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
     fChain->SetBranchAddress("Reco_QQ_mupl_pt_inner", Reco_QQ_mupl_pt_inner, &b_Reco_QQ_mupl_pt_inner);
     fChain->SetBranchAddress("Reco_QQ_mumi_pt_inner", Reco_QQ_mumi_pt_inner, &b_Reco_QQ_mumi_pt_inner);
     fChain->SetBranchAddress("Reco_QQ_mupl_pt_global", Reco_QQ_mupl_pt_global, &b_Reco_QQ_mupl_pt_global);
     fChain->SetBranchAddress("Reco_QQ_mumi_pt_global", Reco_QQ_mumi_pt_global, &b_Reco_QQ_mumi_pt_global);
     fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_inner", Reco_QQ_mupl_ptErr_inner, &b_Reco_QQ_mupl_ptErr_inner);
     fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_inner", Reco_QQ_mumi_ptErr_inner, &b_Reco_QQ_mumi_ptErr_inner);
     fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
     fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
     fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
     fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
     fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
     fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
     fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
     fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
     fChain->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
     fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
     fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
     fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
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
     fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
     fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
     fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
     fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
     }
   
   
   ////// gen branches
   //if (!isAcc){
   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
   fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
   fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
   fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   //}
   //else {
   //fChain->SetBranchAddress("Reco_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   //fChain->SetBranchAddress("Reco_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
   //fChain->SetBranchAddress("Reco_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   //fChain->SetBranchAddress("Reco_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
   //fChain->SetBranchAddress("Reco_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
   //fChain->SetBranchAddress("Reco_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
   //fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
   //fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
   //fChain->SetBranchAddress("Reco_mu_size", &Gen_mu_size, &b_Gen_mu_size);
   //fChain->SetBranchAddress("Reco_mu_type", Gen_mu_type, &b_Gen_mu_type);
   //fChain->SetBranchAddress("Reco_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
   //fChain->SetBranchAddress("Reco_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   //}


   if (!isAcc){   
     if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false);      }
     if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(false); }
     if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(false); }
   }
   
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->GetBranch("Gen_QQ_mupl_4mom")->SetAutoDelete(false); }
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->GetBranch("Gen_QQ_mumi_4mom")->SetAutoDelete(false); }
   fChain->SetBranchStatus("*",0);
   
   if (fChain->GetBranch("runNb"))             { fChain->SetBranchStatus("runNb",1);             }
   if (fChain->GetBranch("Centrality"))        { fChain->SetBranchStatus("Centrality",1);        }
   if (!isAcc)
     {
       if (fChain->GetBranch("Reco_QQ_size"))      { fChain->SetBranchStatus("Reco_QQ_size",1);      }
       if (fChain->GetBranch("Reco_QQ_sign"))      { fChain->SetBranchStatus("Reco_QQ_sign",1);      }
       if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->SetBranchStatus("Reco_QQ_4mom",1);      }
       if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1); }
       if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1); }
       if (fChain->GetBranch("Reco_QQ_ctau"))      { fChain->SetBranchStatus("Reco_QQ_ctau",1);      }
     }
   if (fChain->GetBranch("Gen_QQ_4mom"))      { fChain->SetBranchStatus("Gen_QQ_4mom",1);      }
   if (fChain->GetBranch("Gen_QQ_size"))      { fChain->SetBranchStatus("Gen_QQ_size",1);      }
   if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1); }
   if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1); }

   if (!isAcc){   
     if (fChain->GetBranch("HLTriggers")) fChain->SetBranchStatus("HLTriggers",1);
     if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchStatus("Reco_QQ_trig",1);
     if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchStatus("Reco_QQ_VtxProb",1);
     if (fChain->GetBranch("Reco_QQ_mupl_SelectionType")) fChain->SetBranchStatus("Reco_QQ_mupl_SelectionType",1);
     if (fChain->GetBranch("Reco_QQ_mumi_SelectionType")) fChain->SetBranchStatus("Reco_QQ_mumi_SelectionType",1);
     if (fChain->GetBranch("Reco_QQ_mupl_isGoodMuon")) fChain->SetBranchStatus("Reco_QQ_mupl_isGoodMuon",1);
     if (fChain->GetBranch("Reco_QQ_mumi_isGoodMuon")) fChain->SetBranchStatus("Reco_QQ_mumi_isGoodMuon",1);
     if (fChain->GetBranch("Reco_QQ_mupl_nTrkWMea")) fChain->SetBranchStatus("Reco_QQ_mupl_nTrkWMea",1);
     if (fChain->GetBranch("Reco_QQ_mumi_nTrkWMea")) fChain->SetBranchStatus("Reco_QQ_mumi_nTrkWMea",1);
     if (fChain->GetBranch("Reco_QQ_mupl_nPixWMea")) fChain->SetBranchStatus("Reco_QQ_mupl_nPixWMea",1);
     if (fChain->GetBranch("Reco_QQ_mumi_nPixWMea")) fChain->SetBranchStatus("Reco_QQ_mumi_nPixWMea",1);
     if (fChain->GetBranch("Reco_QQ_mupl_dxy")) fChain->SetBranchStatus("Reco_QQ_mupl_dxy",1);
     if (fChain->GetBranch("Reco_QQ_mumi_dxy")) fChain->SetBranchStatus("Reco_QQ_mumi_dxy",1);
     if (fChain->GetBranch("Reco_QQ_mupl_dz")) fChain->SetBranchStatus("Reco_QQ_mupl_dz",1);
     if (fChain->GetBranch("Reco_QQ_mumi_dz")) fChain->SetBranchStatus("Reco_QQ_mumi_dz",1);
     if (fChain->GetBranch("Reco_QQ_mupl_nTrkHits")) fChain->SetBranchStatus("Reco_QQ_mupl_nTrkHits",1);
     if (fChain->GetBranch("Reco_QQ_mumi_nTrkHits")) fChain->SetBranchStatus("Reco_QQ_mumi_nTrkHits",1);
     if (fChain->GetBranch("Reco_QQ_mupl_normChi2_global")) fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_global",1);
     if (fChain->GetBranch("Reco_QQ_mumi_normChi2_global")) fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_global",1);
     if (fChain->GetBranch("Reco_QQ_mupl_normChi2_inner")) fChain->SetBranchStatus("Reco_QQ_mupl_normChi2_inner",1);
     if (fChain->GetBranch("Reco_QQ_mumi_normChi2_inner")) fChain->SetBranchStatus("Reco_QQ_mumi_normChi2_inner",1);
     if (fChain->GetBranch("Reco_QQ_mupl_TrkMuArb")) fChain->SetBranchStatus("Reco_QQ_mupl_TrkMuArb",1);
     if (fChain->GetBranch("Reco_QQ_mumi_TrkMuArb")) fChain->SetBranchStatus("Reco_QQ_mumi_TrkMuArb",1);
   }
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

Bool_t oniaTree::isGlobalMuonInAccept2015 (TLorentzVector* Muon)
{
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
};

Bool_t oniaTree::areMuonsInAcceptance2015 (Int_t iRecoQQ)
{
  TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
  TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
  return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
};

Bool_t oniaTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
{
  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenQQ);
  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenQQ);
  return (isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi));
};

Bool_t oniaTree::passQualityCuts2015 (Int_t iRecoQQ)
{
  Bool_t cond = true;
  cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
  cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));
  // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);                                                                                                                                               
  cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
  cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
  cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
  cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
  cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);

  cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
  cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));
  // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);                                                                                                                                               
  cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
  cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
  cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
  cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
  cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);
  cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);

  return cond;
};

Double_t oniaTree::deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

Bool_t oniaTree::isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
{
  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);

  bool isMatched(false);
  int iGenMuon(0);
  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
    {
      TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
      TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
      double dRpl = deltaR(GenMuonpl,RecoMuonpl);
      double dRmi = deltaR(GenMuonmi,RecoMuonmi);
      if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  )
	{
	  isMatched = true;
	  matchGen = (TLorentzVector*) Gen_QQ_4mom->At(iGenMuon);
	}
      iGenMuon++;
    }
  return isMatched;
};

Bool_t oniaTree::isMatchedGenDiMuon(int iGenDiMuon, double maxDeltaR)
{
  TLorentzVector* GenMuonpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuonmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuon = (TLorentzVector*) Gen_QQ_4mom->At(iGenDiMuon);

  bool isMatched(false);
  int iRecoMuon(0);
  while ( !isMatched && (iRecoMuon < Reco_QQ_size) )
    {
      TLorentzVector *RecoMuonpl = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRecoMuon);
      TLorentzVector *RecoMuonmi = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRecoMuon);
      double dRpl = deltaR(GenMuonpl,RecoMuonpl);
      double dRmi = deltaR(GenMuonmi,RecoMuonmi);
      if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR) )
	{
	  isMatched = true;
	  matchReco = (TLorentzVector*) Reco_QQ_4mom->At(iRecoMuon);
	  ibest = iRecoMuon;
	}
      iRecoMuon++;
    }

  return isMatched;
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

double oniaTree::readSyst(const char* systfile, double zedmin, double zedmax, double ymin, double ymax) {
  double ans;
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
    if (zmin<zedmin+0.001 && zmin>zedmin-0.001 && zmax<zedmax+0.001 && zmax>zedmax-0.001 && rapmin<ymin+0.001 && rapmin>ymin-0.001 && rapmax<ymax+0.001 && rapmax>ymax-0.001)
      ans = value;
  }
  file.close();
  if (ans==0) cout<<"[WARNING] syst =0";
  return ans;
}


#endif // #ifdef makeAccEff_cxx
