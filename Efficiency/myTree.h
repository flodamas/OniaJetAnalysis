//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 23 15:36:01 2017 by ROOT version 6.02/13
// from TTree myTree/My TTree of dimuons
// found on file: /data_CMS/cms/mnguyen/jPsiJet/HiForestAOD.root
//////////////////////////////////////////////////////////

#ifndef myTree_h
#define myTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class myTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.


   bool          isMc;
   bool          isPr;

   // Declaration of leaf types
   UInt_t          eventNb;
   UInt_t          runNb;
   UInt_t          LS;
   Float_t         zVtx;
   Float_t         nPV;
   Int_t           Centrality;
   Int_t           nTrig;
   Int_t           trigPrescale[15];   //[nTrig]
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
   Int_t           Reco_QQ_type[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_sign[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_4mom;
   TClonesArray    *Reco_QQ_mupl_4mom;
   TClonesArray    *Reco_QQ_mumi_4mom;
   ULong64_t       Reco_QQ_trig[15];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mupl_trig[15];   //[Reco_QQ_size]
   ULong64_t       Reco_QQ_mumi_trig[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_isCowboy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctau3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_ctauErr3D[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_VtxProb[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_dca[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_MassErr[15];   //[Reco_QQ_size]
   TClonesArray    *Reco_QQ_vtx;
   Int_t           Reco_QQ_Ntrk[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_SelectionType[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_SelectionType[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_isGoodMuon[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_isGoodMuon[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_highPurity[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_highPurity[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TrkMuArb[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TrkMuArb[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mupl_TMOneStaTight[15];   //[Reco_QQ_size]
   Bool_t          Reco_QQ_mumi_TMOneStaTight[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nMuValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nMuValHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkHits[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkHits[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_normChi2_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_normChi2_global[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nPixWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nPixWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_nTrkWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_nTrkWMea[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mupl_StationsMatched[15];   //[Reco_QQ_size]
   Int_t           Reco_QQ_mumi_StationsMatched[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxy[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dxyErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dxyErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dz[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dz[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_dzErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_dzErr[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_pt_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_pt_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_inner[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mupl_ptErr_global[15];   //[Reco_QQ_size]
   Float_t         Reco_QQ_mumi_ptErr_global[15];   //[Reco_QQ_size]
   Int_t           Reco_mu_size;
   Int_t           Reco_mu_type[6];   //[Reco_mu_size]
   Int_t           Reco_mu_SelectionType[6];   //[Reco_mu_size]
   Int_t           Reco_mu_charge[6];   //[Reco_mu_size]
   TClonesArray    *Reco_mu_4mom;
   ULong64_t       Reco_mu_trig[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_isGoodMuon[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_highPurity[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TrkMuArb[6];   //[Reco_mu_size]
   Bool_t          Reco_mu_TMOneStaTight[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nMuValHits[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkHits[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_normChi2_global[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nPixWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_nTrkWMea[6];   //[Reco_mu_size]
   Int_t           Reco_mu_StationsMatched[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxy[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dxyErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dz[6];   //[Reco_mu_size]
   Float_t         Reco_mu_dzErr[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_pt_global[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_inner[6];   //[Reco_mu_size]
   Float_t         Reco_mu_ptErr_global[6];   //[Reco_mu_size]


   Int_t           Gen_QQ_size;
   Int_t           Gen_QQ_type[3];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_4mom;
   Int_t           Gen_QQ_momId[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau[3];   //[Gen_QQ_size]
   Float_t         Gen_QQ_ctau3D[3];   //[Gen_QQ_size]
   TClonesArray    *Gen_QQ_mupl_4mom;
   TClonesArray    *Gen_QQ_mumi_4mom;
   Int_t           Gen_mu_size;
   Int_t           Gen_mu_type[18];   //[Gen_mu_size]
   Int_t           Gen_mu_charge[18];   //[Gen_mu_size]
   TClonesArray    *Gen_mu_4mom;


   // Declaration of leaf types from HltTree.h


       Int_t           HLT_HIL1DoubleMu0_v1;
       Int_t           HLT_HIL1DoubleMu0ForPPRef_v1;



   // Declaration of leaf types from HltTree1.h
   Int_t           Onia2MuMuPAT;
   Int_t           ana_step;
   Int_t           pHBHENoiseFilterResultProducer;
   Int_t           HBHENoiseFilterResult;
   Int_t           HBHENoiseFilterResultRun1;
   Int_t           HBHENoiseFilterResultRun2Loose;
   Int_t           HBHENoiseFilterResultRun2Tight;
   Int_t           HBHEIsoNoiseFilterResult;
   Int_t           pPAprimaryVertexFilter;
   Int_t           pBeamScrapingFilter;
   Int_t           pVertexFilterCutG;
   Int_t           pVertexFilterCutGloose;
   Int_t           pVertexFilterCutGtight;
   Int_t           pVertexFilterCutGplus;
   Int_t           pVertexFilterCutE;
   Int_t           pVertexFilterCutEandG;

   // Declaration of leaf types from HiTree.h
   UInt_t          run;
   //ULong64_t       evt;
   UInt_t          lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;

   Int_t           ProcessID;
   // Float_t         pthat;
   Float_t         weight;
   Float_t         alphaQCD;
   Float_t         alphaQED;
   Float_t         qScale;
   Int_t           nMEPartons;
   Int_t           nMEPartonsFiltered;
 //pair<int,int>   *pdfID;
   //Int_t           first;
   //Int_t           second;
 //pair<float,float> *pdfX;
   //Float_t         first;
   //Float_t         second;
 //pair<float,float> *pdfXpdf;
   Float_t         first;
   Float_t         second;
   vector<float>   *ttbar_w;
   vector<int>     *npus;
   vector<float>   *tnpus;

   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiHFplusEta4;
   Float_t         hiHFminusEta4;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixelTracks;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;

   // Declaration of leaf types from t.h
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[56];   //[nref]
   Float_t         jtpt[56];   //[nref]
   Float_t         jteta[56];   //[nref]
   Float_t         jty[56];   //[nref]
   Float_t         jtphi[56];   //[nref]
   Float_t         jtpu[56];   //[nref]
   Float_t         jtm[56];   //[nref]
   Float_t         jtarea[56];   //[nref]
   Float_t         jtPfCHF[56];   //[nref]
   Float_t         jtPfNHF[56];   //[nref]
   Float_t         jtPfCEF[56];   //[nref]
   Float_t         jtPfNEF[56];   //[nref]
   Float_t         jtPfMUF[56];   //[nref]
   Int_t           jtPfCHM[56];   //[nref]
   Int_t           jtPfNHM[56];   //[nref]
   Int_t           jtPfCEM[56];   //[nref]
   Int_t           jtPfNEM[56];   //[nref]
   Int_t           jtPfMUM[56];   //[nref]
   Float_t         jttau1[56];   //[nref]
   Float_t         jttau2[56];   //[nref]
   Float_t         jttau3[56];   //[nref]
   Float_t         discr_jetID_cuts[56];   //[nref]
   Float_t         discr_jetID_bdt[56];   //[nref]
   Float_t         discr_fr01[56];   //[nref]
   Float_t         trackMax[56];   //[nref]
   Float_t         trackSum[56];   //[nref]
   Int_t           trackN[56];   //[nref]
   Float_t         trackHardSum[56];   //[nref]
   Int_t           trackHardN[56];   //[nref]
   Float_t         chargedMax[56];   //[nref]
   Float_t         chargedSum[56];   //[nref]
   Int_t           chargedN[56];   //[nref]
   Float_t         chargedHardSum[56];   //[nref]
   Int_t           chargedHardN[56];   //[nref]
   Float_t         photonMax[56];   //[nref]
   Float_t         photonSum[56];   //[nref]
   Int_t           photonN[56];   //[nref]
   Float_t         photonHardSum[56];   //[nref]
   Int_t           photonHardN[56];   //[nref]
   Float_t         neutralMax[56];   //[nref]
   Float_t         neutralSum[56];   //[nref]
   Int_t           neutralN[56];   //[nref]
   Float_t         hcalSum[56];   //[nref]
   Float_t         ecalSum[56];   //[nref]
   Float_t         eMax[56];   //[nref]
   Float_t         eSum[56];   //[nref]
   Int_t           eN[56];   //[nref]
   Float_t         muMax[56];   //[nref]
   Float_t         muSum[56];   //[nref]
   Int_t           muN[56];   //[nref]
   Float_t         discr_ssvHighEff[56];   //[nref]
   Float_t         discr_ssvHighPur[56];   //[nref]
   Float_t         discr_csvV1[56];   //[nref]
   Float_t         discr_csvV2[56];   //[nref]
   Float_t         discr_muByIp3[56];   //[nref]
   Float_t         discr_muByPt[56];   //[nref]
   Float_t         discr_prob[56];   //[nref]
   Float_t         discr_probb[56];   //[nref]
   Float_t         discr_tcHighEff[56];   //[nref]
   Float_t         discr_tcHighPur[56];   //[nref]
   Float_t         ndiscr_ssvHighEff[56];   //[nref]
   Float_t         ndiscr_ssvHighPur[56];   //[nref]
   Float_t         ndiscr_csvV1[56];   //[nref]
   Float_t         ndiscr_csvV2[56];   //[nref]
   Float_t         ndiscr_muByPt[56];   //[nref]
   Float_t         pdiscr_csvV1[56];   //[nref]
   Float_t         pdiscr_csvV2[56];   //[nref]
   Int_t           nsvtx[56];   //[nref]
   Int_t           svtxntrk[56];   //[nref]
   Float_t         svtxdl[56];   //[nref]
   Float_t         svtxdls[56];   //[nref]
   Float_t         svtxdl2d[56];   //[nref]
   Float_t         svtxdls2d[56];   //[nref]
   Float_t         svtxm[56];   //[nref]
   Float_t         svtxpt[56];   //[nref]
   Float_t         svtxmcorr[56];   //[nref]
   Int_t           nIPtrk[56];   //[nref]
   Int_t           nselIPtrk[56];   //[nref]
   Float_t         mue[56];   //[nref]
   Float_t         mupt[56];   //[nref]
   Float_t         mueta[56];   //[nref]
   Float_t         muphi[56];   //[nref]
   Float_t         mudr[56];   //[nref]
   Float_t         muptrel[56];   //[nref]
   Int_t           muchg[56];   //[nref]



   Int_t           beamId1;
   Int_t           beamId2;
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



   //List of branshes from HltTree.h

   TBranch        *b_HLT_HIL1DoubleMu0_v1;   //!
   TBranch        *b_HLT_HIL1DoubleMu0ForPPRef_v1;   //!


   //List of branshes from HltTree1.h
   TBranch        *b_Onia2MuMuPAT;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_pHBHENoiseFilterResultProducer;   //!
   TBranch        *b_HBHENoiseFilterResult;   //!
   TBranch        *b_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_HBHEIsoNoiseFilterResult;   //!
   TBranch        *b_pPAprimaryVertexFilter;   //!
   TBranch        *b_pBeamScrapingFilter;   //!
   TBranch        *b_pVertexFilterCutG;   //!
   TBranch        *b_pVertexFilterCutGloose;   //!
   TBranch        *b_pVertexFilterCutGtight;   //!
   TBranch        *b_pVertexFilterCutGplus;   //!
   TBranch        *b_pVertexFilterCutE;   //!
   TBranch        *b_pVertexFilterCutEandG;   //!

   // List of branches from HiTree.h
   TBranch        *b_run;   //!
   //   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!

   TBranch        *b_ProcessID;   //!
   //TBranch        *b_pthat;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_alphaQCD;   //!
   TBranch        *b_alphaQED;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_nMEPartons;   //!
   TBranch        *b_nMEPartonsFiltered;   //!
   TBranch        *b_pdfID_first;   //!
   TBranch        *b_pdfID_second;   //!
   TBranch        *b_pdfX_first;   //!
   TBranch        *b_pdfX_second;   //!
   TBranch        *b_pdfXpdf_first;   //!
   TBranch        *b_pdfXpdf_second;   //!
   TBranch        *b_ttbar_w;   //!
   TBranch        *b_npus;   //!
   TBranch        *b_tnpus;   //!


   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiHFplusEta4;   //!
   TBranch        *b_hiHFminusEta4;   //!
   TBranch        *b_hiZDC;   //!
   TBranch        *b_hiZDCplus;   //!
   TBranch        *b_hiZDCminus;   //!
   TBranch        *b_hiHFhit;   //!
   TBranch        *b_hiHFhitPlus;   //!
   TBranch        *b_hiHFhitMinus;   //!
   TBranch        *b_hiET;   //!
   TBranch        *b_hiEE;   //!
   TBranch        *b_hiEB;   //!
   TBranch        *b_hiEEplus;   //!
   TBranch        *b_hiEEminus;   //!
   TBranch        *b_hiNpix;   //!
   TBranch        *b_hiNpixelTracks;   //!
   TBranch        *b_hiNtracks;   //!
   TBranch        *b_hiNtracksPtCut;   //!
   TBranch        *b_hiNtracksEtaCut;   //!
   TBranch        *b_hiNtracksEtaPtCut;   //!

   // List of branches from t.h
   TBranch        *b_evt;   //!
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
   TBranch        *b_discr_ssvHighEff;   //!
   TBranch        *b_discr_ssvHighPur;   //!
   TBranch        *b_discr_csvV1;   //!
   TBranch        *b_discr_csvV2;   //!
   TBranch        *b_discr_muByIp3;   //!
   TBranch        *b_discr_muByPt;   //!
   TBranch        *b_discr_prob;   //!
   TBranch        *b_discr_probb;   //!
   TBranch        *b_discr_tcHighEff;   //!
   TBranch        *b_discr_tcHighPur;   //!
   TBranch        *b_ndiscr_ssvHighEff;   //!
   TBranch        *b_ndiscr_ssvHighPur;   //!
   TBranch        *b_ndiscr_csvV1;   //!
   TBranch        *b_ndiscr_csvV2;   //!
   TBranch        *b_ndiscr_muByPt;   //!
   TBranch        *b_pdiscr_csvV1;   //!
   TBranch        *b_pdiscr_csvV2;   //!
   TBranch        *b_nsvtx;   //!
   TBranch        *b_svtxntrk;   //!
   TBranch        *b_svtxdl;   //!
   TBranch        *b_svtxdls;   //!
   TBranch        *b_svtxdl2d;   //!
   TBranch        *b_svtxdls2d;   //!
   TBranch        *b_svtxm;   //!
   TBranch        *b_svtxpt;   //!
   TBranch        *b_svtxmcorr;   //!
   TBranch        *b_nIPtrk;   //!
   TBranch        *b_nselIPtrk;   //!
   TBranch        *b_mue;   //!
   TBranch        *b_mupt;   //!
   TBranch        *b_mueta;   //!
   TBranch        *b_muphi;   //!
   TBranch        *b_mudr;   //!
   TBranch        *b_muptrel;   //!
   TBranch        *b_muchg;   //!

   TBranch        *b_beamId1;   //!
   TBranch        *b_beamId2;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refm;   //!
   TBranch        *b_refarea;   //!
   TBranch        *b_reftau1;   //!
   TBranch        *b_reftau2;   //!
   TBranch        *b_reftau3;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_subid;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_gentau1;   //!
   TBranch        *b_gentau2;   //!
   TBranch        *b_gentau3;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_genm;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!
   TBranch        *b_gensubid;   //!

   myTree(Bool_t mc = false, Bool_t pr = true);

   virtual ~myTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     EffCalc();
   virtual void     ANEffPlots();
   virtual void     Plot();
   virtual void     ClosureTest();
   virtual void     Unfolding();
   virtual void     JetPtRange();
   virtual void     EffSyst(int jtPtRange=0);
   virtual void     EffStatToy(int nToys=100);
   virtual void     EffStat(int jtPtRange=0);
   virtual void     TnpToy(int min=0, int max=100);
   virtual void     TnpStat(int jtPtRange=0);
   virtual void     EffMisMod(int jtPtRange=0);
   virtual vector<double> readSyst(const char* systfile);
   virtual void     FullEffSyst(int jtPtRange = 0);
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
};

#endif

#ifdef myTree_cxx
myTree::myTree(Bool_t mc = false , Bool_t pr = true) : fChain(0)
{
  //const char* input_file;
  TFile* f(0x0);
  isMc=mc;
  isPr=pr;
    if (isMc)
      {
	if (isPr)
	  f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/v9_ext/merged_HiForestAOD.root");
	else
	  f= (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/v9_ext/merged_HiForestAOD.root");
      }
    else 
  f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data_CMS/cms/mnguyen/jPsiJet/data/v6/merged_HiForestAOD.root");

    if (!f || !f->IsOpen()) 
      {
	if (isMc)
	  {
	    if (isPr)
	      f = new TFile("/data_CMS/cms/mnguyen/jPsiJet/mc/prompt/v9_ext/merged_HiForestAOD.root");
	    else
	      f= new TFile("/data_CMS/cms/mnguyen/jPsiJet/mc/nonprompt/v9_ext/merged_HiForestAOD.root");
	  }
	else 
	  f = new TFile("/data_CMS/cms/mnguyen/jPsiJet/data/v6/merged_HiForestAOD.root");
      }
    if(isMc)
      {
	if(isPr)
	  cout<<endl<<"-----prompt MC-----"<<endl;
	else
	  cout<<endl<<"-----non prompt MC-----"<<endl;
      }
    else
      {
	if (isPr)
	  cout<<endl<<"-----data with pr corrections-----"<<endl;
	else
	  cout<<endl<<"-----data with nonpr corrections-----"<<endl;

      }
      TTree * tree (0x0);
      TTree *hitree (0x0);
      TTree *hlttree(0x0);
      TTree *hlttree1(0x0);
      TTree *tr(0x0);
      TDirectory * dir = (TDirectory*)f->Get("hionia");
      dir->GetObject("myTree",tree);
      if(!tree)
	cout <<"error in the onia tree"<<endl;

      TDirectory * dirhi = (TDirectory*)f->Get("hiEvtAnalyzer");
      dirhi ->GetObject("HiTree",hitree);
 if(!hitree)
   {cout <<"error in the hitree"<<endl;}

      TDirectory * dirhlt = (TDirectory*)f->Get("hltanalysis");
      dirhlt ->GetObject("HltTree",hlttree);
 if(!hlttree)
   {cout <<"error in the hlttree"<<endl;}

      TDirectory * dirhlt1 = (TDirectory*)f->Get("skimanalysis");
      dirhlt1->GetObject("HltTree",hlttree1);
 if(!hlttree1)
   {cout <<"error in the hlttree1"<<endl;}
      
      TDirectory * dirt = (TDirectory*)f->Get("ak4PFJetAnalyzer");
      dirt->GetObject("t",tr);
 if(!tr)
   {cout <<"error in the t"<<endl;}

      tree -> AddFriend(hitree);
      tree -> AddFriend(hlttree);
      tree -> AddFriend(hlttree1);
      tree ->AddFriend(tr);

      Init(tree);
}

myTree::~myTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t myTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t myTree::LoadTree(Long64_t entry)
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

void myTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   fChain->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   fChain->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   fChain->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
   fChain->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
   fChain->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
   fChain->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
   fChain->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
   fChain->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
   fChain->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
   fChain->SetBranchAddress("hiET", &hiET, &b_hiET);
   fChain->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
   fChain->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
   fChain->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
   fChain->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
   fChain->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
   fChain->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
   fChain->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
   fChain->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
   fChain->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
   fChain->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("hiHF", &hiHF, &b_hiHF);
   fChain->SetBranchAddress("hiHFplus", &hiHFplus, &b_hiHFplus);
   fChain->SetBranchAddress("hiHFminus", &hiHFminus, &b_hiHFminus);
   fChain->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4, &b_hiHFplusEta4);
   fChain->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4, &b_hiHFminusEta4);
   fChain->SetBranchAddress("hiZDC", &hiZDC, &b_hiZDC);
   fChain->SetBranchAddress("hiZDCplus", &hiZDCplus, &b_hiZDCplus);
   fChain->SetBranchAddress("hiZDCminus", &hiZDCminus, &b_hiZDCminus);
   fChain->SetBranchAddress("hiHFhit", &hiHFhit, &b_hiHFhit);
   fChain->SetBranchAddress("hiHFhitPlus", &hiHFhitPlus, &b_hiHFhitPlus);
   fChain->SetBranchAddress("hiHFhitMinus", &hiHFhitMinus, &b_hiHFhitMinus);
   fChain->SetBranchAddress("hiET", &hiET, &b_hiET);
   fChain->SetBranchAddress("hiEE", &hiEE, &b_hiEE);
   fChain->SetBranchAddress("hiEB", &hiEB, &b_hiEB);
   fChain->SetBranchAddress("hiEEplus", &hiEEplus, &b_hiEEplus);
   fChain->SetBranchAddress("hiEEminus", &hiEEminus, &b_hiEEminus);
   fChain->SetBranchAddress("hiNpix", &hiNpix, &b_hiNpix);
   fChain->SetBranchAddress("hiNpixelTracks", &hiNpixelTracks, &b_hiNpixelTracks);
   fChain->SetBranchAddress("hiNtracks", &hiNtracks, &b_hiNtracks);
   fChain->SetBranchAddress("hiNtracksPtCut", &hiNtracksPtCut, &b_hiNtracksPtCut);
   fChain->SetBranchAddress("hiNtracksEtaCut", &hiNtracksEtaCut, &b_hiNtracksEtaCut);
   fChain->SetBranchAddress("hiNtracksEtaPtCut", &hiNtracksEtaPtCut, &b_hiNtracksEtaPtCut);

   fChain->SetBranchAddress("Onia2MuMuPAT", &Onia2MuMuPAT, &b_Onia2MuMuPAT);
   fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
   fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
   fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
   fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
   fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
   fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
   fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
   fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
   fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
   fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
   fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
   fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
   fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
   fChain->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
   fChain->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
   fChain->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
   fChain->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
   fChain->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
   fChain->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
   fChain->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
   fChain->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
   fChain->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
   fChain->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
   fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
   fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
   fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
   fChain->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts, &b_discr_jetID_cuts);
   fChain->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt, &b_discr_jetID_bdt);
   fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
   fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
   fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
   fChain->SetBranchAddress("trackN", trackN, &b_trackN);
   fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
   fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
   fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
   fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
   fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("photonN", photonN, &b_photonN);
   fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
   fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
   fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
   fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
   fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
   fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
   fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
   fChain->SetBranchAddress("eMax", eMax, &b_eMax);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("eN", eN, &b_eN);
   fChain->SetBranchAddress("muMax", muMax, &b_muMax);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("muN", muN, &b_muN);
   fChain->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff, &b_discr_ssvHighEff);
   fChain->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur, &b_discr_ssvHighPur);
   fChain->SetBranchAddress("discr_csvV1", discr_csvV1, &b_discr_csvV1);
   fChain->SetBranchAddress("discr_csvV2", discr_csvV2, &b_discr_csvV2);
   fChain->SetBranchAddress("discr_muByIp3", discr_muByIp3, &b_discr_muByIp3);
   fChain->SetBranchAddress("discr_muByPt", discr_muByPt, &b_discr_muByPt);
   fChain->SetBranchAddress("discr_prob", discr_prob, &b_discr_prob);
   fChain->SetBranchAddress("discr_probb", discr_probb, &b_discr_probb);
   fChain->SetBranchAddress("discr_tcHighEff", discr_tcHighEff, &b_discr_tcHighEff);
   fChain->SetBranchAddress("discr_tcHighPur", discr_tcHighPur, &b_discr_tcHighPur);
   fChain->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff, &b_ndiscr_ssvHighEff);
   fChain->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur, &b_ndiscr_ssvHighPur);
   fChain->SetBranchAddress("ndiscr_csvV1", ndiscr_csvV1, &b_ndiscr_csvV1);
   fChain->SetBranchAddress("ndiscr_csvV2", ndiscr_csvV2, &b_ndiscr_csvV2);
   fChain->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt, &b_ndiscr_muByPt);
   fChain->SetBranchAddress("pdiscr_csvV1", pdiscr_csvV1, &b_pdiscr_csvV1);
   fChain->SetBranchAddress("pdiscr_csvV2", pdiscr_csvV2, &b_pdiscr_csvV2);
   fChain->SetBranchAddress("nsvtx", nsvtx, &b_nsvtx);
   fChain->SetBranchAddress("svtxntrk", svtxntrk, &b_svtxntrk);
   fChain->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
   fChain->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
   fChain->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
   fChain->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
   fChain->SetBranchAddress("svtxm", svtxm, &b_svtxm);
   fChain->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
   fChain->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
   fChain->SetBranchAddress("nIPtrk", nIPtrk, &b_nIPtrk);
   fChain->SetBranchAddress("nselIPtrk", nselIPtrk, &b_nselIPtrk);
   fChain->SetBranchAddress("mue", mue, &b_mue);
   fChain->SetBranchAddress("mupt", mupt, &b_mupt);
   fChain->SetBranchAddress("mueta", mueta, &b_mueta);
   fChain->SetBranchAddress("muphi", muphi, &b_muphi);
   fChain->SetBranchAddress("mudr", mudr, &b_mudr);
   fChain->SetBranchAddress("muptrel", muptrel, &b_muptrel);
   fChain->SetBranchAddress("muchg", muchg, &b_muchg);



    if (fChain->GetBranch("HLT_HIL1DoubleMu0_v1")) fChain->SetBranchAddress("HLT_HIL1DoubleMu0_v1", &HLT_HIL1DoubleMu0_v1, &b_HLT_HIL1DoubleMu0_v1);


       if (fChain->GetBranch("HLT_HIL1DoubleMu0ForPPRef_v1")) fChain->SetBranchAddress("HLT_HIL1DoubleMu0ForPPRef_v1", &HLT_HIL1DoubleMu0ForPPRef_v1, &b_HLT_HIL1DoubleMu0ForPPRef_v1);



   ////// gen branches
       if (isMc)
	 {	
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

	   fChain->SetBranchAddress("beamId1", &beamId1, &b_beamId1);
	   fChain->SetBranchAddress("beamId2", &beamId2, &b_beamId2);
	   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
	   fChain->SetBranchAddress("refpt", refpt, &b_refpt);
	   fChain->SetBranchAddress("refeta", refeta, &b_refeta);
	   fChain->SetBranchAddress("refy", refy, &b_refy);
	   fChain->SetBranchAddress("refphi", refphi, &b_refphi);
	   fChain->SetBranchAddress("refm", refm, &b_refm);
	   fChain->SetBranchAddress("refarea", refarea, &b_refarea);
	   fChain->SetBranchAddress("reftau1", reftau1, &b_reftau1);
	   fChain->SetBranchAddress("reftau2", reftau2, &b_reftau2);
	   fChain->SetBranchAddress("reftau3", reftau3, &b_reftau3);
	   fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
	   fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
	   fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
	   fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
	   fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
	   fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
	   fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
	   fChain->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
	   fChain->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
	   fChain->SetBranchAddress("subid", subid, &b_subid);
	   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
	   fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
	   fChain->SetBranchAddress("genpt", genpt, &b_genpt);
	   fChain->SetBranchAddress("geneta", geneta, &b_geneta);
	   fChain->SetBranchAddress("geny", geny, &b_geny);
	   fChain->SetBranchAddress("gentau1", gentau1, &b_gentau1);
	   fChain->SetBranchAddress("gentau2", gentau2, &b_gentau2);
	   fChain->SetBranchAddress("gentau3", gentau3, &b_gentau3);
	   fChain->SetBranchAddress("genphi", genphi, &b_genphi);
	   fChain->SetBranchAddress("genm", genm, &b_genm);
	   fChain->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
	   fChain->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);
	   fChain->SetBranchAddress("gensubid", gensubid, &b_gensubid);

	   fChain->SetBranchAddress("ProcessID", &ProcessID, &b_ProcessID);
	   //fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
	   fChain->SetBranchAddress("weight", &weight, &b_weight);
	   fChain->SetBranchAddress("alphaQCD", &alphaQCD, &b_alphaQCD);
	   fChain->SetBranchAddress("alphaQED", &alphaQED, &b_alphaQED);
	   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
	   fChain->SetBranchAddress("nMEPartons", &nMEPartons, &b_nMEPartons);
	   fChain->SetBranchAddress("nMEPartonsFiltered", &nMEPartonsFiltered, &b_nMEPartonsFiltered);
	   fChain->SetBranchAddress("first", &first, &b_pdfID_first);
	   fChain->SetBranchAddress("second", &second, &b_pdfID_second);
	   //    fChain->SetBranchAddress("first", &first, &b_pdfX_first);
	   //    fChain->SetBranchAddress("second", &second, &b_pdfX_second);
	   //    fChain->SetBranchAddress("first", &first, &b_pdfXpdf_first);
	   //    fChain->SetBranchAddress("second", &second, &b_pdfXpdf_second);
	   fChain->SetBranchAddress("ttbar_w", &ttbar_w, &b_ttbar_w);
	   fChain->SetBranchAddress("npus", &npus, &b_npus);
	   fChain->SetBranchAddress("tnpus", &tnpus, &b_tnpus);
	 }

   Notify();
}

Bool_t myTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void myTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t myTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef myTree_cxx
