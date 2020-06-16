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

double rms(vector<double> v, bool isrelative) {
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

int getCentFromName(TString sName) {
  if (sName.Contains("cent0")) return 0;
  for (int i=2 ; i<200 ;i++)
    if (sName.Contains(Form("cent%d",i)) && !sName.Contains(Form("cent%d",i*10))) return i;
  return 200;
}

void testAccEffStat(bool isPr, bool isPbPb, bool isAcc, int dimension, string caseTag) {

  Double_t ptbinsAna []= {6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 14.0, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 50.0};
  int nBin = sizeof(ptbinsAna)/sizeof(double)-1;

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s_NoBkg_jetR3_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR",isPbPb?"PbPb":"PP"),"READ");
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw; int centr;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  trNom->SetBranchAddress("corr_ptw",&corr_ptw);
  trNom->SetBranchAddress("centr",&centr);

  TEfficiency* tempPrEff = NULL;
  TEfficiency* tempPrAcc = NULL;
  TEfficiency* tempNprEff = NULL;
  TObjArray* tempCentArr = NULL;

  TObjArray* hisArr = new TObjArray(); //hisArr1624->SetOwner(kTRUE);
  TObjArray* zhisArr = new TObjArray(); //hisArr1624->SetOwner(kTRUE);

  TFile *prcorrFile = TFile::Open(Form("FilesAccxEff%s/toyMC/pr%s100Toys_%s%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?(isAcc?"PP":"PbPb"):"PP",(dimension==1)?"_1D":""));
  TObjArray* prcorArr = (TObjArray*) prcorrFile->Get("accEffArray");

  TH1F* histVar = NULL;
  TH2D* zhist2D = NULL;//new TH2D("zhist2D",";z;p_{T,jet}",6,0.064,1,5,10,60);

  TFile* nomFile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str()),"READ");
  TEfficiency* prnomEff = (TEfficiency*) nomFile->Get(Form("hcorr_Jpsi_%s_pr%s_Eff",isPbPb?"PbPb":"PP",(dimension==1)?"_1D":""));
  TEfficiency* prnomAcc = (TEfficiency*) nomFile->Get(Form("hcorr_Jpsi_PP_pr%s_Acc",(dimension==1)?"_1D":""));

  int *centbins = new int[100];
  int ncentbins=0;

  TList* lcorr = nomFile->GetListOfKeys();
  TIter nextCorr(lcorr);

  TObjArray* corrCentNom = new TObjArray();
  corrCentNom->SetOwner(kTRUE);
  if (dimension==3) {
    TObjString* fname(0x0);
    while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
      {
	TEfficiency* h = static_cast<TEfficiency*>(nomFile->FindObjectAny(fname->GetString().Data()));
	
	TString sName(h->GetName());
	if ( sName.Contains("cent") ){
	  corrCentNom->Add(h->Clone());
	  centbins[ncentbins] = getCentFromName(sName);
	  ncentbins++;
	  cout <<"adding "<<sName<<" to the corrction array and adding "<<centbins[ncentbins-1]<<" to the centrality array"<<endl;
	}
      }
  }
  centbins[ncentbins] = 200;
  
  TObjArray* corrCentArr = new TObjArray();
  corrCentArr->SetOwner(kTRUE);
  corrCentArr = (TObjArray*) (prcorrFile->Get("accEffArray_centAll"));
  /*  
  for (int iCent=0; iCent<ncentbins; iCent++) {
    TObjArray* h = static_cast<TObjArray*> (prcorrFile->Get(Form("accEffArray_cent%d%d",centbins[iCent],centbins[iCent+1])));
    h->SetName(Form("accEffArray_cent%d%d",centbins[iCent],centbins[iCent+1]));
    corrCentArr->Add(h->Clone());
  }
  */

  double bf = 1.0;
  double prEff = 1.0;
  double nprEff = 1.0;
  double totEff = 1.0;

  int nentries = trNom->GetEntries();
  //nentries = 10;
  for (int i=0; i<=100; i++)
    {
      histVar = new TH1F (Form("histTot_%d",i), "", nBin, ptbinsAna);
      zhist2D = new TH2D(Form("zhist2D_%d",i),";z;p_{T,jet}",6,0.064,1,5,10,60);

      if (i==0) cout<<"[INFO] Applying the nominal AccxEff"<<endl;
      else cout<<"[INFO] Applying toy "<<i<<"/100"<<endl;
      if (i==0) {
	tempPrEff=prnomEff;
	tempPrAcc=prnomAcc; 
      }
      else {
	if (isAcc){
	  tempPrAcc=(TEfficiency*) prcorArr->At(i-1);
	  tempPrEff=prnomEff;

	}
	else {
	  tempPrAcc=prnomAcc;
	  tempPrEff=(TEfficiency*) prcorArr->At(i-1);
	}
      }

      for (int jentry=0; jentry<nentries; jentry++) {
	trNom->GetEntry(jentry);
	if (abs(jp_rap)>2.4) continue;
	if (jp_pt<6.5 || jp_pt>50) continue;
	if (dimension==1)
	  prEff = (tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_pt)))*(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_pt)));
	else if (dimension==2)
	  prEff = (tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_rap,jp_pt)))*(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt)));
	else {
	  for (int iCent=0; iCent<ncentbins; iCent++) {
	    if (centr<=centbins[iCent+1]) {
	      //centTag = Form("cent%d%d",centbins[i],centbins[i+1]);
	      if (i==0) { 
		//tempCentArr = corrCentNom;
		tempPrEff = (TEfficiency*) corrCentNom->At(iCent);
	      }
	      else {
		//tempCentArr = (TObjArray*) corrCentArr;//->At(iCent);
		tempPrEff = (TEfficiency*) corrCentArr->FindObject(Form("effToy_cent%d%d_%d",centbins[iCent],centbins[iCent+1],(i-1)));
	      }
	      break;
	    }
	  }
	  prEff = (tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_rap,jp_pt)))*(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt)));
	  //cout <<"acc = "<<(tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_rap,jp_pt)))<<", eff = "<<(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt)))<<", for cent = "<<centr<<", pt = "<<jp_pt<<", rap = "<<jp_rap<<", reading efficiency from efficiency: "<<tempPrEff->GetName()<<" centArr = "<<tempCentArr->GetName()<<endl;
	}
	totEff=1.0/prEff;
	totEff = totEff*corr_ptw;
	
	histVar->Fill(jp_pt, totEff);
	if (z>=0 && z<=1)
	  zhist2D->Fill(z, jt_pt,totEff);
      }// end of the tree entries
      hisArr->Add(histVar);
      zhisArr->Add(zhist2D);
    } // end of the variations
  TFile* fsave = new TFile(Form("FilesAccxEff%s/ptHistsWith100%sToys_%s_%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"),"RECREATE");
  hisArr->Write("ptDistArr",TObject::kSingleKey);
  zhisArr->Write("zDistArr",TObject::kSingleKey);
  fsave->Close();
  
  cout<<"[INFO] making csv files"<<endl;
  ofstream fileOut(Form("../Fitter/Systematics/csv/syst_Raa_NJpsi_%s_%s_%sStat%s.csv", isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", isAcc?"Acc":"Eff", caseTag.c_str()));

  fileOut<<"AccxEff stat"<<endl;

  vector<double> vCount;

  vCount.clear();

  TH1F* temp = NULL;
  TH2D* ztemp = NULL;

  for (int j=0; j<nBin; j++)
    {
      vCount.clear();
      for (int i=0; i<=100; i++) {
	temp = (TH1F*) hisArr->At(i);
	vCount.push_back(temp->GetBinContent(temp->FindBin(ptbinsAna[j]+0.001)));
      }
      if (isPbPb)
        fileOut<< "0, 2.4, "<<ptbinsAna[j]<<", "<<ptbinsAna[j+1]<<", 0, 101, 0, 180, "<< rms(vCount,true) << endl;
      else
	fileOut<< "0, 2.4, "<<ptbinsAna[j]<<", "<<ptbinsAna[j+1]<<", 0, 101, 0, 200, "<< rms(vCount,true) << endl;
    }
  fileOut.close();

  ofstream zfileOut(Form("../Fitter/Systematics/csv/syst_midJtPt_NJpsi_%s_%s_%sStat%s.csv", isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", isAcc?"Acc":"Eff", caseTag.c_str()));
  zfileOut<<"AccxEff stat"<<endl;
  vCount.clear();

  int nZBin = 6;

  for (int j=0; j<nZBin; j++)
    {
      vCount.clear();
      for (int i=0; i<=100; i++) {
	ztemp = (TH2D*) zhisArr->At(i);
	vCount.push_back(ztemp->GetBinContent(ztemp->FindBin((0.064+j*0.156)+0.001,35+0.001)));
      }
      if (isPbPb)
        zfileOut<< "0, 2.4, 6.5, 100, "<<(0.064+j*0.156)<<", "<< (0.064+(j+1)*0.156)<<", 0, 180, "<< rms(vCount,true) << endl;
      else
        zfileOut<< "0, 2.4, 6.5, 100, "<<(0.064+j*0.156)<<", "<< (0.064+(j+1)*0.156)<<", 0, 200, "<< rms(vCount,true) << endl;
    }
  zfileOut.close();

  prcorrFile->Close();
  nomFile->Close();
  //nprcorrFile->Close();
  trFile->Close();
}
