#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TBranchElement.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;

void compute(bool doPrompt = true, bool doMid = true){

  string filename = "";
  string outputfile = "";

  if(doPrompt && doMid){
    cout << "prompt mid" << endl;
    filename = "../../unfOutput/matrixOper/matrixOperation_Prompt_Mid_newNominal.root";
    outputfile = "../../unfOutput/matrixOper/systUnc_Promt_Mid_newNominal.root";
  }

  if(doPrompt && !doMid){
    cout << "prompt forward" << endl;
    filename = "../../unfOutput/matrixOper/matrixOperation_Prompt_Fwd_newNominal.root";
    outputfile = "../../unfOutput/matrixOper/systUnc_Promt_Fwd_newNominal.root";
  }

  if(!doPrompt && doMid){
    cout << "nonprompt mid" << endl;
    filename = "../../unfOutput/matrixOper/matrixOperation_NonPrompt_Mid.root";
    outputfile = "../../unfOutput/matrixOper/systUnc_NonPromt_Mid.root";
  }

  if(!doPrompt && !doMid){
    cout << "nonprompt forward" << endl;
    filename = "../../unfOutput/matrixOper/matrixOperation_NonPrompt_Fwd.root";
    outputfile = "../../unfOutput/matrixOper/systUnc_NonPromt_Fwd.root";
  }
  
  
  TFile *file = new TFile(filename.c_str());

  TH1D *h1_zUnf = (TH1D*)file->Get("nominalZUnf");

  if(doMid) h1_zUnf->Rebin(7);
  else h1_zUnf->Rebin(10);
  
  TH1D* h1_zUnfToys[100];

  float sums[7];

  
  sums[0] = 0.;
  sums[1] = 0.;
  sums[2] = 0.;
  sums[3] = 0.;
  sums[4] = 0.;
  sums[5] = 0.;
  sums[6] = 0.;
  
  
  float errs[7];

  
  errs[0] = 0.;
  errs[1]= 0.;
  errs[2]= 0.;
  errs[3]= 0.;
  errs[4]= 0.;
  errs[5]= 0.;
  errs[6]= 0.;
  
  
  
  float relErrs[7];

  relErrs[0] = 0.;
  relErrs[1] = 0.;
  relErrs[2] = 0.;
  relErrs[3] = 0.;
  relErrs[4] = 0.;
  relErrs[5] = 0.;
  relErrs[6] = 0.;
  
  
  int nToys = 100;
  
  for(int i = 0; i < nToys; i++){

    h1_zUnfToys[i] = (TH1D*)file->Get(Form("zUnfSmear_toy%i",i+1));

    if(doMid) h1_zUnfToys[i]->Rebin(7);
    else h1_zUnfToys[i]->Rebin(10);
    
    float nominalBinVal = 0.;
    float toyBinVal = 0.;
    
    for(int ibin = 0; ibin < h1_zUnf->GetNbinsX(); ibin++){

      nominalBinVal = h1_zUnf->GetBinContent(ibin+1);
      toyBinVal = h1_zUnfToys[i]->GetBinContent(ibin+1);

      sums[ibin] += (nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal);

      //cout << "i = " << i << " ibin = " << ibin << endl;
      //cout << "(nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal) = " << (nominalBinVal-toyBinVal)*(nominalBinVal-toyBinVal) << endl;
    }
    
    
  }

  TH1D *h1_zUnf_newSyst = (TH1D*)h1_zUnf->Clone("h1_zUnf_newSyst");
  h1_zUnf_newSyst->Reset();

  int nBins = 5;
  if(doMid) nBins = 7;
  
  for(int ibin = 0; ibin < nBins; ibin++){

    errs[ibin] = TMath::Sqrt(sums[ibin]*1.0/nToys);
    if(h1_zUnf->GetBinContent(ibin+1)>0) relErrs[ibin] = errs[ibin]/(h1_zUnf->GetBinContent(ibin+1));
    else relErrs[ibin] = 0.;

    cout << "z bin = " << ibin+1 << " , err = " << errs[ibin] << " , rel error (in %) = " << relErrs[ibin]*100 << endl;
    h1_zUnf_newSyst->SetBinContent(ibin+1,h1_zUnf->GetBinContent(ibin+1));
    h1_zUnf_newSyst->SetBinError(ibin+1,errs[ibin]);
    
  }

  
  TCanvas * can0 = new TCanvas("can0","can0",600,600);
  h1_zUnf_newSyst->Draw("EP");
  can0->SaveAs("../../plots/trMatrixSystCheck.png");
  

  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  h1_zUnf_newSyst->Write("zUnf_trMatrixSyst");
  outfile->Close();
  
  file->Close();
}

// syst unc from MC transfer matrix stat limitation

void compSystUnc_newBins(){

  //compute(true,true);
  //compute(true,false);

  //compute(false,true);
  compute(false,false);

    
    
}
