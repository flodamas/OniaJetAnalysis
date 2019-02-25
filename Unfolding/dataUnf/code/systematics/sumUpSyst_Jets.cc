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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

using namespace std;

void add(bool doPrompt = true, bool doMid = true){

  string name_unfTrMatrix = "";
  string name_unfSI = "";
  string name_jets = "";
  string name_jets_add = "";
  
  string outputfile = "";

  if(doPrompt && doMid){
    cout <<"prompt mid" << endl;
    name_unfTrMatrix = "../../unfOutput/matrixOper/systUnc_Promt_Mid_newNominal.root";
    name_unfSI = "../../unfOutput/numberSIErrs/SIerr_Prompt_Mid_newNominal.root";
    name_jets = "../../unfOutput/JERJESErrs/JERJESerr_Prompt_Mid_newNominal.root";
    name_jets_add = "../../unfOutput/JERJESErrs/addJESerr_Prompt_Mid_newNominal.root";
    outputfile = "../../unfOutput/totalSystErr/systErrs_Prompt_Mid_Jets_newNominal.root";
  }

  if(doPrompt && !doMid){
    cout <<"prompt fwd"<< endl;
    name_unfTrMatrix = "../../unfOutput/matrixOper/systUnc_Promt_Fwd_newNominal.root";
    name_unfSI = "../../unfOutput/numberSIErrs/SIerr_Prompt_Fwd_newNominal.root";
    name_jets ="../../unfOutput/JERJESErrs/JERJESerr_Prompt_Fwd_newNominal.root";
    name_jets_add = "../../unfOutput/JERJESErrs/addJESerr_Prompt_Fwd_newNominal.root";
    outputfile = "../../unfOutput/totalSystErr/systErrs_Prompt_Fwd_Jets_newNominal.root";
  }

    
  TFile *file_unfTrMatrix = new TFile(name_unfTrMatrix.c_str());
  TFile *file_unfSI = new TFile(name_unfSI.c_str());
  TFile *file_jets = new TFile(name_jets.c_str());
  TFile *file_jets_add = new TFile(name_jets_add.c_str());
  
  TH1D * h_unfTrMatrix = (TH1D*)file_unfTrMatrix->Get("zUnf_trMatrixSyst");
  TH1D * h_unfSI = (TH1D*)file_unfSI->Get("zUnfSI3_newUncert_scaledErr");
  TH1D * h_JES = (TH1D*)file_jets->Get("zUnf_wJESErr");
  TH1D * h_JER = (TH1D*)file_jets->Get("zUnf_wJERErr");
  TH1D * h_JES_add = (TH1D*)file_jets_add->Get("zUnf_wAddJES");
  
  TH1D *h_totalSyst = (TH1D*) h_unfTrMatrix->Clone("h_totalSyst");
  h_totalSyst->Reset();
  
  int zBins = 5;
  if(doMid) zBins = 7;

  float binVal_unfTrMatrix = 0;
  float binErr_unfTrMatrix = 0;

  float binVal_unfSI = 0;
  float binErr_unfSI = 0;

  float binVal_JES = 0;
  float binErr_JES = 0;

  float binVal_JER = 0;
  float binErr_JER = 0;

  float finalBinCont = 0;
  float finalBinErr = 0;

  float binVal_JES_add = 0;
  float binErr_JES_add = 0;
    
  for(int i = 0; i < zBins; i++){


    binVal_unfTrMatrix = h_unfTrMatrix->GetBinContent(i+1);
    binErr_unfTrMatrix = h_unfTrMatrix->GetBinError(i+1);

    binVal_unfSI = h_unfSI->GetBinContent(i+1);
    binErr_unfSI = h_unfSI->GetBinError(i+1);

    binVal_JES = h_JES->GetBinContent(i+1);
    binErr_JES = h_JES->GetBinError(i+1);

    binVal_JES_add = h_JES_add->GetBinContent(i+1);
    binErr_JES_add = h_JES_add->GetBinError(i+1);
    
    
    binVal_JER = h_JER->GetBinContent(i+1);
    binErr_JER = h_JER->GetBinError(i+1);

    finalBinCont = binVal_unfTrMatrix;
    finalBinErr = TMath::Sqrt(binErr_unfTrMatrix*binErr_unfTrMatrix+binErr_unfSI*binErr_unfSI+binErr_JES*binErr_JES+binErr_JER*binErr_JER+binErr_JES_add*binErr_JES_add);

    h_totalSyst->SetBinContent(i+1,finalBinCont);
    h_totalSyst->SetBinError(i+1,finalBinErr);
    

    cout << "z bin = " << i << endl;

    cout << " bin content = " << binVal_unfTrMatrix  << " , " << binVal_unfSI << " , " << binVal_JES << " , " << binVal_JES_add <<   " , " << binVal_JER << endl;
    
    cout << "rel errors (in %) from :" << endl;
    cout << "unf tr matrix = " << binErr_unfTrMatrix*100/binVal_unfTrMatrix << endl;
    cout << "unf # of SI = " << binErr_unfSI*100/binVal_unfSI << endl;
    cout << "JES = " << binErr_JES*100/binVal_JES << endl;
    cout << "JES add = " << binErr_JES_add*100/binVal_JES_add << endl;
    cout << "JER = " << binErr_JER*100/binVal_JER << endl;
    cout << "TOTAL = " << finalBinErr*100/finalBinCont << endl;
    cout << "*****" << endl;
        
  }


  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");

  h_unfTrMatrix->Write("unfTrMatrixSyst");
  h_unfSI->Write("unfSISyst");
  h_JES->Write("JESSyst");
  h_JES_add->Write("JESAddSyst");
  h_JER->Write("JERSyst");
  h_totalSyst->Write("totalSyst");
  
  outfile->Close();

  
}

void sumUpSyst_Jets(){

  //add(true,true);
  add(true,false);
 
}
