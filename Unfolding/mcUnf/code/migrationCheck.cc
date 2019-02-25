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

TH2D* CorrelationHist (const TMatrixT<double> cov, const char* name, const char* title, Double_t lo, Double_t hi)
{
  Int_t nb= cov.GetNrows();
  TH2D* h= new TH2D (name, title, nb, lo, hi, nb, lo, hi);
  h->SetAxisRange (-1.0, 1.0, "Z");
  for(int i=0; i < nb; i++)
    for(int j=0; j < nb; j++) {
      Double_t Viijj= cov(i,i)*cov(j,j);
      if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }
  return h;
}

void compute(bool doPrompt = true, bool doMid = true){

  string filename = "";
  string outputfile = "";

  if(doPrompt && doMid){
    cout <<"prompt mid" << endl; 
    filename = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter.root";
    outputfile = "../unfOutput/quarkoniaSystErr/statErrs_Prompt_Mid.root";
  }

  if(!doPrompt && doMid){
    cout <<"nonprompt mid"<< endl;
    filename = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter.root";
    outputfile = "../unfOutput/quarkoniaSystErr/statErrs_NonPrompt_Mid.root";
  }

    
  TFile *file = new TFile(filename.c_str());

  // get forward transfer matrix
  TH2D *h_fwdMatrix  = (TH2D *)file->Get("trmatrix");
  h_fwdMatrix->Sumw2();

  // 2D histo for its errors
  TH2D *h_fwdMatrixErr  = (TH2D *)h_fwdMatrix->Clone("h_fwdMatrixErr");
  
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get("invtrmat3");
  TMatrixT<double> *covMatrix  = (TMatrixT<double> *)file->Get("covmat3");
  
  //unfolded
  TH2D *h2_UnfData  = (TH2D*)file->Get("hReco_Iter3");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  h_invMatrix->Sumw2();
  
  TH2D *h_invMatrixErr = (TH2D *) h_invMatrix->Clone("h_invMatrixErr");

  TH2D *h_covMatrix = new TH2D(*covMatrix);
  h_covMatrix->Sumw2();  

  TH2D *h_corrHist ;
  h_corrHist = CorrelationHist(*covMatrix, "corr", "Unfolded correlation matrix", 0, 15);
  h_corrHist->Draw("colz");
  
  TH2D *h_invMatrixOut = (TH2D*)h_invMatrix->Clone("h_invMatrixOut");
  h_invMatrixOut->Reset();

  int nBinJet = 3;
  int nBinZ = 5;

  //normalize 2d histogram, such that each measured bin integrates to unity
  for(int jjet = 0; jjet < nBinJet; jjet++){
    for(int jz = 0; jz < nBinZ; jz++){
      int binY = 5*jjet+jz+1;

      float valIntegral = 0;
      float errIntegral = 0;

      float valFwdIntegral = 0;
      float errFwdIntegral = 0;
            
      for(int ijet = 0; ijet < nBinJet; ijet++){
	for(int iz = 0; iz < nBinZ; iz++){
	  int binX = 5*ijet+iz+1;

	  valIntegral += h_invMatrix->GetBinContent(binY,binX);
	  //errIntegral += h_covMatrix->GetBinContent(binY,binX);
	  errIntegral += h_invMatrixErr->GetBinError(binY,binX);
	  
	}
      }

      for(int ijet = 0; ijet < nBinJet; ijet++){
	for(int iz = 0; iz < nBinZ; iz++){
	  int binX = 5*ijet+iz+1;

	  if(valIntegral>0)h_invMatrix->SetBinContent(binY,binX, h_invMatrix->GetBinContent(binY,binX)/valIntegral);
	  //if(errIntegral>0)h_invMatrixErr->SetBinContent(binY,binX, h_covMatrix->GetBinContent(binY,binX)/errIntegral);
	  if(errIntegral>0)h_invMatrixErr->SetBinContent(binY,binX, h_invMatrixErr->GetBinError(binY,binX)/errIntegral);
	  
	}
      }

    }
  }


  //multiply back measured to see if we can get unfolded
  
  TH2D *h2_UnfDataFromMatrix = (TH2D*)h2_UnfData->Clone("h2_UnfDataFromMatrix");
  h2_UnfDataFromMatrix->Reset();
    
  for(int jjet = 0; jjet < nBinJet; jjet++){
    for(int jz = 0; jz < nBinZ; jz++){
      int binY = 5*jjet+jz+1;

      float measVal = h2_MeasData->GetBinContent(jz+1,jjet+1);
      float measErr = h2_MeasData->GetBinError(jz+1,jjet+1);

      double sumErr = 0.;
      
      for(int ijet = 0; ijet < nBinJet; ijet++){
	for(int iz = 0; iz < nBinZ; iz++){
	  int binX = 5*ijet+iz+1;

	  if(binX == binY) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  if(binX < binY) sumErr+=2*h_corrHist->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  //if(binX < binY) sumErr+=2*h_covMatrix->GetBinContent(binY,binX);
	  
	  //if(binX == binY) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  //if(binX == binY) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  //if(binX < binY) sumErr+=2*h_corrHist->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  	    
	  //sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX)*(1+h_corrHist->GetBinContent(binY,binX));
	  //sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  //sumErr+=h_covMatrix->GetBinContent(binY,binX)*h_covMatrix->GetBinContent(binY,binX)*h_corrHist->GetBinContent(binY,binX);
	  //sumErr+=h_covMatrix->GetBinContent(binY,binX)*h_covMatrix->GetBinContent(binY,binX);
	  //}
	}
      }
      
      if(sumErr>0){
	sumErr = sqrt(sumErr);
      }

      for(int ijet = 0; ijet < nBinJet; ijet++){
	for(int iz = 0; iz < nBinZ; iz++){
	  int binX = 5*ijet+iz+1;


	  float valWeight = h_invMatrix->GetBinContent(binY,binX);
	  float errWeight = h_invMatrixErr->GetBinContent(binY,binX);
	  float inVal = h2_UnfDataFromMatrix->GetBinContent(iz+1,ijet+1);
	  float inErr = h2_UnfDataFromMatrix->GetBinError(iz+1,ijet+1);

	  h2_UnfDataFromMatrix->SetBinContent(iz+1,ijet+1,inVal+valWeight*measVal);
	  h_invMatrixOut->SetBinContent(binY,binX,valWeight*measVal);

	  float outErr = 0.;
	  if(sumErr>0)outErr=errWeight*measErr/sumErr;
	  h2_UnfDataFromMatrix->SetBinError(iz+1,ijet+1,sqrt(outErr*outErr+inErr*inErr));
	  h_invMatrixOut->SetBinError(binY,binX,outErr);

	  /*	  
	  if(iz==nBinZ-1 && ijet==1) {
	    cout <<"contribution to the bin with iz = " << iz << " and ijet = " << ijet << endl;
	    cout << "from the measured bin with jz = " << jz << " and jjet = " << jjet << endl;
	    cout << " >> measVal = " << measVal << " and measErr = " << measErr << " sumErr = " << sumErr << endl;
	    cout <<  ">> valWeight = " << valWeight << " valWeight*measVal = " << valWeight*measVal <<  " , inErr = " << inErr << " , outErr = " << outErr <<  " and final error = " << sqrt(outErr*outErr+inErr*inErr) << endl;
	    cout << " **********************" << endl;
	  }
	  */
	  
	  
	  
	}
      }

    }
  }

  
}

void migrationCheck(){

  compute(true,true);
  //compute(true,false);

 
}

