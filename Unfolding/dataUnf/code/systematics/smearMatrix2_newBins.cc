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

void operate(bool doPrompt = true, bool doMid = true, int toyNumber = 1){

  string filename = "";
  string outputfile = "";

  if(doPrompt && doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter_49z15ptBins7zMeasBins_newNominal_nominal.root";
    outputfile = "../../unfOutput/matrixOper/matrixOperation_Prompt_Mid_newNominal.root";
  }

  if(doPrompt && !doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter_50z15ptBins_newNominal_nominal.root";
    outputfile = "../../unfOutput/matrixOper/matrixOperation_Prompt_Fwd_newNominal.root";
  }

  if(!doPrompt && doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter_49z15ptBins7zMeasBins_nominal.root";
    outputfile = "../../unfOutput/matrixOper/matrixOperation_NonPrompt_Mid.root";
  }

  if(!doPrompt && !doMid){
    filename = "../../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter_50z15ptBins_nominal.root";
    outputfile = "../../unfOutput/matrixOper/matrixOperation_NonPrompt_Fwd.root";
  }
  
   
  TFile *file = new TFile(filename.c_str());

  // get forward transfer matrix
  TH2D *h_fwdMatrix  = (TH2D *)file->Get("trmatrix");
  
  // 2D histo for its errors
  TH2D *h_fwdMatrixErr  = (TH2D *)h_fwdMatrix->Clone("h_fwdMatrixErr");
  
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get("invtrmat3");

  TH2D *h2_UnfData  = (TH2D*)file->Get("hReco_Iter3;1");
  TH2D *h2_MeasData  = (TH2D*)file->Get("fh2MeasData;1");

  TH1D *h1_zUnf = (TH1D*)h2_UnfData->ProjectionX("h1_zUnf",6,10);
  
  //  TH1D *h1_zUnf = (TH1D*)file->Get("hMUnf_1_Iter3;1");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  h_invMatrix->Sumw2();
  
  TH2D *h_invMatrixErr = (TH2D *) h_invMatrix->Clone("h_invMatrixErr");

  // 4D inv matrix, norm so its proj gives the truth
  TH2D *h_invMatrixOut = (TH2D*)h_invMatrix->Clone("h_invMatrixOut");
  h_invMatrixOut->Reset();

  int nBinJet_reco = 3;
  //int nBinZ_reco = 5;

  
  int nBinJet_gen = 15;
  //int nBinZ_gen = 50;
  
  
  int nBinZ_gen = 50;
  if(doMid) nBinZ_gen = 49;
  int nBinZ_reco = 5;
  if(doMid) nBinZ_reco = 7;

  
  cout << "here" << endl;
  
  //normalize 2d histogram, such that each measured bin integrates to unity
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){

      int binY = 5*jjet+jz+1;
      if(doMid) binY = 7*jjet+jz+1;
      
      float valIntegral = 0;
      float errIntegral = 0;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;
	  
	  valIntegral += h_invMatrix->GetBinContent(binY,binX);
	  errIntegral += h_invMatrixErr->GetBinError(binY,binX);
	  
	}
      }

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;
	  
	  if(valIntegral>0)h_invMatrix->SetBinContent(binY,binX, h_invMatrix->GetBinContent(binY,binX)/valIntegral);
	  if(errIntegral>0)h_invMatrixErr->SetBinContent(binY,binX, h_invMatrixErr->GetBinError(binY,binX)/errIntegral);

	}
      }

    }
  }

  //multiply back measured to see if we can get unfolded
  
  TH2D *h2_UnfDataFromMatrix = (TH2D*)h2_UnfData->Clone("h2_UnfDataFromMatrix");
  h2_UnfDataFromMatrix->Reset();
    
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = 5*jjet+jz+1;
      if(doMid) binY = 7*jjet+jz+1;

      float measVal = h2_MeasData->GetBinContent(jz+1,jjet+1);
      float measErr = h2_MeasData->GetBinError(jz+1,jjet+1);

      double sumErr = 0.;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;
	  
	  if(h_invMatrixErr->GetBinContent(binY,binX)>0) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	}
      }
      if(sumErr>0)sumErr = sqrt(sumErr);

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;
	  
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
	}
      }

    }
  }


  //fix the errors in h_invMatrixOut
  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = 5*jjet+jz+1;
      if(doMid) binY = 7*jjet+jz+1;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;
	  
	  float fwdMatixVal = h_fwdMatrix->GetBinContent(binY,binX);
	  float fwdMatrixErr = h_fwdMatrix->GetBinError(binY,binX);

	  if(fwdMatixVal>0)cout << "fwd matrix relative error = " << fwdMatrixErr/fwdMatixVal << endl;
	  if(fwdMatrixErr/fwdMatixVal>0.99)cout <<"fwd rel error is 100% in the bin with jjet = " << jjet << " , jz = " << jz << " , ijet = " << ijet << " , iz = " << iz << endl; 
	  
	  float invMatrixVal = h_invMatrixOut->GetBinContent(binY,binX);
	  cout << "invMatrixVal  = " << invMatrixVal << endl;
	  float invMatrixErr = 0;
	  if(fwdMatixVal>0) invMatrixErr = fwdMatrixErr*invMatrixVal/fwdMatixVal;

	  if(fwdMatixVal>0)cout << "inv matrix relative error = " << invMatrixErr/invMatrixVal << endl;
	  if(invMatrixErr/invMatrixVal>0.99)cout<<"fwd rel error is 100% in the bin with jjet = " << jjet << " , jz = " << jz << " , ijet = " << ijet << " , iz = " <<iz << endl;

	  h_invMatrixOut->SetBinError(binY,binX,invMatrixErr);  

	}
      }

    }
  }

  //project the nominal, take the unfolded nominal jet pt bin, loop over the unfolded z bins, sum over all measured bins.

  TH1D *h_zUnfNomBin = (TH1D *)h1_zUnf->Clone("h_zUnfNomBin");
  h_zUnfNomBin->Reset();
  float zBinsVals[nBinZ_gen];  
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    //int ijet = 1;
    zBinsVals[iz] = 0;
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      int binX = 50*ijet+iz+1;
      if(doMid) binX = 49*ijet+iz+1;
          
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  int binY = 5*jjet+jz+1;
	  if(doMid) binY = 7*jjet+jz+1;
	  
	  if(ijet>4 && ijet<10) {
	    zBinsVals[iz] += h_invMatrixOut->GetBinContent(binY,binX);
	    cout << "bin content to add = " << h_invMatrixOut->GetBinContent(binY,binX) << endl;
	  }
	  
	}
      }
    }
    
    cout << "iz+1 = " << iz+1 << " bin content = " << zBinsVals[iz] << endl;
    h_zUnfNomBin->SetBinContent(iz+1,zBinsVals[iz]);
    
  }

  //compare the unfolded result from the one obtained with the matrix projection
  
  for(int iz = 0; iz < nBinZ_gen; iz++){
    cout << "bin = " << iz << endl;
    cout << "projected bin content = "<< h_zUnfNomBin->GetBinContent(iz+1) << endl;
    cout << "unfolded bin content = "<< h1_zUnf->GetBinContent(iz+1) << endl;
  }

  /*
  TCanvas *can_nominal = new TCanvas("can_nominal","can_nominal",1200,600);
  can_nominal->Divide(2,1);
  can_nominal->cd(1);
  h1_zUnf->GetXaxis()->SetTitle("unfolded result");
  h1_zUnf->Draw();
  can_nominal->cd(2);
  h_zUnfNomBin->GetXaxis()->SetTitle("projection from the matrix");
  h_zUnfNomBin->Draw("EP");
  can_nominal->SaveAs("../plots/nominalBinCheck.png");
  */
  
  //smear h_invMatrixOut -> ok
  TH2D *h_invMatrixOutSmear = (TH2D*)h_invMatrixOut->Clone("h_invMatrixOutSmear");
  h_invMatrixOutSmear->Reset();
    
  TRandom3 * myRand = new TRandom3();
  myRand->SetSeed(0);

  for(int jjet = 0; jjet < nBinJet_reco; jjet++){
    for(int jz = 0; jz < nBinZ_reco; jz++){
      int binY = 5*jjet+jz+1;
      if(doMid) binY = 7*jjet+jz+1;
      
      float oldVal = 0;
      float errVal = 0;
      float smearVal = 0;

      for(int ijet = 0; ijet < nBinJet_gen; ijet++){
	for(int iz = 0; iz < nBinZ_gen; iz++){
	  int binX = 50*ijet+iz+1;
	  if(doMid) binX = 49*ijet+iz+1;

	  oldVal = h_invMatrixOut->GetBinContent(binY,binX);
	  errVal = h_invMatrixOut->GetBinError(binY,binX);
	  smearVal = myRand->Gaus(oldVal,errVal);

	  if(smearVal<0) {
	    smearVal = 0;
	  }
	  
	  h_invMatrixOutSmear->SetBinContent(binY,binX,smearVal);
	  
	}
      }

    }

  }

  //project the correct axis
  
  TH1D *h_zUnfNomBinSmear = (TH1D *)h1_zUnf->Clone("h_zUnfNomBinSmear");
  h_zUnfNomBinSmear->Reset();
  float zBinsValsSmear[nBinZ_gen];

  for(int iz = 0; iz < nBinZ_gen; iz++){
    zBinsValsSmear[iz] = 0;
    // int ijet = 1; // not 1 jet pt bin, but 5,6,7,8,9
    for(int ijet = 0; ijet < nBinJet_gen; ijet++){
      
      int binX = 50*ijet+iz+1;
      if(doMid) binX = 49*ijet+iz+1;

      //      zBinsValsSmear[iz] = 0;
      
      for(int jjet = 0; jjet < nBinJet_reco; jjet++){
	for(int jz = 0; jz < nBinZ_reco; jz++){
	  int binY = 5*jjet+jz+1;
	  if(doMid)  binY = 7*jjet+jz+1;
	  // do it only for the unf jt pt 25-35 GeV
	  if(ijet>4 && ijet<10)zBinsValsSmear[iz] += h_invMatrixOutSmear->GetBinContent(binY,binX);
	}
      }
      
    }
    
    h_zUnfNomBinSmear->SetBinContent(iz+1,zBinsValsSmear[iz]);
    
  }
  

  h_zUnfNomBinSmear->Draw("EP");
    
  
  TFile *outfile = new TFile(outputfile.c_str(),"UPDATE");
  if(toyNumber == 1) h1_zUnf->Write("nominalZUnf");
  h_zUnfNomBinSmear->Write(Form("zUnfSmear_toy%i",toyNumber));
  outfile->Close();
  
}

void smearMatrix2_newBins(){

  for(int i = 1; i < 101; i++){
    operate(true,true,i);
    operate(true,false,i);
    operate(false,true,i);
    operate(false,false,i);
  }
}
