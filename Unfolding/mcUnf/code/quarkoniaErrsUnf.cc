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

void compute(bool doPrompt = true, bool doMid = true){

  string filename = "";
  string dataSystErr = "";
  string outputfile = "";

  if(doPrompt && doMid){
    cout <<"prompt mid" << endl; 
    filename = "../unfOutput/step3/UnfoldedDistributions_Prompt_Mid_8iter.root";
    dataSystErr = "../unfOutput/quarkoniaSystErr/meas_data_prompt_mid_systErrs.root";
    outputfile = "../unfOutput/quarkoniaSystErr/systErrs_Prompt_Mid.root";
  }

  if(doPrompt && !doMid){
    cout <<"prompt fwd"<< endl;
    filename = "../unfOutput/step3/UnfoldedDistributions_Prompt_Fwd_8iter.root";
    dataSystErr = "../unfOutput/quarkoniaSystErr/meas_data_prompt_fwd_systErrs.root";
    outputfile = "../unfOutput/quarkoniaSystErr/systErrs_Prompt_Fwd.root";
  }

  if(!doPrompt && doMid){
    cout <<"nonprompt mid"<< endl;
    filename = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Mid_8iter.root";
    dataSystErr = "../unfOutput/quarkoniaSystErr/meas_data_nonprompt_mid_systErrs.root";
    outputfile = "../unfOutput/quarkoniaSystErr/systErrs_NonPrompt_Mid.root";
  }

  if(!doPrompt && !doMid){
    cout <<"nonprompt fwd"<< endl;
    filename = "../unfOutput/step3/UnfoldedDistributions_NonPrompt_Fwd_8iter.root";
    dataSystErr = "../unfOutput/quarkoniaSystErr/meas_data_nonprompt_fwd_systErrs.root";
    outputfile = "../unfOutput/quarkoniaSystErr/systErrs_NonPrompt_Fwd.root";
  }
    
  TFile *file = new TFile(filename.c_str());
  TFile *fileDataSyst = new TFile(dataSystErr.c_str());

  // get forward transfer matrix
  TH2D *h_fwdMatrix  = (TH2D *)file->Get("trmatrix");
  h_fwdMatrix->Sumw2();

  // 2D histo for its errors
  TH2D *h_fwdMatrixErr  = (TH2D *)h_fwdMatrix->Clone("h_fwdMatrixErr");
  
  // get inverse transfer matrix
  TMatrixT<double> *invMatrix  = (TMatrixT<double> *)file->Get("invtrmat3");

  //unfolded
  TH2D *h2_UnfData  = (TH2D*)file->Get("hReco_Iter3");

  //measured with syst uncertanties 
  TH2D *h2_MeasData  = (TH2D*)fileDataSyst->Get("h_Meas");
  
  TH2D *h_invMatrix = new TH2D(*invMatrix);
  h_invMatrix->Sumw2();
  
  TH2D *h_invMatrixErr = (TH2D *) h_invMatrix->Clone("h_invMatrixErr");

  //why?
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
	  errIntegral += h_invMatrixErr->GetBinError(binY,binX);

	  valFwdIntegral += h_fwdMatrix->GetBinContent(binY,binX);
	  errFwdIntegral += h_fwdMatrixErr->GetBinError(binY,binX);
	  
	  
	}
      }

      for(int ijet = 0; ijet < nBinJet; ijet++){
	for(int iz = 0; iz < nBinZ; iz++){
	  int binX = 5*ijet+iz+1;

	  if(valFwdIntegral>0)h_fwdMatrix->SetBinContent(binY,binX,h_fwdMatrix->GetBinContent(binY,binX)/valFwdIntegral);
	  if(errFwdIntegral>0)h_fwdMatrixErr->SetBinContent(binY,binX,h_fwdMatrixErr->GetBinError(binY,binX)/errFwdIntegral);
	    
	  if(valIntegral>0)h_invMatrix->SetBinContent(binY,binX, h_invMatrix->GetBinContent(binY,binX)/valIntegral);
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
	  if(h_invMatrixErr->GetBinContent(binY,binX)>0) sumErr+=h_invMatrixErr->GetBinContent(binY,binX)*h_invMatrixErr->GetBinContent(binY,binX);
	  /*
	  if(iz==nBinZ-1 && ijet==1) {
	    cout <<"contribution to the bin with iz = " << iz << " and ijet = " << ijet << "from the measured bin with jz = " << jz << " and jjet = " << jjet << endl;
	    cout << "measVal = " << measVal << " and measErr = " << measErr << " sumErr = " << sumErr << endl;
	  }
	  */
	  
	}
      }
      if(sumErr>0)sumErr = sqrt(sumErr);
      //if(jjet==0&&jz==1)cout<< " sumErr "<<sumErr<<endl;

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
	  //if(ijet==2&&iz==0)cout<<" inErr "<<inErr<<" outErr "<<outErr<<" errWeight "<<errWeight<<" measErr "<<measErr<<" product "<<errWeight*measErr/sumErr<<endl;
	  h2_UnfDataFromMatrix->SetBinError(iz+1,ijet+1,sqrt(outErr*outErr+inErr*inErr));
	  h_invMatrixOut->SetBinError(binY,binX,outErr);

	  if(iz==nBinZ-1 && ijet==1) {
	    cout <<"contribution to the bin with iz = " << iz << " and ijet = " << ijet << endl;
	    cout << "from the measured bin with jz = " << jz << " and jjet = " << jjet << endl;
	    cout << " >> measVal = " << measVal << " and measErr = " << measErr << " sumErr = " << sumErr << endl;
	    cout <<  ">> valWeight = " << valWeight << " valWeight*measVal = " << valWeight*measVal <<  " , inErr = " << inErr << " , outErr = " << outErr <<  " and final error = " << sqrt(outErr*outErr+inErr*inErr) << endl;
	    cout << " **********************" << endl;
	  }
	  
	  
	}
      }

    }
  }

  TH1D *h_quarkoniaErrUnf =  h2_UnfDataFromMatrix->ProjectionX("h_quarkoniaErrUnf",2,2);
  h_quarkoniaErrUnf->Draw("EP");

    
  TH1D *h_quarkoniaErrBeforeUnf = (TH1D*)fileDataSyst->Get("zMeasCentBin");

  /*
  TCanvas * can0 = new TCanvas("can0","can0",1200,600);
  can0->Divide(2,1);
  can0->cd(1);
  h_quarkoniaErrUnf->GetXaxis()->SetTitle("after errors are proparated through the unfolding");
  h_quarkoniaErrUnf->Draw("EP");
  
  can0->cd(2);
  h_quarkoniaErrBeforeUnf->GetXaxis()->SetTitle("before errors are proparated through the unfolding");
  h_quarkoniaErrBeforeUnf->Draw("EP");
  
  if(doPrompt && doMid) can0->SaveAs("../plots/quarkoniaSystErrorsTest_prompt_mid.png");
  if(doPrompt && !doMid) can0->SaveAs("../plots/quarkoniaSystErrorsTest_prompt_fwd.png");
  if(!doPrompt && doMid) can0->SaveAs("../plots/quarkoniaSystErrorsTest_nonprompt_mid.png");
  if(!doPrompt && !doMid) can0->SaveAs("../plots/quarkoniaSystErrorsTest_nonprompt_fwd.png");
  
  
  for(int iz = 0; iz < 5; iz++){
    //    cout << " before unf propagation err =  " << h_quarkoniaErrBeforeUnf->GetBinError(iz+1) << " , after err = " << h_quarkoniaErrUnf->GetBinError(iz+1) << endl;
    cout << "iz = " << iz << endl;
    float relErr1 = h_quarkoniaErrBeforeUnf->GetBinError(iz+1)/h_quarkoniaErrBeforeUnf->GetBinContent(iz+1);
    float relErr2 = h_quarkoniaErrUnf->GetBinError(iz+1)/h_quarkoniaErrUnf->GetBinContent(iz+1); 
    cout <<"rel error before (in %) = " << relErr1*100 << " after (in %) = " << relErr2*100 << endl;
    cout << "content before = " << h_quarkoniaErrBeforeUnf->GetBinContent(iz+1) << " after = " << h_quarkoniaErrUnf->GetBinContent(iz+1) << " , ratio = " << h_quarkoniaErrBeforeUnf->GetBinContent(iz+1)/h_quarkoniaErrUnf->GetBinContent(iz+1)  <<   endl;
    cout << "error before = " << h_quarkoniaErrBeforeUnf->GetBinError(iz+1) << " after = " << h_quarkoniaErrUnf->GetBinError(iz+1) << ", ratio = " << h_quarkoniaErrBeforeUnf->GetBinError(iz+1)/h_quarkoniaErrUnf->GetBinError(iz+1) << endl;
    cout << "******************" << endl;
  }
  */

  for(int iz = 0; iz < 5; iz++){

    float relErr1 = 0 ;
    if(h_quarkoniaErrBeforeUnf->GetBinContent(iz+1)>0) relErr1 = h_quarkoniaErrBeforeUnf->GetBinError(iz+1)/h_quarkoniaErrBeforeUnf->GetBinContent(iz+1);
    float relErr2 = 0;
    if(h_quarkoniaErrUnf->GetBinContent(iz+1)>0) relErr2 = h_quarkoniaErrUnf->GetBinError(iz+1)/h_quarkoniaErrUnf->GetBinContent(iz+1);
    
    cout << "iz = " << iz << endl;
    cout << "rel error (in %) : before unf = " << relErr1*100 << " , after stat errors propagation = " << relErr2*100 << endl;
    cout << "******************" << endl;
  }
  
  cout<< " printing ratio of input and output "<<endl;

  cout << "unfolded from the matrix multiplication :" << endl; 
  for(int ijet = 0; ijet < nBinJet; ijet++){
    for(int iz = 0; iz < nBinZ; iz++){
      cout<<" ijet "<<ijet<<" iz "<<iz<<" bin content = "<< h2_UnfDataFromMatrix->GetBinContent(iz+1,ijet+1)<<" +- "<<h2_UnfDataFromMatrix->GetBinError(iz+1,ijet+1)<<endl;
    }
  }

  cout << "original unfolded distribution with stat errors only:" << endl;
  for(int ijet = 0; ijet < nBinJet; ijet++){
    for(int iz = 0; iz < nBinZ; iz++){
      cout<<" ijet "<<ijet<<" iz "<<iz<<" bin content = "<<h2_UnfData->GetBinContent(iz+1,ijet+1)<<" +- "<<h2_UnfData->GetBinError(iz+1,ijet+1)<<endl;
    }
  }

  /*
  TFile *outfile = new TFile(outputfile.c_str(),"RECREATE");
  
  h_quarkoniaErrBeforeUnf->Write("zBeforeUnf");  
  h_quarkoniaErrUnf->Write("zUnf_wQuarkoniaSyst");
  
  outfile->Close();
  */
  
}

void quarkoniaErrsUnf(){

  compute(true,true);
  
  //  compute(true,false);
  //  compute(false,true);
  //  compute(false,false);
 
}
