#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif

using namespace std;

void create(bool doPrompt = true, bool doMid = true, Int_t stepNumber = 1, double SF = 1.1){
  
  //#ifdef __CINT__
  //gSystem->Load("/home/ikucher/newRooUnfoldVersion/RooUnfold/libRooUnfold");
  //#endif
  
  string inputName = "";
  string outputName = "";
  string partOfOutput = "response";

  string SF_name = "";
  if(SF == 1.1) SF_name = "_nominal";
  if(SF == 1.2) SF_name = "_up";
  if(SF == 1.0) SF_name= "_down";

  cout << "step # =" << stepNumber << endl;
    
  if(doPrompt && doMid) {
    inputName = Form("dataUnfNewMidBins/unfInput/step%i/unfolding_4D_prompt_midRapidity_49z15ptBins7zMeasBins%s.root",stepNumber,SF_name.c_str());
  }
  if(doPrompt && !doMid) {
    inputName = Form("dataUnfNewMidBins/unfInput/step%i/unfolding_4D_prompt_fwdRapidity_50z15ptBins%s.root",stepNumber,SF_name.c_str());
  }
  if(!doPrompt && doMid) {
    inputName = Form("dataUnfNewMidBins/unfInput/step%i/unfolding_4D_nonprompt_midRapidity_49z15ptBins7zMeasBins%s.root",stepNumber,SF_name.c_str());
  }
  if(!doPrompt && !doMid) {
    inputName = Form("dataUnfNewMidBins/unfInput/step%i/unfolding_4D_nonprompt_fwdRapidity_50z15ptBins%s.root",stepNumber,SF_name.c_str());
  }

  outputName = inputName;
  outputName.replace(33,9,partOfOutput);

  cout << "outputName = " << outputName << endl;
  
  TFile *f = new TFile(inputName.c_str());
  f->ls();

  //Get response

  string thnSparseName = "";
  //take normalized tr matrix 
  thnSparseName = "hs_newJetPtNorm;1";

  
  THnSparseF *hn = static_cast<THnSparseF*>(f->Get(thnSparseName.c_str()));
  //hn->Sumw2();
  
  TH1F * h_z_gen = (TH1F*)f->Get("h_z_gen;1");
  TH1F * hGenZJetPtCentBin = (TH1F*)f->Get("hGenZJetPtCentBin;1");
  TH1F * hGenZJetPtLowBin = (TH1F*)f->Get("hGenZJetPtLowBin;1");
  TH1F * hGenZJetPtHighBin = (TH1F*)f->Get("hGenZJetPtHighBin;1");

  TH2D * h_trMatrix = (TH2D*)f->Get("h_trMatrix;1");
  
  const Int_t ndim = 4;
  Int_t dim[ndim];
  for(Int_t i = 0; i<ndim; i++)
    dim[i] = i;

  Int_t nDim = hn->GetNdimensions();
  cout <<"nDim = " << nDim << endl;
  
  Int_t iPtTrue   = 0;
  Int_t iZTrue  = 1;
  Int_t iPtDet  = 2;
  Int_t iZDet = 3;

  TH2D *fh2Smear = dynamic_cast<TH2D*>(hn->Projection(2,3,"E"));
  TH2D *fh2Prior = dynamic_cast<TH2D*>(hn->Projection(0,1,"E"));

  Int_t nBinPt[2] = {3,15};
  Double_t ptmin[2] = {15.0,15.0};
  Double_t ptmax[2] = {45.0,45.0};

  /*
  Int_t nBinZ[2] = {5,50};
  Double_t mmin[2] = {0.,0.};
  Double_t mmax[2] = {1.0,1.0};
  */

  int n_zBins = 5;
  if(doMid) n_zBins = 7;

  double z_min = 0.;
  if(doMid) z_min = 0.02;
  double z_max = 1.0;

  int n_zGenBins = 50;
  if(doMid) n_zGenBins = 49;

  Int_t nBinZ[2] = {n_zBins,n_zGenBins};
  Double_t mmin[2] = {z_min,z_min};
  Double_t mmax[2] = {z_max,z_max};
    
  //dimensions of measured axis
  TH2D *fh2RespDimM = new TH2D("fh2RespDimM","fh2RespDimM",nBinZ[0],mmin[0],mmax[0],nBinPt[0],ptmin[0],ptmax[0]);
  //dimensions of true axis
  TH2D *fh2RespDimT = new TH2D("fh2RespDimT","fh2RespDimT",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);  
  //feed-out of response
  TH2D *fh2Miss     = new TH2D("fh2Miss","fh2Miss",nBinZ[1],mmin[1],mmax[1],nBinPt[1],ptmin[1],ptmax[1]);
  cout << "fh2Smear->GetEntries() " << fh2Smear->GetEntries() << endl;

  //fill detector-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimM->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimM->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimM->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Smear->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Smear->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimM->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimM->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimM->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Smear->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Smear->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Smear->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      //cout << "ix = " << ix << " , iy = " << iy << " jxlow = " << jxlow << " jxup = " << jxup << " jylow = " << jylow << " jyup = " << jyup << endl;
      //cout << "con = " << con << endl;
      fh2RespDimM->SetBinContent(ix,iy,con);
      fh2RespDimM->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimM");

  //fill particle-level distribution
  for(Int_t ix = 1; ix<=fh2RespDimT->GetNbinsX(); ix++) {
    Double_t xlow = fh2RespDimT->GetXaxis()->GetBinLowEdge(ix);
    Double_t xup = fh2RespDimT->GetXaxis()->GetBinUpEdge(ix);
    Int_t jxlow = fh2Prior->GetXaxis()->FindBin(xlow+0.000001);
    Int_t jxup = fh2Prior->GetXaxis()->FindBin(xup-0.000001);
    for(Int_t iy = 1; iy<=fh2RespDimT->GetNbinsY(); iy++) {
      Double_t ylow = fh2RespDimT->GetYaxis()->GetBinLowEdge(iy);
      Double_t yup = fh2RespDimT->GetYaxis()->GetBinUpEdge(iy);
      Int_t jylow = fh2Prior->GetYaxis()->FindBin(ylow+0.000001);
      Int_t jyup = fh2Prior->GetYaxis()->FindBin(yup-0.000001);

      Double_t err = 0.;
      Double_t con = fh2Prior->IntegralAndError(jxlow,jxup,jylow,jyup,err);
      //cout << "ix = " << ix << " , iy = " << iy << endl;
      //cout << "con = " << con << endl;
      fh2RespDimT->SetBinContent(ix,iy,con);
      fh2RespDimT->SetBinError(ix,iy,err);
    }
  }
  Printf("Created fh2RespDimT");

  //response object for RooUnfold
  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(fh2RespDimM,fh2RespDimT);

  //Fill RooUnfoldResponse object
  
  Int_t* coord = new Int_t[nDim];
  //  Int_t nbin = fhnSparseReduced->GetNbins();
  Int_t nbin = hn->GetNbins();
  
  //  cout << "nbin = " << nbin << endl;

  /*
  cout << " fhnSparseReduced->GetAxis(0)->GetXmin() = " << fhnSparseReduced->GetAxis(0)->GetXmin() << endl;
  cout << " fhnSparseReduced->GetAxis(1)->GetXmin() = "<< fhnSparseReduced->GetAxis(1)->GetXmin() << endl;
  cout << " fhnSparseReduced->GetAxis(2)->GetXmin() = "<< fhnSparseReduced->GetAxis(2)->GetXmin() << endl;
  cout << " fhnSparseReduced->GetAxis(3)->GetXmin() = "<< fhnSparseReduced->GetAxis(3)->GetXmin() << endl;

  cout << " fhnSparseReduced->GetAxis(0)->GetXmax() = "<< fhnSparseReduced->GetAxis(0)->GetXmax() << endl;
  cout << " fhnSparseReduced->GetAxis(1)->GetXmax() = "<< fhnSparseReduced->GetAxis(1)->GetXmax() << endl;
  cout << " fhnSparseReduced->GetAxis(2)->GetXmax() = "<< fhnSparseReduced->GetAxis(2)->GetXmax() << endl;
  cout << " fhnSparseReduced->GetAxis(3)->GetXmax() = "<< fhnSparseReduced->GetAxis(3)->GetXmax() << endl;
  */
  
  for(Int_t bin=0; bin<nbin; bin++) {

    //cout << "bin = " << bin << endl;

    /*
    Double_t w = fhnSparseReduced->GetBinContent(bin,coord);
    Double_t pttrue = fhnSparseReduced->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t ztrue = fhnSparseReduced->GetAxis(1)->GetBinCenter(coord[1]);
    Double_t ptdet = fhnSparseReduced->GetAxis(2)->GetBinCenter(coord[2]);
    Double_t zdet = fhnSparseReduced->GetAxis(3)->GetBinCenter(coord[3]);
    */
    
    Double_t w = hn->GetBinContent(bin,coord);
    Double_t pttrue = hn->GetAxis(0)->GetBinCenter(coord[0]);
    Double_t ztrue = hn->GetAxis(1)->GetBinCenter(coord[1]);
    Double_t ptdet = hn->GetAxis(2)->GetBinCenter(coord[2]);
    Double_t zdet = hn->GetAxis(3)->GetBinCenter(coord[3]);
    
    
    /*
    cout << "all : " << endl;
    
    cout << "pttrue = " << pttrue << endl;
    cout << "ztrue = " << ztrue << endl;
    cout << "ptdet = " << ptdet << endl;
    cout << "zdet = "<< zdet<< endl;
    cout << "w = " <<w << endl;
    */
    
    if(zdet>=mmin[0] && zdet<=mmax[0]
       && ztrue>=mmin[1] && ztrue<=mmax[1]
       && ptdet>=ptmin[0] && ptdet<=ptmax[0]
       && pttrue>=ptmin[1] && pttrue<=ptmax[1]
       ){
         fResponse->Fill(zdet,ptdet,ztrue,pttrue,w);
    } 
    else {

      /*
      cout << "failed conditions : " << endl;
      
      cout << "bin = " << bin << endl;
      
      cout << "coord[0] = " << coord[0] << endl;
      cout << "coord[1] = " << coord[1] << endl;
      cout << "coord[2] = " << coord[2] << endl;
      cout << "coord[3] = " << coord[3] << endl;
      
      cout << "pttrue = " << pttrue << endl;
      cout << "ztrue = " << ztrue << endl;
      cout << "ptdet = " << ptdet << endl;
      cout << "zdet = " << zdet << endl;
      cout << "w = " << w << endl;

      cout << "******************" << endl;
      */
      
      fResponse->Miss(ztrue,pttrue,w);
      fh2Miss->Fill(ztrue,pttrue,w);
    }
  }

  delete [] coord;

  //Write response + 2D histos to file
  TFile *fout = new TFile(outputName.c_str(),"RECREATE");
  hn->Write("fhn");
  fResponse->Write("resp");
  fh2Smear->Write("fh2Smear");
  fh2Prior->Write("fh2Prior");
  fh2RespDimM->Write();
  fh2RespDimT->Write();
  fh2Miss->Write();

  h_trMatrix->Write();
  h_z_gen->Write();
  hGenZJetPtCentBin->Write();
  hGenZJetPtLowBin->Write();
  hGenZJetPtHighBin->Write();
    
  fout->Write();
  fout->Close();
}

void createRooUnfoldResponse(Int_t step = 1, double SF = 1.1){

  create(true,true,step,SF);
  create(true,false,step,SF);

  create(false,true,step,SF);
  create(false,false,step,SF);
  
}

  
