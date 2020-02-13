// header with input
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
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"

using namespace std;

///////////////////////////////////////////////
//           general parameters              //
///////////////////////////////////////////////
bool sameSample = false;
bool mc2015 = false;
bool flatPrior = true;
int nIter = 3;
int nSIter = 99;
int nSIter_pp = 6;

double normPP = 3.4058383e-09;
double normPbPb = 1.5292e-08;

///////////////////////////////////////////////
//              systematics                  //
///////////////////////////////////////////////
int centShift = 0; // 0 is nominal, -1 for systDown and +1 for systUp

int JESsyst = 0; // 0 is nominal, -1 for systDown and +1 for systUp
int JERsyst = 0; // 0 is nominal, -1 for systDown and +1 for systUp
double SF = 1.1; // 1.1 is nominal, 1.0 for systDown, 1.2 for systUp

string systTag = "";

///////////////////////////////////////////////
//             jet parameters                //
///////////////////////////////////////////////
float jetR = 0.3;
//float jetR = 0.4;
float min_z = 0.064; 
float max_z = 1.0;

float min_jetpt = 10.;
float max_jetpt = 60.;

int nBinZ_gen = 48;
int nBinZ_reco = 6;
  
int nBinJet_gen = 25;
int nBinJet_reco = 5;

float jetPt_gen_binWidth = (max_jetpt-min_jetpt)*1.0/nBinJet_gen;
float z_gen_binWidth = (max_z-min_z)*1.0/nBinZ_gen;
float jetPt_reco_binWidth = (max_jetpt-min_jetpt)*1.0/nBinJet_reco;
float z_reco_binWidth = (max_z-min_z)*1.0/nBinZ_reco;

float midLowerPt = (max_jetpt+min_jetpt-jetPt_reco_binWidth)*0.5;
float midUpperPt = (max_jetpt+min_jetpt+jetPt_reco_binWidth)*0.5;

int midLowerId = ( (int) (nBinJet_reco/2) )*(nBinJet_gen/nBinJet_reco);

float unfStart = 0.22; // where we separate underflow from measured

///////////////////////////////////////////////
//            jpsi parameters                //
///////////////////////////////////////////////
float min_jp_pt = 6.5;
float max_jp_pt = 100.;//100. or 35.
float max_jp_eta = 2.4;//2.4 or 1.6

///////////////////////////////////////////////
//            style parameters               //
///////////////////////////////////////////////
TColor *pal = new TColor();
// good for primary marker colors
Int_t kmagenta = pal->GetColor(124,  0,124);
Int_t kviolet  = pal->GetColor( 72,  0,190);
Int_t kblue    = pal->GetColor(  9,  0,200);
Int_t kazure   = pal->GetColor(  0, 48, 97);
Int_t kcyan    = pal->GetColor(  0, 83, 98);
Int_t kteal    = pal->GetColor(  0, 92, 46);
Int_t kgreen   = pal->GetColor( 15, 85, 15);
Int_t kspring  = pal->GetColor( 75, 97, 53);
Int_t kyellow  = pal->GetColor(117,118,  0);
Int_t korange  = pal->GetColor(101, 42,  0);
Int_t kred     = pal->GetColor(190,  0,  3);
Int_t kpink    = pal->GetColor(180, 35,145);
// good for systematic band fill
Int_t kmagentaLight = pal->GetColor(215,165,215);
Int_t kvioletLight  = pal->GetColor(200,160,255);
Int_t kblueLight    = pal->GetColor(178,185,254);
Int_t kazureLight   = pal->GetColor(153,195,225);
Int_t kcyanLight    = pal->GetColor(140,209,224);
Int_t ktealLight    = pal->GetColor( 92,217,141);
Int_t kgreenLight   = pal->GetColor(135,222,135);
Int_t kspringLight  = pal->GetColor(151,207,116);
Int_t kyellowLight  = pal->GetColor(225,225,100);
Int_t korangeLight  = pal->GetColor(255,168,104);
Int_t kredLight     = pal->GetColor(253,169,179);
Int_t kpinkLight    = pal->GetColor(255,192,224);

Int_t col[] = {kBlack, kblue, kGreen+2, kCyan+2, kred, kGray+1};
int markerStyle[] = {kFullSquare,kFullCircle,kOpenTriangleUp,kOpenTriangleDown,kOpenCross,kOpenCrossX};
int markerSize[] = {1,1,1,1,1,1};
int lineStyle[] = {1, 1, 7, 7, 8, 7};
int lineWidth[] = {1,1,1,1,1,1};

///////////////////////////////////////////////
//               functions                   //
///////////////////////////////////////////////
void printInput() {
  cout << "############### Unfolding input ###############"<<endl;
  cout << "### "<<min_z<<" < z < "<<max_z<<endl;
  cout << "### "<<min_jetpt<<" < jetpt < "<<max_jetpt<<endl;
  cout << "### Reco bins: "<<nBinZ_reco<<" z bins and "<<nBinJet_reco<<" jtpt bins"<<endl;
  cout << "### Gen bins: "<<nBinZ_gen<<" z bins and "<<nBinJet_gen<<" jtpt bins"<<endl;
  cout << "### jetPt_gen_binWidth = "<<jetPt_gen_binWidth<<endl;
  cout << "### z_gen_binWidth = "<<z_gen_binWidth<<endl;
  cout << "### jetPt_reco_binWidth = "<<jetPt_reco_binWidth<<endl;
  cout << "### z_reco_binWidth = "<<z_reco_binWidth<<endl;
  cout << "### midLowerPt = "<<midLowerPt<<endl;
  cout << "### midUpperPt = "<<midUpperPt<<endl;
  cout << "### midLowerId = "<<midLowerId<<endl;
  cout << "### unfStart = "<<unfStart<<endl;
  cout << "###############################################"<<endl;
}

bool setSystTag () {
  systTag = "";

  if (JESsyst!=0 && (JERsyst!=0 || SF!=1.1 || centShift!=0) ) { cout <<"[ERROR] please check your systematic options"<<endl; return false;}
  if (JERsyst!=0 && (JESsyst!=0 || SF!=1.1 || centShift!=0) ) { cout <<"[ERROR] please check your systematic options"<<endl; return false;}
  if (SF!=1.1 && (JERsyst!=0 || JESsyst!=0 || centShift!=0) ) { cout <<"[ERROR] please check your systematic options"<<endl; return false;}
  if (centShift!=0 && (JERsyst!=0 || JESsyst!=0 || SF!=1.1) ) { cout <<"[ERROR] please check your systematic options"<<endl; return false;}
  
  if (JESsyst == 1) systTag = systTag + "_JESSystUp";
  else if (JESsyst == -1) systTag = systTag + "_JESSystDown";
  
  if (JERsyst == 1) systTag = systTag + "_JERSystUp";
  else if (JERsyst == -1) systTag = systTag + "_JERSystDown";
  
  if (SF == 1.2) systTag = systTag + "_SFSystUp";
  else if (SF == 1.0) systTag = systTag + "_SFSystDown";
  
  if (systTag == "") systTag = "_nominal";
  
  if (centShift == 1) systTag = systTag + "_centShiftSystUp";
  else if (centShift == -1) systTag = systTag + "_centShiftSystDown";
  
  return true;
}
