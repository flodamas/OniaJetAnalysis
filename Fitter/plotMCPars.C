#define plotMCPars_cxx

#include "plotMCPars.h"

#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TLine.h"
#include "TMath.h"

using namespace std;

const char* fName = "";

struct paramList {
  double alpha_Jpsi_val;
  double alpha_Jpsi_err;
  double f_Jpsi_val;
  double f_Jpsi_err;
  double n_Jpsi_val;
  double n_Jpsi_err;
  double sigma1_Jpsi_val;
  double sigma1_Jpsi_err;
  double rSigma21_Jpsi_val;
  double rSigma21_Jpsi_err;
};

Double_t getWaverage(TGraphErrors* graph)
{
  Double_t num = 0;
  Double_t denom(0);
  int npoints = graph->GetN();
  Double_t x(0);
  Double_t y(0);
  Double_t err(0);
  for (int i = 0 ; i < npoints ; i++)
  {
    graph->GetPoint(i,x,y);
    err = graph->GetErrorY(i);
    
    num += y/err;
    denom += 1./err;
  }
  num /= denom;
  return num;
}

void plotParams(map< string, map<anabin, paramList> > paramFullMap, const char* fileName, const char* fName, bool isRBin, bool isPtBin, bool isCentBin, bool wMean);

void plotMCPars::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L plotMCPars.C+
  //      Root > plotMCPars t("inputFile",arguments)
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
//  gROOT->SetBatch(true);
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  const int nFunct = 2;
  const char* fNameList[nFunct] = {"DoubleCrystalBall","GaussianAndCrystalBall"};
  
  for (int i = 0 ; i < nFunct ; i++)
  {
    fName = fNameList[i];
    
    map< string, map<anabin, paramList> > paramFullMap;
    map<anabin, paramList> paramMapPP;
    map<anabin, paramList> paramMapPbPb;
    paramFullMap["PP"] = paramMapPP;
    paramFullMap["PbPb"] = paramMapPbPb;
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      map<anabin, paramList> paramMap = paramFullMap[collSystem];
      
      paramList list;
      
      anabin thebin(zmin, zmax, ymin, ymax, ptmin, ptmax, centmin, centmax);
      if (paramMap.count(thebin)) list = paramMap[thebin];
      
      if (strcmp(jpsiName,fName)) continue;
      
      list.alpha_Jpsi_val = alpha_Jpsi_val;
      list.alpha_Jpsi_err = alpha_Jpsi_err;
      list.f_Jpsi_val = f_Jpsi_val;
      list.f_Jpsi_err = f_Jpsi_err;
      list.n_Jpsi_val = n_Jpsi_val;
      list.n_Jpsi_err = n_Jpsi_err;
      list.sigma1_Jpsi_val = sigma1_Jpsi_val;
      list.sigma1_Jpsi_err = sigma1_Jpsi_err;
      list.rSigma21_Jpsi_val = rSigma21_Jpsi_val;
      list.rSigma21_Jpsi_err = rSigma21_Jpsi_err;
      
      paramMap[thebin] = list;
      paramFullMap[collSystem] = paramMap;
    }
    
    plotParams(paramFullMap,saveName.Data(),fName,isRBinning, isPBinning, isCBinning,wMean);
  }
}

void plotParams(map< string, map<anabin, paramList> > paramFullMap, const char* fileName, const char* fName, bool isRBin, bool isPtBin, bool isCBin, bool wMean)
{
  ///////////////////////
  // Create and fill graphs
  ///////////////////////
  
  TGraphErrors* gPPInt_alpha = new TGraphErrors();
  TGraphErrors* gPPPt_alpha = new TGraphErrors();
  TGraphErrors* gPPRap_alpha = new TGraphErrors();
  TGraphErrors* gPbPbInt_alpha = new TGraphErrors();
  TGraphErrors* gPbPbPt_alpha = new TGraphErrors();
  TGraphErrors* gPbPbCent_alpha = new TGraphErrors();
  TGraphErrors* gPbPbRap_alpha = new TGraphErrors();

  TGraphErrors* gPPInt_f = new TGraphErrors();
  TGraphErrors* gPPPt_f = new TGraphErrors();
  TGraphErrors* gPPRap_f = new TGraphErrors();
  TGraphErrors* gPbPbInt_f = new TGraphErrors();
  TGraphErrors* gPbPbPt_f = new TGraphErrors();
  TGraphErrors* gPbPbCent_f = new TGraphErrors();
  TGraphErrors* gPbPbRap_f = new TGraphErrors();

  TGraphErrors* gPPInt_n = new TGraphErrors();
  TGraphErrors* gPPPt_n = new TGraphErrors();
  TGraphErrors* gPPRap_n = new TGraphErrors();
  TGraphErrors* gPbPbInt_n = new TGraphErrors();
  TGraphErrors* gPbPbPt_n = new TGraphErrors();
  TGraphErrors* gPbPbCent_n = new TGraphErrors();
  TGraphErrors* gPbPbRap_n = new TGraphErrors();

  TGraphErrors* gPPInt_sigma = new TGraphErrors();
  TGraphErrors* gPPPt_sigma = new TGraphErrors();
  TGraphErrors* gPPRap_sigma = new TGraphErrors();
  TGraphErrors* gPbPbInt_sigma = new TGraphErrors();
  TGraphErrors* gPbPbPt_sigma = new TGraphErrors();
  TGraphErrors* gPbPbCent_sigma = new TGraphErrors();
  TGraphErrors* gPbPbRap_sigma = new TGraphErrors();

  TGraphErrors* gPPInt_rsigma = new TGraphErrors();
  TGraphErrors* gPPPt_rsigma = new TGraphErrors();
  TGraphErrors* gPPRap_rsigma = new TGraphErrors();
  TGraphErrors* gPbPbInt_rsigma = new TGraphErrors();
  TGraphErrors* gPbPbPt_rsigma = new TGraphErrors();
  TGraphErrors* gPbPbCent_rsigma = new TGraphErrors();
  TGraphErrors* gPbPbRap_rsigma = new TGraphErrors();

  int sysN = 2;
  const char* sysname[2] = {"PP", "PbPb"};
  
  TString rapRange("");
  TString ptRange("");
  TString centRange("");
  double rapMin(0.);
  double rapMax(0.);
  
  for ( int i = 0 ; i < sysN ; i++ )
  {
    int centbin = 0;
    int pTbin = 0;
    int Rapbin = 0;
    
    map<anabin, paramList> paramMap = paramFullMap[sysname[i]];
   
    map<anabin, paramList>::const_iterator itm;
    for (itm=paramMap.begin(); itm!=paramMap.end(); itm++)
    {
      anabin thebin = itm->first;
      paramList list = itm->second;
      
      if ((!strcmp(sysname[i],"PbPb")) && (isPtBin && !isCBin)) centRange = Form("%d-%d %%",(int)(thebin.centbin().low()/2.),(int)(thebin.centbin().high()/2.));
      
      if ( ((thebin.ptbin().low() == 6.5) || (thebin.ptbin().low() == 3.0)) && ((thebin.ptbin().high() == 30.0) || (thebin.ptbin().high() == 50.0) || (thebin.ptbin().high() == 6.5))) //pt integrated
      {
        // Note that there is a pT bin 3.0-6.5 that we count as integrated since we study the cent dependence in that bin. Be careful if in an eta range the pT bin 3-6.5 also exists
        if ( ((thebin.centbin().low() == 0.) && (thebin.centbin().high() == 200.)) || ((thebin.centbin().low() == -1.) && (thebin.centbin().high() == 1000.)) || (!isRBin && isPtBin && !isCBin) ) //pt-cent integrated
        {
          rapMin = thebin.rapbin().low();
          rapMax = thebin.rapbin().high();
          
          if (!isRBin)
          {
            rapRange = Form("%.1f<|y|<%.1f",thebin.rapbin().low(),thebin.rapbin().high());
            ptRange = Form("%.1f<p_{T}<%.1f",thebin.ptbin().low(),thebin.ptbin().high());
          }
          else
          {
            ptRange = Form("%.1f<p_{T}<%.1f",thebin.ptbin().low(),thebin.ptbin().high());
          }
          
          if (isRBin)
          {
            if ((thebin.rapbin().low() == 0) && (thebin.rapbin().high() > 2.39)) // We need the fully integrated bin (pt-cent-rap) for this to work
            {
              if (!strcmp(sysname[i],"PP"))
              {
                gPPInt_alpha->SetPoint(0,0.,list.alpha_Jpsi_val);
                gPPInt_alpha->SetPointError(0,0.,list.alpha_Jpsi_err);
                
                gPPInt_f->SetPoint(0,0.,list.f_Jpsi_val);
                gPPInt_f->SetPointError(0,0.,list.f_Jpsi_err);
                
                gPPInt_n->SetPoint(0,0.,list.n_Jpsi_val);
                gPPInt_n->SetPointError(0,0.,list.n_Jpsi_err);
                
                gPPInt_sigma->SetPoint(0,0.,list.sigma1_Jpsi_val);
                gPPInt_sigma->SetPointError(0,0.,list.sigma1_Jpsi_err);
                
                gPPInt_rsigma->SetPoint(0,0.,list.rSigma21_Jpsi_val);
                gPPInt_rsigma->SetPointError(0,0.,list.rSigma21_Jpsi_err);
              }
              if (!strcmp(sysname[i],"PbPb"))
              {
                gPbPbInt_alpha->SetPoint(0,0.,list.alpha_Jpsi_val);
                gPbPbInt_alpha->SetPointError(0,0.,list.alpha_Jpsi_err);
                
                gPbPbInt_f->SetPoint(0,0.,list.f_Jpsi_val);
                gPbPbInt_f->SetPointError(0,0.,list.f_Jpsi_err);
                
                gPbPbInt_n->SetPoint(0,0.,list.n_Jpsi_val);
                gPbPbInt_n->SetPointError(0,0.,list.n_Jpsi_err);
                
                gPbPbInt_sigma->SetPoint(0,0.,list.sigma1_Jpsi_val);
                gPbPbInt_sigma->SetPointError(0,0.,list.sigma1_Jpsi_err);
                
                gPbPbInt_rsigma->SetPoint(0,0.,list.rSigma21_Jpsi_val);
                gPbPbInt_rsigma->SetPointError(0,0.,list.rSigma21_Jpsi_err);
              }
              
            }
            else
            {
              double bin = (thebin.rapbin().high() + thebin.rapbin().low())/2.;
              double binSize = thebin.rapbin().high() - bin;
              
              if (!strcmp(sysname[i],"PP"))
              {
                gPPRap_alpha->SetPoint(Rapbin,bin,list.alpha_Jpsi_val);
                gPPRap_alpha->SetPointError(Rapbin,binSize,list.alpha_Jpsi_err);
                
                gPPRap_f->SetPoint(Rapbin,bin,list.f_Jpsi_val);
                gPPRap_f->SetPointError(Rapbin,binSize,list.f_Jpsi_err);
                
                gPPRap_n->SetPoint(Rapbin,bin,list.n_Jpsi_val);
                gPPRap_n->SetPointError(Rapbin,binSize,list.n_Jpsi_err);
                
                gPPRap_sigma->SetPoint(Rapbin,bin,list.sigma1_Jpsi_val);
                gPPRap_sigma->SetPointError(Rapbin,binSize,list.sigma1_Jpsi_err);
                
                gPPRap_rsigma->SetPoint(Rapbin,bin,list.rSigma21_Jpsi_val);
                gPPRap_rsigma->SetPointError(Rapbin,binSize,list.rSigma21_Jpsi_err);
              }
              if (!strcmp(sysname[i],"PbPb"))
              {
                gPbPbRap_alpha->SetPoint(Rapbin,bin,list.alpha_Jpsi_val);
                gPbPbRap_alpha->SetPointError(Rapbin,binSize,list.alpha_Jpsi_err);
                
                gPbPbRap_f->SetPoint(Rapbin,bin,list.f_Jpsi_val);
                gPbPbRap_f->SetPointError(Rapbin,binSize,list.f_Jpsi_err);
                
                gPbPbRap_n->SetPoint(Rapbin,bin,list.n_Jpsi_val);
                gPbPbRap_n->SetPointError(Rapbin,binSize,list.n_Jpsi_err);
                
                gPbPbRap_sigma->SetPoint(Rapbin,bin,list.sigma1_Jpsi_val);
                gPbPbRap_sigma->SetPointError(Rapbin,binSize,list.sigma1_Jpsi_err);
                
                gPbPbRap_rsigma->SetPoint(Rapbin,bin,list.rSigma21_Jpsi_val);
                gPbPbRap_rsigma->SetPointError(Rapbin,binSize,list.rSigma21_Jpsi_err);
              }
              
              Rapbin++;
            }
          }
          else
          {
            if (!strcmp(sysname[i],"PP"))
            {
              gPPInt_alpha->SetPoint(0,0.,list.alpha_Jpsi_val);
              gPPInt_alpha->SetPointError(0,0.,list.alpha_Jpsi_err);
              
              gPPInt_f->SetPoint(0,0.,list.f_Jpsi_val);
              gPPInt_f->SetPointError(0,0.,list.f_Jpsi_err);
              
              gPPInt_n->SetPoint(0,0.,list.n_Jpsi_val);
              gPPInt_n->SetPointError(0,0.,list.n_Jpsi_err);
              
              gPPInt_sigma->SetPoint(0,0.,list.sigma1_Jpsi_val);
              gPPInt_sigma->SetPointError(0,0.,list.sigma1_Jpsi_err);
              
              gPPInt_rsigma->SetPoint(0,0.,list.rSigma21_Jpsi_val);
              gPPInt_rsigma->SetPointError(0,0.,list.rSigma21_Jpsi_err);
            }
            if (!strcmp(sysname[i],"PbPb"))
            {
              gPbPbInt_alpha->SetPoint(0,0.,list.alpha_Jpsi_val);
              gPbPbInt_alpha->SetPointError(0,0.,list.alpha_Jpsi_err);
              
              gPbPbInt_f->SetPoint(0,0.,list.f_Jpsi_val);
              gPbPbInt_f->SetPointError(0,0.,list.f_Jpsi_err);
              
              gPbPbInt_n->SetPoint(0,0.,list.n_Jpsi_val);
              gPbPbInt_n->SetPointError(0,0.,list.n_Jpsi_err);
              
              gPbPbInt_sigma->SetPoint(0,0.,list.sigma1_Jpsi_val);
              gPbPbInt_sigma->SetPointError(0,0.,list.sigma1_Jpsi_err);
              
              gPbPbInt_rsigma->SetPoint(0,0.,list.rSigma21_Jpsi_val);
              gPbPbInt_rsigma->SetPointError(0,0.,list.rSigma21_Jpsi_err);
            }
          }
          
        }
        else
        {
          double bin = (thebin.centbin().high() + thebin.centbin().low())/2.;
          double binSize = thebin.centbin().high() - bin;
          
          gPbPbCent_alpha->SetPoint(centbin,bin,list.alpha_Jpsi_val);
          gPbPbCent_alpha->SetPointError(centbin,binSize,list.alpha_Jpsi_err);
          
          gPbPbCent_f->SetPoint(centbin,bin,list.f_Jpsi_val);
          gPbPbCent_f->SetPointError(centbin,binSize,list.f_Jpsi_err);
          
          gPbPbCent_n->SetPoint(centbin,bin,list.n_Jpsi_val);
          gPbPbCent_n->SetPointError(centbin,binSize,list.n_Jpsi_err);
          
          gPbPbCent_sigma->SetPoint(centbin,bin,list.sigma1_Jpsi_val);
          gPbPbCent_sigma->SetPointError(centbin,binSize,list.sigma1_Jpsi_err);
          
          gPbPbCent_rsigma->SetPoint(centbin,bin,list.rSigma21_Jpsi_val);
          gPbPbCent_rsigma->SetPointError(centbin,binSize,list.rSigma21_Jpsi_err);
          
          centbin++;
        }
      }
      else
      {
        double bin = (thebin.ptbin().high() + thebin.ptbin().low())/2.;
        double binSize = thebin.ptbin().high() - bin;
        
        if (!strcmp(sysname[i],"PP"))
        {
          gPPPt_alpha->SetPoint(pTbin,bin,list.alpha_Jpsi_val);
          gPPPt_alpha->SetPointError(pTbin,binSize,list.alpha_Jpsi_err);
          
          gPPPt_f->SetPoint(pTbin,bin,list.f_Jpsi_val);
          gPPPt_f->SetPointError(pTbin,binSize,list.f_Jpsi_err);
          
          gPPPt_n->SetPoint(pTbin,bin,list.n_Jpsi_val);
          gPPPt_n->SetPointError(pTbin,binSize,list.n_Jpsi_err);
          
          gPPPt_sigma->SetPoint(pTbin,bin,list.sigma1_Jpsi_val);
          gPPPt_sigma->SetPointError(pTbin,binSize,list.sigma1_Jpsi_err);
          
          gPPPt_rsigma->SetPoint(pTbin,bin,list.rSigma21_Jpsi_val);
          gPPPt_rsigma->SetPointError(pTbin,binSize,list.rSigma21_Jpsi_err);
        }
        if (!strcmp(sysname[i],"PbPb"))
        {
          gPbPbPt_alpha->SetPoint(pTbin,bin,list.alpha_Jpsi_val);
          gPbPbPt_alpha->SetPointError(pTbin,binSize,list.alpha_Jpsi_err);
          
          gPbPbPt_f->SetPoint(pTbin,bin,list.f_Jpsi_val);
          gPbPbPt_f->SetPointError(pTbin,binSize,list.f_Jpsi_err);
          
          gPbPbPt_n->SetPoint(pTbin,bin,list.n_Jpsi_val);
          gPbPbPt_n->SetPointError(pTbin,binSize,list.n_Jpsi_err);
          
          gPbPbPt_sigma->SetPoint(pTbin,bin,list.sigma1_Jpsi_val);
          gPbPbPt_sigma->SetPointError(pTbin,binSize,list.sigma1_Jpsi_err);
          
          gPbPbPt_rsigma->SetPoint(pTbin,bin,list.rSigma21_Jpsi_val);
          gPbPbPt_rsigma->SetPointError(pTbin,binSize,list.rSigma21_Jpsi_err);
        }

        pTbin++;
      }

    }
    
  }
  
  ///////////////////////
  // Plot and save canvases
  ///////////////////////
  TString sFileName(fileName);
  TString sPath(fileName);
  sPath.Remove(sPath.Last('/'),sPath.Sizeof());

  TObjArray* aSave = new TObjArray();
  aSave->SetOwner(true);
  
  double avPbPb(0.);
  double avPP(0.);
  double maxPbPb(0);
  double maxPP(0);
  double minPbPb(0);
  double minPP(0);
  double maxTot(0.);
  double minTot(0.);
  if (isRBin)
  {
    ///////////////////////////////
    /////////////Rap dep////////////
    ///////////////////////////////
    
    TCanvas* c1 = new TCanvas(Form("alphaVsRap_%s",fName),"canvas",81,98,700,504);
    c1->Range(-5.547577,-0.004024593,31.93896,0.03073326);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.1479885);
    c1->SetRightMargin(0.05172414);
    c1->SetTopMargin(0.08421053);
    c1->SetBottomMargin(0.1157895);
    c1->SetFrameBorderMode(0);
    c1->SetFrameBorderMode(0);
    
    maxPbPb = TMath::MaxElement(gPbPbRap_alpha->GetN(),gPbPbRap_alpha->GetY());
    minPbPb = TMath::MinElement(gPbPbRap_alpha->GetN(),gPbPbRap_alpha->GetY());
    
    minPP = TMath::MinElement(gPPRap_alpha->GetN(),gPPRap_alpha->GetY());
    maxPP = TMath::MaxElement(gPPRap_alpha->GetN(),gPPRap_alpha->GetY());
    
    maxTot = max(maxPbPb, maxPP);
    maxTot = maxPbPb*1.2;
    minTot = min(minPbPb, minPP);
    minTot = minPbPb - minPbPb*0.2;
    
    TH1* hYdummya = new TH1I("hYdummya","#alpha vs. |y|;|y|;#alpha",10, 0.,2.4);
    hYdummya->GetYaxis()->SetRangeUser(minTot,maxTot);
    hYdummya->GetYaxis()->SetLabelSize(0.05);
    hYdummya->GetYaxis()->SetTitleSize(0.06);
    hYdummya->GetXaxis()->SetLabelSize(0.05);
    hYdummya->GetXaxis()->SetTitleSize(0.06);
    hYdummya->GetXaxis()->SetTitleOffset(0.8);
    hYdummya->SetStats(0);
    hYdummya->Draw();
    
    TLegend* l1 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
    
//    gPbPbInt_alpha->SetMarkerStyle(21);
//    gPbPbInt_alpha->SetMarkerColor(2);
//    gPbPbInt_alpha->SetLineColor(2);
//    l1->AddEntry(gPbPbInt_alpha,"int. PbPb","lp");
//    gPbPbInt_alpha->Draw("psame");
//    
//    gPPInt_alpha->SetMarkerStyle(21);
//    gPPInt_alpha->SetMarkerColor(4);
//    gPPInt_alpha->SetLineColor(4);
//    l1->AddEntry(gPPInt_alpha,"int. pp","lp");
//    gPPInt_alpha->Draw("psame");
    
    gPbPbRap_alpha->SetMarkerStyle(20);
    gPbPbRap_alpha->SetMarkerColor(2);
    gPbPbRap_alpha->SetLineColor(2);
    l1->AddEntry(gPbPbRap_alpha,"diff. PbPb","lp");
    gPbPbRap_alpha->Draw("psame");
    
    gPPRap_alpha->SetMarkerStyle(20);
    gPPRap_alpha->SetMarkerColor(4);
    gPPRap_alpha->SetLineColor(4);
    l1->AddEntry(gPPRap_alpha,"diff. pp","lp");
    gPPRap_alpha->Draw("psame");
    
    avPbPb = gPbPbRap_alpha->GetMean(2);
    if (wMean) avPbPb = getWaverage(gPbPbRap_alpha);
    TLine* avPbPb1 = new TLine(hYdummya->GetXaxis()->GetXmin(),avPbPb,hYdummya->GetXaxis()->GetXmax(),avPbPb);
    avPbPb1->SetLineColor(2);
    avPbPb1->SetLineStyle(2);
    avPbPb1->SetLineWidth(3);
    l1->AddEntry(avPbPb1,"avg. PbPb","l");
    avPbPb1->Draw("same");
    
    avPP = gPPRap_alpha->GetMean(2);
    if (wMean) avPP = getWaverage(gPPRap_alpha);
    TLine* avPP1 = new TLine(hYdummya->GetXaxis()->GetXmin(),avPP,hYdummya->GetXaxis()->GetXmax(),avPP);
    avPP1->SetLineColor(4);
    avPP1->SetLineStyle(2);
    avPP1->SetLineWidth(3);
    l1->AddEntry(avPP1,"avg. pp","l");
    avPP1->Draw("same");
    
    l1->Draw("same");
    
    TLatex *  panel = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
    panel->SetNDC();
    panel->SetTextFont(42);
    panel->SetTextSize(0.05);
    panel->SetLineWidth(2);
    panel->Draw("same");
    
    TLatex *  text = new TLatex(0.4913793,0.2494759,Form("%s ; 0-100%%",ptRange.Data()));
    text->SetNDC();
    text->SetTextFont(42);
    text->SetTextSize(0.06708595);
    text->SetLineWidth(2);
    text->Draw("same");
    
    TLatex *  text2 = new TLatex(0.4913793,0.17,fName);
    text2->SetNDC();
    text2->SetTextFont(42);
    text2->SetTextSize(0.06708595);
    text2->SetLineWidth(2);
    text2->Draw("same");
    
    c1->SaveAs(Form("%s/alphaVsRap_%s.pdf",sPath.Data(),fName));
    aSave->Add(c1);
    
    cout << "Alpha values for PP rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    int nbins = gPPRap_alpha->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPPRap_alpha->GetY()[j] << endl;
    }
    
    cout << "Alpha values for PbPb rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    nbins = gPbPbRap_alpha->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPbPbRap_alpha->GetY()[j] << endl;
    }
    ///////////////////////////////
    ///////////////////////////////
    
    TCanvas* c2 = new TCanvas(Form("fVsRap_%s",fName),"canvas",81,98,700,504);
    c2->Range(-5.547577,-0.004024593,31.93896,0.03073326);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    c2->SetLeftMargin(0.1479885);
    c2->SetRightMargin(0.05172414);
    c2->SetTopMargin(0.08421053);
    c2->SetBottomMargin(0.1157895);
    c2->SetFrameBorderMode(0);
    c2->SetFrameBorderMode(0);
    
    maxPbPb = TMath::MaxElement(gPbPbRap_f->GetN(),gPbPbRap_f->GetY());
    minPbPb = TMath::MinElement(gPbPbRap_f->GetN(),gPbPbRap_f->GetY());
    
    minPP = TMath::MinElement(gPPRap_f->GetN(),gPPRap_f->GetY());
    maxPP = TMath::MaxElement(gPPRap_f->GetN(),gPPRap_f->GetY());
    
    maxTot = max(maxPbPb, maxPP);
    maxTot = maxPbPb*1.2;
    minTot = min(minPbPb, minPP);
    minTot = minPbPb - minPbPb*0.2;
    
    TH1* hYdummyf = new TH1I("hYdummyf","f vs. |y|;|y|;f",10, 0.,2.4);
    hYdummyf->GetYaxis()->SetRangeUser(minTot,maxTot);
    hYdummyf->GetYaxis()->SetLabelSize(0.05);
    hYdummyf->GetYaxis()->SetTitleSize(0.06);
    hYdummyf->GetXaxis()->SetLabelSize(0.05);
    hYdummyf->GetXaxis()->SetTitleSize(0.06);
    hYdummyf->GetXaxis()->SetTitleOffset(0.8);
    hYdummyf->SetStats(0);
    hYdummyf->Draw();
    
    TLegend* l2 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
    
//    gPbPbInt_f->SetMarkerStyle(21);
//    gPbPbInt_f->SetMarkerColor(2);
//    gPbPbInt_f->SetLineColor(2);
//    l2->AddEntry(gPbPbInt_f,"int. PbPb","lp");
//    gPbPbInt_f->Draw("psame");
//    
//    gPPInt_f->SetMarkerStyle(21);
//    gPPInt_f->SetMarkerColor(4);
//    gPPInt_f->SetLineColor(4);
//    l2->AddEntry(gPPInt_f,"int. pp","lp");
//    gPPInt_f->Draw("psame");
    
    gPbPbRap_f->SetMarkerStyle(20);
    gPbPbRap_f->SetMarkerColor(2);
    gPbPbRap_f->SetLineColor(2);
    l2->AddEntry(gPbPbRap_f,"diff. PbPb","lp");
    gPbPbRap_f->Draw("psame");
    
    gPPRap_f->SetMarkerStyle(20);
    gPPRap_f->SetMarkerColor(4);
    gPPRap_f->SetLineColor(4);
    l2->AddEntry(gPPRap_f,"diff. pp","lp");
    gPPRap_f->Draw("psame");
    
    avPbPb = gPbPbRap_f->GetMean(2);
    if (wMean) avPbPb = getWaverage(gPbPbRap_f);
    TLine* avPbPb2 = new TLine(hYdummyf->GetXaxis()->GetXmin(),avPbPb,hYdummyf->GetXaxis()->GetXmax(),avPbPb);
    avPbPb2->SetLineColor(2);
    avPbPb2->SetLineStyle(2);
    avPbPb2->SetLineWidth(3);
    l2->AddEntry(avPbPb2,"avg. PbPb","l");
    avPbPb2->Draw("same");
    
    avPP = gPPRap_f->GetMean(2);
    if (wMean) avPP = getWaverage(gPPRap_f);
    TLine* avPP2 = new TLine(hYdummyf->GetXaxis()->GetXmin(),avPP,hYdummyf->GetXaxis()->GetXmax(),avPP);
    avPP2->SetLineColor(4);
    avPP2->SetLineStyle(2);
    avPP2->SetLineWidth(3);
    l2->AddEntry(avPP2,"avg. pp","l");
    avPP2->Draw("same");
    
    l2->Draw("same");
    
    TLatex *  panel2 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
    panel2->SetNDC();
    panel2->SetTextFont(42);
    panel2->SetTextSize(0.05);
    panel2->SetLineWidth(2);
    panel2->Draw("same");
    
    TLatex *  textf = new TLatex(0.4913793,0.2494759,Form("%s ; 0-100%%",ptRange.Data()));
    textf->SetNDC();
    textf->SetTextFont(42);
    textf->SetTextSize(0.06708595);
    textf->SetLineWidth(2);
    textf->Draw("same");
    
    TLatex *  text2f = new TLatex(0.4913793,0.17,fName);
    text2f->SetNDC();
    text2f->SetTextFont(42);
    text2f->SetTextSize(0.06708595);
    text2f->SetLineWidth(2);
    text2f->Draw("same");
    
    c2->SaveAs(Form("%s/fVsRap_%s.pdf",sPath.Data(),fName));
    aSave->Add(c2);
    
    ///////////////////////////////
    ///////////////////////////////
    
    TCanvas* c3 = new TCanvas(Form("nVsRap_%s",fName),"canvas",81,98,700,504);
    c3->Range(-5.547577,-0.004024593,31.93896,0.03073326);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetBorderSize(2);
    c3->SetLeftMargin(0.1479885);
    c3->SetRightMargin(0.05172414);
    c3->SetTopMargin(0.08421053);
    c3->SetBottomMargin(0.1157895);
    c3->SetFrameBorderMode(0);
    c3->SetFrameBorderMode(0);
    
    maxPbPb = TMath::MaxElement(gPbPbRap_n->GetN(),gPbPbRap_n->GetY());
    minPbPb = TMath::MinElement(gPbPbRap_n->GetN(),gPbPbRap_n->GetY());
    
    minPP = TMath::MinElement(gPPRap_n->GetN(),gPPRap_n->GetY());
    maxPP = TMath::MaxElement(gPPRap_n->GetN(),gPPRap_n->GetY());
    
    maxTot = max(maxPbPb, maxPP);
    maxTot = maxPbPb*1.2;
    minTot = min(minPbPb, minPP);
    minTot = minPbPb - minPbPb*0.2;
    
    TH1* hYdummyn = new TH1I("hYdummyn","n vs. |y|;|y|;n",10, 0.,2.4);
    hYdummyn->GetYaxis()->SetRangeUser(minTot,maxTot);
    hYdummyn->GetYaxis()->SetLabelSize(0.05);
    hYdummyn->GetYaxis()->SetTitleSize(0.06);
    hYdummyn->GetXaxis()->SetLabelSize(0.05);
    hYdummyn->GetXaxis()->SetTitleSize(0.06);
    hYdummyn->GetXaxis()->SetTitleOffset(0.8);
    hYdummyn->SetStats(0);
    hYdummyn->Draw();
    
    TLegend* l3 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
    
//    gPbPbInt_n->SetMarkerStyle(21);
//    gPbPbInt_n->SetMarkerColor(2);
//    gPbPbInt_n->SetLineColor(2);
//    l3->AddEntry(gPbPbInt_n,"int. PbPb","lp");
//    gPbPbInt_n->Draw("psame");
//    
//    gPPInt_n->SetMarkerStyle(21);
//    gPPInt_n->SetMarkerColor(4);
//    gPPInt_n->SetLineColor(4);
//    l3->AddEntry(gPPInt_n,"int. pp","lp");
//    gPPInt_n->Draw("psame");
    
    gPbPbRap_n->SetMarkerStyle(20);
    gPbPbRap_n->SetMarkerColor(2);
    gPbPbRap_n->SetLineColor(2);
    l3->AddEntry(gPbPbRap_n,"diff. PbPb","lp");
    gPbPbRap_n->Draw("psame");
    
    gPPRap_n->SetMarkerStyle(20);
    gPPRap_n->SetMarkerColor(4);
    gPPRap_n->SetLineColor(4);
    l3->AddEntry(gPPRap_n,"diff. pp","lp");
    gPPRap_n->Draw("psame");
    
    avPbPb = gPbPbRap_n->GetMean(2);
    if (wMean) avPbPb = getWaverage(gPbPbRap_n);
    TLine* avPbPb3 = new TLine(hYdummyn->GetXaxis()->GetXmin(),avPbPb,hYdummyn->GetXaxis()->GetXmax(),avPbPb);
    avPbPb3->SetLineColor(2);
    avPbPb3->SetLineStyle(2);
    avPbPb3->SetLineWidth(3);
    l3->AddEntry(avPbPb3,"avg. PbPb","l");
    avPbPb3->Draw("same");
    
    avPP = gPPRap_n->GetMean(2);
    if (wMean) avPP = getWaverage(gPPRap_n);
    TLine* avPP3 = new TLine(hYdummyn->GetXaxis()->GetXmin(),avPP,hYdummyn->GetXaxis()->GetXmax(),avPP);
    avPP3->SetLineColor(4);
    avPP3->SetLineStyle(2);
    avPP3->SetLineWidth(3);
    l3->AddEntry(avPP3,"avg. pp","l");
    avPP3->Draw("same");
    
    l3->Draw("same");
    
    TLatex *  panel3 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
    panel3->SetNDC();
    panel3->SetTextFont(42);
    panel3->SetTextSize(0.05);
    panel3->SetLineWidth(2);
    panel3->Draw("same");
    
    TLatex *  textn = new TLatex(0.4913793,0.2494759,Form("%s ; 0-100%%",ptRange.Data()));
    textn->SetNDC();
    textn->SetTextFont(42);
    textn->SetTextSize(0.06708595);
    textn->SetLineWidth(2);
    textn->Draw("same");
    
    TLatex *  text2n = new TLatex(0.4913793,0.17,fName);
    text2n->SetNDC();
    text2n->SetTextFont(42);
    text2n->SetTextSize(0.06708595);
    text2n->SetLineWidth(2);
    text2n->Draw("same");
    
    c3->SaveAs(Form("%s/nVsRap_%s.pdf",sPath.Data(),fName));
    aSave->Add(c3);
    
    cout << "n values for PP rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    nbins = gPPRap_n->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPPRap_n->GetY()[j] << endl;
    }
    
    cout << "n values for PbPb rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    nbins = gPbPbRap_n->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPbPbRap_n->GetY()[j] << endl;
    }
    
    ///////////////////////////////
    ///////////////////////////////
    
    TCanvas* c4 = new TCanvas(Form("sVsRap_%s",fName),"canvas",81,98,700,504);
    c4->Range(-5.547577,-0.004024593,31.93896,0.03073326);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetBorderSize(2);
    c4->SetLeftMargin(0.1479885);
    c4->SetRightMargin(0.05172414);
    c4->SetTopMargin(0.08421053);
    c4->SetBottomMargin(0.1157895);
    c4->SetFrameBorderMode(0);
    c4->SetFrameBorderMode(0);
    
    maxPbPb = TMath::MaxElement(gPbPbRap_sigma->GetN(),gPbPbRap_sigma->GetY());
    minPbPb = TMath::MinElement(gPbPbRap_sigma->GetN(),gPbPbRap_sigma->GetY());
    
    minPP = TMath::MinElement(gPPRap_sigma->GetN(),gPPRap_sigma->GetY());
    maxPP = TMath::MaxElement(gPPRap_sigma->GetN(),gPPRap_sigma->GetY());
    
    maxTot = max(maxPbPb, maxPP);
    maxTot = maxPbPb*1.2;
    minTot = min(minPbPb, minPP);
    minTot = minPbPb - minPbPb*0.2;
    
    TH1* hYdummys = new TH1I("hYdummys","#sigma vs. |y|;|y|;#sigma_{1} (MeV/c^{2})",10, 0.,2.4);
    hYdummys->GetYaxis()->SetRangeUser(minTot,maxTot);
    hYdummys->GetYaxis()->SetLabelSize(0.05);
    hYdummys->GetYaxis()->SetTitleSize(0.06);
    hYdummys->GetXaxis()->SetLabelSize(0.05);
    hYdummys->GetXaxis()->SetTitleSize(0.06);
    hYdummys->GetXaxis()->SetTitleOffset(0.8);
    hYdummys->SetStats(0);
    hYdummys->Draw();
    
    TLegend* l4 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
    
//    gPbPbInt_sigma->SetMarkerStyle(21);
//    gPbPbInt_sigma->SetMarkerColor(2);
//    gPbPbInt_sigma->SetLineColor(2);
//    l4->AddEntry(gPbPbInt_sigma,"int. PbPb","lp");
//    gPbPbInt_sigma->Draw("psame");
//    
//    gPPInt_sigma->SetMarkerStyle(21);
//    gPPInt_sigma->SetMarkerColor(4);
//    gPPInt_sigma->SetLineColor(4);
//    l4->AddEntry(gPPInt_sigma,"int. pp","lp");
//    gPPInt_sigma->Draw("psame");
    
    gPbPbRap_sigma->SetMarkerStyle(20);
    gPbPbRap_sigma->SetMarkerColor(2);
    gPbPbRap_sigma->SetLineColor(2);
    l4->AddEntry(gPbPbRap_sigma,"diff. PbPb","lp");
    gPbPbRap_sigma->Draw("psame");
    
    gPPRap_sigma->SetMarkerStyle(20);
    gPPRap_sigma->SetMarkerColor(4);
    gPPRap_sigma->SetLineColor(4);
    l4->AddEntry(gPPRap_sigma,"diff. pp","lp");
    gPPRap_sigma->Draw("psame");
    
    avPbPb = gPbPbRap_sigma->GetMean(2);
    if (wMean) avPbPb = getWaverage(gPbPbRap_sigma);
    TLine* avPbPb4 = new TLine(hYdummys->GetXaxis()->GetXmin(),avPbPb,hYdummys->GetXaxis()->GetXmax(),avPbPb);
    avPbPb4->SetLineColor(2);
    avPbPb4->SetLineStyle(2);
    avPbPb4->SetLineWidth(3);
    l4->AddEntry(avPbPb4,"avg. PbPb","l");
    avPbPb4->Draw("same");
    
    avPP = gPPRap_sigma->GetMean(2);
    if (wMean) avPP = getWaverage(gPPRap_sigma);
    TLine* avPP4 = new TLine(hYdummys->GetXaxis()->GetXmin(),avPP,hYdummys->GetXaxis()->GetXmax(),avPP);
    avPP4->SetLineColor(4);
    avPP4->SetLineStyle(2);
    avPP4->SetLineWidth(3);
    l4->AddEntry(avPP4,"avg. pp","l");
    avPP4->Draw("same");
    
    l4->Draw("same");
    
    TLatex *  panel4 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %2.1f GeV/c^{2} (%2.1f-%2.1f)}{av. PbPb =%2.1f GeV/c^{2} (%2.1f-%2.1f)}",avPP*1000,minPP*1000,maxPP*1000,avPbPb*1000,minPbPb*1000,maxPbPb*1000));
    panel4->SetNDC();
    panel4->SetTextFont(42);
    panel4->SetTextSize(0.05);
    panel4->SetLineWidth(2);
    panel4->Draw("same");
    
    TLatex *  texts = new TLatex(0.4913793,0.2494759,Form("%s ; 0-100%%",ptRange.Data()));
    texts->SetNDC();
    texts->SetTextFont(42);
    texts->SetTextSize(0.06708595);
    texts->SetLineWidth(2);
    texts->Draw("same");
    
    TLatex *  text2s = new TLatex(0.4913793,0.17,fName);
    text2s->SetNDC();
    text2s->SetTextFont(42);
    text2s->SetTextSize(0.06708595);
    text2s->SetLineWidth(2);
    text2s->Draw("same");
    
    c4->SaveAs(Form("%s/sVsRap_%s.pdf",sPath.Data(),fName));
    aSave->Add(c4);
    
    ///////////////////////////////
    ///////////////////////////////
    
    TCanvas* c5 = new TCanvas(Form("rsVsRap_%s",fName),"canvas",81,98,700,504);
    c5->Range(-5.547577,-0.004024593,31.93896,0.03073326);
    c5->SetFillColor(0);
    c5->SetBorderMode(0);
    c5->SetBorderSize(2);
    c5->SetLeftMargin(0.1479885);
    c5->SetRightMargin(0.05172414);
    c5->SetTopMargin(0.08421053);
    c5->SetBottomMargin(0.1157895);
    c5->SetFrameBorderMode(0);
    c5->SetFrameBorderMode(0);
    
    maxPbPb = TMath::MaxElement(gPbPbRap_rsigma->GetN(),gPbPbRap_rsigma->GetY());
    minPbPb = TMath::MinElement(gPbPbRap_rsigma->GetN(),gPbPbRap_rsigma->GetY());
    
    minPP = TMath::MinElement(gPPRap_rsigma->GetN(),gPPRap_rsigma->GetY());
    maxPP = TMath::MaxElement(gPPRap_rsigma->GetN(),gPPRap_rsigma->GetY());
    
    maxTot = max(maxPbPb, maxPP);
    maxTot = maxPbPb*1.2;
    minTot = min(minPbPb, minPP);
    minTot = minPbPb - minPbPb*0.2;
    
    TH1* hYdummyrs = new TH1I("hpTdummyrs","#sigma_{2}/#sigma_{1} vs. |y|;|y|;#sigma_{2}/#sigma_{1}",10, 0.,2.4);
    hYdummyrs->GetYaxis()->SetRangeUser(minTot,maxTot);
    hYdummyrs->GetYaxis()->SetLabelSize(0.05);
    hYdummyrs->GetYaxis()->SetTitleSize(0.06);
    hYdummyrs->GetXaxis()->SetLabelSize(0.05);
    hYdummyrs->GetXaxis()->SetTitleSize(0.06);
    hYdummyrs->GetXaxis()->SetTitleOffset(0.8);
    hYdummyrs->SetStats(0);
    hYdummyrs->Draw();
    
    TLegend* l5 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
    
//    gPbPbInt_rsigma->SetMarkerStyle(21);
//    gPbPbInt_rsigma->SetMarkerColor(2);
//    gPbPbInt_rsigma->SetLineColor(2);
//    l5->AddEntry(gPbPbInt_rsigma,"int. PbPb","lp");
//    gPbPbInt_rsigma->Draw("psame");
//    
//    gPPInt_rsigma->SetMarkerStyle(21);
//    gPPInt_rsigma->SetMarkerColor(4);
//    gPPInt_rsigma->SetLineColor(4);
//    l5->AddEntry(gPPInt_rsigma,"int. pp","lp");
//    gPPInt_rsigma->Draw("psame");
    
    gPbPbRap_rsigma->SetMarkerStyle(20);
    gPbPbRap_rsigma->SetMarkerColor(2);
    gPbPbRap_rsigma->SetLineColor(2);
    l5->AddEntry(gPbPbRap_rsigma,"diff. PbPb","lp");
    gPbPbRap_rsigma->Draw("psame");
    
    gPPRap_rsigma->SetMarkerStyle(20);
    gPPRap_rsigma->SetMarkerColor(4);
    gPPRap_rsigma->SetLineColor(4);
    l5->AddEntry(gPPRap_rsigma,"diff. pp","lp");
    gPPRap_rsigma->Draw("psame");
    
    avPbPb = gPbPbRap_rsigma->GetMean(2);
    if (wMean) avPbPb = getWaverage(gPbPbRap_rsigma);
    TLine* avPbPb5 = new TLine(hYdummyrs->GetXaxis()->GetXmin(),avPbPb,hYdummyrs->GetXaxis()->GetXmax(),avPbPb);
    avPbPb5->SetLineColor(2);
    avPbPb5->SetLineStyle(2);
    avPbPb5->SetLineWidth(3);
    l5->AddEntry(avPbPb5,"avg. PbPb","l");
    avPbPb5->Draw("same");
    
    avPP = gPPRap_rsigma->GetMean(2);
    if (wMean) avPP = getWaverage(gPPRap_rsigma);
    TLine* avPP5 = new TLine(hYdummyrs->GetXaxis()->GetXmin(),avPP,hYdummyrs->GetXaxis()->GetXmax(),avPP);
    avPP5->SetLineColor(4);
    avPP5->SetLineStyle(2);
    avPP5->SetLineWidth(3);
    l5->AddEntry(avPP5,"avg. pp","l");
    avPP5->Draw("same");
    
    l5->Draw("same");
    
    TLatex *  panel5 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
    panel5->SetNDC();
    panel5->SetTextFont(42);
    panel5->SetTextSize(0.05);
    panel5->SetLineWidth(2);
    panel5->Draw("same");
    
    TLatex *  textrs = new TLatex(0.4913793,0.2494759,Form("%s ; 0-100%%",ptRange.Data()));
    textrs->SetNDC();
    textrs->SetTextFont(42);
    textrs->SetTextSize(0.06708595);
    textrs->SetLineWidth(2);
    textrs->Draw("same");
    
    TLatex *  text2rs = new TLatex(0.4913793,0.17,fName);
    text2rs->SetNDC();
    text2rs->SetTextFont(42);
    text2rs->SetTextSize(0.06708595);
    text2rs->SetLineWidth(2);
    text2rs->Draw("same");
    
    c5->SaveAs(Form("%s/rsVsRap_%s.pdf",sPath.Data(),fName));
    aSave->Add(c5);
    
    cout << "rsigma values for PP rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    nbins = gPPRap_rsigma->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPPRap_rsigma->GetY()[j] << endl;
    }
    
    cout << "rsigma values for PbPb rapidity bins in " << Form("%s ; 0-100%%",ptRange.Data()) << " :" << endl;
    nbins = gPbPbRap_rsigma->GetN();
    for (int j = 0 ; j < nbins ; j++)
    {
      cout << gPbPbRap_rsigma->GetY()[j] << endl;
    }
    
  }
  else
  {
    ///////////////////////////////
    /////////////pT dep////////////
    ///////////////////////////////
    if (isPtBin)
    {
      TCanvas* c1 = new TCanvas(Form("alphaVsPt_%s",fName),"canvas",81,98,700,504);
      c1->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c1->SetFillColor(0);
      c1->SetBorderMode(0);
      c1->SetBorderSize(2);
      c1->SetLeftMargin(0.1479885);
      c1->SetRightMargin(0.05172414);
      c1->SetTopMargin(0.08421053);
      c1->SetBottomMargin(0.1157895);
      c1->SetFrameBorderMode(0);
      c1->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbPt_alpha->GetN(),gPbPbPt_alpha->GetY()),TMath::MaxElement(gPbPbInt_alpha->GetN(),gPbPbInt_alpha->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbPt_alpha->GetN(),gPbPbPt_alpha->GetY()),TMath::MinElement(gPbPbInt_alpha->GetN(),gPbPbInt_alpha->GetY()));
      
      maxPP = max(TMath::MaxElement(gPPPt_alpha->GetN(),gPPPt_alpha->GetY()),TMath::MaxElement(gPPInt_alpha->GetN(),gPPInt_alpha->GetY()));
      minPP = min(TMath::MinElement(gPPPt_alpha->GetN(),gPPPt_alpha->GetY()),TMath::MinElement(gPPInt_alpha->GetN(),gPPInt_alpha->GetY()));
      
      maxTot = max(maxPbPb,maxPP);
      maxTot = maxTot*1.2;
      minTot = min(minPbPb,minPP);
      minTot = minTot - minTot*0.2;
      
      TH1* hpTdummya = new TH1I("hpTdummya","#alpha vs. p_{T};p_{T}(GeV/c);#alpha",50, 0.,50.);
      hpTdummya->GetYaxis()->SetRangeUser(minTot,maxTot);
      hpTdummya->GetYaxis()->SetLabelSize(0.05);
      hpTdummya->GetYaxis()->SetTitleSize(0.06);
      hpTdummya->GetXaxis()->SetLabelSize(0.05);
      hpTdummya->GetXaxis()->SetTitleSize(0.06);
      hpTdummya->GetXaxis()->SetTitleOffset(0.8);
      hpTdummya->SetStats(0);
      hpTdummya->Draw();
      
      TLegend* l1 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_alpha->SetMarkerStyle(21);
      gPbPbInt_alpha->SetMarkerColor(2);
      gPbPbInt_alpha->SetLineColor(2);
      l1->AddEntry(gPbPbInt_alpha,"int. PbPb","lp");
      gPbPbInt_alpha->Draw("psame");
      
      gPPInt_alpha->SetMarkerStyle(21);
      gPPInt_alpha->SetMarkerColor(4);
      gPPInt_alpha->SetLineColor(4);
      l1->AddEntry(gPPInt_alpha,"int. pp","lp");
      gPPInt_alpha->Draw("psame");
      
      gPbPbPt_alpha->SetMarkerStyle(20);
      gPbPbPt_alpha->SetMarkerColor(2);
      gPbPbPt_alpha->SetLineColor(2);
      l1->AddEntry(gPbPbPt_alpha,"diff. PbPb","lp");
      gPbPbPt_alpha->Draw("psame");
      
      gPPPt_alpha->SetMarkerStyle(20);
      gPPPt_alpha->SetMarkerColor(4);
      gPPPt_alpha->SetLineColor(4);
      l1->AddEntry(gPPPt_alpha,"diff. pp","lp");
      gPPPt_alpha->Draw("psame");
      
      avPbPb = gPbPbPt_alpha->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbPt_alpha);
      TLine* avPbPb1 = new TLine(hpTdummya->GetXaxis()->GetXmin(),avPbPb,hpTdummya->GetXaxis()->GetXmax(),avPbPb);
      avPbPb1->SetLineColor(2);
      avPbPb1->SetLineStyle(2);
      avPbPb1->SetLineWidth(3);
      l1->AddEntry(avPbPb1,"avg. PbPb","l");
      avPbPb1->Draw("same");
      
      avPP = gPPPt_alpha->GetMean(2);
      if (wMean) avPP = getWaverage(gPPPt_alpha);
      TLine* avPP1 = new TLine(hpTdummya->GetXaxis()->GetXmin(),avPP,hpTdummya->GetXaxis()->GetXmax(),avPP);
      avPP1->SetLineColor(4);
      avPP1->SetLineStyle(2);
      avPP1->SetLineWidth(3);
      l1->AddEntry(avPP1,"avg. pp","l");
      avPP1->Draw("same");
      
      l1->Draw("same");
      
      TLatex *  panel = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
      panel->SetNDC();
      panel->SetTextFont(42);
      panel->SetTextSize(0.05);
      panel->SetLineWidth(2);
      panel->Draw("same");
      
      TLatex *  text = new TLatex(0.4913793,0.2494759,isCBin ?  Form("%s ; 0-100%%",rapRange.Data()) : Form("%s ; %s",rapRange.Data(),centRange.Data()));
      text->SetNDC();
      text->SetTextFont(42);
      text->SetTextSize(0.06708595);
      text->SetLineWidth(2);
      text->Draw("same");
      
      TLatex *  text2 = new TLatex(0.4913793,0.17,fName);
      text2->SetNDC();
      text2->SetTextFont(42);
      text2->SetTextSize(0.06708595);
      text2->SetLineWidth(2);
      text2->Draw("same");
      
      c1->SaveAs(Form("%s/alphaVsPt_%s.pdf",sPath.Data(),fName));
      aSave->Add(c1);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c2 = new TCanvas(Form("fVsPt_%s",fName),"canvas",81,98,700,504);
      c2->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c2->SetFillColor(0);
      c2->SetBorderMode(0);
      c2->SetBorderSize(2);
      c2->SetLeftMargin(0.1479885);
      c2->SetRightMargin(0.05172414);
      c2->SetTopMargin(0.08421053);
      c2->SetBottomMargin(0.1157895);
      c2->SetFrameBorderMode(0);
      c2->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbPt_f->GetN(),gPbPbPt_f->GetY()),TMath::MaxElement(gPbPbInt_f->GetN(),gPbPbInt_f->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbPt_f->GetN(),gPbPbPt_f->GetY()),TMath::MinElement(gPbPbInt_f->GetN(),gPbPbInt_f->GetY()));
      
      maxPP = max(TMath::MaxElement(gPPPt_f->GetN(),gPPPt_f->GetY()),TMath::MaxElement(gPPInt_f->GetN(),gPPInt_f->GetY()));
      minPP = min(TMath::MinElement(gPPPt_f->GetN(),gPPPt_f->GetY()),TMath::MinElement(gPPInt_f->GetN(),gPPInt_f->GetY()));
      
      maxTot = max(maxPbPb,maxPP);
      maxTot = maxTot*1.2;
      minTot = min(minPbPb,minPP);
      minTot = minTot - minTot*0.2;
      
      TH1* hpTdummyf = new TH1I("hpTdummyf","f vs. p_{T};p_{T}(GeV/c);f",50, 0.,50.);
      hpTdummyf->GetYaxis()->SetRangeUser(minTot,maxTot);
      hpTdummyf->GetYaxis()->SetLabelSize(0.05);
      hpTdummyf->GetYaxis()->SetTitleSize(0.06);
      hpTdummyf->GetXaxis()->SetLabelSize(0.05);
      hpTdummyf->GetXaxis()->SetTitleSize(0.06);
      hpTdummyf->GetXaxis()->SetTitleOffset(0.8);
      hpTdummyf->SetStats(0);
      hpTdummyf->Draw();
      
      TLegend* l2 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_f->SetMarkerStyle(21);
      gPbPbInt_f->SetMarkerColor(2);
      gPbPbInt_f->SetLineColor(2);
      l2->AddEntry(gPbPbInt_f,"int. PbPb","lp");
      gPbPbInt_f->Draw("psame");
      
      gPPInt_f->SetMarkerStyle(21);
      gPPInt_f->SetMarkerColor(4);
      gPPInt_f->SetLineColor(4);
      l2->AddEntry(gPPInt_f,"int. pp","lp");
      gPPInt_f->Draw("psame");
      
      gPbPbPt_f->SetMarkerStyle(20);
      gPbPbPt_f->SetMarkerColor(2);
      gPbPbPt_f->SetLineColor(2);
      l2->AddEntry(gPbPbPt_f,"diff. PbPb","lp");
      gPbPbPt_f->Draw("psame");
      
      gPPPt_f->SetMarkerStyle(20);
      gPPPt_f->SetMarkerColor(4);
      gPPPt_f->SetLineColor(4);
      l2->AddEntry(gPPPt_f,"diff. pp","lp");
      gPPPt_f->Draw("psame");
      
      avPbPb = gPbPbPt_f->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbPt_f);
      TLine* avPbPb2 = new TLine(hpTdummyf->GetXaxis()->GetXmin(),avPbPb,hpTdummyf->GetXaxis()->GetXmax(),avPbPb);
      avPbPb2->SetLineColor(2);
      avPbPb2->SetLineStyle(2);
      avPbPb2->SetLineWidth(3);
      l2->AddEntry(avPbPb2,"avg. PbPb","l");
      avPbPb2->Draw("same");
      
      avPP = gPPPt_f->GetMean(2);
      if (wMean) avPP = getWaverage(gPPPt_f);
      TLine* avPP2 = new TLine(hpTdummyf->GetXaxis()->GetXmin(),avPP,hpTdummyf->GetXaxis()->GetXmax(),avPP);
      avPP2->SetLineColor(4);
      avPP2->SetLineStyle(2);
      avPP2->SetLineWidth(3);
      l2->AddEntry(avPP2,"avg. pp","l");
      avPP2->Draw("same");
      
      l2->Draw("same");
      
      TLatex *  panel2 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
      panel2->SetNDC();
      panel2->SetTextFont(42);
      panel2->SetTextSize(0.05);
      panel2->SetLineWidth(2);
      panel2->Draw("same");
      
      TLatex *  textf = new TLatex(0.4913793,0.2494759,isCBin ? Form("%s ; 0-100%%",rapRange.Data()) : Form("%s ; %s",rapRange.Data(),centRange.Data()));
      textf->SetNDC();
      textf->SetTextFont(42);
      textf->SetTextSize(0.06708595);
      textf->SetLineWidth(2);
      textf->Draw("same");
      
      TLatex *  text2f = new TLatex(0.4913793,0.17,fName);
      text2f->SetNDC();
      text2f->SetTextFont(42);
      text2f->SetTextSize(0.06708595);
      text2f->SetLineWidth(2);
      text2f->Draw("same");
      
      c2->SaveAs(Form("%s/fVsPt_%s.pdf",sPath.Data(),fName));
      aSave->Add(c2);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c3 = new TCanvas(Form("nVsPt_%s",fName),"canvas",81,98,700,504);
      c3->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c3->SetFillColor(0);
      c3->SetBorderMode(0);
      c3->SetBorderSize(2);
      c3->SetLeftMargin(0.1479885);
      c3->SetRightMargin(0.05172414);
      c3->SetTopMargin(0.08421053);
      c3->SetBottomMargin(0.1157895);
      c3->SetFrameBorderMode(0);
      c3->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbPt_n->GetN(),gPbPbPt_n->GetY()),TMath::MaxElement(gPbPbInt_n->GetN(),gPbPbInt_n->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbPt_n->GetN(),gPbPbPt_n->GetY()),TMath::MinElement(gPbPbInt_n->GetN(),gPbPbInt_n->GetY()));
      
      maxPP = max(TMath::MaxElement(gPPPt_n->GetN(),gPPPt_n->GetY()),TMath::MaxElement(gPPInt_n->GetN(),gPPInt_n->GetY()));
      minPP = min(TMath::MinElement(gPPPt_n->GetN(),gPPPt_n->GetY()),TMath::MinElement(gPPInt_n->GetN(),gPPInt_n->GetY()));
      
      maxTot = max(maxPbPb,maxPP);
      maxTot = maxTot*1.2;
      minTot = min(minPbPb,minPP);
      minTot = minTot - minTot*0.2;
      
      TH1* hpTdummyn = new TH1I("hpTdummyn","n vs. p_{T};p_{T}(GeV/c);n",50, 0.,50.);
      hpTdummyn->GetYaxis()->SetRangeUser(minTot,maxTot);
      hpTdummyn->GetYaxis()->SetLabelSize(0.05);
      hpTdummyn->GetYaxis()->SetTitleSize(0.06);
      hpTdummyn->GetXaxis()->SetLabelSize(0.05);
      hpTdummyn->GetXaxis()->SetTitleSize(0.06);
      hpTdummyn->GetXaxis()->SetTitleOffset(0.8);
      hpTdummyn->SetStats(0);
      hpTdummyn->Draw();
      
      TLegend* l3 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_n->SetMarkerStyle(21);
      gPbPbInt_n->SetMarkerColor(2);
      gPbPbInt_n->SetLineColor(2);
      l3->AddEntry(gPbPbInt_n,"int. PbPb","lp");
      gPbPbInt_n->Draw("psame");
      
      gPPInt_n->SetMarkerStyle(21);
      gPPInt_n->SetMarkerColor(4);
      gPPInt_n->SetLineColor(4);
      l3->AddEntry(gPPInt_n,"int. pp","lp");
      gPPInt_n->Draw("psame");
      
      gPbPbPt_n->SetMarkerStyle(20);
      gPbPbPt_n->SetMarkerColor(2);
      gPbPbPt_n->SetLineColor(2);
      l3->AddEntry(gPbPbPt_n,"diff. PbPb","lp");
      gPbPbPt_n->Draw("psame");
      
      gPPPt_n->SetMarkerStyle(20);
      gPPPt_n->SetMarkerColor(4);
      gPPPt_n->SetLineColor(4);
      l3->AddEntry(gPPPt_n,"diff. pp","lp");
      gPPPt_n->Draw("psame");
      
      avPbPb = gPbPbPt_n->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbPt_n);
      TLine* avPbPb3 = new TLine(hpTdummyn->GetXaxis()->GetXmin(),avPbPb,hpTdummyn->GetXaxis()->GetXmax(),avPbPb);
      avPbPb3->SetLineColor(2);
      avPbPb3->SetLineStyle(2);
      avPbPb3->SetLineWidth(3);
      l3->AddEntry(avPbPb3,"avg. PbPb","l");
      avPbPb3->Draw("same");
      
      avPP = gPPPt_n->GetMean(2);
      if (wMean) avPP = getWaverage(gPPPt_n);
      TLine* avPP3 = new TLine(hpTdummyn->GetXaxis()->GetXmin(),avPP,hpTdummyn->GetXaxis()->GetXmax(),avPP);
      avPP3->SetLineColor(4);
      avPP3->SetLineStyle(2);
      avPP3->SetLineWidth(3);
      l3->AddEntry(avPP3,"avg. pp","l");
      avPP3->Draw("same");
      
      l3->Draw("same");
      
      TLatex *  panel3 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
      panel3->SetNDC();
      panel3->SetTextFont(42);
      panel3->SetTextSize(0.05);
      panel3->SetLineWidth(2);
      panel3->Draw("same");
      
      TLatex *  textn = new TLatex(0.4913793,0.2494759,isCBin ? Form("%s ; 0-100%%",rapRange.Data()) : Form("%s ; %s",rapRange.Data(),centRange.Data()));
      textn->SetNDC();
      textn->SetTextFont(42);
      textn->SetTextSize(0.06708595);
      textn->SetLineWidth(2);
      textn->Draw("same");
      
      TLatex *  text2n = new TLatex(0.4913793,0.17,fName);
      text2n->SetNDC();
      text2n->SetTextFont(42);
      text2n->SetTextSize(0.06708595);
      text2n->SetLineWidth(2);
      text2n->Draw("same");
      
      c3->SaveAs(Form("%s/nVsPt_%s.pdf",sPath.Data(),fName));
      aSave->Add(c3);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c4 = new TCanvas(Form("sVsPt_%s",fName),"canvas",81,98,700,504);
      c4->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c4->SetFillColor(0);
      c4->SetBorderMode(0);
      c4->SetBorderSize(2);
      c4->SetLeftMargin(0.1479885);
      c4->SetRightMargin(0.05172414);
      c4->SetTopMargin(0.08421053);
      c4->SetBottomMargin(0.1157895);
      c4->SetFrameBorderMode(0);
      c4->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbPt_sigma->GetN(),gPbPbPt_sigma->GetY()),TMath::MaxElement(gPbPbInt_sigma->GetN(),gPbPbInt_sigma->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbPt_sigma->GetN(),gPbPbPt_sigma->GetY()),TMath::MinElement(gPbPbInt_sigma->GetN(),gPbPbInt_sigma->GetY()));
      
      maxPP = max(TMath::MaxElement(gPPPt_sigma->GetN(),gPPPt_sigma->GetY()),TMath::MaxElement(gPPInt_sigma->GetN(),gPPInt_sigma->GetY()));
      minPP = min(TMath::MinElement(gPPPt_sigma->GetN(),gPPPt_sigma->GetY()),TMath::MinElement(gPPInt_sigma->GetN(),gPPInt_sigma->GetY()));
      
      maxTot = max(maxPbPb,maxPP);
      maxTot = maxTot*1.2;
      minTot = min(minPbPb,minPP);
      minTot = minTot - minTot*0.2;
      
      TH1* hpTdummys = new TH1I("hpTdummys","#sigma vs. p_{T};p_{T}(GeV/c);#sigma_{1} (MeV/c^{2})",50, 0.,50.);
      hpTdummys->GetYaxis()->SetRangeUser(minTot,maxTot);
      hpTdummys->GetYaxis()->SetLabelSize(0.05);
      hpTdummys->GetYaxis()->SetTitleSize(0.06);
      hpTdummys->GetXaxis()->SetLabelSize(0.05);
      hpTdummys->GetXaxis()->SetTitleSize(0.06);
      hpTdummys->GetXaxis()->SetTitleOffset(0.8);
      hpTdummys->SetStats(0);
      hpTdummys->Draw();
      
      TLegend* l4 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_sigma->SetMarkerStyle(21);
      gPbPbInt_sigma->SetMarkerColor(2);
      gPbPbInt_sigma->SetLineColor(2);
      l4->AddEntry(gPbPbInt_sigma,"int. PbPb","lp");
      gPbPbInt_sigma->Draw("psame");
      
      gPPInt_sigma->SetMarkerStyle(21);
      gPPInt_sigma->SetMarkerColor(4);
      gPPInt_sigma->SetLineColor(4);
      l4->AddEntry(gPPInt_sigma,"int. pp","lp");
      gPPInt_sigma->Draw("psame");
      
      gPbPbPt_sigma->SetMarkerStyle(20);
      gPbPbPt_sigma->SetMarkerColor(2);
      gPbPbPt_sigma->SetLineColor(2);
      l4->AddEntry(gPbPbPt_sigma,"diff. PbPb","lp");
      gPbPbPt_sigma->Draw("psame");
      
      gPPPt_sigma->SetMarkerStyle(20);
      gPPPt_sigma->SetMarkerColor(4);
      gPPPt_sigma->SetLineColor(4);
      l4->AddEntry(gPPPt_sigma,"diff. pp","lp");
      gPPPt_sigma->Draw("psame");
      
      avPbPb = gPbPbPt_sigma->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbPt_sigma);
      TLine* avPbPb4 = new TLine(hpTdummys->GetXaxis()->GetXmin(),avPbPb,hpTdummys->GetXaxis()->GetXmax(),avPbPb);
      avPbPb4->SetLineColor(2);
      avPbPb4->SetLineStyle(2);
      avPbPb4->SetLineWidth(3);
      l4->AddEntry(avPbPb4,"avg. PbPb","l");
      avPbPb4->Draw("same");
      
      avPP = gPPPt_sigma->GetMean(2);
      if (wMean) avPP = getWaverage(gPPPt_sigma);
      TLine* avPP4 = new TLine(hpTdummys->GetXaxis()->GetXmin(),avPP,hpTdummys->GetXaxis()->GetXmax(),avPP);
      avPP4->SetLineColor(4);
      avPP4->SetLineStyle(2);
      avPP4->SetLineWidth(3);
      l4->AddEntry(avPP4,"avg. pp","l");
      avPP4->Draw("same");
      
      l4->Draw("same");
      
      TLatex *  panel4 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %2.1f GeV/c^{2} (%2.1f-%2.1f)}{av. PbPb =%2.1f GeV/c^{2} (%2.1f-%2.1f)}",avPP*1000,minPP*1000,maxPP*1000,avPbPb*1000,minPbPb*1000,maxPbPb*1000));
      panel4->SetNDC();
      panel4->SetTextFont(42);
      panel4->SetTextSize(0.05);
      panel4->SetLineWidth(2);
      panel4->Draw("same");
      
      TLatex *  texts = new TLatex(0.4913793,0.2494759,isCBin ? Form("%s ; 0-100%%",rapRange.Data()) : Form("%s ; %s",rapRange.Data(),centRange.Data()));
      texts->SetNDC();
      texts->SetTextFont(42);
      texts->SetTextSize(0.06708595);
      texts->SetLineWidth(2);
      texts->Draw("same");
      
      TLatex *  text2s = new TLatex(0.4913793,0.17,fName);
      text2s->SetNDC();
      text2s->SetTextFont(42);
      text2s->SetTextSize(0.06708595);
      text2s->SetLineWidth(2);
      text2s->Draw("same");
      
      c4->SaveAs(Form("%s/sVsPt_%s.pdf",sPath.Data(),fName));
      aSave->Add(c4);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c5 = new TCanvas(Form("rsVsPt_%s",fName),"canvas",81,98,700,504);
      c5->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c5->SetFillColor(0);
      c5->SetBorderMode(0);
      c5->SetBorderSize(2);
      c5->SetLeftMargin(0.1479885);
      c5->SetRightMargin(0.05172414);
      c5->SetTopMargin(0.08421053);
      c5->SetBottomMargin(0.1157895);
      c5->SetFrameBorderMode(0);
      c5->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbPt_rsigma->GetN(),gPbPbPt_rsigma->GetY()),TMath::MaxElement(gPbPbInt_rsigma->GetN(),gPbPbInt_rsigma->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbPt_rsigma->GetN(),gPbPbPt_rsigma->GetY()),TMath::MinElement(gPbPbInt_rsigma->GetN(),gPbPbInt_rsigma->GetY()));
      
      maxPP = max(TMath::MaxElement(gPPPt_rsigma->GetN(),gPPPt_rsigma->GetY()),TMath::MaxElement(gPPInt_rsigma->GetN(),gPPInt_rsigma->GetY()));
      minPP = min(TMath::MinElement(gPPPt_rsigma->GetN(),gPPPt_rsigma->GetY()),TMath::MinElement(gPPInt_rsigma->GetN(),gPPInt_rsigma->GetY()));
      
      maxTot = max(maxPbPb,maxPP);
      maxTot = maxTot*1.2;
      minTot = min(minPbPb,minPP);
      minTot = minTot - minTot*0.2;
      
      TH1* hpTdummyrs = new TH1I("hpTdummyrs","#sigma_{2}/#sigma_{1} vs. p_{T};p_{T}(GeV/c);#sigma_{2}/#sigma_{1}",50, 0.,50.);
      hpTdummyrs->GetYaxis()->SetRangeUser(minTot,maxTot);
      hpTdummyrs->GetYaxis()->SetLabelSize(0.05);
      hpTdummyrs->GetYaxis()->SetTitleSize(0.06);
      hpTdummyrs->GetXaxis()->SetLabelSize(0.05);
      hpTdummyrs->GetXaxis()->SetTitleSize(0.06);
      hpTdummyrs->GetXaxis()->SetTitleOffset(0.8);
      hpTdummyrs->SetStats(0);
      hpTdummyrs->Draw();
      
      TLegend* l5 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_rsigma->SetMarkerStyle(21);
      gPbPbInt_rsigma->SetMarkerColor(2);
      gPbPbInt_rsigma->SetLineColor(2);
      l5->AddEntry(gPbPbInt_rsigma,"int. PbPb","lp");
      gPbPbInt_rsigma->Draw("psame");
      
      gPPInt_rsigma->SetMarkerStyle(21);
      gPPInt_rsigma->SetMarkerColor(4);
      gPPInt_rsigma->SetLineColor(4);
      l5->AddEntry(gPPInt_rsigma,"int. pp","lp");
      gPPInt_rsigma->Draw("psame");
      
      gPbPbPt_rsigma->SetMarkerStyle(20);
      gPbPbPt_rsigma->SetMarkerColor(2);
      gPbPbPt_rsigma->SetLineColor(2);
      l5->AddEntry(gPbPbPt_rsigma,"diff. PbPb","lp");
      gPbPbPt_rsigma->Draw("psame");
      
      gPPPt_rsigma->SetMarkerStyle(20);
      gPPPt_rsigma->SetMarkerColor(4);
      gPPPt_rsigma->SetLineColor(4);
      l5->AddEntry(gPPPt_rsigma,"diff. pp","lp");
      gPPPt_rsigma->Draw("psame");
      
      avPbPb = gPbPbPt_rsigma->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbPt_rsigma);
      TLine* avPbPb5 = new TLine(hpTdummyrs->GetXaxis()->GetXmin(),avPbPb,hpTdummyrs->GetXaxis()->GetXmax(),avPbPb);
      avPbPb5->SetLineColor(2);
      avPbPb5->SetLineStyle(2);
      avPbPb5->SetLineWidth(3);
      l5->AddEntry(avPbPb5,"avg. PbPb","l");
      avPbPb5->Draw("same");
      
      avPP = gPPPt_rsigma->GetMean(2);
      if (wMean) avPP = getWaverage(gPPPt_rsigma);
      TLine* avPP5 = new TLine(hpTdummyrs->GetXaxis()->GetXmin(),avPP,hpTdummyrs->GetXaxis()->GetXmax(),avPP);
      avPP5->SetLineColor(4);
      avPP5->SetLineStyle(2);
      avPP5->SetLineWidth(3);
      l5->AddEntry(avPP4,"avg. pp","l");
      avPP5->Draw("same");
      
      l5->Draw("same");
      
      TLatex *  panel5 = new TLatex(0.2,0.8,Form("#splitline{av. pp = %0.2f (%0.2f-%0.2f)}{av. PbPb =%0.2f (%0.2f-%0.2f)}",avPP,minPP,maxPP,avPbPb,minPbPb,maxPbPb));
      panel5->SetNDC();
      panel5->SetTextFont(42);
      panel5->SetTextSize(0.05);
      panel5->SetLineWidth(2);
      panel5->Draw("same");
      
      TLatex *  textrs = new TLatex(0.4913793,0.2494759,isCBin ? Form("%s ; 0-100%%",rapRange.Data()) : Form("%s ; %s",rapRange.Data(),centRange.Data()));
      textrs->SetNDC();
      textrs->SetTextFont(42);
      textrs->SetTextSize(0.06708595);
      textrs->SetLineWidth(2);
      textrs->Draw("same");
      
      TLatex *  text2rs = new TLatex(0.4913793,0.17,fName);
      text2rs->SetNDC();
      text2rs->SetTextFont(42);
      text2rs->SetTextSize(0.06708595);
      text2rs->SetLineWidth(2);
      text2rs->Draw("same");
      
      c5->SaveAs(Form("%s/rsVsPt_%s.pdf",sPath.Data(),fName));
      aSave->Add(c5);
    }
    ///////////////////////////////
    /////////////Cent dep////////////
    ///////////////////////////////
    if (isCBin)
    {
      TCanvas* c6 = new TCanvas(Form("alphaVsCent_%s",fName),"canvas",81,98,700,504);
      c6->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c6->SetFillColor(0);
      c6->SetBorderMode(0);
      c6->SetBorderSize(2);
      c6->SetLeftMargin(0.1479885);
      c6->SetRightMargin(0.05172414);
      c6->SetTopMargin(0.08421053);
      c6->SetBottomMargin(0.1157895);
      c6->SetFrameBorderMode(0);
      c6->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbCent_alpha->GetN(),gPbPbCent_alpha->GetY()),TMath::MaxElement(gPbPbInt_alpha->GetN(),gPbPbInt_alpha->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbCent_alpha->GetN(),gPbPbCent_alpha->GetY()),TMath::MinElement(gPbPbInt_alpha->GetN(),gPbPbInt_alpha->GetY()));
      
      maxTot = maxPbPb*1.2;
      minTot = minPbPb - minPbPb*0.2;
      
      TH1* hCentdummya = new TH1I("hCentdummya","#alpha vs. cent.;cent.;#alpha",200, 0.,200);
      hCentdummya->GetYaxis()->SetRangeUser(minTot,maxTot);
      hCentdummya->GetYaxis()->SetLabelSize(0.05);
      hCentdummya->GetYaxis()->SetTitleSize(0.06);
      hCentdummya->GetXaxis()->SetLabelSize(0.05);
      hCentdummya->GetXaxis()->SetTitleSize(0.06);
      hCentdummya->GetXaxis()->SetTitleOffset(0.8);
      hCentdummya->SetStats(0);
      hCentdummya->Draw();
      
      TLegend* l6 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_alpha->SetMarkerStyle(21);
      gPbPbInt_alpha->SetMarkerColor(2);
      gPbPbInt_alpha->SetLineColor(2);
      l6->AddEntry(gPbPbInt_alpha,"int. PbPb","lp");
      gPbPbInt_alpha->Draw("psame");
      
      gPPInt_alpha->SetMarkerStyle(21);
      gPPInt_alpha->SetMarkerColor(4);
      gPPInt_alpha->SetLineColor(4);
      l6->AddEntry(gPPInt_alpha,"int. pp","lp");
      gPPInt_alpha->Draw("psame");
      
      gPbPbCent_alpha->SetMarkerStyle(20);
      gPbPbCent_alpha->SetMarkerColor(2);
      gPbPbCent_alpha->SetLineColor(2);
      l6->AddEntry(gPbPbCent_alpha,"diff. PbPb","lp");
      gPbPbCent_alpha->Draw("psame");
      
      avPbPb = gPbPbCent_alpha->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbCent_alpha);
      TLine* avPbPb6 = new TLine(hCentdummya->GetXaxis()->GetXmin(),avPbPb,hCentdummya->GetXaxis()->GetXmax(),avPbPb);
      avPbPb6->SetLineColor(2);
      avPbPb6->SetLineStyle(2);
      avPbPb6->SetLineWidth(3);
      l6->AddEntry(avPbPb6,"avg. PbPb","l");
      avPbPb6->Draw("same");
      
      l6->Draw("same");
      
      TLatex *  panel6 = new TLatex(0.2,0.8,Form("av. PbPb =%0.2f (%0.2f-%0.2f)",avPbPb,minPbPb,maxPbPb));
      panel6->SetNDC();
      panel6->SetTextFont(42);
      panel6->SetTextSize(0.05);
      panel6->SetLineWidth(2);
      panel6->Draw("same");
      
      TLatex *  textCent = new TLatex(0.4126074,0.2494759,Form("%s ; %s",rapRange.Data(),ptRange.Data()));
      textCent->SetNDC();
      textCent->SetTextFont(42);
      textCent->SetTextSize(0.06708595);
      textCent->SetLineWidth(2);
      textCent->Draw("same");
      
      TLatex *  text2Cent = new TLatex(0.4126074,0.17,fName);
      text2Cent->SetNDC();
      text2Cent->SetTextFont(42);
      text2Cent->SetTextSize(0.06708595);
      text2Cent->SetLineWidth(2);
      text2Cent->Draw("same");
      
      c6->SaveAs(Form("%s/alphaVsCent_%s.pdf",sPath.Data(),fName));
      aSave->Add(c6);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c7 = new TCanvas(Form("fVsCent_%s",fName),"canvas",81,98,700,504);
      c7->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c7->SetFillColor(0);
      c7->SetBorderMode(0);
      c7->SetBorderSize(2);
      c7->SetLeftMargin(0.1479885);
      c7->SetRightMargin(0.05172414);
      c7->SetTopMargin(0.08421053);
      c7->SetBottomMargin(0.1157895);
      c7->SetFrameBorderMode(0);
      c7->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbCent_f->GetN(),gPbPbCent_f->GetY()),TMath::MaxElement(gPbPbInt_f->GetN(),gPbPbInt_f->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbCent_f->GetN(),gPbPbCent_f->GetY()),TMath::MinElement(gPbPbInt_f->GetN(),gPbPbInt_f->GetY()));
      
      maxTot = maxPbPb*1.2;
      minTot = minPbPb - minPbPb*0.2;
      
      TH1* hCentdummyf = new TH1I("hCentdummyf","f vs. cent.;cent.;f",200, 0.,200);
      hCentdummyf->GetYaxis()->SetRangeUser(minTot,maxTot);
      hCentdummyf->GetYaxis()->SetLabelSize(0.05);
      hCentdummyf->GetYaxis()->SetTitleSize(0.06);
      hCentdummyf->GetXaxis()->SetLabelSize(0.05);
      hCentdummyf->GetXaxis()->SetTitleSize(0.06);
      hCentdummyf->GetXaxis()->SetTitleOffset(0.8);
      hCentdummyf->SetStats(0);
      hCentdummyf->Draw();
      
      TLegend* l7 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_f->SetMarkerStyle(21);
      gPbPbInt_f->SetMarkerColor(2);
      gPbPbInt_f->SetLineColor(2);
      l7->AddEntry(gPbPbInt_f,"int. PbPb","lp");
      gPbPbInt_f->Draw("psame");
      
      gPPInt_f->SetMarkerStyle(21);
      gPPInt_f->SetMarkerColor(4);
      gPPInt_f->SetLineColor(4);
      l7->AddEntry(gPPInt_f,"int. pp","lp");
      gPPInt_f->Draw("psame");
      
      gPbPbCent_f->SetMarkerStyle(20);
      gPbPbCent_f->SetMarkerColor(2);
      gPbPbCent_f->SetLineColor(2);
      l7->AddEntry(gPbPbCent_f,"diff. PbPb","lp");
      gPbPbCent_f->Draw("psame");
      
      avPbPb = gPbPbCent_f->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbCent_f);
      TLine* avPbPb7 = new TLine(hCentdummyf->GetXaxis()->GetXmin(),avPbPb,hCentdummyf->GetXaxis()->GetXmax(),avPbPb);
      avPbPb7->SetLineColor(2);
      avPbPb7->SetLineStyle(2);
      avPbPb7->SetLineWidth(3);
      l7->AddEntry(avPbPb7,"avg. PbPb","l");
      avPbPb7->Draw("same");
      
      l7->Draw("same");
      
      TLatex *  panel7 = new TLatex(0.2,0.8,Form("av. PbPb =%0.2f (%0.2f-%0.2f)",avPbPb,minPbPb,maxPbPb));
      panel7->SetNDC();
      panel7->SetTextFont(42);
      panel7->SetTextSize(0.05);
      panel7->SetLineWidth(2);
      panel7->Draw("same");
      
      TLatex *  textCentf = new TLatex(0.4126074,0.2494759,Form("%s ; %s",rapRange.Data(),ptRange.Data()));
      textCentf->SetNDC();
      textCentf->SetTextFont(42);
      textCentf->SetTextSize(0.06708595);
      textCentf->SetLineWidth(2);
      textCentf->Draw("same");
      
      TLatex *  text2Centf = new TLatex(0.4126074,0.17,fName);
      text2Centf->SetNDC();
      text2Centf->SetTextFont(42);
      text2Centf->SetTextSize(0.06708595);
      text2Centf->SetLineWidth(2);
      text2Centf->Draw("same");
      
      c7->SaveAs(Form("%s/fVsCent_%s.pdf",sPath.Data(),fName));
      aSave->Add(c7);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c8 = new TCanvas(Form("nVsCent_%s",fName),"canvas",81,98,700,504);
      c8->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c8->SetFillColor(0);
      c8->SetBorderMode(0);
      c8->SetBorderSize(2);
      c8->SetLeftMargin(0.1479885);
      c8->SetRightMargin(0.05172414);
      c8->SetTopMargin(0.08421053);
      c8->SetBottomMargin(0.1157895);
      c8->SetFrameBorderMode(0);
      c8->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbCent_n->GetN(),gPbPbCent_n->GetY()),TMath::MaxElement(gPbPbInt_n->GetN(),gPbPbInt_n->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbCent_n->GetN(),gPbPbCent_n->GetY()),TMath::MinElement(gPbPbInt_n->GetN(),gPbPbInt_n->GetY()));
      
      maxTot = maxPbPb*1.2;
      minTot = minPbPb - minPbPb*0.2;
      
      TH1* hCentdummyn = new TH1I("hCentdummyn","n vs. cent.;cent.;n",200, 0.,200);
      hCentdummyn->GetYaxis()->SetRangeUser(minTot,maxTot);
      hCentdummyn->GetYaxis()->SetLabelSize(0.05);
      hCentdummyn->GetYaxis()->SetTitleSize(0.06);
      hCentdummyn->GetXaxis()->SetLabelSize(0.05);
      hCentdummyn->GetXaxis()->SetTitleSize(0.06);
      hCentdummyn->GetXaxis()->SetTitleOffset(0.8);
      hCentdummyn->SetStats(0);
      hCentdummyn->Draw();
      
      TLegend* l8 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_n->SetMarkerStyle(21);
      gPbPbInt_n->SetMarkerColor(2);
      gPbPbInt_n->SetLineColor(2);
      l8->AddEntry(gPbPbInt_n,"int. PbPb","lp");
      gPbPbInt_n->Draw("psame");
      
      gPPInt_n->SetMarkerStyle(21);
      gPPInt_n->SetMarkerColor(4);
      gPPInt_n->SetLineColor(4);
      l8->AddEntry(gPPInt_n,"int. pp","lp");
      gPPInt_n->Draw("psame");
      
      gPbPbCent_n->SetMarkerStyle(20);
      gPbPbCent_n->SetMarkerColor(2);
      gPbPbCent_n->SetLineColor(2);
      l8->AddEntry(gPbPbCent_n,"diff. PbPb","lp");
      gPbPbCent_n->Draw("psame");
      
      avPbPb = gPbPbCent_n->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbCent_n);
      TLine* avPbPb8 = new TLine(hCentdummyn->GetXaxis()->GetXmin(),avPbPb,hCentdummyn->GetXaxis()->GetXmax(),avPbPb);
      avPbPb8->SetLineColor(2);
      avPbPb8->SetLineStyle(2);
      avPbPb8->SetLineWidth(3);
      l8->AddEntry(avPbPb8,"avg. PbPb","l");
      avPbPb8->Draw("same");
      
      l8->Draw("same");
      
      TLatex *  panel8 = new TLatex(0.2,0.8,Form("av. PbPb =%0.2f (%0.2f-%0.2f)",avPbPb,minPbPb,maxPbPb));
      panel8->SetNDC();
      panel8->SetTextFont(42);
      panel8->SetTextSize(0.05);
      panel8->SetLineWidth(2);
      panel8->Draw("same");
      
      TLatex *  textCentn = new TLatex(0.4126074,0.2494759,Form("%s ; %s",rapRange.Data(),ptRange.Data()));
      textCentn->SetNDC();
      textCentn->SetTextFont(42);
      textCentn->SetTextSize(0.06708595);
      textCentn->SetLineWidth(2);
      textCentn->Draw("same");
      
      TLatex *  text2Centn = new TLatex(0.4126074,0.17,fName);
      text2Centn->SetNDC();
      text2Centn->SetTextFont(42);
      text2Centn->SetTextSize(0.06708595);
      text2Centn->SetLineWidth(2);
      text2Centn->Draw("same");
      
      c8->SaveAs(Form("%s/nVsCent_%s.pdf",sPath.Data(),fName));
      aSave->Add(c8);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c9 = new TCanvas(Form("sVsCent_%s",fName),"canvas",81,98,700,504);
      c9->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c9->SetFillColor(0);
      c9->SetBorderMode(0);
      c9->SetBorderSize(2);
      c9->SetLeftMargin(0.1479885);
      c9->SetRightMargin(0.05172414);
      c9->SetTopMargin(0.08421053);
      c9->SetBottomMargin(0.1157895);
      c9->SetFrameBorderMode(0);
      c9->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbCent_sigma->GetN(),gPbPbCent_sigma->GetY()),TMath::MaxElement(gPbPbInt_sigma->GetN(),gPbPbInt_sigma->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbCent_sigma->GetN(),gPbPbCent_sigma->GetY()),TMath::MinElement(gPbPbInt_sigma->GetN(),gPbPbInt_sigma->GetY()));
      
      maxTot = maxPbPb*1.2;
      minTot = minPbPb - minPbPb*0.2;
      
      TH1* hCentdummys = new TH1I("hCentdummys","#sigma vs. cent.;cent.;#sigma_{1} (MeV/c^{2})",200, 0.,200);
      hCentdummys->GetYaxis()->SetRangeUser(minTot,maxTot);
      hCentdummys->GetYaxis()->SetLabelSize(0.05);
      hCentdummys->GetYaxis()->SetTitleSize(0.06);
      hCentdummys->GetXaxis()->SetLabelSize(0.05);
      hCentdummys->GetXaxis()->SetTitleSize(0.06);
      hCentdummys->GetXaxis()->SetTitleOffset(0.8);
      hCentdummys->SetStats(0);
      hCentdummys->Draw();
      
      TLegend* l9 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_sigma->SetMarkerStyle(21);
      gPbPbInt_sigma->SetMarkerColor(2);
      gPbPbInt_sigma->SetLineColor(2);
      l9->AddEntry(gPbPbInt_sigma,"int. PbPb","lp");
      gPbPbInt_sigma->Draw("psame");
      
      gPPInt_sigma->SetMarkerStyle(21);
      gPPInt_sigma->SetMarkerColor(4);
      gPPInt_sigma->SetLineColor(4);
      l9->AddEntry(gPPInt_sigma,"int. pp","lp");
      gPPInt_sigma->Draw("psame");
      
      gPbPbCent_sigma->SetMarkerStyle(20);
      gPbPbCent_sigma->SetMarkerColor(2);
      gPbPbCent_sigma->SetLineColor(2);
      l9->AddEntry(gPbPbCent_sigma,"diff. PbPb","lp");
      gPbPbCent_sigma->Draw("psame");
      
      avPbPb = gPbPbCent_sigma->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbCent_sigma);
      TLine* avPbPb9 = new TLine(hCentdummys->GetXaxis()->GetXmin(),avPbPb,hCentdummys->GetXaxis()->GetXmax(),avPbPb);
      avPbPb9->SetLineColor(2);
      avPbPb9->SetLineStyle(2);
      avPbPb9->SetLineWidth(3);
      l9->AddEntry(avPbPb9,"avg. PbPb","l");
      avPbPb9->Draw("same");
      
      l9->Draw("same");
      
      TLatex *  panel9 = new TLatex(0.2,0.8,Form("av. PbPb =%2.1f GeV/c^{2} (%2.1f-%2.1f)",avPbPb*1000,minPbPb*1000,maxPbPb*1000));
      panel9->SetNDC();
      panel9->SetTextFont(42);
      panel9->SetTextSize(0.05);
      panel9->SetLineWidth(2);
      panel9->Draw("same");
      
      TLatex *  textCents = new TLatex(0.4126074,0.2494759,Form("%s ; %s",rapRange.Data(),ptRange.Data()));
      textCents->SetNDC();
      textCents->SetTextFont(42);
      textCents->SetTextSize(0.06708595);
      textCents->SetLineWidth(2);
      textCents->Draw("same");
      
      TLatex *  text2Cents = new TLatex(0.4126074,0.17,fName);
      text2Cents->SetNDC();
      text2Cents->SetTextFont(42);
      text2Cents->SetTextSize(0.06708595);
      text2Cents->SetLineWidth(2);
      text2Cents->Draw("same");
      
      c9->SaveAs(Form("%s/sVsCent_%s.pdf",sPath.Data(),fName));
      aSave->Add(c9);
      
      ///////////////////////////////
      ///////////////////////////////
      
      TCanvas* c10 = new TCanvas(Form("rsVsCent_%s",fName),"canvas",81,98,700,504);
      c10->Range(-5.547577,-0.004024593,31.93896,0.03073326);
      c10->SetFillColor(0);
      c10->SetBorderMode(0);
      c10->SetBorderSize(2);
      c10->SetLeftMargin(0.1479885);
      c10->SetRightMargin(0.05172414);
      c10->SetTopMargin(0.08421053);
      c10->SetBottomMargin(0.1157895);
      c10->SetFrameBorderMode(0);
      c10->SetFrameBorderMode(0);
      
      maxPbPb = max(TMath::MaxElement(gPbPbCent_rsigma->GetN(),gPbPbCent_rsigma->GetY()),TMath::MaxElement(gPbPbInt_rsigma->GetN(),gPbPbInt_rsigma->GetY()));
      minPbPb = min(TMath::MinElement(gPbPbCent_rsigma->GetN(),gPbPbCent_rsigma->GetY()),TMath::MinElement(gPbPbInt_rsigma->GetN(),gPbPbInt_rsigma->GetY()));
      
      maxTot = maxPbPb*1.2;
      minTot = minPbPb - minPbPb*0.2;
      
      TH1* hCentdummyrs = new TH1I("hCentdummyrs","#sigma_{2}/#sigma_{1} vs. cent.;cent.;#sigma_{2}/#sigma_{1}",200, 0.,200);
      hCentdummyrs->GetYaxis()->SetRangeUser(minTot,maxTot);
      hCentdummyrs->GetYaxis()->SetLabelSize(0.05);
      hCentdummyrs->GetYaxis()->SetTitleSize(0.06);
      hCentdummyrs->GetXaxis()->SetLabelSize(0.05);
      hCentdummyrs->GetXaxis()->SetTitleSize(0.06);
      hCentdummyrs->GetXaxis()->SetTitleOffset(0.8);
      hCentdummyrs->SetStats(0);
      hCentdummyrs->Draw();
      
      TLegend* l10 = new TLegend(0.7363897,0.6924686,0.9355301,0.8933054);
      
      gPbPbInt_rsigma->SetMarkerStyle(21);
      gPbPbInt_rsigma->SetMarkerColor(2);
      gPbPbInt_rsigma->SetLineColor(2);
      l10->AddEntry(gPbPbInt_rsigma,"int. PbPb","lp");
      gPbPbInt_rsigma->Draw("psame");
      
      gPPInt_rsigma->SetMarkerStyle(21);
      gPPInt_rsigma->SetMarkerColor(4);
      gPPInt_rsigma->SetLineColor(4);
      l10->AddEntry(gPPInt_rsigma,"int. pp","lp");
      gPPInt_rsigma->Draw("psame");
      
      gPbPbCent_rsigma->SetMarkerStyle(20);
      gPbPbCent_rsigma->SetMarkerColor(2);
      gPbPbCent_rsigma->SetLineColor(2);
      l10->AddEntry(gPbPbPt_rsigma,"diff. PbPb","lp");
      gPbPbCent_rsigma->Draw("psame");
      
      avPbPb = gPbPbCent_rsigma->GetMean(2);
      if (wMean) avPbPb = getWaverage(gPbPbCent_rsigma);
      TLine* avPbPb10 = new TLine(hCentdummyrs->GetXaxis()->GetXmin(),avPbPb,hCentdummyrs->GetXaxis()->GetXmax(),avPbPb);
      avPbPb10->SetLineColor(2);
      avPbPb10->SetLineStyle(2);
      avPbPb10->SetLineWidth(3);
      l10->AddEntry(avPbPb10,"avg. PbPb","l");
      avPbPb10->Draw("same");
      
      l10->Draw("same");
      
      TLatex *  panel10 = new TLatex(0.2,0.8,Form("av. PbPb =%0.2f (%0.2f-%0.2f)",avPbPb,minPbPb,maxPbPb));
      panel10->SetNDC();
      panel10->SetTextFont(42);
      panel10->SetTextSize(0.05);
      panel10->SetLineWidth(2);
      panel10->Draw("same");
      
      TLatex *  textCentrs = new TLatex(0.4126074,0.2494759,Form("%s ; %s",rapRange.Data(),ptRange.Data()));
      textCentrs->SetNDC();
      textCentrs->SetTextFont(42);
      textCentrs->SetTextSize(0.06708595);
      textCentrs->SetLineWidth(2);
      textCentrs->Draw("same");
      
      TLatex *  text2Centrs = new TLatex(0.4126074,0.17,fName);
      text2Centrs->SetNDC();
      text2Centrs->SetTextFont(42);
      text2Centrs->SetTextSize(0.06708595);
      text2Centrs->SetLineWidth(2);
      text2Centrs->Draw("same");
      
      c10->SaveAs(Form("%s/rsVsCent_%s.pdf",sPath.Data(),fName));
      aSave->Add(c10);
    }
  }
  
  TFile* fsave = new TFile(Form("%s_%s.root",sFileName.Data(),fName),"RECREATE");
  aSave->Write("mcParamsEvolution");
  fsave->Close();
  
  cout << "Closed " << Form("%s_%s.root",sFileName.Data(),fName) << endl;
}
