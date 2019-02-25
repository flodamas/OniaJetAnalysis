#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "Riostream.h"
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>


void MixAccEff() 
{
  TFile* prFile = TFile::Open("Utilities/pr_correction_AccEff.root");
  TFile* nprFile = TFile::Open("Utilities/npr_correction_AccEff.root");
  cout<<"[INFO] Getting prompt AccEff"<<endl;
  TEfficiency* prEff = (TEfficiency*) prFile->Get("hcorr_Jpsi_PP");
  cout<<"[INFO] Getting nonprompt AccEff"<<endl;
  TEfficiency* nprEff = (TEfficiency*) nprFile->Get("hcorr_Jpsi_PP");
  prFile->Close();
  nprFile->Close();
  cout<<"[INFO] Saving both AccEff"<<endl;
  TFile* fsave = new TFile("../Fitter/Input/correction_AccEff.root","RECREATE");
  prEff->SetName("hcorr_Jpsi_PP_pr");
  nprEff->SetName("hcorr_Jpsi_PP_npr");
  prEff->Write("hcorr_Jpsi_PP_pr");
  nprEff->Write("hcorr_Jpsi_PP_npr");
  fsave->Close();

}
void BfracExtraction()
{
  Double_t ptbins [] = {3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35};
  int nptbins = ((sizeof(ptbins)/sizeof(double))-1);
  Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
  int netabins = ((sizeof(etabins)/sizeof(double))-1);
  Double_t bfracbins [] = {3, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35,50};
  int nbfracbins = ((sizeof(bfracbins)/sizeof(double))-1); 

  /////////fit the b-fraction
  TF1  *fitfun0 = new TF1("fitfun0","exp([0]+[1]*pow(x,1)+[2]*pow(x,2)+[3]*pow(x,3))", 3, 35);
  TF1  *fitfun1 = new TF1("fitfun1","[0]+[1]*pow(x,0.5)+[2]*x+[3]*pow(x,1.5)", 3, 50);
  TF1  *fitfun2 = new TF1("fitfun2","exp([0]+[1]*pow(x,1)+[2]*pow(x,2)+[3]*pow(x,3))", 3, 50);
  TF1  *fitfun3 = new TF1("fitfun3","pol4", 3, 50);
  TF1  *fitfun4 = new TF1("fitfun4","[0]*pow(1-[1]/x,[2])", 3, 50);
  //fitfun4->SetParLimits(0, 0, 1);
  //fitfun4->SetParLimits(1, 3, 7);
  TH1F *bhist = new TH1F ("bhist", "b fracion", nbfracbins, bfracbins);
  gStyle->SetOptStat(0);

  ifstream in;
  in.open("Utilities/bfraction.dat");

  Float_t x,y,z,w;
  Int_t nlines = 0;

  while (1)
    {
      in >> x >> y >> z >> w;
      if (!in.good()) break;
      //if (nlines < 5) printf("x=%8f, z=%8f\n",x,z);
      bhist->SetBinContent(bhist->FindBin((x+y)*0.5), z);
      bhist->SetBinError(bhist->FindBin((x+y)*0.5), w);
      nlines++;
    }

  in.close();

  bhist->Fit("fitfun0","R");
  bhist->Fit("fitfun1","R");
  bhist->Fit("fitfun2","R");
  //bhist->Fit("fitfun3","R");
  //bhist->Fit("fitfun4");
  TCanvas *c = new TCanvas("c","",1000,1000);
  bhist->Draw();
  bhist->SetMinimum(0);
  bhist->SetMaximum(0.7);
  fitfun0->SetLineColor(kRed);
  fitfun1->SetLineColor(kBlue);
  fitfun2->SetLineColor(kGreen);
  fitfun3->SetLineColor(kMagenta);
  fitfun4->SetLineColor(kBlue);
  TLegend *leg = new TLegend(0.15, 0.75, 0.25, 0.85);
  leg->AddEntry(fitfun0, "exp(3d ord pol)", "l");
  leg->AddEntry(fitfun1, "Francois's pol", "l");
  leg->AddEntry(fitfun2, "exp(3d ord pol)", "l");
  //leg->AddEntry(fitfun3, "pol4", "l");
  //leg->AddEntry(fitfun4, "pol4", "l");
  leg->SetBorderSize(0);
  TLine* lim = new TLine(35, 0, 35, 0.7);
  lim->SetLineColor(kBlack);
  lim->SetLineStyle(2); 
  fitfun0->Draw("same");
  fitfun1->Draw("same");
  fitfun2->Draw("same");
  //fitfun3->Draw("same");
  //fitfun4->Draw("same");
  leg->Draw("same");
  lim->Draw("same");

  cout<<"[INFO] the function to use is: f(x) = exp(["<<fitfun0->GetParameter(0)<<"]+["<<fitfun0->GetParameter(1)<<"]*pow(x,1)+["<<fitfun0->GetParameter(2)<<"]*pow(x,2)+["<<fitfun0->GetParameter(3)<<"]*pow(x,3))"<<endl;
  c->SaveAs("Utilities/bFracFit.pdf");
}/// end of mixAccEff function
