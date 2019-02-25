///// this is a macro to check the consistecy of the histograms for the efficiency
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;
using namespace  RooFit;

void checkHist(const char*file = "") 
{
  TFile * histFile = TFile::Open(file);
  TH2F * pass = (TH2F*) histFile->Get("hcorr_his_num");
  TH2F * tot = (TH2F*) histFile->Get("hcorr_his_deno");
  TEfficiency* eff = (TEfficiency*) histFile->Get("hcorr_Jpsi_PP");
  int nbins = (pass->GetNbinsX()+2)*(pass->GetNbinsY()+2);
  cout<< "pass xbins = " << pass->GetNbinsX() << " ybins = " << pass->GetNbinsY() << "; tot xbins = " << tot->GetNbinsX() << " ybins = " << tot->GetNbinsY() << endl;
  for (int i = 0; i < pass->GetNbinsX()+2; i++)
    for (int j = 0; j < pass->GetNbinsY()+2; j++)
      //if (pass->GetBinContent(i, j) > tot->GetBinContent(i, j))
      cout<< "In the bin with y = "<< pass->GetXaxis()->GetBinCenter(i) << " and pt = " << pass->GetYaxis()->GetBinCenter(j) << " (pass = " << pass->GetBinContent(i, j) << "+-"<< pass->GetBinError(i, j)<<") ; (tot = " << tot->GetBinContent(i, j) << "+-"<< tot->GetBinError(i, j)<< ") ; (eff = " << eff->GetEfficiency(eff->FindFixBin(pass->GetXaxis()->GetBinCenter(i), pass->GetYaxis()->GetBinCenter(j)))<< "+- errL = " << eff->GetEfficiencyErrorLow(eff->FindFixBin(pass->GetXaxis()->GetBinCenter(i), pass->GetYaxis()->GetBinCenter(j)))<<"+- errH = "<< eff->GetEfficiencyErrorUp(eff->FindFixBin(pass->GetXaxis()->GetBinCenter(i), pass->GetYaxis()->GetBinCenter(j)))<<")"<<endl;
}
