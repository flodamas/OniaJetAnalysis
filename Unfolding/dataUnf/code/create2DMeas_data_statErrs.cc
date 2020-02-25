#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void create2DMeas_data_statErrs(bool doPrompt = true, bool doPbPb = false){

  gSystem->mkdir(Form("%s/dataUnf/plots",unfPath.c_str()));
  gSystem->mkdir(Form("%s/dataUnf/data_results",unfPath.c_str()));
    
  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string filename5 = "";
  
  filename1 = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_lowerJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowerJtPt_rap_jetR3_statErr.root";
  filename2 = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_lowJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowJtPt_rap_jetR3_statErr.root";
  filename3 = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_midJtPt/ctauMass/DATA/fitsPars/unfoldingInput_midJtPt_rap_jetR3_statErr.root";
  filename4 = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_highJtPt/ctauMass/DATA/fitsPars/unfoldingInput_highJtPt_rap_jetR3_statErr.root";
  filename5 = "/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_higherJtPt/ctauMass/DATA/fitsPars/unfoldingInput_higherJtPt_rap_jetR3_statErr.root";
  //filename1 = "../../../../Output/FromLizardo/DataFits_lowerJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowerJtPt_rap_jetR3_statErr.root";
  //filename2 = "../../../../Output/FromLizardo/DataFits_lowJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowJtPt_rap_jetR3_statErr.root";
  //filename3 = "../../../../Output/FromLizardo/DataFits_midJtPt/ctauMass/DATA/fitsPars/unfoldingInput_midJtPt_rap_jetR3_statErr.root";
  //filename4 = "../../../../Output/FromLizardo/DataFits_highJtPt/ctauMass/DATA/fitsPars/unfoldingInput_highJtPt_rap_jetR3_statErr.root";
  //filename5 = "../../../../Output/FromLizardo/DataFits_higherJtPt/ctauMass/DATA/fitsPars/unfoldingInput_higherJtPt_rap_jetR3_statErr.root";
  
  TFile *file1 = new TFile(filename1.c_str(),"UPDATE");
  TFile *file2 = new TFile(filename2.c_str(),"UPDATE");
  TFile *file3 = new TFile(filename3.c_str(),"UPDATE");
  TFile *file4 = new TFile(filename4.c_str(),"UPDATE");
  TFile *file5 = new TFile(filename5.c_str(),"UPDATE");

  TH1F *h_lowerBin;
  TH1F *h_lowBin;
  TH1F *h_centBin;
  TH1F *h_highBin;
  TH1F *h_higherBin;


  h_lowerBin = (TH1F*)file1->Get(Form("%sHist_%s_lowerJtPt_rap_statErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP"));
  h_lowBin = (TH1F*)file2->Get(Form("%sHist_%s_lowJtPt_rap_statErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP"));
  h_centBin = (TH1F*)file3->Get(Form("%sHist_%s_midJtPt_rap_statErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP"));
  h_highBin = (TH1F*)file4->Get(Form("%sHist_%s_highJtPt_rap_statErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP"));
  h_higherBin = (TH1F*)file5->Get(Form("%sHist_%s_higherJtPt_rap_statErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP"));
  
  TH2D *h_Meas = new TH2D("h_Meas","h_Meas",nBinZ_reco,min_z,max_z,nBinJet_reco,min_jetpt,max_jetpt);
  h_Meas->Sumw2();
  
  int nbins = h_lowBin->GetNbinsX();

  float lowerBinVals[nbins];
  float lowBinVals[nbins];
  float centBinVals[nbins];
  float highBinVals[nbins];
  float higherBinVals[nbins];

  float lowerBinErrs[nbins];
  float lowBinErrs[nbins];
  float centBinErrs[nbins];
  float highBinErrs[nbins];
  float higherBinErrs[nbins];
  
  
  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_lowBin : %i ", i) << h_lowBin->GetBinContent(i) << endl;
    lowerBinVals[i-1] = h_lowerBin->GetBinContent(i);
    lowerBinErrs[i-1] = h_lowerBin->GetBinError(i);
  }

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_lowBin : %i ", i) << h_lowBin->GetBinContent(i) << endl;
    lowBinVals[i-1] = h_lowBin->GetBinContent(i);
    lowBinErrs[i-1] = h_lowBin->GetBinError(i);
  }

  cout << "***" << endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_centBin : %i ", i) << h_centBin->GetBinContent(i) << endl;
    centBinVals[i-1] = h_centBin->GetBinContent(i);
    centBinErrs[i-1] = h_centBin->GetBinError(i);
  }

  cout << "***"<< endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_highBin : %i ", i) << h_highBin->GetBinContent(i) << endl;
    highBinVals[i-1] = h_highBin->GetBinContent(i);
    highBinErrs[i-1] = h_highBin->GetBinError(i);
  }


  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_highBin : %i ", i) << h_highBin->GetBinContent(i) << endl;
    higherBinVals[i-1] = h_higherBin->GetBinContent(i);
    higherBinErrs[i-1] = h_higherBin->GetBinError(i);
  }

  cout << "***"<< endl;

  float content = 0.0;
  float errVal = 0.0;
  
  for(int icount=0; icount < nBinJet_reco; icount++){
    for(int jcount=0; jcount < nBinZ_reco; jcount++){

      if(icount == 0) {
	content = lowerBinVals[jcount];
	errVal = lowerBinErrs[jcount];
      }
      if(icount == 1) {
	content = lowBinVals[jcount];
	errVal = lowBinErrs[jcount];
      }

      if(icount == 2) {
	content = centBinVals[jcount];
	errVal = centBinErrs[jcount];
      }

      if(icount == 3) {
	content = highBinVals[jcount];
	errVal = highBinErrs[jcount];
      }
      if(icount == 4) {
	content = higherBinVals[jcount];
	errVal = higherBinErrs[jcount];
      }

      cout << "jet bin = " << icount << " z bin = " << jcount << endl;
      cout << "content = " << content << endl;
      
      h_Meas->SetBinContent(jcount+1,icount+1,content);
      h_Meas->SetBinError(jcount+1,icount+1,errVal);
    }
  }

  TCanvas * can0 = new TCanvas ("can0","can0",1500,500);
  can0->Divide(5,1);

  can0->cd(1);
  h_lowerBin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+0*jetPt_reco_binWidth), (int)(min_jetpt+1*jetPt_reco_binWidth)));
  h_lowerBin->Draw();

  can0->cd(2);
  h_lowBin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+1*jetPt_reco_binWidth), (int)(min_jetpt+2*jetPt_reco_binWidth)));
  h_lowBin->Draw();

  can0->cd(3);
  h_centBin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+2*jetPt_reco_binWidth), (int)(min_jetpt+3*jetPt_reco_binWidth)));
  h_centBin->Draw();

  can0->cd(4);
  h_highBin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+3*jetPt_reco_binWidth), (int)(min_jetpt+4*jetPt_reco_binWidth)));
  h_highBin->Draw();

  can0->cd(5);
  h_higherBin->SetTitle(Form("%d < jet p_{T} < %d",(int) (min_jetpt+4*jetPt_reco_binWidth), (int)(min_jetpt+5*jetPt_reco_binWidth)));
  h_higherBin->Draw();

  can0->SaveAs(Form("../plots/data_%s_z_distr_%s_statErr.pdf",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt"));

  TCanvas * can1 = new TCanvas ("can1","can1",600,600);
  can1->SetRightMargin(0.2);
  
  h_Meas->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"PP"));
  
  h_Meas->SetStats(0);
  h_Meas->GetXaxis()->SetTitle("z");
  h_Meas->GetYaxis()->SetTitle("jet p_{T}");
  h_Meas->GetYaxis()->SetTitleOffset(1.2);
  h_Meas->Draw("TEXTcolz");
  can1->SaveAs(Form("../plots/data_%s_meas_%s_statErr.pdf",doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt"));

  
  string outputfile = "";
  outputfile = Form("%s/dataUnf/data_results/meas_%s_data_%s_statErrs.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt");
  
  TFile *file_data_meas = new TFile(outputfile.c_str(),"RECREATE");

  h_Meas->Write();
  h_lowerBin->Write("zMeasLowerBin");
  h_lowBin->Write("zMeasLowBin");
  h_centBin->Write("zMeasCentBin");
  h_highBin->Write("zMeasHighBin");
  h_higherBin->Write("zMeasHigherBin");

}
