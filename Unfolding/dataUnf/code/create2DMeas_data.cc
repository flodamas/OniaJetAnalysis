#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
#endif

void create2DMeas_data(bool doPrompt = true, bool doPbPb = false){
  bool statErr = !systErr;
  //if ((doCent&&doPeri) ||!doPbPb) {doCent = false; doPeri = false;}

  gSystem->mkdir(Form("%s/dataUnf/plots",unfPath.c_str()));
  gSystem->mkdir(Form("%s/dataUnf/data_results",unfPath.c_str()));
    
  string filename0 = "";
  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename4 = "";
  string filename5 = "";
  
  filename0 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_lowestJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowestJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  filename1 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_lowerJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowerJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  filename2 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_lowJtPt/ctauMass/DATA/fitsPars/unfoldingInput_lowJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  filename3 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_midJtPt/ctauMass/DATA/fitsPars/unfoldingInput_midJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  filename4 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_highJtPt/ctauMass/DATA/fitsPars/unfoldingInput_highJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  filename5 = Form("/home/llr/cms/diab/JpsiInJetsPbPb/Fitter/Output/DataFits/DataFits_higherJtPt/ctauMass/DATA/fitsPars/unfoldingInput_higherJtPt_rap_jetR3_%sErr.root",statErr?"stat":"syst");
  
  TFile *file0 = new TFile(filename0.c_str(),"READ");
  TFile *file1 = new TFile(filename1.c_str(),"READ");
  TFile *file2 = new TFile(filename2.c_str(),"READ");
  TFile *file3 = new TFile(filename3.c_str(),"READ");
  TFile *file4 = new TFile(filename4.c_str(),"READ");
  TFile *file5 = new TFile(filename5.c_str(),"READ");

  TH1F *h_lowestBin;
  TH1F *h_lowerBin;
  TH1F *h_lowBin;
  TH1F *h_centBin;
  TH1F *h_highBin;
  TH1F *h_higherBin;


  h_lowestBin = (TH1F*)file0->Get(Form("%sHist_%s_lowestJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  h_lowerBin = (TH1F*)file1->Get(Form("%sHist_%s_lowerJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  h_lowBin = (TH1F*)file2->Get(Form("%sHist_%s_lowJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  h_centBin = (TH1F*)file3->Get(Form("%sHist_%s_midJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  h_highBin = (TH1F*)file4->Get(Form("%sHist_%s_highJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  h_higherBin = (TH1F*)file5->Get(Form("%sHist_%s_higherJtPt_rap%s_%sErr",doPrompt?"pr":"npr",doPbPb?"PbPb":"PP",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst"));
  
  TH2D *h_Meas = new TH2D("h_Meas","h_Meas",nBinZ_reco,min_z,max_z,nBinJet_reco,min_jetpt,max_jetpt);
  h_Meas->Sumw2();
  
  int nbins = h_lowBin->GetNbinsX();

  float lowestBinVals[nbins];
  float lowerBinVals[nbins];
  float lowBinVals[nbins];
  float centBinVals[nbins];
  float highBinVals[nbins];
  float higherBinVals[nbins];

  float lowestBinErrs[nbins];
  float lowerBinErrs[nbins];
  float lowBinErrs[nbins];
  float centBinErrs[nbins];
  float highBinErrs[nbins];
  float higherBinErrs[nbins];
  
  double totIntegral = 0;
  double totErr = 0;
  for(int i = 1; i < nbins+1; i++){
    //if (doPbPb && (h_lowestBin->GetBinCenter(i)<0.65)) {lowestBinVals[i-1] = 0; lowestBinErrs[i-1]=0; continue;}
    //if ((h_lowestBin->GetBinCenter(i)<0.65)) {lowestBinVals[i-1] = 0; lowestBinErrs[i-1]=0; continue;}
    cout << Form("h_lowestBin : %i ", i) << h_lowestBin->GetBinContent(i) << endl;
    lowestBinVals[i-1] = h_lowestBin->GetBinContent(i);
    lowestBinErrs[i-1] = h_lowestBin->GetBinError(i);
    totIntegral=totIntegral+lowestBinVals[i-1];
    totErr=totErr+lowestBinErrs[i-1]*lowestBinErrs[i-1];
  }
  
  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_lowerBin : %i ", i) << h_lowerBin->GetBinContent(i) << endl;
    lowerBinVals[i-1] = h_lowerBin->GetBinContent(i);
    lowerBinErrs[i-1] = h_lowerBin->GetBinError(i);
    totIntegral=totIntegral+lowerBinVals[i-1];
    totErr=totErr+lowerBinErrs[i-1]*lowerBinErrs[i-1];
  }

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_lowBin : %i ", i) << h_lowBin->GetBinContent(i) << endl;
    lowBinVals[i-1] = h_lowBin->GetBinContent(i);
    lowBinErrs[i-1] = h_lowBin->GetBinError(i);
    totIntegral=totIntegral+lowBinVals[i-1];
    totErr=totErr+lowBinErrs[i-1]*lowBinErrs[i-1];
  }

  cout << "***" << endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_centBin : %i ", i) << h_centBin->GetBinContent(i) << endl;
    centBinVals[i-1] = h_centBin->GetBinContent(i);
    centBinErrs[i-1] = h_centBin->GetBinError(i);
    totIntegral=totIntegral+centBinVals[i-1];
    totErr=totErr+centBinErrs[i-1]*centBinErrs[i-1];

  }

  cout << "***"<< endl;

  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_highBin : %i ", i) << h_highBin->GetBinContent(i) << endl;
    highBinVals[i-1] = h_highBin->GetBinContent(i);
    highBinErrs[i-1] = h_highBin->GetBinError(i);
    totIntegral=totIntegral+highBinVals[i-1];
    totErr=totErr+highBinErrs[i-1]*highBinErrs[i-1];
  }


  for(int i = 1; i < nbins+1; i++){
    cout << Form("h_higherBin : %i ", i) << h_higherBin->GetBinContent(i) << endl;
    higherBinVals[i-1] = h_higherBin->GetBinContent(i);
    higherBinErrs[i-1] = h_higherBin->GetBinError(i);
    totIntegral=totIntegral+higherBinVals[i-1];
    totErr=totErr+higherBinErrs[i-1]*higherBinErrs[i-1];

    if(noProbFit && doPbPb && (h_higherBin->GetBinCenter(i)<0.22) && doCent) {
      cout <<"the bin error is too big:"<<higherBinErrs[i-1]/higherBinVals[i-1]<<", setting it to:0.1"<<endl; 
      higherBinErrs[i-1]=0.1*higherBinVals[i-1];
    }
    else if(doPbPb && (h_higherBin->GetBinCenter(i)<0.22) && doCent)
      cout <<"the bin error is:"<<higherBinErrs[i-1]/higherBinVals[i-1]<<endl;
  }
  totErr = sqrt(totErr);

  cout << "***"<< endl;

  float content = 0.0;
  float errVal = 0.0;
  
  for(int icount=0; icount < nBinJet_reco; icount++){
    for(int jcount=0; jcount < nBinZ_reco; jcount++){

      if(icount == 0) {
	content = lowestBinVals[jcount];
	errVal = lowestBinErrs[jcount];
      }

      if(icount == 1) {
	content = lowerBinVals[jcount];
	errVal = lowerBinErrs[jcount];
      }
      if(icount == 2) {
	content = lowBinVals[jcount];
	errVal = lowBinErrs[jcount];
      }

      if(icount == 3) {
	content = centBinVals[jcount];
	errVal = centBinErrs[jcount];
      }

      if(icount == 4) {
	content = highBinVals[jcount];
	errVal = highBinErrs[jcount];
      }
      if(icount == 5) {
	content = higherBinVals[jcount];
	errVal = higherBinErrs[jcount];
      }

      //cout << "jet bin = " << icount << " z bin = " << jcount << endl;
      //cout << "content = " << content << endl;
      
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

  can0->SaveAs(Form("%s/dataUnf/plots/data_%s_z_distr_%s_%sErr.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",statErr?"stat":"syst"));

  TCanvas * can1 = new TCanvas ("can1","can1",600,600);
  can1->SetRightMargin(0.2);
  
  h_Meas->SetTitle(Form("%s %s",doPrompt?"prompt":"nonprompt",doPbPb?"PbPb":"PP"));
  
  h_Meas->SetStats(0);
  h_Meas->GetXaxis()->SetTitle("z");
  h_Meas->GetYaxis()->SetTitle("jet p_{T}");
  h_Meas->GetYaxis()->SetTitleOffset(1.2);
  h_Meas->Draw("TEXTcolz");
  can1->SaveAs(Form("%s/dataUnf/plots/data_%s_meas_%s_%sErr.pdf",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",statErr?"stat":"syst"));

  
  string outputfile = "";
  outputfile = Form("%s/dataUnf/data_results/meas_%s_data_%s%s_%sErrs.root",unfPath.c_str(),doPbPb?"PbPb":"PP",doPrompt?"prompt":"nonprompt",doCent?"_centBin":doPeri?"_periBin":"",statErr?"stat":"syst");
  
  TFile *file_data_meas = new TFile(outputfile.c_str(),"RECREATE");

  h_Meas->Write();
  h_lowestBin->Write("zMeasLowestBin");
  h_lowerBin->Write("zMeasLowerBin");
  h_lowBin->Write("zMeasLowBin");
  h_centBin->Write("zMeasCentBin");
  h_highBin->Write("zMeasHighBin");
  h_higherBin->Write("zMeasHigherBin");
  file_data_meas->Close();

  cout <<"for "<<(doPbPb?"PbPb":"PP")<<" ,total number of events = "<<totIntegral<<", err = "<<totErr<<endl;
  if (!doPbPb) {
    totIntegral = totIntegral*normPP;
    totErr = totErr*normPP;
  }
  else if(!doCent && !doPeri) {
    totIntegral = totIntegral*normPbPb;
    totErr = totErr*normPbPb;
  }
  cout <<"for "<<(doPbPb?"PbPb":"PP")<<" ,Normalized total number of events = "<<totIntegral<<", err = "<<totErr<<endl;
}
