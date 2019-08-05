#include "makeAccEff.h" 
//macro to get all AccEff Syst

//Double_t ptbins []= {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 100.0};
Double_t ptbins []= {3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.5, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.25, 9.5, 9.75, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 19., 20.0, 25.0, 30.0, 50, 100.0};
Double_t ybins []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


float min_z = 0.064;
float max_z = 1.0;

float min_jetpt = 10.;
float max_jetpt = 60.;

int nBinZ = 6;

float z_reco_binWidth = (max_z-min_z)*1.0/nBinZ;

void oniaTree::AccEffMisMod(string caseLabel) {
} // end of the function


void oniaTree::AccEffStatToy(int nToys) {
  // Randomise TEfficiency 100 times. Output will be a TObjArray with 100 TH2
  TFile* effFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ"); 
  TEfficiency *eff = (TEfficiency*) effFile->Get(Form("hcorr_Jpsi_%s_%s",isPbPb?"PbPb":"PP",isPr?"pr":"npr"));
  
  TObjArray *arrEffs = new TObjArray(); // Array to store efficiencies
  arrEffs->SetOwner(true);

  TRandom* rnd = new TRandom3(); // For randomisation

  for (int i=0 ; i < nToys ; i++)
    {
      if (i%10==0) cout <<"[INFO] toy "<<i<<"/"<<nToys<<endl;
      TH2* histoTot = (TH2*) eff->GetTotalHistogram()->Clone(Form("hTotToy_%d",i)); // Get corresponding histo and number of bins
      histoTot->Sumw2();
      TH2* histoPass = (TH2*) histoTot->Clone(Form("hPassToy_%d",i)); // Get corresponding histo and number of bins
      histoPass->Sumw2();
      int nBinsX = histoTot->GetNbinsX();
      int nBinsY = histoTot->GetNbinsY();
  
      for (int j = 1 ; j <= nBinsX ; j++)
	{
	  for (int k = 1 ; k <= nBinsY ; k++)
	    {
	      int bin = histoTot->GetBin(j,k);      
	      double effVal = eff->GetEfficiency(bin);
	      int ntot = histoTot->GetBinContent(bin);
	      int newPass = rnd->Binomial(ntot,effVal);
	      histoPass->SetBinContent(histoPass->GetBin(j,k),newPass);
	    }
	}
      TEfficiency* newHeff = (TEfficiency*) eff->Clone(Form("effToy_%d",i)); // Histo to store new 2D eff
      newHeff->SetPassedHistogram(*histoPass, "f");
      newHeff->SetTotalHistogram(*histoTot, "f");
      arrEffs->Add(newHeff);
    }
  effFile->Close();
  gSystem->mkdir("FilesAccxEff/toyMC");
  TFile* fsave = new TFile(Form("FilesAccxEff/toyMC/%sAccEff%dToys_%s.root",isPr?"pr":"npr", nToys,isPbPb?"PbPb":"PP"),"RECREATE");
  arrEffs->Write("accEffArray", TObject::kSingleKey);
  fsave->Close();
}

void oniaTree::AccEffStat(string caseLabel) {
  int jtPtmin=10;
  int jtPtmax=60;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=30; jtPtmax=40;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=20; jtPtmax=30;}
  else if (caseLabel.find("lowerJtPt")!=std::string::npos) {jtPtmin=10; jtPtmax=20;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=40; jtPtmax=50;}
  else if (caseLabel.find("higherJtPt")!=std::string::npos) {jtPtmin=50; jtPtmax=60;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}


  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s_NoBkg_jetR3_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR",isPbPb?"PbPb":"PP"),"READ");
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  trNom->SetBranchAddress("corr_ptw",&corr_ptw);

  TEfficiency* tempPrEff = NULL;
  TEfficiency* tempNprEff = NULL;

  TObjArray* hisArr = new TObjArray(); //hisArr1624->SetOwner(kTRUE);

  TFile *prcorrFile = TFile::Open(Form("FilesAccxEff/toyMC/prAccEff100Toys_%s.root",isPbPb?"PbPb":"PP"));
  TObjArray* prcorArr = (TObjArray*) prcorrFile->Get("accEffArray");

  TH1F* histVar = NULL;

  TFile* nomFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ");
  TEfficiency* prnomEff = (TEfficiency*) nomFile->Get(Form("hcorr_Jpsi_%s_pr",isPbPb?"PbPb":"PP"));

  double bf = 1.0;
  double prEff = 1.0;
  double nprEff = 1.0;
  double totEff = 1.0;

  int nentries = trNom->GetEntries();
  for (int i=0; i<=100; i++)
    {
	histVar = new TH1F (Form("histTot_%d",i), "", nBinZ,  min_z, max_z);

      if (i==0) cout<<"[INFO] Applying the nominal AccxEff"<<endl;
      else cout<<"[INFO] Applying toy "<<i<<"/100"<<endl;
      if (i==0) {
	tempPrEff=prnomEff;
	//tempNprEff=nprnomEff; 
      }
      else {
	tempPrEff=(TEfficiency*) prcorArr->At(i-1);
      }

      for (int jentry=0; jentry<nentries; jentry++) {
	trNom->GetEntry(jentry);
	if (abs(jp_rap)>2.4) continue;
	if (jp_pt<6.5|| jp_pt>100) continue;
	if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;

	prEff = tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt));
	totEff=1.0/prEff;
	totEff = totEff*corr_ptw;

	  histVar->Fill(z, totEff);
      }// end of the tree entries
      hisArr->Add(histVar);
    } // end of the variations
  TFile* fsave = new TFile(Form("FilesAccEff/zHistsWith100AccEffToys_%s_%s.root",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"),"RECREATE");
  hisArr->Write();
  fsave->Close();

  cout<<"[INFO] making csv files"<<endl;
  ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_AccEffStat.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP"));

  fileOut<<"AccxEff stat"<<endl;

  vector<double> vCount;

  vCount.clear();

  TH1F* temp = NULL;

  for (int j=0; j<nBinZ; j++)
    {
      vCount.clear();
      for (int i=0; i<=100; i++) {
	temp = (TH1F*) hisArr->At(i);
	vCount.push_back(temp->GetBinContent(temp->FindBin(j*z_reco_binWidth+min_z+0.001)));
      }
      if (isPbPb)
        fileOut<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 180, "<< rms(vCount,true) << endl;
      else
        fileOut<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 200, "<< rms(vCount,true) << endl;
    }
  fileOut.close();
  prcorrFile->Close();
  nomFile->Close();
  //nprcorrFile->Close();
  trFile->Close();
}


void oniaTree::TnpSyst(string caseLabel) {

  int npt2D = sizeof(ptbins)/sizeof(double)-1;
  int ny2D = sizeof(ybins)/sizeof(double)-1;

  int jtPtmin=10;
  int jtPtmax=60;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=30; jtPtmax=40;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=20; jtPtmax=30;}
  else if (caseLabel.find("lowerJtPt")!=std::string::npos) {jtPtmin=10; jtPtmax=20;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=40; jtPtmax=50;}
  else if (caseLabel.find("higherJtPt")!=std::string::npos) {jtPtmin=50; jtPtmax=60;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s_NoBkg_jetR3_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR",isPbPb?"PbPb":"PP"),"READ"); 
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  trNom->SetBranchAddress("corr_ptw",&corr_ptw);

  cout<<"[INFO] Loading the corrections"<<endl;
  //TFile* prAccFile = TFile::Open(Form("FilesAccxEff/Acc/prAccHists_%s.root",isPbPb?"PbPb":"PP"),"READ"); 
  TFile* prAccFile = TFile::Open(Form("FilesAccxEff/Acc/prAccHists_%s.root","PP"),"READ"); //use pp acceptance for pp and PbPb
  TFile* prEffFile = TFile::Open(Form("FilesAccxEff/Eff/prEffHists_%s.root",isPbPb?"PbPb":"PP"),"READ");

  double prEff = 1.0;
  double totEff = 1.0;

  string corrName_pp [] = 
    {
      "nominal", //nominal
      "muidtrg_plus1sig_syst",
      "muidtrg_minus1sig_syst",
      "muidtrg_plus1sig_stat",
      "muidtrg_minus1sig_stat",
      "glb_plus1sig_syst",
      "glb_minus1sig_syst",
      "glb_plus1sig_stat",
      "glb_minus1sig_stat"
    };
  string corrName_pbpb [] = {
    "nominal",
    "muid_plus1sig_syst",
    "muid_minus1sig_syst",
    "muid_plus1sig_stat",
    "muid_minus1sig_stat",
    "trg_plus1sig_syst",
    "trg_minus1sig_syst",
    "trg_plus1sig_stat",
    "trg_minus1sig_stat",
    "trk_plus1sig_syst",
    "trk_minus1sig_syst",
    "trk_plus1sig_stat",
    "trk_minus1sig_stat"
  };

  string systName_pp [] = {
      "muidtrgSyst",
      "muidtrgStat",
      "glbSyst",
      "glbStat"
    };

  string systName_pbpb [] = {
    "muidSyst",
    "muidStat",
    "trgSyst",
    "trgStat",
    "trkSyst",
    "trkStat"
  };

  int nCorr= sizeof(corrName_pp)/sizeof(corrName_pp[0]);
  int nSyst= sizeof(systName_pp)/sizeof(systName_pp[0]);
  if (isPbPb) {
    nCorr = sizeof(corrName_pbpb)/sizeof(corrName_pbpb[0]);
    nSyst=sizeof(systName_pbpb)/sizeof(systName_pbpb[0]);
  }

  cout <<"nCorr = "<<nCorr<<endl;

  string *systName = new string [nSyst];//new string[sizeof(systName_pp)/sizeof(systName_pp[0])];
  string *corrName = new string [nCorr];//new string[sizeof(corrName_pp)/sizeof(corrName_pp[0])];
  
  cout<<"systName[0]="<<systName[0]<<endl;

  for (int i=0;i<nSyst;i++){
    if (isPbPb)
      systName[i]=systName_pbpb[i];
    else
      systName[i]=systName_pp[i];
    cout<<"systName = "<<systName[i]<<endl;
  }

  for (int i=0;i<nCorr;i++){
    if (isPbPb)
      corrName[i]=corrName_pbpb[i];
    else
      corrName[i]=corrName_pp[i];
    cout<<"corrName = "<<corrName[i]<<endl;
  }
  

  TObjArray *countArr = new TObjArray(15); countArr->SetOwner(true);
  TEfficiency *prcorrTemp = NULL; TEfficiency *nprcorrTemp = NULL; TH1F* countHist = NULL;

  vector<double> vCount;

  TH2F *prAccNum = (TH2F*) prAccFile->Get("hnum_2d_nominal");
  TH2F *prAccDen = (TH2F*) prAccFile->Get("hdeno_2d");
  TH2F *prEffNum = NULL;//(TH2F*) prEffFile->Get("hnum_nominal");
  TH2F *prEffDen = (TH2F*) prEffFile->Get("hdeno_pty");
  prEffDen->Multiply(prAccDen);

  if (!prAccNum||!prAccDen||!prEffNum||!prEffDen) cout<<"histogram not found!";
  cout<<"[INFO] all the hists are loaded"<<endl;

  for(int i=0; i<nCorr; i++) {
    prEffNum = (TH2F*) prEffFile->Get(Form("hnum_%s", corrName[i].c_str()));
    prEffNum->Multiply(prAccNum);

    cout<<"syst i = "<<i<<endl;
    prcorrTemp = new TEfficiency(Form("prCorr_%s", corrName[i].c_str()), "AccxEff(y,pt); y; pt; eff", ny2D, ybins, npt2D, ptbins);
    prcorrTemp->SetStatisticOption(TEfficiency::kBBayesian);
    prcorrTemp->SetPassedHistogram(*prEffNum,"f");
    prcorrTemp->SetTotalHistogram(*prEffDen,"f");

    countHist = new TH1F (Form("his_%s",corrName[i].c_str()), Form("tnp %s",corrName[i].c_str()), nBinZ, min_z, max_z);

    //filling the histograms
    cout<<"[INFO] Filling the histograms"<<endl;
    int nentries = trNom->GetEntries();
    for (int jentry=0; jentry<nentries; jentry++) {
      if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<endl;
      trNom->GetEntry(jentry);
      
      if (abs(jp_rap)>2.4) continue;
      if (jp_pt<6.5|| jp_pt>100) continue;
      if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;

      prEff = prcorrTemp->GetEfficiency(prcorrTemp->FindFixBin(jp_rap,jp_pt));
      totEff = 1.0/prEff;
      totEff = totEff*corr_ptw;
      countHist->Fill(z, totEff);
    }
    countArr->Add(countHist);
  } 

  TFile *fsave = new TFile(Form("FilesAccxEff/tnpSystHists_%s_%s.root",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"),"RECREATE");
  countArr->Write();
  fsave->Close();
  cout<<"[INFO] Getting the ratios and filling the csv files"<<endl;
  
    ////// filling the different syst files
  for (int i = 0; i<nSyst ; i++)
    { 
      ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt",isPbPb?"PbPb":"PP",systName[i].c_str()));

      fileOut << "tnp " << systName[i] << endl;
      vCount.clear();

      for (int j=0; j<nBinZ; j++){
	vCount.clear();
	countHist = (TH1F*) countArr->At(0);
	double binX = countHist->FindBin(j*z_reco_binWidth+min_z+0.001);
	vCount.push_back((double)(countHist->GetBinContent(binX)));
	
	countHist = (TH1F*) countArr->At(2*i+1);
	vCount.push_back((double)(countHist->GetBinContent(binX)));
	
	countHist = (TH1F*) countArr->At(2*i+2);
	vCount.push_back((double)(countHist->GetBinContent(binX)));
	
	if (isPbPb)
	  fileOut<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 180, "<< maxdiff(vCount, 1) << endl;
	else
	  fileOut<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 200, "<< maxdiff(vCount, 1) << endl;
      }
      fileOut.close();
    }
  
  ofstream fileAll(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_fullAccEff.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt",isPbPb?"PbPb":"PP"));
  fileAll<<"AccxEff"<<endl;
  double val1 [] = {0,0,0,0,0,0,0,0};
  for (int j=0; j<nBinZ; j++){
    double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "AccEffStat"), min_z+j*z_reco_binWidth, min_z+(j+1)*z_reco_binWidth, 0, 2.4);
    cout <<"for AccEffStat v1 = "<<v1<<endl;
    double syst=v1*v1;
    
    for (int i=0; i<nSyst; i++){

      v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", systName[i].c_str()), min_z+j*z_reco_binWidth, min_z+(j+1)*z_reco_binWidth, 0, 2.4);
      cout <<"for "<< systName[i] <<" v1 = "<<v1<<endl;
      syst = syst + v1*v1;
    }
    if (isPbPb)
      fileAll<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 180, "<< sqrt(syst) << endl;
    else
      fileAll<< "0, 2.4, 6.5, 100, "<<min_z+j*z_reco_binWidth<<", "<< min_z+(j+1)*z_reco_binWidth<<", 0, 200, "<< sqrt(syst) << endl;
  }
  fileAll.close();
}


void oniaTree::TnpToy(int min, int max) {
}

void oniaTree::TnpStat(string caseLabel) {

}


void oniaTree::AccEffSyst_all() {

  //string systName [] = {"midJtPt_016","midJtPt_1624", "lowJtPt_016","lowJtPt_1624", "highJtPt_016","highJtPt_1624", "NoJets_total"};
  string systName [] = {"midJtPt", "lowJtPt", "highJtPt", "lowerJtPt","higherJtPt"};
  int systN = sizeof(systName)/sizeof(systName[0]);

  //cout<<"--------[INFO] Sarting 100 toys for AccxEff stats----------"<<endl;
  //AccEffStatToy(100);

  for (int i=0; i<systN; i++) {
    cout<<"--------[INFO] Sarting AccEffStat for "<<systName[i]<<"----------"<<endl;
    AccEffStat(systName[i]);
    cout<<"--------[INFO] Sarting TnpSyst for "<<systName[i]<<"----------"<<endl;
    TnpSyst(systName[i]);
    //cout<<"--------[INFO] Sarting TnpStat for "<<systName[i]<<"----------"<<endl;
    //TnpStat(systName[i]);
    //cout<<"--------[INFO] Sarting AccEffMisMod for "<<systName[i]<<"----------"<<endl;
    //AccEffMisMod(systName[i]);
    cout<<"--------[INFO] All systematics done for "<<systName[i]<<"----------"<<endl;
  }
}
