//macro to get all AccEff Syst
Double_t ptbins []= {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50};
Double_t ybins []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

void oniaTree::AccEffMisMod(string caseLabel = "") {
  double totBins [] = {0, 1.6, 2.4};

  int jtPtmin;
  int jtPtmax;
  bool isMid = false;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=25; jtPtmax=35;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=15; jtPtmax=25;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=35; jtPtmax=45;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}

  if (caseLabel.find("016")!=std::string::npos) isMid = true;


  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR"),"READ"); //always apply it on non prompt
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  if (isPr) trNom->SetBranchAddress("corr_ptw",&corr_ptw);

  TFile *corrFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ");
  TEfficiency* prHist = (TEfficiency*) corrFile->Get("hcorr_Jpsi_PP_pr");
  TEfficiency* nprHist = (TEfficiency*) corrFile->Get("hcorr_Jpsi_PP_npr");

  TH1F* countHist = NULL;
  TH1F* nomHist = NULL;

  if (jtPtmin==0) {
    countHist = new TH1F ("countHist", "", 2, totBins);
    nomHist = new TH1F ("nomHist", "", 2, totBins);
  }
  else if (isMid){
    countHist = new TH1F ("countHist", "", 7, 0.02, 1);
    nomHist = new TH1F ("nomHist", "", 7, 0.02, 1);
  }
  else {
    countHist = new TH1F ("countHist", "", 5, 0, 1);
    nomHist = new TH1F ("nomHist", "", 5, 0, 1);
  }
  TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
  double prEff = 1.0; double nprEff = 1.0; double bf = 1.0; double varEff = 1.0; double totEff = 1.0;

  cout<<"[INFO] Processing the entries"<<endl;
  int nentries = trNom->GetEntries();
  for (int jentry=0; jentry<nentries; jentry++) {
    trNom->GetEntry(jentry);

    if (abs(jp_rap)>2.4) continue;
    if (jp_pt<3 || jp_pt>35) continue;
    if (isMid && jtPtmin>10 && abs(jp_rap)>1.6) continue;
    if (!isMid && jtPtmin>10 && abs(jp_rap)<=1.6) continue; 
    if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;

    bf = bfrac->Eval(jp_pt);
    prEff = prHist->GetEfficiency(prHist->FindFixBin(jp_rap,jp_pt));
    nprEff = nprHist->GetEfficiency(nprHist->FindFixBin(jp_rap,jp_pt));
    totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
    if (isPr)
      totEff = totEff*corr_ptw;

    if (jtPtmin > 10) {
      nomHist->Fill(z, totEff);
      if (isPr) countHist->Fill(z, corr_ptw*1.0/prEff);
      else countHist->Fill(z, 1.0/nprEff);
    }
    else {
      nomHist->Fill(jp_rap, totEff);
      if (isPr) countHist->Fill(jp_rap, corr_ptw*1.0/prEff);
      else countHist->Fill(jp_rap, 1.0/nprEff);
    }
  }// end of the tree entries
  cout<<"[INFO] filling csv files"<<endl;
  ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_AccEffMisMod.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt"));
  fileOut<<"AccxEff mismodeling"<<endl;
  double val, var;
  if (jtPtmin > 10) {
    int jmax =5;
    if (isMid) jmax =6;
    for (int j=1; j<jmax; j++)
      {
	if (isMid){
	  val=nomHist->GetBinContent(nomHist->FindBin(j*0.14+0.2));
	  var=countHist->GetBinContent(countHist->FindBin(j*0.14+0.2));
	  fileOut<< "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14<<", 0, 200, "<< 1.*(val-var)/val << endl;
	}
	else 
	  {
	    val=nomHist->GetBinContent(nomHist->FindBin(j*0.2+0.1));
	    var=countHist->GetBinContent(countHist->FindBin(j*0.2+0.1));
	    fileOut<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< 1.*(val-var)/val << endl;
	  }
      }
  }
  else {
      val=nomHist->GetBinContent(nomHist->FindBin(0));
      var=countHist->GetBinContent(countHist->FindBin(0));
      fileOut<< "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< 1.*(val-var)/val << endl;
      val=nomHist->GetBinContent(nomHist->FindBin(2));
      var=countHist->GetBinContent(countHist->FindBin(2));
      fileOut<< "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< 1.*(val-var)/val << endl;
    }
  fileOut.close();
  //corrFile->Close();
  trFile->Close();
} // end of the function


void oniaTree::AccEffStatToy(int nToys) {
  // Randomise TEfficiency 100 times. Output will be a TObjArray with 100 TH2
  TFile* effFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ"); 
  TEfficiency *eff = (TEfficiency*) effFile->Get(Form("hcorr_Jpsi_PP_%s",isPr?"pr":"npr"));
  
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
  TFile* fsave = new TFile(Form("FilesAccxEff/toyMC/%sAccEff%dToys.root",isPr?"pr":"npr", nToys),"RECREATE");
  arrEffs->Write("accEffArray", TObject::kSingleKey);
  fsave->Close();
}

void oniaTree::AccEffStat(string caseLabel = "") {
  double totBins [] = {0, 1.6, 2.4};

  int jtPtmin;
  int jtPtmax;
  bool isMid = false;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=25; jtPtmax=35;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=15; jtPtmax=25;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=35; jtPtmax=45;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}

  if (caseLabel.find("016")!=std::string::npos) isMid = true;

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR"),"READ");
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  if (isPr) trNom->SetBranchAddress("corr_ptw",&corr_ptw);

  TEfficiency* tempPrEff = NULL;
  TEfficiency* tempNprEff = NULL;

  TObjArray* hisArr = new TObjArray(); //hisArr1624->SetOwner(kTRUE);

  TFile *prcorrFile = TFile::Open("FilesAccxEff/toyMC/prAccEff100Toys.root");
  TFile *nprcorrFile = TFile::Open("FilesAccxEff/toyMC/nprAccEff100Toys.root");
  TObjArray* prcorArr = (TObjArray*) prcorrFile->Get("accEffArray");
  TObjArray* nprcorArr = (TObjArray*) nprcorrFile->Get("accEffArray");

  TH1F* histVar = NULL;

  TFile* nomFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ");
  TEfficiency* prnomEff = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP_pr");
  TEfficiency* nprnomEff = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP_npr");
  TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
  double bf = 1.0;
  double prEff = 1.0;
  double nprEff = 1.0;
  double totEff = 1.0;

  int nentries = trNom->GetEntries();
  for (int i=0; i<=100; i++)
    {
      if (jtPtmin==0)
	histVar = new TH1F (Form("histTot_%d",i), "", 2, totBins);
      else if (isMid)
	histVar = new TH1F (Form("histTot_%d",i), "", 7, 0.02, 1);
      else 
	histVar = new TH1F (Form("histTot_%d",i), "", 5, 0, 1);

      if (i==0) cout<<"[INFO] Applying the nominal AccxEff"<<endl;
      else cout<<"[INFO] Applying toy "<<i<<"/100"<<endl;
      if (i==0) {
	tempPrEff=prnomEff;
	tempNprEff=nprnomEff; 
      }
      else {
	tempPrEff=(TEfficiency*) prcorArr->At(i-1);
	tempNprEff=(TEfficiency*) nprcorArr->At(i-1);
      }

      for (int jentry=0; jentry<nentries; jentry++) {
	trNom->GetEntry(jentry);
	if (abs(jp_rap)>2.4) continue;
	if (jp_pt<3|| jp_pt>35) continue;
	if (isMid && jtPtmin>10 && abs(jp_rap)>1.6) continue;
	if (!isMid && jtPtmin>10 && abs(jp_rap)<=1.6) continue; 
	if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;

	bf = bfrac->Eval(jp_pt);
	prEff = tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt));
	nprEff = tempNprEff->GetEfficiency(tempNprEff->FindFixBin(jp_rap,jp_pt));
	totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
	if (isPr) totEff = totEff*corr_ptw;

	if (jtPtmin > 10)
	  histVar->Fill(z, totEff);
	else 
	  histVar->Fill(abs(jp_rap), totEff);
      }// end of the tree entries
      hisArr->Add(histVar);
    } // end of the variations

  cout<<"[INFO] making csv files"<<endl;
  ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_AccEffStat.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt"));

  fileOut<<"AccxEff stat"<<endl;

  vector<double> vCount;

  vCount.clear();

  TH1F* temp = NULL;

  int jmax =5;
  if (isMid) jmax = 6;
  if (jtPtmin > 10) {
    for (int j=1; j<jmax; j++)
      {
	vCount.clear();
	for (int i=0; i<=100; i++) {
	  temp = (TH1F*) hisArr->At(i);
          vCount.push_back(temp->GetBinContent(temp->FindBin(j*0.2+0.1)));
	}
	if (isMid)
	  fileOut<< "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< rms(vCount,true) << endl;
	else 
	  fileOut<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< j*0.2+0.2 <<", 0, 200, "<< rms(vCount,true) << endl;
      }
  }
  else {
    vCount.clear();
    for (int i=0; i<=100; i++) {
      temp = (TH1F*) hisArr->At(i);
      vCount.push_back(temp->GetBinContent(temp->FindBin(1)));
    }
      fileOut<< "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<<  rms(vCount,true)<< endl;
      vCount.clear();
      for (int i=0; i<=100; i++) {
	temp = (TH1F*) hisArr->At(i);
	vCount.push_back(temp->GetBinContent(temp->FindBin(2)));
      }
      fileOut<< "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< rms(vCount,true) << endl;
    }
  
  fileOut.close();
  prcorrFile->Close();
  nomFile->Close();
  nprcorrFile->Close();
  trFile->Close();
}


void oniaTree::TnpSyst(string caseLabel = "") {
  double totBins [] = {0, 1.6, 2.4};

  int npt2D = sizeof(ptbins)/sizeof(double)-1;
  int ny2D = sizeof(ybins)/sizeof(double)-1;

  int jtPtmin;
  int jtPtmax;
  bool isMid = false;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=25; jtPtmax=35;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=15; jtPtmax=25;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=35; jtPtmax=45;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}

  if (caseLabel.find("016")!=std::string::npos) isMid = true;

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR"),"READ"); 
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  if (isPr) trNom->SetBranchAddress("corr_ptw",&corr_ptw);

  cout<<"[INFO] Loading the corrections"<<endl;
  TFile* prAccFile = TFile::Open("FilesAccxEff/Acc/prAccHists.root","READ"); 
  TFile* nprAccFile = TFile::Open("FilesAccxEff/Acc/nprAccHists.root","READ");
  TFile* prEffFile = TFile::Open("FilesAccxEff/Eff/prEffHists.root","READ");
  TFile* nprEffFile = TFile::Open("FilesAccxEff/Eff/nprEffHists.root","READ");

  TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
  double bf = 1.0;
  double prEff = 1.0;
  double nprEff = 1.0;
  double totEff = 1.0;

  string corrName [] = 
    {
      "nominal", //nominal
      "binned",
      "trg_plus1sig",
      "trg_minus1sig",
      "muid_sta",
      "muid",
      "muid_plus1sig",
      "muid_minus1sig",
      "sta",
      "sta_plus1sig",
      "sta_minus1sig",
      "trk_plus1sig",
      "trk_minus1sig"
    };

  string systName [] =
    {                                          
      "binned",
      "trgSyst",
      "muidSyst",
      "staSyst",
      "trkSyst"
  };
  //TObjArray *prcorrHis = new TObjArray(15);
  //TObjArray *nprcorrHis = new TObjArray(15);
  TObjArray *countArr = new TObjArray(15);
  TEfficiency *prcorrTemp = NULL; TEfficiency *nprcorrTemp = NULL; TH1F* countHist = NULL;

  vector<double> vCount;

  TH2F *prAccNum = (TH2F*) prAccFile->Get("hnum_2d_nominal");
  TH2F *prAccDen = (TH2F*) prAccFile->Get("hdeno_2d");
  TH2F *nprAccNum = (TH2F*) nprAccFile->Get("hnum_2d_nominal");
  TH2F *nprAccDen = (TH2F*) nprAccFile->Get("hdeno_2d");
  TH2F *prEffNum = (TH2F*) prEffFile->Get("hnum_2d_nominal");
  TH2F *prEffDen = (TH2F*) prEffFile->Get("hdeno_2d");
  TH2F *nprEffNum = (TH2F*) nprEffFile->Get("hnum_2d_nominal");
  TH2F *nprEffDen = (TH2F*) nprEffFile->Get("hdeno_2d");

  for(int i=0; i<13; i++) {
    prEffNum = (TH2F*) prEffFile->Get(Form("hnum_2d_%s", corrName[i].c_str()));
    nprEffNum = (TH2F*) nprEffFile->Get(Form("hnum_2d_%s", corrName[i].c_str()));
    prEffDen = (TH2F*) prEffFile->Get("hdeno_2d");
    nprEffDen = (TH2F*) nprEffFile->Get("hdeno_2d");

    prEffNum->Multiply(prAccNum);
    prEffDen->Multiply(prAccDen);
    nprEffNum->Multiply(nprAccNum);
    nprEffDen->Multiply(nprAccDen);

    prcorrTemp = new TEfficiency(Form("prCorr_%s", corrName[i].c_str()), "AccxEff(y,pt); y; pt; eff", ny2D, ybins, npt2D, ptbins);
    prcorrTemp->SetStatisticOption(TEfficiency::kBBayesian);
    prcorrTemp->SetPassedHistogram(*prEffNum,"f");
    prcorrTemp->SetTotalHistogram(*prEffDen,"f");
    //prcorrHis->Add(corrTemp);

    nprcorrTemp = new TEfficiency(Form("nprCorr_%s", corrName[i].c_str()), "AccxEff(y,pt); y; pt; eff", ny2D, ybins, npt2D, ptbins);
    nprcorrTemp->SetStatisticOption(TEfficiency::kBBayesian);
    nprcorrTemp->SetPassedHistogram(*nprEffNum,"f");
    nprcorrTemp->SetTotalHistogram(*nprEffDen,"f");
    //nprcorrHis->Add(corrTemp);

    if (jtPtmin > 10)
      if (isMid)
	countHist = new TH1F (Form("his_%s",corrName[i].c_str()), Form("tnp %s",corrName[i].c_str()), 7, 0.02, 1);
      else
	countHist = new TH1F (Form("his_%s",corrName[i].c_str()), Form("tnp %s",corrName[i].c_str()), 5, 0, 1);
    else 
      countHist = new TH1F (Form("his_%s",corrName[i].c_str()), Form("tnp %s",corrName[i].c_str()), 2, totBins);

    //filling the histograms
    cout<<"[INFO] Filling the histograms"<<endl;
    int nentries = trNom->GetEntries();
    for (int jentry=0; jentry<nentries; jentry++) {
      if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<endl;
      trNom->GetEntry(jentry);
      
      if (abs(jp_rap)>2.4) continue;
      if (jp_pt<3|| jp_pt>35) continue;
      if (isMid && jtPtmin>10 && abs(jp_rap)>1.6) continue;
      if (!isMid && jtPtmin>10 && abs(jp_rap)<=1.6) continue; 
      if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;
      
      bf = bfrac->Eval(jp_pt);
      prEff = prcorrTemp->GetEfficiency(prcorrTemp->FindFixBin(jp_rap,jp_pt));
      nprEff = nprcorrTemp->GetEfficiency(nprcorrTemp->FindFixBin(jp_rap,jp_pt));
      totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
      if (isPr) totEff = totEff*corr_ptw;
      
      if (jtPtmin > 10)
	countHist->Fill(z, totEff);
      else 
	countHist->Fill(abs(jp_rap), totEff);
    }
    countArr->Add(countHist);
  }

 
    cout<<"[INFO] Getting the ratios and filling the csv files"<<endl;

    ////// filling the different syst files
    for (int i = 0; i<5 ; i++)
      { 
	ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_tnp%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt",systName[i].c_str()));
	string systLabel [] =
	  {
	    "binned",
	    "trg Syst",
	    "muid Syst",
	    "sta Syst",
	    "trk Syst"
	  };
	
	fileOut << "tnp " << systLabel[i] << endl;
	vCount.clear();
	int jmax;
	if (jtPtmin>10) 
	  if (isMid)
	    jmax =6;
	  else 
	    jmax =5;

	else jmax =3;
 
	for (int j=1; j<jmax; j++){
	  vCount.clear();
	  countHist = (TH1F*) countArr->At(0);
	  double binX = 0;
	  if (jtPtmin>10)
	    if (isMid)
	    binX =countHist->FindBin(j*0.14+0.20);
	    else 
	      binX =countHist->FindBin(j*0.2+0.1);
	  else
	    binX = countHist->FindBin(j);

	    vCount.push_back((double)(countHist->GetBinContent(binX)));
	    vCount.push_back((double)(countHist->GetBinContent(binX)));

	  if (i==0) { /////binned
	    countHist = (TH1F*) countArr->At(1);
	    vCount.push_back((double)(countHist->GetBinContent(binX)));
	  }
	  if (i==1) { /////trgSyst
	    countHist = (TH1F*) countArr->At(2);
            vCount.push_back((double)(countHist->GetBinContent(binX)));

	    countHist = (TH1F*) countArr->At(3);
            vCount.push_back((double)(countHist->GetBinContent(binX)));
	  }
	  if (i==2) { ////muidSyst
	    countHist = (TH1F*) countArr->At(6);
            vCount.push_back((double)(countHist->GetBinContent(binX)));

	    countHist = (TH1F*) countArr->At(7);
            vCount.push_back((double)(countHist->GetBinContent(binX)));
	  }
	  if (i==3) { //////staSyst
	    countHist = (TH1F*) countArr->At(9);
            vCount.push_back((double)(countHist->GetBinContent(binX)));

	    countHist = (TH1F*) countArr->At(10);
            vCount.push_back((double)(countHist->GetBinContent(binX)));
	  }
	  if (i==4) { ///////trkSyst
	    countHist = (TH1F*) countArr->At(11);
            vCount.push_back((double)(countHist->GetBinContent(binX)));

	    countHist = (TH1F*) countArr->At(12);
            vCount.push_back((double)(countHist->GetBinContent(binX)));
	  }
	  if (jtPtmin>10)
	    if (isMid)
	      fileOut << "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< maxdiff(vCount, 1) << endl;
	    else
	      fileOut << "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< maxdiff(vCount, 1) << endl;
	  else 
	    if (j==1)
	      fileOut << "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< maxdiff(vCount, 1) << endl;
	    else
	      fileOut << "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< maxdiff(vCount, 1) << endl;
	}
	fileOut.close();
      }
}


void oniaTree::TnpToy(int min, int max) {

  int npt2D = sizeof(ptbins)/sizeof(double)-1;
  int ny2D = sizeof(ybins)/sizeof(double)-1;

  TFile* AccFile = TFile::Open(Form("FilesAccxEff/Acc/%sAccHists.root",isPr?"pr":"npr"),"READ");
  TFile* EffFile = TFile::Open(Form("FilesAccxEff/Eff/%sEffHists.root",isPr?"pr":"npr"),"READ");

  TH2F *AccNum = (TH2F*) AccFile->Get("hnum_2d_nominal");
  TH2F *AccDen = (TH2F*) AccFile->Get("hdeno_2d");
  TH2F *EffNum = (TH2F*) EffFile->Get("hnum_2d_nominal");
  TH2F *EffDen = (TH2F*) EffFile->Get("hdeno_2d");
  
  TObjArray *trgEff = new TObjArray(); trgEff->SetOwner(kTRUE);
  TObjArray *muidEff = new TObjArray(); muidEff->SetOwner(kTRUE);
  TObjArray *staEff = new TObjArray(); staEff->SetOwner(kTRUE);
  
  TEfficiency* trgTemp = NULL;
  TEfficiency* muidTemp = NULL;
  TEfficiency* staTemp = NULL;
  
  TH2F* trgNum = NULL;
  TH2F* muidNum = NULL;
  TH2F* staNum = NULL;
  
  for (int a=min; a<=max; a++)
    {
      trgNum = new TH2F (Form("trgNum%d",a), "", ny2D, ybins, npt2D, ptbins); trgNum->Sumw2();
      muidNum = new TH2F (Form("muidNum%d",a), "", ny2D, ybins, npt2D, ptbins); muidNum->Sumw2();
      staNum = new TH2F (Form("staNum%d",a), "", ny2D, ybins, npt2D, ptbins);staNum->Sumw2();

      TH2F *EffDen = (TH2F*) EffFile->Get("hdeno_2d");

      Long64_t nentries =fChain->GetEntries();
      //nentries = 2000000;
      
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<" in toy "<<a<<" ("<<min<<"-"<<max<<")"<<endl;
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m=GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt<3 || jpsi_pt>50) continue;
	      if (abs(jpsi_rap)>=2.4) continue;
	      if (!areGenMuonsInAcceptance2015(iQQ)) continue;

	      ibest = -1;
	      if (!isMatchedGenDiMuon(iQQ)) continue;//fill the num with the matching reco                                                                                                                  
	      if (ibest < 0) continue;
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(ibest);
	      TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(ibest);
	      TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(ibest);
	      if (!areMuonsInAcceptance2015(ibest)) continue;
	      if (!passQualityCuts2015(ibest)) continue;
	      if (!isTriggerMatch(ibest, triggerIndex_PP)) continue;
	      if (Reco_QQ_sign[ibest]!=0) continue;
	      if (RecoQQ4mom->Pt()<3 || RecoQQ4mom->Pt()>50) continue;
	      if (abs(RecoQQ4mom->Rapidity())>=2.4) continue;
	      if (RecoQQ4mom->M()<2.6 || RecoQQ4mom->M()>3.5) continue;

	      tnp_weight=1.0;

	      //trg 
	      tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
		tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	      trgNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);
	      
	      //muid
	      tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
		tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
		tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	      muidNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);
	      
	      //sta
	      tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
		tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),a+1) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),a+1) *
		tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
	      staNum->Fill(jpsi_rap, jpsi_pt, tnp_weight);
	    } //end of genQQ   
	}//end of events

      staNum->Multiply(AccNum);
      trgNum->Multiply(AccNum);
      muidNum->Multiply(AccNum);
      EffDen->Multiply(AccDen);

      trgTemp = new TEfficiency(Form("trgTemp%d",a), "trgToy AccxEff(y,pt); y; pt; AccxEff", ny2D, ybins, npt2D, ptbins);
      trgTemp->SetStatisticOption(TEfficiency::kBBayesian);
      trgTemp->SetPassedHistogram(*trgNum,"f");
      trgTemp->SetTotalHistogram(*EffDen,"f");
      trgTemp->SetName(Form("trgToy%d",a));
      trgEff->Add(trgTemp);
      
      muidTemp = new TEfficiency(Form("muidTemp%d",a), "muidToy AccxEff(y,pt); y; pt; AccxEff", ny2D, ybins, npt2D, ptbins);
      muidTemp->SetStatisticOption(TEfficiency::kBBayesian);
      muidTemp->SetPassedHistogram(*muidNum,"f");
      muidTemp->SetTotalHistogram(*EffDen,"f");
      muidTemp->SetName(Form("muidToy%d",a));	
      muidEff->Add(muidTemp);
      
      staTemp = new TEfficiency(Form("staTemp%d",a), "staToy AccxEff(y,pt); y; pt; AccxEff", ny2D, ybins, npt2D, ptbins);
      staTemp->SetStatisticOption(TEfficiency::kBBayesian);
      staTemp->SetPassedHistogram(*staNum,"f");
      staTemp->SetTotalHistogram(*EffDen,"f");
      staTemp->SetName(Form("staToy%d",a));
      staEff->Add(staTemp);
    }//endl of 100 toy loop
  gSystem->mkdir("FilesAccxEff/toyMC");
  TFile *fsave = new TFile (Form("FilesAccxEff/toyMC/%stoys%d%d.root", (isPr)?"pr":"npr", min, max),"RECREATE");
  trgEff->Write(Form("trg%dtoys", max-min+1), TObject::kSingleKey);
  muidEff->Write(Form("muid%dtoys", max-min+1), TObject::kSingleKey);
  staEff->Write(Form("sta%dtoys", max-min+1), TObject::kSingleKey);
  fsave->Close();
}

void oniaTree::TnpStat(string caseLabel = "") {
  double totBins [] = {0, 1.6, 2.4};

  int jtPtmin;
  int jtPtmax;
  bool isMid = false;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtmin=25; jtPtmax=35;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtmin=15; jtPtmax=25;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtmin=35; jtPtmax=45;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtmin=0; jtPtmax=1000;}

  if (caseLabel.find("016")!=std::string::npos) isMid = true;

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_PP_NoBkg_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR"),"READ");
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  if (isPr) trNom->SetBranchAddress("corr_ptw",&corr_ptw);
  TH1F* temp1 = NULL;
  TH1F* temp2 = NULL;

  TObjArray* trgArr = new TObjArray(); trgArr->SetOwner(kTRUE);
  TObjArray* muidArr = new TObjArray(); muidArr->SetOwner(kTRUE);
  TObjArray* staArr = new TObjArray(); staArr->SetOwner(kTRUE);
  

  TF1 *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
  double prEff = 1.0; double nprEff = 1.0; double bf = 1.0;
  double totEff=1.0;

  cout<<"[INFO] Creating arrays of histograms."<<endl;
  for (int i =0; i<100; i++)
    {
      if (jtPtmin>10) {
	if (isMid){
	  temp1 = new TH1F (Form("trg_%d",i), Form("trg %d",i), 7, 0.02, 1);
	  trgArr->Add(temp1);
	  temp1 = new TH1F (Form("muid_%d",i), Form("muid %d",i), 7, 0.02, 1);
	  muidArr->Add(temp1);
	  temp1 = new TH1F (Form("sta_%d",i), Form("sta %d",i), 7, 0.02, 1);
	  staArr->Add(temp1);
	}
	else {
          temp1 = new TH1F (Form("trg_%d",i), Form("trg %d",i), 5, 0, 1);
          trgArr->Add(temp1);
          temp1 = new TH1F (Form("muid_%d",i), Form("muid %d",i), 5, 0, 1);
          muidArr->Add(temp1);
          temp1 = new TH1F (Form("sta_%d",i), Form("sta %d",i), 5, 0, 1);
          staArr->Add(temp1);
	}
      }
      else {
	temp1 = new TH1F (Form("trg_%d",i), Form("trg %d",i), 2, totBins);
        trgArr->Add(temp1);
        temp1 = new TH1F (Form("muid_%d",i), Form("muid %d",i), 2, totBins);
        muidArr->Add(temp1);
        temp1 = new TH1F (Form("sta_%d",i), Form("sta %d",i), 2, totBins);
        staArr->Add(temp1);
      }
    }

  cout<<"[INFO] Loading the nominal corrections"<<endl;
  TFile* nomFile = TFile::Open("../Fitter/Input/correction_AccEff.root","READ"); 
  TEfficiency* prCorNom = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP_pr");
  TEfficiency* nprCorNom = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP_npr");
  TEfficiency* prTemp3 = NULL;
  TEfficiency* prTemp2 = NULL;
  TEfficiency* prTemp1 = NULL;
  TEfficiency* nprTemp3 = NULL;
  TEfficiency* nprTemp2 = NULL;
  TEfficiency* nprTemp1 = NULL;

  TH1F* nomHist = NULL;
  if (jtPtmin>10) 
    if (isMid) nomHist = new TH1F ("nomHist", "", 7, 0.02, 1);
    else nomHist = new TH1F ("nomHist", "", 5, 0, 1);
  else nomHist = new TH1F ("nomHist", "", 2, totBins);

  TFile* prtoyFile = NULL;
  TObjArray* prArr1 = NULL;
  TObjArray* prArr2 = NULL;
  TObjArray* prArr3 = NULL;
  TFile* nprtoyFile = NULL;
  TObjArray* nprArr1 = NULL;
  TObjArray* nprArr2 = NULL;
  TObjArray* nprArr3 = NULL;

  cout<<"[INFO] starting to process the tree."<<endl;

  int nentries = trNom->GetEntries();
  for (int i=0; i<10; i++)
    {
      prtoyFile = TFile::Open(Form("FilesAccxEff/toyMC/prtoys%d%d.root", i*10, i*10+9));
      nprtoyFile = TFile::Open(Form("FilesAccxEff/toyMC/nprtoys%d%d.root", i*10, i*10+9)); 
      prArr1 = (TObjArray*) prtoyFile->Get("trg10toys");
      prArr2 = (TObjArray*) prtoyFile->Get("muid10toys");
      prArr3 = (TObjArray*) prtoyFile->Get("sta10toys");
      nprArr1 = (TObjArray*) nprtoyFile->Get("trg10toys");
      nprArr2 = (TObjArray*) nprtoyFile->Get("muid10toys");
      nprArr3 = (TObjArray*) nprtoyFile->Get("sta10toys");
      for (int j=0; j<10; j++)
	{
	  prTemp1 = (TEfficiency*) prArr1->At(j);
	  prTemp2 = (TEfficiency*) prArr2->At(j);
	  prTemp3 = (TEfficiency*) prArr3->At(j);
	  nprTemp1 = (TEfficiency*) nprArr1->At(j);
	  nprTemp2 = (TEfficiency*) nprArr2->At(j);
	  nprTemp3 = (TEfficiency*) nprArr3->At(j);
	  
	  for (int jentry=0; jentry<nentries; jentry++) {
	    if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<" in toy "<<i*10+j<<endl;
	    trNom->GetEntry(jentry);
	    
	    if (abs(jp_rap)>2.4) continue;
	    if (jp_pt<3|| jp_pt>35) continue;
	    if (isMid && jtPtmin>10 && abs(jp_rap)>1.6) continue;
	    if (!isMid && jtPtmin>10 && abs(jp_rap)<=1.6) continue; 
	    if (jt_pt<jtPtmin || jt_pt>jtPtmax) continue;
	    
	    bf = bfrac->Eval(jp_pt);
	    if (i==0 && j==0) {
	      prEff = prCorNom->GetEfficiency(prCorNom->FindFixBin(jp_rap,jp_pt));
	      nprEff = nprCorNom->GetEfficiency(nprCorNom->FindFixBin(jp_rap,jp_pt));
	      totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
	      if (isPr) totEff = totEff*corr_ptw;
	      if (jtPtmin>10)
		nomHist->Fill(z, totEff);
	      else 
		nomHist->Fill(jp_rap, totEff);
	    }

	    temp1 = (TH1F*) trgArr->At(i*10+j);
	    prEff = prTemp1->GetEfficiency(prTemp1->FindFixBin(jp_rap,jp_pt));
	    nprEff = nprTemp1->GetEfficiency(nprTemp1->FindFixBin(jp_rap,jp_pt));
	    totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
	    if (isPr) totEff = totEff*corr_ptw;
	    if (jtPtmin>10)
	      temp1->Fill(z,totEff);
	    else 
	      temp1->Fill(jp_rap, totEff);
	    temp1 = (TH1F*) muidArr->At(i*10+j);
	    prEff = prTemp2->GetEfficiency(prTemp2->FindFixBin(jp_rap,jp_pt));
	    nprEff = nprTemp2->GetEfficiency(nprTemp2->FindFixBin(jp_rap,jp_pt));
	    totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
	    if (isPr) totEff = totEff*corr_ptw;
	    if (jtPtmin>10)
              temp1->Fill(z,totEff);
            else
              temp1->Fill(jp_rap, totEff);

	    temp1 = (TH1F*) staArr->At(i*10+j);
	    prEff = prTemp3->GetEfficiency(prTemp3->FindFixBin(jp_rap,jp_pt));
	    nprEff = nprTemp3->GetEfficiency(nprTemp3->FindFixBin(jp_rap,jp_pt));
	    totEff = 1.0/(bf*nprEff + (1-bf)*prEff);
	    if (isPr) totEff = totEff*corr_ptw;
	    if (jtPtmin>10)
              temp1->Fill(z,totEff);
            else
              temp1->Fill(jp_rap, totEff);
	  }// end of tree entries loop
	}// end of tobjarrays inside every file
      prtoyFile->Close();
      nprtoyFile->Close();
    }// end of files
  ////// filling the different syst files  
  ofstream file_trg(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_tnptrgStat.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt"));
  ofstream file_muid(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_tnpmuidStat.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt"));
  ofstream file_sta(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_tnpstaStat.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt"));

  file_trg << "tnp trg stat" << endl;
  file_muid << "tnp muid stat" << endl;
  file_sta << "tnp sta stat" << endl;

  vector<double> v_trg;
  vector<double> v_muid;
  vector<double> v_sta;

  v_trg.clear();
  v_muid.clear();
  v_sta.clear();

  int jmax;
  if (jtPtmin>10) 
    if (isMid) jmax =6;
    else jmax =5;
  else jmax =3;
  
  for (int j=1; j<jmax; j++){
    v_trg.clear();
    v_muid.clear();
    v_sta.clear();
    double binX = 0;
    if (jtPtmin>10)
      if (isMid)
	binX =nomHist->FindBin(j*0.14+0.2);
      else 
      binX =nomHist->FindBin(j*0.2+0.1);
    else
      binX = nomHist->FindBin(j);
    
    v_trg.push_back(nomHist->GetBinContent(binX));
    v_muid.push_back(nomHist->GetBinContent(binX));
    v_sta.push_back(nomHist->GetBinContent(binX));

    for (int i = 0; i<100; i++)
      {
	temp1 = (TH1F*) trgArr->At(i);
	v_trg.push_back(temp1->GetBinContent(binX));
	temp1 = (TH1F*) muidArr->At(i);
	v_muid.push_back(temp1->GetBinContent(binX));
	temp1 = (TH1F*) staArr->At(i);
	v_sta.push_back(temp1->GetBinContent(binX));
      }
    if (jtPtmin>10 && isMid) {
      file_trg<< "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< rms(v_trg,true) << endl;
      file_muid<< "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< rms(v_muid, true) << endl;
      file_sta<< "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< rms(v_sta, true) << endl;
    }
    else if (jtPtmin>10 && !isMid) {
      file_trg<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< rms(v_trg,true) << endl;
      file_muid<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< rms(v_muid, true) << endl;
      file_sta<< "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< rms(v_sta, true) << endl;
    }
    else {
      if (j==1){
      file_trg<< "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< rms(v_trg,true) << endl;
      file_muid<< "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< rms(v_muid, true) << endl;
      file_sta<< "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< rms(v_sta, true) << endl;
      }
      else {
	file_trg<< "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< rms(v_trg,true) << endl;
	file_muid<< "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< rms(v_muid, true) << endl;
	file_sta<< "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< rms(v_sta, true) << endl;
      }
    }
  }// end of z-range loop
  file_trg.close();
  file_muid.close();
  file_sta.close();
} //end of fct TnpStat


void oniaTree::FullAccEffSyst(string caseLabel) {
  int jtPtRange;
  bool isMid = false;

  if (caseLabel.find("midJtPt")!=std::string::npos) {jtPtRange = 0;}
  else if (caseLabel.find("lowJtPt")!=std::string::npos) {jtPtRange = -1;}
  else if (caseLabel.find("highJtPt")!=std::string::npos) {jtPtRange = 1;}
  else if (caseLabel.find("NoJets")!=std::string::npos) {jtPtRange = 100;}

  if (caseLabel.find("016")!=std::string::npos) isMid = true;

  string systName [] = {"tnpbinned", "tnptrgSyst", "tnpmuidSyst", "tnpstaSyst", "tnptrkSyst", "tnptrgStat", "tnpmuidStat", "tnpstaStat", "AccEffStat", "AccEffMisMod"};

  ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_fullAccEff.csv",caseLabel.c_str(),isPr?"prompt":"nonprompt"));

  fileOut<<"AccxEff"<<endl;
  double val1 [] = {0,0,0,0,0,0};
  int jmax;
  if (jtPtRange == 100) jmax = 2;
  else {
    if (isMid) jmax = 5;
    else jmax = 4;
  }

  if (jtPtRange != 100 && isMid)
    for (int j=1; j<6; j++){
      double syst=0;
      for (int i=0; i<(sizeof(systName)/sizeof(systName[0])); i++){
	double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", systName[i].c_str()), 0.16+j*0.14, 0.16+(j+1)*0.14, 0, 1.6);
	syst = syst + v1*v1;
      }
      fileOut << "0, 1.6, 6.5, 35, "<<0.16+j*0.14<<", "<< 0.16+(j+1)*0.14 <<", 0, 200, "<< sqrt(syst) << endl;
    }
  else if (jtPtRange != 100 && !isMid)
    for (int j=1; j<5; j++){
      double syst=0;
      for (int i=0; i<(sizeof(systName)/sizeof(systName[0])); i++){
        double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", systName[i].c_str()), j*0.2, (j+1)*0.2, 1.6, 2.4);
        syst = syst + v1*v1;
      }
      fileOut << "1.6, 2.4, 3, 35, "<<j*0.2<<", "<< (j+1)*0.2 <<", 0, 200, "<< sqrt(syst) << endl;
    }
  else {
    double syst=0;
    for (int i=0; i<(sizeof(systName)/sizeof(systName[0])); i++){
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", systName[i].c_str()), 0, 101, 0, 1.6);
      syst = syst + v1*v1;
    }
    fileOut << "0, 1.6, 6.5, 35, 0, 101, 0, 200, "<< sqrt(syst) << endl;
    syst=0;
    for (int i=0; i<(sizeof(systName)/sizeof(systName[0])); i++){
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_PP_%s.csv", caseLabel.c_str(), isPr?"prompt":"nonprompt", systName[i].c_str()), 0, 101, 1.6, 2.4);
      syst = syst + v1*v1;
    }
    fileOut << "1.6, 2.4, 3, 35, 0, 101, 0, 200, "<< sqrt(syst) << endl;
  }
  fileOut.close();
}

void oniaTree::AccEffSyst_all() {

  //string systName [] = {"midJtPt_016","midJtPt_1624", "lowJtPt_016","lowJtPt_1624", "highJtPt_016","highJtPt_1624", "NoJets_total"};
  string systName [] = {"midJtPt_016", "lowJtPt_016", "highJtPt_016"};
  int systN = sizeof(systName)/sizeof(systName[0]);

  for (int i=0; i<systN; i++) {
    cout<<"--------[INFO] Sarting TnpSyst for "<<systName[i]<<"----------"<<endl;
    //TnpSyst(systName[i]);
    cout<<"--------[INFO] Sarting AccEffStat for "<<systName[i]<<"----------"<<endl;
    AccEffStat(systName[i]);
    cout<<"--------[INFO] Sarting TnpStat for "<<systName[i]<<"----------"<<endl;
    //TnpStat(systName[i]);
    cout<<"--------[INFO] Sarting AccEffMisMod for "<<systName[i]<<"----------"<<endl;
    //AccEffMisMod(systName[i]);
    cout<<"--------[INFO] Sarting FullAccEffSyst for "<<systName[i]<<"----------"<<endl;
    FullAccEffSyst(systName[i]);
    cout<<"--------[INFO] All systematics done for "<<systName[i]<<"----------"<<endl;
  }
}
