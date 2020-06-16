#include "makeAccEff.h" 
//macro to get all AccEff Syst
float min_z = 0.064;
float max_z = 1.0;

float min_jetpt = 10.;
float max_jetpt = 60.;

int nBinZ = 6;

float z_reco_binWidth = (max_z-min_z)*1.0/nBinZ;

Double_t ptbins_01btw6590 [] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.25, 9.5, 9.75, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 19., 20.0, 25.0, 30.0, 50, 100.0};
Double_t ptbins_02btw6590 [] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 50, 100.0};
Double_t ptbins_025btw6580 []= {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 100.0};
Double_t ptbins_005t75_01t90_25t100_100t200 [] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.55, 6.6, 6.65, 6.7, 6.75, 6.8, 6.85, 6.9, 6.95, 7.0, 7.05, 7.1, 7.15, 7.2, 7.25, 7.3, 7.35, 7.4, 7.45, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.25, 9.5, 9.75, 10.0, 10, 11, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 25.0, 30.0, 50, 100.0};

Double_t ybins_24EvenBins []= {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
Double_t ybins_absEta_12EvenBins []= {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

Double_t ybins_1bin1010_15Bins []= {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
Double_t ybins_absEta_1bin010_8Bins []= {0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

Double_t ybins_48EvenBins []= {-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};

int centbins_10t40_20t100 [] = {0, 10, 20, 30, 40, 60, 80, 100, 140, 200};
int centbins_5t20_10t40_20t100 [] = {0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 200}; 

Double_t zbins_ref [] = {0.064, 0.22, 0.376, 0.532, 0.688, 0.844, 1.};
Double_t jtptbins_ref [] = {6.5, 10, 20, 30, 40, 50, 60};
string jtptTags [] = {"lowestJtPt","lowerJtPt","lowJtPt","midJtPt","highJtPt","higherJtPt"};

void oniaTree::setBins(string caseTag) {
  ////////for pt
  if (caseTag.find("pt_FinerThenCoarser")!=std::string::npos) {
    nptbins = sizeof(ptbins_005t75_01t90_25t100_100t200)/sizeof(double);
    for (int i=0; i<nptbins; i++)
      ptbins[i] = ptbins_005t75_01t90_25t100_100t200[i];
  }
  else if (caseTag.find("pt_SizeDoubled")!=std::string::npos){
    nptbins = sizeof(ptbins_02btw6590)/sizeof(double);
    for (int i=0; i<nptbins; i++)
      ptbins[i] = ptbins_02btw6590[i];
  }
  else {
    nptbins = sizeof(ptbins_01btw6590)/sizeof(double);
    for (int i=0; i<nptbins; i++)
      ptbins[i] = ptbins_01btw6590[i];
  }

  nptbins = nptbins-1;
  /////for rapidity
  if (caseTag.find("absEta_12EvenBins")!=std::string::npos){
    nybins = sizeof(ybins_absEta_12EvenBins)/sizeof(double);
    for (int i=0; i<nybins; i++)
      ybins[i] = ybins_absEta_12EvenBins[i];
  }
  else if (caseTag.find("absEta_1bin010_8Bins")!=std::string::npos) {
    nybins = sizeof(ybins_absEta_1bin010_8Bins)/sizeof(double);
    for (int i=0; i<nybins; i++)
      ybins[i] = ybins_absEta_1bin010_8Bins[i];
  }
  else if (caseTag.find("rap_1bin1010_15Bins")!=std::string::npos) {
    nybins = sizeof(ybins_1bin1010_15Bins)/sizeof(double);
    for (int i=0; i<nybins; i++)
      ybins[i] = ybins_1bin1010_15Bins[i];
  } 
  else if (caseTag.find("rap_48EvenBins")!=std::string::npos) {
    //cout <<"using rap_48EvenBins"<<endl;
    nybins = sizeof(ybins_48EvenBins)/sizeof(double);
    for (int i=0; i<nybins; i++)
      ybins[i] = ybins_48EvenBins[i];
  }
  else {
    nybins = sizeof(ybins_24EvenBins)/sizeof(double);
    for (int i=0; i<nybins; i++)
      ybins[i] = ybins_24EvenBins[i];
  }

  nybins = nybins-1;
  /////for centrality
  if (caseTag.find("cent_5t20")!=std::string::npos){
    ncentbins = sizeof(centbins_5t20_10t40_20t100)/sizeof(int);
    for (int i=0; i<ncentbins; i++)
      centbins[i] = centbins_5t20_10t40_20t100[i];
  }
  else {
    ncentbins = sizeof(centbins_10t40_20t100)/sizeof(int);
    for (int i=0; i<ncentbins; i++)
      centbins[i] = centbins_10t40_20t100[i];
  }
  ncentbins = ncentbins-1;
  cout <<"ncentbins = "<<ncentbins<<endl;
  ////for z
  nzbins = sizeof(zbins_ref)/sizeof(double);
  for (int i=0; i<nzbins; i++)
    zbins[i] = zbins_ref[i];
  nzbins = nzbins-1;
  ////for jet pt 
  njtptbins = sizeof(jtptbins_ref)/sizeof(double);
  for (int i=0; i<njtptbins; i++)
    jtptbins[i] = jtptbins_ref[i];
  njtptbins = njtptbins-1;
}

void oniaTree::AccEffMisMod(string caseTag) {
  setBins(caseTag);
  TFile* fNom = TFile::Open(Form("FilesAccxEff%s/ptHistsWith100%sToys_%s_%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"));
  TObjArray* zhisArrNom = (TObjArray*) fNom->Get("zDistArr");
  TH2D* zHistNom = (TH2D*) zhisArrNom->At(0);

  TObjArray* zhisArrNomCent = NULL;
  TH2D* zHistNomCent = NULL;
  TObjArray* zhisArrNomPeri = NULL;
  TH2D* zHistNomPeri = NULL;
    if (isPbPb) {
      zhisArrNomCent = (TObjArray*) fNom->Get("zDistArrCent");
      zHistNomCent = (TH2D*) zhisArrNomCent->At(0);
      zhisArrNomPeri = (TObjArray*) fNom->Get("zDistArrPeri");
      zHistNomPeri = (TH2D*) zhisArrNomPeri->At(0);
    }

  caseTag.replace(caseTag.find("_NoWeights"),10,"");

  TFile* fSyst = TFile::Open(Form("FilesAccxEff%s/ptHistsWith100%sToys_%s_%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"));
  TObjArray* zhisArrSyst = (TObjArray*) fSyst->Get("zDistArr");
  TH2D* zHistSyst = (TH2D*) zhisArrSyst->At(0);

  TObjArray* zhisArrSystCent = NULL;
  TH2D* zHistSystCent = NULL;
  TObjArray* zhisArrSystPeri = NULL;
  TH2D* zHistSystPeri = NULL;
  if (isPbPb) {
    zhisArrSystCent = (TObjArray*) fSyst->Get("zDistArrCent");
    zHistSystCent = (TH2D*) zhisArrSystCent->At(0);
    zhisArrSystPeri = (TObjArray*) fSyst->Get("zDistArrPeri"); 
    zHistSystPeri = (TH2D*) zhisArrSystPeri->At(0);
  }

  for (int k=0; k<njtptbins; k++) {
    ofstream zfileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_EffMisMod.csv", jtptTags[k].c_str(),isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP"));
    zfileOut<<"Eff mismod"<<endl;

    for (int j=0; j<nzbins; j++) {
      double nNom = zHistNom->GetBinContent(zHistNom->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));
      double nSyst = zHistSyst->GetBinContent(zHistSyst->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));


      //cout <<"for Jetpt = "<<jtptTags[k]<<", z = "<<zbins[j]<<", nNom = "<<nNom<<", nSyst = "<<nSyst<<", err = "<<fabs(nNom - nSyst)*1./nNom<<endl;      
      if (isPbPb) {
	zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 180, "<< fabs(nNom - nSyst)*1./nNom << endl;
	nNom = zHistNomCent->GetBinContent(zHistNomCent->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));
	nSyst = zHistSystCent->GetBinContent(zHistSystCent->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));
	zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, "<<centCut<<", "<< fabs(nNom - nSyst)*1./nNom << endl;
	nNom = zHistNomPeri->GetBinContent(zHistNomPeri->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));
	nSyst = zHistSystPeri->GetBinContent(zHistSystPeri->FindBin(zbins[j]+0.001,jtptbins[k]+0.001));
	zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", "<<centCut<<", 180, "<< fabs(nNom - nSyst)*1./nNom << endl;
      }
      else
	zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 200, "<<  fabs(nNom - nSyst)*1./nNom<< endl;
    }
    zfileOut.close();
  }
  
} // end of the function


void oniaTree::AccEffStatToy_Acc(int nToys, string caseTag) {
  // Randomise TEfficiency 100 times. Output will be a TObjArray with 100 TH2
  setBins(caseTag);
  TFile* effFile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str())); 
  TEfficiency *eff = (TEfficiency*) effFile->Get(Form("hcorr_Jpsi_%s_%s_Acc",isPbPb?"PbPb":"PP",isPr?"pr":"npr"));
  
  TObjArray *arrEffs = new TObjArray(); // Array to store efficiencies
  arrEffs->SetOwner(true);

  TRandom* rnd = new TRandom3(); // For randomisation

  for (int i=0 ; i < nToys ; i++)
    {
      //if (i%10==0) 
      cout <<"[INFO] toy "<<i<<"/"<<nToys<<endl;
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
	      double ntot = histoTot->GetBinContent(bin);
	      double newPass = rnd->Binomial(ntot,effVal);
	      histoPass->SetBinContent(histoPass->GetBin(j,k),newPass);
	    }
	}
      TEfficiency* newHeff = (TEfficiency*) eff->Clone(Form("effToy_%d",i)); // Histo to store new 2D eff
      newHeff->SetPassedHistogram(*histoPass, "f");
      newHeff->SetTotalHistogram(*histoTot, "f");
      arrEffs->Add(newHeff);
    }
  effFile->Close();
  gSystem->mkdir(Form("FilesAccxEff%s/toyMC",caseTag.c_str()));
  TFile* fsave = new TFile(Form("FilesAccxEff%s/toyMC/%sAcc%dToys_%s.root",caseTag.c_str(),isPr?"pr":"npr", nToys,isPbPb?"PbPb":"PP"),"RECREATE");
  arrEffs->Write("accEffArray", TObject::kSingleKey);
  fsave->Close();
}

void oniaTree::AccEffStatToy_Eff(int nToys,string caseTag) {
  // Randomise TEfficiency 100 times. Output will be a TObjArray with 100 TH2
  TFile* effFile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str())); 
  TEfficiency *eff = (TEfficiency*) effFile->Get(Form("hcorr_Jpsi_%s_%s_Eff",isPbPb?"PbPb":"PP",isPr?"pr":"npr"));
  
  TObjArray *arrEffs = new TObjArray(); // Array to store efficiencies
  arrEffs->SetOwner(true);

  TRandom* rnd = new TRandom3(); // For randomisation

  for (int i=0 ; i < nToys ; i++)
    {
      //if (i%10==0) 
      cout <<"[INFO] toy "<<i<<"/"<<nToys<<endl;
      TH2* histoTot = (TH2*) eff->GetTotalHistogram()->Clone(Form("hTotToy_%d",i)); // Get corresponding histo and number of bins
      histoTot->Sumw2();
      TH2* histoPass = (TH2*) histoTot->Clone(Form("hPassToy_%d",i)); // Get corresponding histo and number of bins
      histoPass->Sumw2();
      int nBinsX = histoTot->GetNbinsX();
      int nBinsY = histoTot->GetNbinsY();
      
      for (int j = 1 ; j <= nBinsX ; j++)
	{
	  //cout <<"start bin "<<j<<endl;
	  for (int k = 1 ; k <= nBinsY ; k++)
	    {
	      int bin = histoTot->GetBin(j,k);      
	      double effVal = eff->GetEfficiency(bin);
	      double ntot = histoTot->GetBinContent(bin);
	      double newPass = rnd->Binomial(ntot,effVal);
	      histoPass->SetBinContent(histoPass->GetBin(j,k),newPass);
	    }
	}
      TEfficiency* newHeff = (TEfficiency*) eff->Clone(Form("effToy_%d",i)); // Histo to store new 2D eff
      newHeff->SetPassedHistogram(*histoPass, "f");
      newHeff->SetTotalHistogram(*histoTot, "f");
      arrEffs->Add(newHeff);
    }

  gSystem->mkdir(Form("FilesAccxEff%s/toyMC",caseTag.c_str()));
  TFile* fsave = new TFile(Form("FilesAccxEff%s/toyMC/%sEff%dToys_%s.root",caseTag.c_str(),isPr?"pr":"npr", nToys,isPbPb?"PbPb":"PP"),"RECREATE");
  eff->Write("nominal");
  arrEffs->Write("accEffArray", TObject::kSingleKey);
  
  
  if (caseTag.find("centBins")!=std::string::npos && isPbPb) {
    setBins(caseTag);
    TObjArray *arrEffsCent = new TObjArray();
    arrEffsCent->SetOwner(true);
    
    for (int iCent = 0; iCent<ncentbins; iCent++) {
      cout <<"[INFO] centrality[i] = "<<iCent<<" = cent_"<<centbins[iCent]<<centbins[iCent+1]<<endl;
      TEfficiency* effTemp = (TEfficiency*) effFile->Get(Form("hcorr_Jpsi_%s_%s_Eff_cent%d%d",isPbPb?"PbPb":"PP",isPr?"pr":"npr",centbins[iCent],centbins[iCent+1]));
      cout <<"Reading from efficiency "<<effTemp->GetName()<<endl;
      TH2* histoTot = (TH2*) effTemp->GetTotalHistogram()->Clone(Form("hTotToy_cent%d%d",centbins[iCent],centbins[iCent+1])); // Get corresponding histo and number of bins
      for (int i=0 ; i < nToys ; i++)
	{
	  cout <<"[INFO] toy "<<i<<"/"<<nToys<<endl;
	  histoTot->Sumw2();
	  TH2* histoPass = (TH2*) histoTot->Clone(Form("hPassToy_cent%d%d_%d",centbins[iCent],centbins[iCent+1],i)); // Get corresponding histo and number of bins
	  histoPass->Sumw2();
	  int nBinsX = histoTot->GetNbinsX();
	  int nBinsY = histoTot->GetNbinsY();
	  for (int j = 1 ; j <= nBinsX ; j++)
	    {
	      for (int k = 1 ; k <= nBinsY ; k++)
		{
		  int bin = histoTot->GetBin(j,k);      
		  double effVal = effTemp->GetEfficiency(bin);
		  double ntot = histoTot->GetBinContent(bin);
		  double newPass = rnd->Binomial(ntot,effVal);
		  histoPass->SetBinContent(histoPass->GetBin(j,k),newPass);
		  //cout <<"bin = ["<<j<<","<<k<<"]"<<", effVal = "<<effVal<<", ntot = "<<ntot<<", oldEff = "<<effVal<<", newEff = "<<newPass/ntot<<endl;
		}
	    }
	  TEfficiency* newHeff = (TEfficiency*) effTemp->Clone(Form("effToy_cent%d%d_%d",centbins[iCent],centbins[iCent+1],i)); // Histo to store new 2D eff
	  newHeff->SetPassedHistogram(*histoPass, "f");
	  newHeff->SetTotalHistogram(*histoTot, "f");
	  arrEffsCent->Add(newHeff);
	}
      effTemp->Write(Form("nominal_cent%d%d",centbins[iCent],centbins[iCent+1]));
    }
    arrEffsCent->Write(Form("accEffArray_centAll"), TObject::kSingleKey);
  }
  
  //effFile->Close();
  fsave->Close();
}

void oniaTree::AccEffStat(string caseTag) {
  setBins(caseTag);
  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s_NoBkg_jetR3_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR",isPbPb?"PbPb":"PP"),"READ");
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw; int centr;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  trNom->SetBranchAddress("corr_ptw",&corr_ptw);
  trNom->SetBranchAddress("centr",&centr);

  TEfficiency* tempPrEff = NULL;
  TEfficiency* tempPrAcc = NULL;
  TEfficiency* tempNprEff = NULL;
  TObjArray* tempCentArr = NULL;

  TObjArray* zhisArr = new TObjArray(); //hisArr1624->SetOwner(kTRUE);
  TObjArray* zhisArrCent = new TObjArray();
  TObjArray* zhisArrPeri = new TObjArray();

  TFile *prcorrFile = TFile::Open(Form("FilesAccxEff%s/toyMC/pr%s100Toys_%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?(isAcc?"PP":"PbPb"):"PP"));
  TObjArray* prcorArr = (TObjArray*) prcorrFile->Get("accEffArray");

  //TH1F* histVar = NULL;
  TH2D* zhist2D = NULL;//new TH2D("zhist2D",";z;p_{T,jet}",6,0.064,1,5,10,60);
  TH2D* zhist2DCent = NULL;
  TH2D* zhist2DPeri = NULL;

  TFile* nomFile = TFile::Open(Form("../Fitter/Input/correction_AccEff_centMaps%s.root",caseTag.c_str()),"READ");
  TEfficiency* prnomEff = (TEfficiency*) nomFile->Get(Form("hcorr_Jpsi_%s_pr_Eff",isPbPb?"PbPb":"PP"));
  TEfficiency* prnomAcc = (TEfficiency*) nomFile->Get("hcorr_Jpsi_PP_pr_Acc");

  TList* lcorr = nomFile->GetListOfKeys();
  TIter nextCorr(lcorr);

  TObjArray* corrCentNom = new TObjArray();
  corrCentNom->SetOwner(kTRUE);
  if (isPbPb) {
    TObjString* fname(0x0);
    while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
      {
	TEfficiency* h = static_cast<TEfficiency*>(nomFile->FindObjectAny(fname->GetString().Data()));
	
	TString sName(h->GetName());
	if ( sName.Contains("cent") ){
	  corrCentNom->Add(h->Clone());
	  cout <<"adding "<<sName<<" to the corrction array"<<endl;
	}
      }
  }
  
  TObjArray* corrCentArr = new TObjArray();
  corrCentArr->SetOwner(kTRUE);
  corrCentArr = (TObjArray*) (prcorrFile->Get("accEffArray_centAll"));

  double prEff = 1.0;
  double totEff = 1.0;

  int nentries = trNom->GetEntries();
  int nToy =100;
  if ( ! (caseTag.find("_NoWeights")!=std::string::npos) ) nToy =1;
  cout <<"nToy = "<<nToy<<endl;
  for (int i=0; i<=nToy; i++) ///////100 
    {
      zhist2D = new TH2D(Form("zhist2D_%d",i),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);
      zhist2DCent = new TH2D(Form("zhist2DCent_%d",i),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);
      zhist2DPeri = new TH2D(Form("zhist2DPeri_%d",i),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);

      if (i==0) cout<<"[INFO] Applying the nominal AccxEff"<<endl;
      else cout<<"[INFO] Applying toy "<<i<<"/100"<<endl;
      if (i==0) {
	tempPrEff=prnomEff;
	tempPrAcc=prnomAcc; 
      }
      else {
	if (isAcc){
	  tempPrAcc=(TEfficiency*) prcorArr->At(i-1);
	  tempPrEff=prnomEff;
	}
	else {
	  tempPrAcc=prnomAcc;
	  tempPrEff=(TEfficiency*) prcorArr->At(i-1);
	}
      }
      //cout <<"got the correction"<<endl;
      for (int jentry=0; jentry<nentries; jentry++) {
	trNom->GetEntry(jentry);
	if (abs(jp_rap)>2.4) continue;
	if (jp_pt<6.5 ||jp_pt>100) continue;
	if (!isPbPb)
	  prEff = (tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_rap,jp_pt)))*(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt)));
	else {
	  for (int iCent=0; iCent<ncentbins; iCent++) {
	    if (centr<=centbins[iCent+1]) {
	      if (i==0 || isAcc)
		tempPrEff = (TEfficiency*) corrCentNom->At(iCent);
	      else 
		tempPrEff = (TEfficiency*) corrCentArr->FindObject(Form("effToy_cent%d%d_%d",centbins[iCent],centbins[iCent+1],(i-1)));
	      break;
	    }
	  }
	  prEff = (tempPrAcc->GetEfficiency(tempPrAcc->FindFixBin(jp_rap,jp_pt)))*(tempPrEff->GetEfficiency(tempPrEff->FindFixBin(jp_rap,jp_pt)));
	}
	if (prEff>0)
	  totEff=1.0/prEff;
	else totEff=1.0;
	totEff = totEff*corr_ptw;	
	if (z>=0 && z<=1) {
	  zhist2D->Fill(z, jt_pt,totEff);
	  if (isPbPb) {
	    if (centr<centCut)
	      zhist2DCent->Fill(z, jt_pt,totEff);
	    else 
	      zhist2DPeri->Fill(z, jt_pt,totEff);
	  }
	}
      }// end of the tree entries
      //cout <<"done with applying the corrections"<<endl;
      zhisArr->Add(zhist2D);
      if (isPbPb) {
	zhisArrCent->Add(zhist2DCent);
	zhisArrPeri->Add(zhist2DPeri);
      }
    } // end of the variations
  TFile* fsave = new TFile(Form("FilesAccxEff%s/ptHistsWith100%sToys_%s_%s.root",caseTag.c_str(),isAcc?"Acc":"Eff",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"),"RECREATE");
  zhisArr->Write("zDistArr",TObject::kSingleKey);
  if (isPbPb) {
    zhisArrCent->Write("zDistArrCent",TObject::kSingleKey);
    zhisArrPeri->Write("zDistArrPeri",TObject::kSingleKey);
  }
  fsave->Close();
  
  if ( ! (caseTag.find("_NoWeights")!=std::string::npos) ) return;  

  cout<<"[INFO] making csv files"<<endl;

  for (int k=0; k<njtptbins; k++) {
    cout <<"file name" <<Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%sStat.csv", jtptTags[k].c_str(),isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", isAcc?"Acc":"Eff")<<endl;
    ofstream zfileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%sStat.csv", jtptTags[k].c_str(),isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", isAcc?"Acc":"Eff"));
    if (isAcc) zfileOut<<"Acc stat"<<endl;
    else zfileOut<<"Eff stat"<<endl;

    vector<double> vCount;
    vCount.clear();

    vector<double> vCountCent;
    vCountCent.clear();

    vector<double> vCountPeri;
    vCountPeri.clear();

    TH2D* ztemp = NULL;
    
    for (int j=0; j<nzbins; j++)
      {
	vCount.clear();
	vCountCent.clear();
	vCountPeri.clear();
	for (int i=0; i<=100; i++) { // 0 to 1 
	  ztemp = (TH2D*) zhisArr->At(i);
	  vCount.push_back(ztemp->GetBinContent(ztemp->FindBin(zbins[j]+0.001,jtptbins[k]+0.001)));
	  if (isPbPb) {
	    ztemp = (TH2D*) zhisArrCent->At(i);
	    vCountCent.push_back(ztemp->GetBinContent(ztemp->FindBin(zbins[j]+0.001,jtptbins[k]+0.001)));
	    ztemp = (TH2D*) zhisArrPeri->At(i);
	    vCountPeri.push_back(ztemp->GetBinContent(ztemp->FindBin(zbins[j]+0.001,jtptbins[k]+0.001)));
	  }
	}
	if (isPbPb) {
	  zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 180, "<< rms(vCount,true) << endl;
	  zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, "<<centCut<<", "<< rms(vCountCent,true) << endl;
	  zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", "<<centCut<<", 180, "<< rms(vCountPeri,true) << endl;
	}
	else
	  zfileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 200, "<< rms(vCount,true) << endl;
      }
    zfileOut.close();
  }
  prcorrFile->Close();
  nomFile->Close();
  trFile->Close();
}


void oniaTree::TnpSyst(string caseTag) {
  setBins(caseTag);

  cout<<"[INFO] Loading the tree"<<endl;
  TFile* trFile = TFile::Open(Form("../Fitter/TreesForUnfolding/tree_%s_%s_NoBkg_jetR3_AccEff_JEC.root",isPr?"MCJPSIPR":"MCJPSINOPR",isPbPb?"PbPb":"PP"),"READ"); 
  TTree* trNom = (TTree*) trFile->Get("treeForUnfolding");
  float z; float jp_pt; float jp_rap; float jp_mass; float jt_pt; float jt_rap; float corr_ptw; int centr;
  trNom->SetBranchAddress("z",&z);
  trNom->SetBranchAddress("jp_pt",&jp_pt);
  trNom->SetBranchAddress("jp_rap",&jp_rap);
  trNom->SetBranchAddress("jp_mass",&jp_mass);
  trNom->SetBranchAddress("jt_pt",&jt_pt);
  trNom->SetBranchAddress("jt_rap",&jt_rap);
  trNom->SetBranchAddress("corr_ptw",&corr_ptw);
  trNom->SetBranchAddress("centr",&centr);

  cout<<"[INFO] Loading the corrections"<<endl;
  //TFile* prAccFile = TFile::Open(Form("FilesAccxEff/Acc/prAccHists_%s.root",isPbPb?"PbPb":"PP"),"READ"); 
  TFile* prAccFile = TFile::Open(Form("FilesAccxEff%s/Acc/prAccHists_%s.root",caseTag.c_str(),"PP"),"READ"); //use pp acceptance for pp and PbPb
  TFile* prEffFile = TFile::Open(Form("FilesAccxEff%s/Eff/prEffHists_%s.root",caseTag.c_str(),isPbPb?"PbPb":"PP"),"READ");

  double prEff = 1.0;
  double prAcc = 1.0;
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
    "trk_minus1sig_stat",
    "tag_syst",
    "tag_syst"//temporary
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
    "trkStat",
    "tagSyst"
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
  
  
  TObjArray *countArr = new TObjArray(15); countArr->SetOwner(true); TH2D* countHist = NULL;
  TObjArray *countArrCent = new TObjArray(15); countArrCent->SetOwner(true); TH2D* countHistCent = NULL;
  TObjArray *countArrPeri = new TObjArray(15); countArrPeri->SetOwner(true); TH2D* countHistPeri = NULL;
  TEfficiency *prEffTemp = NULL; TEfficiency *nprEffTemp = NULL; 
  TEfficiency *prAccTemp = NULL; TEfficiency *nprAccTemp = NULL;

  vector<double> vCount;
  vector<double> vCountCent;
  vector<double> vCountPeri;

  TH2F *prAccNum = (TH2F*) prAccFile->Get("hnum_2d_nominal");
  TH2F *prAccDen = (TH2F*) prAccFile->Get("hdeno_2d");
  TH2F *prEffNum = NULL;//(TH2F*) prEffFile->Get("hnum_nominal");
  TH2F *prEffDen = (TH2F*) prEffFile->Get("hdeno_pty");
  //prEffDen->Multiply(prAccDen);

  if (!prAccNum||!prAccDen||!prEffNum||!prEffDen) cout<<"histogram not found!";
  cout<<"[INFO] all the hists are loaded"<<endl;

  
  for(int i=0; i<nCorr; i++) {
    prEffNum = (TH2F*) prEffFile->Get(Form("hnum_%s", corrName[i].c_str()));
    //prEffNum->Multiply(prAccNum);

    cout<<"syst i = "<<i<<endl;
    prEffTemp = new TEfficiency(Form("prEffCorr_%s", corrName[i].c_str()), "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
    prEffTemp->SetStatisticOption(TEfficiency::kBBayesian);
    prEffTemp->SetPassedHistogram(*prEffNum,"f");
    prEffTemp->SetTotalHistogram(*prEffDen,"f");

    prAccTemp = new TEfficiency(Form("prAccCorr_%s", corrName[i].c_str()), "AccxEff(y,pt); y; pt; eff", nybins, ybins, nptbins, ptbins);
    prAccTemp->SetStatisticOption(TEfficiency::kBBayesian);
    prAccTemp->SetPassedHistogram(*prAccNum,"f");
    prAccTemp->SetTotalHistogram(*prAccDen,"f");

    //countHist = new TH1F (Form("his_%s",corrName[i].c_str()), Form("tnp %s",corrName[i].c_str()), nBinZ, min_z, max_z);
    countHist = new TH2D(Form("zhist2D_%s",corrName[i].c_str()),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);
    countHistCent = new TH2D(Form("zhist2DCent_%s",corrName[i].c_str()),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);
    countHistPeri = new TH2D(Form("zhist2DPeri_%s",corrName[i].c_str()),";z;p_{T,jet}",nzbins,zbins,njtptbins,jtptbins);
    //filling the histograms
    cout<<"[INFO] Filling the histograms"<<endl;
    int nentries = trNom->GetEntries();
    for (int jentry=0; jentry<nentries; jentry++) {
      if (jentry%1000000==0) cout<<"[INFO] Processing entry "<<jentry<<"/"<<nentries<<endl;
      trNom->GetEntry(jentry);
      
      if (abs(jp_rap)>2.4) continue;
      if (jp_pt<6.5|| jp_pt>100) continue;

      prEff = prEffTemp->GetEfficiency(prEffTemp->FindFixBin(jp_rap,jp_pt));
      prAcc = prAccTemp->GetEfficiency(prAccTemp->FindFixBin(jp_rap,jp_pt));
      totEff = 1.0/(prEff*prAcc);
      totEff = totEff*corr_ptw;
      countHist->Fill(z, jt_pt, totEff);
      if (isPbPb) {
	if (centr<centCut) countHistCent->Fill(z, jt_pt, totEff);
	else countHistPeri->Fill(z, jt_pt, totEff);
      }
    }
    countArr->Add(countHist);
    if (isPbPb) {
      countArrCent->Add(countHistCent);
      countArrPeri->Add(countHistPeri);
    }
  } 
  
  TFile *fsave = new TFile(Form("FilesAccxEff/tnpSystHists_%s_%s.root",isPbPb?"PbPb":"PP",isPr?"prompt":"nonprompt"),"RECREATE");
  countArr->Write();
  countArrCent->Write();
  countArrPeri->Write();
  fsave->Close();
  cout<<"[INFO] Getting the ratios and filling the csv files"<<endl;
  
  ////// filling the different syst files
  for (int k=0; k<njtptbins; k++) {
    for (int i = 0; i<nSyst ; i++) { 
      ofstream fileOut(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt",isPbPb?"PbPb":"PP",systName[i].c_str()));
      fileOut << "tnp " << systName[i] << endl;
      vCount.clear();
      vCountCent.clear();
      vCountPeri.clear();

      for (int j=0; j<nzbins; j++){
	vCount.clear();
	vCountCent.clear();
	vCountPeri.clear();

	countHist = (TH2D*) countArr->At(0);
	double binX = countHist->FindBin(zbins[j]+0.001,jtptbins[k]+0.001);
	vCount.push_back((double)(countHist->GetBinContent(binX)));
	
	countHist = (TH2D*) countArr->At(2*i+1);
	vCount.push_back((double)(countHist->GetBinContent(binX)));
	
	countHist = (TH2D*) countArr->At(2*i+2);
	vCount.push_back((double)(countHist->GetBinContent(binX)));

	if (isPbPb) {
	  countHistCent = (TH2D*) countArrCent->At(0);
	  vCountCent.push_back((double)(countHistCent->GetBinContent(binX)));
	  countHistCent = (TH2D*) countArrCent->At(2*i+1);
	  vCountCent.push_back((double)(countHistCent->GetBinContent(binX)));	
	  countHistCent = (TH2D*) countArrCent->At(2*i+2);
	  vCountCent.push_back((double)(countHistCent->GetBinContent(binX)));

	  countHistPeri = (TH2D*) countArrPeri->At(0);
	  vCountPeri.push_back((double)(countHistPeri->GetBinContent(binX)));
	  countHistPeri = (TH2D*) countArrPeri->At(2*i+1);
	  vCountPeri.push_back((double)(countHistPeri->GetBinContent(binX)));	
	  countHistPeri = (TH2D*) countArrPeri->At(2*i+2);
	  vCountPeri.push_back((double)(countHistPeri->GetBinContent(binX)));
	}

	
	if (isPbPb) {
	  fileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 180, "<< maxdiff(vCount, 1) << endl;
	  fileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, "<<centCut<<", "<< maxdiff(vCountCent, 1) << endl;
	  fileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", "<<centCut<<", 180, "<< maxdiff(vCountPeri, 1) << endl;
	}
	else
	  fileOut<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 200, "<< maxdiff(vCount, 1) << endl;
      }
      fileOut.close();
    }
  

    int centMin = 0;
    int centMax=200;
    if (isPbPb) centMax=180;
    ofstream fileAll(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_fullAccEff.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt",isPbPb?"PbPb":"PP"));
    fileAll<<"AccxEff"<<endl;
    double val1 [] = {0,0,0,0,0,0,0,0};
    for (int j=0; j<nzbins; j++){
      double v1 = 0;
      double v1Cent = 0;
      double v1Peri = 0;
      double syst = 0;
      double systCent = 0;
      double systPeri = 0;
      v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffStat"), zbins[j], zbins[j+1], 0, 2.4,centMin, centMax);
      if (isPbPb) {
	v1Cent = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffStat"), zbins[j], zbins[j+1], 0, 2.4,centMin, centCut);
	v1Peri = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffStat"), zbins[j], zbins[j+1], 0, 2.4,centCut, centMax);
      }
      cout <<"for EffStat v1 = "<<v1<<endl;
      syst=v1*v1;
      systCent = v1Cent*v1Cent;
      systPeri=v1Peri*v1Peri;

      v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "AccStat"), zbins[j], zbins[j+1], 0, 2.4, centMin, centMax);
      syst = syst + v1*v1;
      if (isPbPb) {
	v1Cent = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "AccStat"), zbins[j], zbins[j+1], 0, 2.4, centMin, centCut);
	systCent = systCent + v1Cent*v1Cent;
	
	v1Peri = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "AccStat"), zbins[j], zbins[j+1], 0, 2.4, centCut, centMax);
	systPeri = systPeri + v1Peri*v1Peri;
      }

      v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffMisMod"), zbins[j], zbins[j+1], 0, 2.4, centMin, centMax);
      syst = syst + v1*v1;
      if (isPbPb) {
      v1Cent = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffMisMod"), zbins[j], zbins[j+1], 0, 2.4, centMin, centCut);
      systCent = systCent + v1Cent*v1Cent;

      v1Cent = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", "EffMisMod"), zbins[j], zbins[j+1], 0, 2.4, centCut, centMax);
      systPeri = systPeri + v1Peri*v1Peri;
      }
      for (int i=0; i<nSyst; i++){
	v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", systName[i].c_str()), zbins[j], zbins[j+1], 0, 2.4, centMin, centMax);
	syst = syst + v1*v1;
	if (isPbPb) {
	v1Cent = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", systName[i].c_str()), zbins[j], zbins[j+1], 0, 2.4, centMin, centCut);
	systCent = systCent + v1Cent*v1Cent;
	v1Peri = readSyst(Form("../Fitter/Systematics/csv/syst_%s_NJpsi_%s_%s_tnp%s.csv", jtptTags[k].c_str(), isPr?"prompt":"nonprompt", isPbPb?"PbPb":"PP", systName[i].c_str()), zbins[j], zbins[j+1], 0, 2.4, centCut, centMax);
	systPeri = systPeri + v1Peri*v1Peri;
	}
      }
      if (isPbPb) {
	fileAll<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 180, "<< sqrt(syst) << endl;
	fileAll<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, "<<centCut<<", "<< sqrt(systCent) << endl;
	fileAll<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", "<<centCut<<", 180, "<< sqrt(systPeri) << endl;
      }
      else
	fileAll<< "0, 2.4, 6.5, 100, "<<zbins[j]<<", "<< zbins[j+1]<<", 0, 200, "<< sqrt(syst) << endl;
    }
    fileAll.close();
  }
}


void oniaTree::TnpToy(int min, int max) {
}

void oniaTree::TnpStat(string caseLabel) {

}


void oniaTree::AccEffSyst_all() {
  //AccEffStat("_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
  //isAcc = true;
  //AccEffStat("_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
  //isAcc = false;
  //AccEffMisMod("_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
  TnpSyst("_pt_SizeDoubled_rap_1bin1010_15Bins_centBins_NoWeights");
}
