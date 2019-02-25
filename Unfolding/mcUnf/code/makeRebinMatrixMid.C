void makeRebinMatrixMid(){

  // find the transfer matrix to take a measurement of 50 bins in z and 15 bins in pT to one with 5 bins in z and 3 bins in pT.  

  // number of fine and coarse bins
  int nfpt = 15;
  int ncpt = 3;
  
  int nfz = 49;
  int ncz = 7;

  // bin width in coarse bins
  double wcpt = 10.;
  double wcz = 0.14;

  // bin width in narrow bins
  double wfpt = wcpt*ncpt/nfpt;
  double wfz = wcz*ncz/nfz;

  // overall bin boundaries
  double lopt = 15.;
  double loz = 0.02;
  double hipt = lopt + wcpt*ncpt;
  double hiz = loz + wcz*ncz;
  
    
  TH2D *h2dpt = new TH2D("h2dpt","h2dpt",nfpt,lopt,hipt,ncpt,lopt,hipt);
  TH2D *h2dz = new TH2D("h2dz","h2dz",nfz,loz,hiz,ncz,loz,hiz);

  int fDim = 4.;  
  Int_t* bins_sparce = new Int_t[fDim];
  Double_t *xmin_sparce = new Double_t[fDim];
  Double_t *xmax_sparce = new Double_t[fDim];

  bins_sparce[0] = 15;
  xmin_sparce[0] = 15.0;
  xmax_sparce[0] = 45.0;

  bins_sparce[1] = 49;
  xmin_sparce[1] = 0.02;
  xmax_sparce[1] = 1.0;

  bins_sparce[2] = 3;
  xmin_sparce[2] = 15.0;
  xmax_sparce[2] = 45.0;

  bins_sparce[3] = 7;
  xmin_sparce[3] = 0.02;
  xmax_sparce[3] = 1.0;

  //initial not normalized 4D
  THnSparseF * fSparse = new THnSparseF("hs", "hs", fDim, bins_sparce, xmin_sparce, xmax_sparce);
  fSparse->Sumw2();
  fSparse->CalculateErrors();

  Int_t* bins_sparce2 = new Int_t[fDim];
  Double_t *xmin_sparce2 = new Double_t[fDim];
  Double_t *xmax_sparce2 = new Double_t[fDim];

  bins_sparce2[2] = 15;
  xmin_sparce2[2] = 15.0;
  xmax_sparce2[2] = 45.0;

  bins_sparce2[3] = 49;
  xmin_sparce2[3] = 0.02;
  xmax_sparce2[3] = 1.0;

  bins_sparce2[0] = 3;
  xmin_sparce2[0] = 15.0;
  xmax_sparce2[0] = 45.0;

  bins_sparce2[1] = 7;
  xmin_sparce2[1] = 0.02;
  xmax_sparce2[1] = 1.0;
  
  THnSparseF * fSparseInv = new THnSparseF("hsInv", "hsInv", fDim, bins_sparce2, xmin_sparce2, xmax_sparce2);
  fSparseInv->Sumw2();
  fSparseInv->CalculateErrors();
  

  for(int ipt = 0; ipt < ncpt; ipt++){    
    for(int iz = 0; iz < ncz; iz++){

      float izMid = iz*wcz+wcz/2+loz;

      //cout << "izMid = " << izMid << endl;
      
      float ijetMid = 15.+ipt*wcpt+wcpt/2;

      
      int jzlo = iz*nfz/ncz+loz;
      int jzhi = (iz+1)*nfz/ncz+loz;
      
      int jptlo = ipt*nfpt/ncpt;
      int jpthi = (ipt+1)*nfpt/ncpt;
     
      for(int jpt = jptlo; jpt < jpthi; jpt++){	
	for(int jz = jzlo; jz < jzhi; jz++){
	  
	  h2dpt->SetBinContent(jpt+1,ipt+1,1.);
	  h2dz->SetBinContent(jz+1,iz+1,1.);

	  float jzMid = jz*wfz+wfz/2+loz;

	  //cout << "jzMid = " << jzMid << endl;
	  
	  float jjetMid = 15.+jpt*wfpt+wfpt/2;
	  	  
	  const double x[4] = {jjetMid,jzMid,ijetMid,izMid};
	  int bin = fSparse->GetBin(x);
	  fSparse->SetBinContent(bin, 1.);

	  const double x2[4] = {ijetMid,izMid,jjetMid,jzMid};
	  int bin2 = fSparseInv->GetBin(x2);
	  fSparseInv->SetBinContent(bin2, 1.);
	  
	  
	  //cout << "bin = " << bin << " izMid = " << izMid << " jzMid = " << jzMid << endl;
	}
      }
      
    }
  }

  /*
  TCanvas *c1=new TCanvas("c1","c1",600,600);
  h2dpt->Draw("zcol");
  TCanvas *c2=new TCanvas("c2","c2",600,600);
  h2dz->Draw("zcol");
  */
  
    
  TCanvas *c1_4Dpr=new TCanvas("c1_4Dpr","c1_4Dpr",600,600);
  fSparseInv->Projection(0,2)->Draw("zcol");
  TCanvas *c2_4Dpr=new TCanvas("c2_4Dpr","c2_4Dpr",600,600);
  fSparseInv->Projection(1,3)->Draw("zcol");

  TFile *fout = new TFile("diag4DMatrixInvMid.root","RECREATE");
  fSparse->Write("diag4DMatrix");  
  fSparseInv->Write("diag4DMatrixInv");

}
