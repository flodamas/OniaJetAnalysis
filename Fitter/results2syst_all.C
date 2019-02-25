#ifndef results2syst_all_C
#define results2syst_all_C

#include <iostream>

#include "TString.h"

#include "results2syst.C"

using namespace std;

void results2syst_all(bool testChi2 = true, bool isMid = true, int jtPtRange = 0)
{
  bool is16004 = false;
  TString app("");
  if (isMid) app = "016";
  else app = "1624";

  TString jtapp("");
  if (jtPtRange==-1) jtapp = "lowJtPt";
  else if (jtPtRange==1) jtapp = "highJtPt";  
  else if (jtPtRange==0) jtapp = "midJtPt";
  else if (jtPtRange==100) {jtapp = "NoJets"; app = "total";}
  else {cout<<"[ERROR] Please check the initial settings!"<<endl; return;}

  TString workdir_nominal(Form("DataFits/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  
  bool isXS = false;
  if (jtPtRange==100) isXS = true;
  TString workdir_sys("");
  TString sysFileName("");
  
  // Background systematics
  //NJpsi
  workdir_sys = Form("%s,%s,%s,%s,%s",workdir_nominal.Data(),Form("DataFits_pval25/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()),Form("DataFits_pval100/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()),Form("DataFits_bkgExp/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()),Form("DataFits_mass2634/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  cout << "<<<<<< Computing invariant mass background systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_massBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PbPb",true,"N_Jpsi",testChi2,isMid);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_massBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PP",true,"N_Jpsi",testChi2, isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_massBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PbPb",true,"b_Jpsi",testChi2,isMid);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_massBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_massBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_massBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_massBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_massBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass bkg.",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  
  
  
  // Signal systematics
  //NJpsi
  workdir_sys = Form("%s,%s,%s",workdir_nominal.Data(),Form("DataFits_constrained/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()),Form("DataFits_CBGauss/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  cout << "<<<<<< Computing invariant mass signal systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_massSig.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PbPb",true,"N_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_massSig.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PP",true,"N_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_massSig.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PbPb",true,"b_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_massSig.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_massSig.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_massSig.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_massSig.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_massSig.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"inv. mass sig.",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  
  
  // CtauErr systematics
  //NJpsi
  workdir_sys = Form("%s,%s",workdir_nominal.Data(),Form("DataFits_ctauErrTot/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  cout << "<<<<<< Computing ctau error systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_ctauErr.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PbPb",true,"N_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_ctauErr.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PP",true,"N_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_ctauErr.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PbPb",true,"b_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_ctauErr.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_ctauErr.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_ctauErr.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_ctauErr.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_ctauErr.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau error",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  
  
  
  // CtauTrue systematics
  //NJpsi
  workdir_sys = Form("%s,%s",workdir_nominal.Data(),Form("DataFits_ctauReco/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  cout << "<<<<<< Computing ctau true systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_ctauTrue.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PbPb",true,"N_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_ctauTrue.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PP",true,"N_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_ctauTrue.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PbPb",true,"b_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_ctauTrue.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_ctauTrue.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_ctauTrue.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_ctauTrue.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_ctauTrue.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau true",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  
  
  
  // CtauRes systematics
  //NJpsi
  workdir_sys = Form("%s,%s",workdir_nominal.Data(),Form("DataFits_ctauResPromptMC/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data())/*,Form("DataFits_ctauResNonPromptMC/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data())*/);
  //workdir_sys = Form("%s,%s",workdir_nominal.Data(),Form("DataFits_%s_2D_2CB_polBkg_ctauResPromptMC",jtapp.Data(),app.Data())); // Remove 2D fits with NP MC resolution from the systematics
  cout << "<<<<<< Computing ctau resolution systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_ctauRes.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PbPb",true,"N_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_ctauRes.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PP",true,"N_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_ctauRes.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PbPb",true,"b_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_ctauRes.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_ctauRes.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_ctauRes.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_ctauRes.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_ctauRes.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau res.",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  
  
  
  // CtauBkg systematics
  //NJpsi
  workdir_sys = Form("%s,%s",workdir_nominal.Data(),Form("DataFits_ctauBkgTemplate/DataFits_%s/DataFits_%s",jtapp.Data(),app.Data()));
  cout << "<<<<<< Computing ctau background systematics out of : " << workdir_sys.Data() << " >>>>>>" << endl;
  cout << endl;
  //cout << "For N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_PbPb_ctauBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PbPb",true,"N_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_PP_ctauBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PP",true,"N_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //b fraction
  //cout << "For b fraction in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_BJpsi_PbPb_ctauBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PbPb",true,"b_Jpsi",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For b fraction in pp: " << endl;
  sysFileName = Form("syst_%s_%s_BJpsi_PP_ctauBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PP",true,"b_Jpsi",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_prompt
  //cout << "For prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_prompt_PbPb_ctauBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PbPb",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_prompt_PP_ctauBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PP",true,"N_Jpsi_prompt",testChi2,isMid, isXS);
  cout << endl;
  
  //NJpsi_nonprompt
  //cout << "For non prompt N_Jpsi in Pbpb: " << endl;
  //sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PbPb_ctauBkg.csv",jtapp.Data(),app.Data());
  //results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PbPb",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  //cout << endl;
  cout << "For non prompt N_Jpsi in pp: " << endl;
  sysFileName = Form("syst_%s_%s_NJpsi_nonprompt_PP_ctauBkg.csv",jtapp.Data(),app.Data());
  results2syst(workdir_sys.Data(),sysFileName.Data(),"ctau bkg.",1,"PP",true,"N_Jpsi_nonprompt",testChi2,isMid, isXS);
  cout << endl;
  
  return;
}
#endif
