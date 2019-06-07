#ifndef results2tree_all_C
#define results2tree_all_C

#include <iostream>

#include "TString.h"

#include "results2tree.C"

using namespace std;

void results2tree_all()
{ 
  TString workdir("");
  ///////////////////////////////
  //        Nominal            //
  ///////////////////////////////
  workdir = "DataFits/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1); 
  // total XS
  workdir = Form("DataFits/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  ///////////////////////////////
  //   Background variations   //
  ///////////////////////////////
  workdir = "DataFits_bkgExp/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval25/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval100/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //low and high Jt Pt
  workdir = "DataFits_bkgExp/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval25/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval100/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_bkgExp/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval25/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_pval100/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //total XS
  workdir = Form("DataFits_bkgExp/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = Form("DataFits_pval25/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = Form("DataFits_pval100/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = Form("DataFits_mass2634/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  ///////////////////////////////
  //    Signal variations      //
  ///////////////////////////////
  workdir = "DataFits_constrained/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_CBGauss/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  // low and high JtPt
  workdir = "DataFits_constrained/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_CBGauss/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_constrained/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = "DataFits_CBGauss/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //total XS
  workdir = Form("DataFits_constrained/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  workdir = Form("DataFits_CBGauss/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  ///////////////////////////////
  //    CtauErr variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauErrTot/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //low and high JtPt
  workdir = "DataFits_ctauErrTot/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_ctauErrTot/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //total XS
  workdir = Form("DataFits_ctauErrTot/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  ///////////////////////////////
  //   CtauTrue variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauReco/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //low and JtPt
  workdir = "DataFits_ctauReco/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_ctauReco/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //totalXS
  workdir = Form("DataFits_ctauReco/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  
  ///////////////////////////////
  //   CtauRes variations      //
  ///////////////////////////////
  workdir = "DataFits_ctauResPromptMC/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //low and high JtPt
  workdir = "DataFits_ctauResPromptMC/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_ctauResPromptMC/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //total XS
  workdir = Form("DataFits_ctauResPromptMC/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  ///////////////////////////////
  //   CtauBkg variations      //
  ///////////////////////////////
  workdir = "DataFits_ctauBkgTemplate/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //low and high JtPt
  workdir = "DataFits_ctauBkgTemplate/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  workdir = "DataFits_ctauBkgTemplate/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  //total XS
  workdir = Form("DataFits_ctauBkgTemplate/DataFits_NoJets/DataFits_total");
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.4, "AccEff", 1);

  return;
}

#endif
