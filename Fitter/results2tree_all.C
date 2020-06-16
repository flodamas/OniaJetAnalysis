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
  /*
  workdir = "DataFits/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  
  ///////////////////////////////
  //   Background variations   //
  ///////////////////////////////
  workdir = "DataFits_ExpChev/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_ExpChev/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ExpChev/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ExpChev/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ExpChev/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_ExpChev/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  workdir = "DataFits_mass2634/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_mass2634/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_mass2634/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_mass2634/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  ///////////////////////////////
  //    Signal variations      //
  ///////////////////////////////

  workdir = "DataFits_freeAlpha/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_freeAlpha/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeAlpha/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeAlpha/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeAlpha/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_freeAlpha/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  workdir = "DataFits_freeN/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_freeN/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeN/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeN/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeN/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_freeN/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  workdir = "DataFits_freeRSigma/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_freeRSigma/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeRSigma/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeRSigma/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_freeRSigma/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_freeRSigma/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  workdir = "DataFits_CBGaus/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_CBGaus/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_CBGaus/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_CBGaus/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_CBGaus/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_CBGaus/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 

  */
  ///////////////////////////////
  //    CtauErr variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauErrTot/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_ctauErrTot/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauErrTot/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauErrTot/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauErrTot/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_ctauErrTot/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  
  ///////////////////////////////
  //   CtauRes variations      //
  ///////////////////////////////
  workdir = "DataFits_ctauResPrMC/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_ctauResPrMC/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauResPrMC/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauResPrMC/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauResPrMC/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_ctauResPrMC/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  
  ///////////////////////////////
  //   CtauTrue variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauRecoTemp/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_ctauRecoTemp/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauRecoTemp/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauRecoTemp/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauRecoTemp/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_ctauRecoTemp/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  
  ///////////////////////////////
  //   CtauBkg variations      //
  ///////////////////////////////
  workdir = "DataFits_ctauBkgTemp/DataFits_midJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits_ctauBkgTemp/DataFits_lowJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauBkgTemp/DataFits_lowerJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauBkgTemp/DataFits_lowestJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1);
  workdir = "DataFits_ctauBkgTemp/DataFits_highJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  workdir = "DataFits_ctauBkgTemp/DataFits_higherJtPt";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, 0.3, "AccEff", 1); 
  
  return;
}
  
#endif
