#include "makeAccEff.C"

makeAllAccEffMaps(bool doPP = true) {
  if (doPP)
    makeAllAccEffMaps_pp();
  else 
    makeAllAccEffMaps_pbpb();
}

makeAllAccEffMaps_pp() {
  oniaTree t_pp_Eff(0,1,0);
  oniaTree t_pp_Acc(0,1,1);
  t_pp_Eff.EffCalc("_pt_SizeDoubled_centBins_NoWeights");
  t_pp_Acc.AccCalc("_pt_SizeDoubled_centBins_NoWeights");
}

makeAllAccEffMaps_pbpb() {
  oniaTree t_pbpb_Eff(1,1,0);
  t_pbpb_Eff.EffCalc();
  t_pbpb_Eff.AccEffCalc();
}
