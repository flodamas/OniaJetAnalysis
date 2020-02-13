cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -l -q -b  prepareInputsNewPrNominal.cc+\($1\)

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/
root -l -q -b  dataUnf/code/createRooUnfoldResponseNewPrNominal.cxx+\($1\)
root -l -q -b  dataUnf/code/unfoldStepNewPrNominal.cxx+\($1\)
root -l -q -b  dataUnf/code/unfoldStepNewPrNominalDiag.cxx+\($1\)
