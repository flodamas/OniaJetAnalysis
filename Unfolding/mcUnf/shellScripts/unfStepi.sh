cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/mcUnf/code
root -l -q -b  prepareInputs.cc+\($1\)

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/
root -l -q -b  mcUnf/code/createRooUnfoldResponse.cxx+\($1\)
root -l -q -b  mcUnf/code/unfoldStep.cxx+\($1\)

