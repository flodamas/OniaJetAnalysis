cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -b -l -q  ~/LoadUnfoldinLib.C+ prepareInputsNewPrNominal.cc+\($1\)

cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding
root -b -l -q  ~/LoadUnfoldinLib.C+ dataUnf/code/createRooUnfoldResponseNewPrNominal.cxx+\($1\)
root -b -l -q  ~/LoadUnfoldinLib.C+ dataUnf/code/unfoldStepNewPrNominal.cxx+\($1\)
root -b -l -q  ~/LoadUnfoldinLib.C+ dataUnf/code/unfoldStepNewPrNominalDiag.cxx+\($1\)
