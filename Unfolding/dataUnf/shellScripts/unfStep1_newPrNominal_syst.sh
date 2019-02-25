cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNewPrNominalSyst.cc+\(1\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNewPrNominalSyst.cxx\(1\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNewPrNominalSyst.cxx\(1\)

