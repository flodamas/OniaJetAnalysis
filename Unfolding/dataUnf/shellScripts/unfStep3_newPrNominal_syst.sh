cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNewPrNominalSyst.cc+\(3\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNewPrNominalSyst.cxx\(3\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNewPrNominalSyst.cxx\(3\)

