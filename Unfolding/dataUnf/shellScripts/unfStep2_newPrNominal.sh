cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNewPrNominal.cc+\(2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNewPrNominal.cxx\(2\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNewPrNominal.cxx\(2\)

