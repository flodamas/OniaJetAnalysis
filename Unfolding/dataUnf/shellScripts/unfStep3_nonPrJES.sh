cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNonPrJES.cc+\(3\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNonPrJES.cxx\(3\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNonPrJES.cxx\(3\)

