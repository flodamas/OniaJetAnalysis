cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNonPrJES.cc+\(4\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNonPrJES.cxx\(4\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNonPrJES.cxx\(4\)

