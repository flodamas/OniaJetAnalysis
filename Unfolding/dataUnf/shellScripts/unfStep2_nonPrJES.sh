cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsNonPrJES.cc+\(2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseNonPrJES.cxx\(2\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepNonPrJES.cxx\(2\)

