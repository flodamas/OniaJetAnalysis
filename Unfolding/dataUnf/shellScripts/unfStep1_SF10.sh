cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(1,1.0\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponse.cxx\(1,1.0\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStep.cxx\(1,1.0\)

