cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(3,1.0\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponse.cxx\(3,1.0\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStep.cxx\(3,1.0\)

