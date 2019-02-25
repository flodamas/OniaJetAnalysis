cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(2,1.2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponse.cxx\(2,1.2\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStep.cxx\(2,1.2\)

