cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(3,1.2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponse.cxx\(3,1.2\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStep.cxx\(3,1.2\)

