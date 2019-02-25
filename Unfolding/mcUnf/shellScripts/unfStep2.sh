cd /home/ikucher/newRooUnfoldVersion/RooUnfold/mcUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  mcUnfNewMidBins/code/createRooUnfoldResponse.cxx\(2\)
root -l -q -b  mcUnfNewMidBins/code/unfoldStep.cxx\(2\)

