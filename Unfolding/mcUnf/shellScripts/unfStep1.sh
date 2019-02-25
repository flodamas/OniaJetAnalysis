cd /home/ikucher/newRooUnfoldVersion/RooUnfold/mcUnfNewMidBins/code
root -l -q -b  prepareInputs.cc+\(1\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  mcUnfNewMidBins/code/createRooUnfoldResponse.cxx\(1\)
root -l -q -b  mcUnfNewMidBins/code/unfoldStep.cxx\(1\)

