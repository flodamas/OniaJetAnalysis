cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsSyst.cc+\(1\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseSyst.cxx\(1\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepSyst.cxx\(1\)

