cd /home/ikucher/newRooUnfoldVersion/RooUnfold/dataUnfNewMidBins/code
root -l -q -b  prepareInputsSyst.cc+\(2\)

cd /home/ikucher/newRooUnfoldVersion/RooUnfold
root -l -q -b  dataUnfNewMidBins/code/createRooUnfoldResponseSyst.cxx\(2\)
root -l -q -b  dataUnfNewMidBins/code/unfoldStepSyst.cxx\(2\)

