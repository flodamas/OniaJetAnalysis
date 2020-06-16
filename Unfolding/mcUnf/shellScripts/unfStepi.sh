cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding/mcUnf/code
root -l -q -b  ~/LoadUnfoldinLib.C+ prepareInputs.cc+\($1\)

cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding
root -l -q -b  ~/LoadUnfoldinLib.C+ mcUnf/code/createRooUnfoldResponse.cxx+\($1\)
root -l -q -b  ~/LoadUnfoldinLib.C+ mcUnf/code/unfoldStep.cxx+\($1\)

