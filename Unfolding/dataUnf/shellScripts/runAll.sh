#!/bin/sh

cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -b -l -q ~/LoadUnfoldinLib.C+ create2DMeas_data.cc+\(1,1\)
root -b -l -q ~/LoadUnfoldinLib.C+ create2DMeas_data.cc+\(1,0\)
root -b -l -q ~/LoadUnfoldinLib.C+ makeRebinMatrix.C+
root -b -l -q ~/LoadUnfoldinLib.C+ createRooUnfoldResponseDiag.cxx+

cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding/dataUnf/shellScripts
for i in {1..2}
do
    echo "Step $i : "
    sh unfStepi.sh $i
    echo "Done!"
done

cd /home/llr/cms/diab/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -b -l -q ~/LoadUnfoldinLib.C+ PlotRatios_DataUnfolded_afterDiag.cc+
root -b -l -q ~/LoadUnfoldinLib.C+ finalResults_statErrs.cc+
root -l -q -b plotFinalResults.cc+
echo "Finished!"
