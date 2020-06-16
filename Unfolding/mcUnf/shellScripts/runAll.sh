#!/bin/sh

for i in {0..30}
do
    echo "Step $i : "
    sh unfStepi.sh $i
    echo "Done!"
done

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/mcUnf/code
root -l -q -b ~/LoadUnfoldinLib.C+ PlotRatios_MCUnfoldedTruth.cc+

echo "Finished!"
