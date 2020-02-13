#!/bin/sh

for i in {1..99}
do
    echo "Step $i : "
    sh unfStepi.sh $i
    echo "Done!"
done

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/mcUnf/code
root -l -q -b PlotRatios_MCUnfoldedTruth.cc+

echo "Finished!"
