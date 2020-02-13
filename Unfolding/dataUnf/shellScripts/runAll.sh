#!/bin/sh

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -l -q -b create2DMeas_data_statErrs.cc+\(1,1\)
root -l -q -b create2DMeas_data_statErrs.cc+\(1,0\)

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/shellScripts
for i in {1..99}
do
    echo "Step $i : "
    sh unfStepi.sh $i
    echo "Done!"
done

cd /Users/diab/Phd_LLR/JpsiJetAnalysisPbPb2019/JpsiInJetsPbPb/Unfolding/dataUnf/code
root -l -q -b PlotRatios_DataUnfolded_afterDiag.cc+
root -l -q -b finalResults_statErrs.cc+
#root -l -q -b plotFinalResults.cc+
echo "Finished!"
