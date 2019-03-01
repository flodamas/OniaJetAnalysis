# Jpsi in jets 2019 fitter
Repository for code regarding Jpsi in jet in 5.02TeV data.

## Fitter
* fitter.C: master file
 * The fitter can be run with the following command:  root -l -b -q fitter.C+'("name of work dir")' 
* clean_wd.sh: source this script to clean the working directory (RooDatasets, plots)
* Input/createInput.py: automatically create the input (configuration) files
* results2tree.C: parse an output directory and put the results into a TTree
 * results2tree("name of work dir","var1,var2,var3")
 * NB: the names of the variables ("var1", etc. in the example above) should be the actual name of the variable in the workspace, without "\_PP" or "\_PbPb". e.g. "sigma1\_Jpsi"
* plotNJJ.C: plot the results for fragmentation function (uses the trees from results2tree.C as input)

## Efficiency
* makeAccEff.C: make the histograms (acceptance and efficiency) from the onia trees

## Unfolding
* dataUnf for DATA
* mcUnf for MC