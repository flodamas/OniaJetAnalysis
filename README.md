# Jpsi in jets 2019 fitter
Repository for code regarding Jpsi in jet in 5.02TeV data.

## Fitter
* fitter.C: master file
  - The fitter can be run with the following command:  root -l -b -q fitter.C+'("name of work dir",.....)' -- Please check the settings!!
  - The fitter reads Input/name_of_work_dir/InputTrees.txt for the input trees and Input/name_of_work_dir/InitialParam_\*.csv for the initial parameters
  - if the input directory has many subdirectories it is enough to give the name of the main work directory
* clean_wd.sh: source this script to clean the working directory (RooDatasets, plots)
* results2tree.C: parse an output directory and put the results into a TTree
  - results2tree("name of work dir","var1,var2,var3")
  - NB: the names of the variables ("var1", etc. in the example above) should be the actual name of the variable in the workspace, without "\_PP" or "\_PbPb". e.g. "sigma1\_Jpsi"
  - results2tree_all.C runs results2tree for all directories
* results2syst.C: creates .csv files for the different systematic uncertainties using the trees created by results2tree as input. The files are located in Fitter/Systematics/csv/
  - can be run for all systematics useing results2syst_all.C
* plotNJJ.C: plot the results for fragmentation function (uses the trees from results2tree.C as input)
* plotFitParams.C: can be used to plot different parameters
  - most important functions: getUnfoldingInput and getUnfoldingInput_all that give the histograms used as input for the data unfolding (uses the trees from results2tree.C as input)

## Efficiency
* makeAccEff.C: make the histograms (acceptance and efficiency) from the onia trees
 * For a general case run the following in a root session:
 ```
   .L makeAccEff.C+
   oniaTree t(..,..,..) -- make the class onia tree, the settings correspond to prompt/nonprompt, pp/PbPb, Acc/Eff
   t.EffCalc("caseTag") -- caseTag corresponds to the output directory and is also used for different settings in the code, t should be an efficiency tree (previous step)
   t.AccCalc("caseTag") -- t should be an acceptance tree (previous step)
   AccEffCalc("caseTag")
  ```
 * To get the exact same settings as the jpsi in jet analysis (HIN-19-007), one can simply run runAllAccEff.sh that has the correc settings
 * For the systematic uncertainties (after making the oniaTree t):
   * t.AccEffStatToy_Eff(n, caseTag) creates n toys for the efficiency of caseTag
   * t.AccEffStatToy_Acc(n, caseTag) creates n toys for the acceptance of caseTag
   * t.AccEffStat(caseTag) calculates the statistical uncertainty on the acceptance and efficiency, creates a csv file with the uncertainty in Fitter/Systematics/csv/
   * t.AccEffMisMod(caseTag) calculated the uncertainty for the efficiency uncertainty. For this to work, one must run t.AccEffStatToy_Eff(n, caseTag) with n=100 for the nominal case and another case with the same settings without the "\_NoWeights" in the caseTag.
   * t.TnpSyst(caseTag) calculates the uncertainties of the Tag&probe scale factors (also adds in quadratire the rest of the uncertainties and produces the total uncertainty .csv file)
 
## Unfolding
* unfolding MC
  * Unfolding/mcUnf/shellScripts/runAll.sh runs all macros  
  * in Unfolding/mcUnf/code:
    * inputParams.h contains all parameters and settings
    * unfoldAllSteps.C runs all steps (contains different functions for different settings)
    * prepareInputs.cc creates the 4 dimentional THnSparseF with the correct cuts and normalization (e.g. flat prior)
    * createRooUnfoldResponse.cxx converts the THnSparseF to a transfer matrix
    * unfoldStep.cxx does the unfolding (1SI)
    * PlotRatios_MCUnfoldedTruth.cc plot the closure tests
    
* unfolding data
  * runAll.sh in Unfolding/dataUnf/shellScripts/runAll.sh runs all macros
  * in Unfolding/dataUnf/code:
    * inputParams.h contains all parameters and settings
    * unfoldAllSteps.C runs all steps (contains different functions for different settings, also conatins functions for the systematic uncertainties)
    * makeRebinMatrix.C and createRooUnfoldResponseDiag.cxx create an empty matrix with the unfolding bins that would be used in the unfolding step
    * create2DMeas_data.cc takes the 1D histograms from Fitter/plotFitParams.C and combines them into 2D histogram for the unfolding
    * prepareInputsNewPrNominal.cc creates the 4 dimentional THnSparseF with the correct cuts and normalization (e.g. flat prior)
    * createRooUnfoldResponseNewPrNominal.cxx converts the THnSparseF to a transfer matrix
    * unfoldStepNewPrNominal.cxx and unfoldStepNewPrNominalDiag.cxx do the unfolding (1SI) with the correct error treatement
    * PlotRatios_DataUnfolded_afterDiag.cc plot the comparison between the SI
    * finalResults_statErrs.cc makes the root files that wold be used for the nominal and systematic results
    * smearTransferMatrix.cc makes the 100 toys of the transfer matrix 
    * trStatSystUnc.cc calculates the statistical uncertainty on the transfer matrix
    * plotFinalResults.cc and plotCentResults.cc plot the final results
