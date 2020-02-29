//#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include "inputParams.h"
//#endif

#include "prepareInputs.cc"
#include "createRooUnfoldResponse.cxx"
#include "unfoldStep.cxx"
#include "PlotRatios_MCUnfoldedTruth.cc"


void unfoldAllSteps_Def(bool setSameSample = true, bool setMc2015 = false, bool setFlatPrior=true, int SetCentShift = 0, bool setSmearMeas = false, bool setDataDist = false, int nIt = 3, int nSI = 99, int nSI_pp = 6) {
  sameSample = setSameSample;
  mc2015 = setMc2015;
  flatPrior = setFlatPrior;
  nIter = nIt;
  nSIter = nSI;
  nSIter_pp = nSI_pp;
  centShift = SetCentShift;
  smearMeas = setSmearMeas;
  dataDist = setDataDist;
  
  for (int i=1;i<=nSIter;i++)
    {
      prepareInputs(i);
      createRooUnfoldResponse(i);
      unfoldStep(i);
      std::cout<<"[INFO] SI "<<i<<" done"<<endl;
    }  
  std::cout<<"[INFO] plotting done"<<endl;
  PlotRatios_MCUnfoldedTruth();
}


//void unfoldAllSteps_withSmearing(int SetCentShift = 0, bool setDataDist = false, int nIt = 3, int nSI = 99, int nSI_pp = 6, int istart=1, int istop=20) {
void unfoldAllSteps(int SetCentShift = 0, bool setDataDist = false, int nIt = 3, int nSI = 99, int nSI_pp = 6, int istart=1, int istop=20) {
  sameSample = true;
  mc2015 = false;
  flatPrior = true;
  smearMeas = true;
  centShift = SetCentShift;
  dataDist = setDataDist;
  nIter = nIt;
  nSIter = nSI;
  nSIter_pp = nSI_pp;
  
  for (int i=istart;i<=istop;i++)
  {
    prepareInputs(i);
    createRooUnfoldResponse(i);
    unfoldStep(i);
    std::cout<<"[INFO] SI "<<i<<" done"<<endl;
  }  
  std::cout<<"[INFO] plotting done"<<endl;
  PlotRatios_MCUnfoldedTruth();
}
