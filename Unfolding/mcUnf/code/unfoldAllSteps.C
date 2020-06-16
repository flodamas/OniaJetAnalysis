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
void unfoldAllSteps(int SetCentShift = 0, bool setDataDist = false, int nIt = 7, int nSI = 25, int nSI_pp = 0, int istart=11, int istop=20) {
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
  if (istop>nSI)
    PlotRatios_MCUnfoldedTruth();
  std::cout<<"[INFO] plotting done"<<endl;
}


void unfoldAllSteps_something(int nIt = 3, int nSI = 25, int nSI_pp = 6, int istart=1, int istop=50) {
  sameSample = false;
  mc2015 = false;
  flatPrior = false;//true;
  smearMeas = false;
  centShift = 0;//SetCentShift;
  dataDist = false;//setDataDist;
  nIter = nIt;
  nSIter = nSI;
  nSIter_pp = nSI_pp;
  
  max_jetpt_real = 70.;
  nBinJet_gen = 32;
  
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

void unfoldAllSteps_default(int istart=1, int istop=4, bool setCent=false, bool setPeri=false) {
  doCent=setCent;
  doPeri=setPeri;
  for (int i=istart;i<=istop;i++) {
    prepareInputs(i);
    createRooUnfoldResponse(i);
    unfoldStep(i);
    std::cout<<"[INFO] SI "<<i<<" done"<<endl;
  }
  std::cout<<"[INFO] plotting done"<<endl;
  PlotRatios_MCUnfoldedTruth();
}
