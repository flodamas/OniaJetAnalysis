#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/monster.h"
#include "Systematics/syst.h"
//#include "getMCResults.C"

#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TArrow.h"
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"

using namespace std;

////////////////
// PARAMETERS //
////////////////

//#ifndef poiname_check
//#define poiname_check
//const char* poiname       = "N_Jpsi";
//#endif
//const char* ylabel        = "N_{JJ}";

bool  dopp              = false;
bool  doPbPb            = true;
bool  doprompt          = false;  // prompt Jpsi
bool  dononprompt       = true;  // nonprompt Jpsi
bool  plotMid           = true;
bool  plotFwd           = true;
int   jtPtRange         = 0;
bool  applyEff          = true;
bool  applyAcc          = true;
bool  doLogPt           = false;
bool  includeEffSyst    = false;
bool  excludeNonFitSyst = false;
bool  plotFwdMid        = false;
bool  isPreliminary     = true;
bool  plotUnfolded      = true;
bool  underflowOff      = true;
bool  mcON              = true;
bool  bffON             = false;
bool  ctauCut           = false;
string nameTag_base     = "";    // can put here e.g. "_prompt", "_nonprompt", ...

const bool useNcoll = true; // false -> use TAA / NMB, true -> use Ncoll / lumiPbPb

double histMax = 0;
//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void setOptions(bool adoPbPb, bool adopp, bool adoprompt, bool adononprompt, bool aplotMid, bool aplotFwd, bool aexcludeNonFitSyst, string anameTag_base="", int ajtPtRange=0, bool aplotUnfolded = true, bool aunderflowOff = true, bool amcON = true, bool abffON = false, bool actauCut = false);
void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high);

void plotNJJ(vector<anabin> thecats, string xaxis, string workDirName);

void centrality2npart(TGraphAsymmErrors* tg, bool issyst=false, double xshift=0.);
int color(int i);
int markerstyle(int i);
string nameTag;

class njj_input {
public:
  double npp;
  double dnpp_stat;
  double systpp;
  double naa;
  double dnaa_stat;
  double systaa;
  double taa;
  double lumipp;
  double lumiaa;
  double ncoll;
  double bfracpp;
  double dbfracpp;
  double systbfracpp;
  double bfracaa;
  double dbfracaa;
  double systbfracaa;
  double correlaa;
  double correlpp;
  syst   statpp;
};

map<anabin, njj_input > readResults(const char* resultsFile);
void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color);

/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////
void plotZ(string workDirName) {
  string xaxis = "z";
  vector<anabin> theCats;
  if (plotMid){
    if (plotFwd)
      theCats.push_back(anabin(0.0, 1.0, 0, 2.4, 6.5, 100, 0, 200));
    else
      theCats.push_back(anabin(0.0, 1.0, 0, 1.6, 6.5, 100, 0, 200));
  }
  else if (plotFwd)
    theCats.push_back(anabin(0.0, 1.0, 1.6, 2.4, 3.0, 100, 0, 200));

  plotNJJ(theCats,xaxis,workDirName);
};

void plotAll(string workDirName) {
  plotZ(workDirName);
};

void doAllplots() {
  //doPbPb, dopp, doprompt, dononprompt, plotMid, plotFwd, excludeNonFitSyst, nameTag_base="", jtPtRange, plotUnfolded, underflowOff, mcON, bffON, ctauCut

  ///////////////////////////////////////////////////////
  //        pp, prompt, mid and fwd, ctauCut           //
  ///////////////////////////////////////////////////////
  setOptions(false, true, true, false, true, true, false, "", 0, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_midJtPt");

  setOptions(false, true, true, false, true, true, false, "", 1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_highJtPt");

  setOptions(false, true, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowJtPt");

  setOptions(false, true, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowerJtPt");

  ///////////////////////////////////////////////////////
  //       PbPb, prompt, mid and fwd, ctauCut          //
  ///////////////////////////////////////////////////////
  setOptions(true, false, true, false, true, true, false, "", 0, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_midJtPt");

  setOptions(true, false, true, false, true, true, false, "", 1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_highJtPt");

  setOptions(true, false, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowJtPt");

  setOptions(true, false, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowerJtPt");

  ///////////////////////////////////////////////////////
  //      PbPb/pp, prompt, mid and fwd, ctauCut        //
  ///////////////////////////////////////////////////////
  setOptions(true, true, true, false, true, true, false, "", 0, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_midJtPt");

  setOptions(true, true, true, false, true, true, false, "", 1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_highJtPt");

  setOptions(true, true, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowJtPt");

  setOptions(true, true, true, false, true, true, false, "", -1, false, false, false, false, true);
  printOptions();
  plotAll("DataFits_ctauCut/DataFits_lowerJtPt");

};

/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plotNJJ(vector<anabin> thecats, string xaxis, string outputDir) {
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent" or "z"
  // outputDir is the directory to save the plots
  
  anabin binToSkip(0.3, 0.44, 0, 1.6, 6.5, 35, 0, 200);

  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  
  TString poi("NJpsi");
  if (doprompt) poi = "NJpsi_prompt";
  if (dononprompt) poi = "NJpsi_nonprompt";

  TString sTag(Form("%s_%s", (jtPtRange==-1)?"lowJtPt":(jtPtRange==1)?"highJtPt":"midJtPt", (plotMid)?(plotFwd?"_024":"_016"):"_1624"));
  

  TFile *f(0x0);
  if (ctauCut) f = new TFile(treeFileName(outputDir.c_str(),"DATA","","mass"));
  else f = new TFile(treeFileName(outputDir.c_str(),"DATA","","ctauMass"));
  if (!f || !f->IsOpen()) {
    if (ctauCut) {
      results2tree(outputDir.c_str(),"DATA","","mass");
      f = new TFile(treeFileName(outputDir.c_str(),"DATA","","mass"));
    }
    else {
      results2tree(outputDir.c_str(),"DATA","","ctauMass");
      f = new TFile(treeFileName(outputDir.c_str(),"DATA","","ctauMass"));
    }
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  //else cout << "[INFO] tree found!"<<endl;

  TFile* funf(0x0);
  TH1F* unfHist(0x0);
  TFile* funfSyst(0x0);
  TH1F* unfSyst(0x0);

  if (plotUnfolded)
    {
      funf = TFile::Open(Form("Output/unfData/results/unfResult_%s_%s_statErr.root",(doprompt)?"prompt":"nonprompt", (plotMid)?"mid":"fwd"));
      unfHist = (TH1F*) funf->Get("zUnf");
      funfSyst = TFile::Open(Form("Output/unfData/results/systErrs_%s_%s_Total.root",(doprompt)?"Prompt":"NonPrompt", (plotMid)?"Mid":"Fwd"));
      unfSyst = (TH1F*) funfSyst->Get("totalSyst");
    }
  map<anabin, njj_input> theVars_inputs;
  
  map<anabin, syst> syst_PP = readSyst_all("PP",poi.Data(),sTag.Data(),includeEffSyst,false);
  map<anabin, syst> syst_PbPb = readSyst_all("PbPb",poi.Data(),sTag.Data(),includeEffSyst,false);
  //map<anabin, syst> syst_taa_low = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_low.csv",sTag.Data()),excludeNonFitSyst);
  //map<anabin, syst> syst_taa_high = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_high.csv",sTag.Data()),excludeNonFitSyst);
  //map<anabin, syst> syst_Nmb = readSyst(Form("Systematics/csv/syst_%s_PbPb_Nmb.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> syst_lumipp = readSyst(Form("Systematics/csv/syst_%s_PP_lumi.csv",sTag.Data()),excludeNonFitSyst);
  map<anabin, syst> stat_PP; // for PP statistics
  map<anabin, syst> syst_glb_low; // for the boxes at 1
  map<anabin, syst> syst_glb_high; // for the boxes at 1
  
  vector<double> x, ex, y, ey;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float eff, acc, lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL,bfrac_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  if (ctauCut){
    tr->SetBranchAddress("N_Jpsi_val",&val);
    tr->SetBranchAddress("N_Jpsi_errL",&errL);
    tr->SetBranchAddress("N_Jpsi_errH",&errH);
  } else {
    tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
    tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  }
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    anabin thebin(zmin, zmax, ymin, ymax, ptmin, ptmax, centmin, centmax);
    bool ispp = (TString(collSystem)=="PP");
    if (ispp) {
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp_stat = errL;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_errL;
      theVars_inputs[thebin].systpp = syst_PP[thebin].value;
      syst thestat_PP; thestat_PP.name = "stat_PP"; thestat_PP.value = errL/val;
      stat_PP[thebin] = thestat_PP;
      theVars_inputs[thebin].statpp = thestat_PP;
      theVars_inputs[thebin].lumipp = lumi;
      theVars_inputs[thebin].correlpp = correl;
    } else {
      theVars_inputs[thebin].naa = val;
      theVars_inputs[thebin].dnaa_stat = errL;
      theVars_inputs[thebin].bfracaa = bfrac;
      theVars_inputs[thebin].dbfracaa = bfrac_errL;
      theVars_inputs[thebin].systaa = syst_PbPb[thebin].value;
      theVars_inputs[thebin].lumiaa = lumi;
      theVars_inputs[thebin].taa = taa;
      theVars_inputs[thebin].ncoll = ncoll;
      theVars_inputs[thebin].correlaa = correl;

    }
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned;
  map<anabin, vector<double> > theVarsBinned_stat;
  map<anabin, vector<double> > theVarsBinned_syst_low;
  map<anabin, vector<double> > theVarsBinned_syst_high;
  map<anabin, TGraphAsymmErrors* > theGraphs;
  map<anabin, TGraphAsymmErrors* > theGraphs_syst;
  
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned[*it] = vector<double>();
  }
  cout <<"before the Vars_inputs loop"<<endl;
  for (map<anabin, njj_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    cout<<"[INFO] initializing the input "<<endl;
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    njj_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    if (thebinOrig == binToSkip && underflowOff) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    njj_input spp = theVars_inputs[thebinPP];
    cout << "Result: " << endl;
    thebinOrig.print();
    cout <<  spp.npp << " +- " << spp.dnpp_stat << " ; " << spp.bfracpp << endl;
    //if (spp.npp <= 0) continue;
    //if ((doprompt || dononprompt) && spp.bfracpp<=0) continue;

    if (doPbPb && s.naa <= 0) continue;
    if (dopp && spp.npp <= 0) continue;

    theBins[thebin].push_back(it->first);

    double normfactorpp = 1., normfactoraa = 1.;

    normfactorpp = 1./spp.lumipp; cout <<"normfactorpp = 1/lumi = 1/" << spp.lumipp << " = "<<normfactorpp<<endl;
    if (useNcoll) {
      normfactoraa = 1./s.lumiaa;
      normfactoraa *= 1./(208.*208.*(HI::findNcollAverage(it->first.centbin().low(),it->first.centbin().high())/HI::findNcollAverage(0,200)));
    } else {
      double myNmb = NMB * s.lumiaa / lumipbpb_ABCD;
      normfactoraa = 1./(myNmb*s.taa*1e-3); // the 1e-3 factor is because taa is in mb-1 while lumis are in mub-1
    }
    normfactoraa *= 200./(it->first.centbin().high()-it->first.centbin().low());

    cout <<"normfactoraa = "<<normfactoraa<<" ; lumiaa = "<<s.lumiaa<< endl;
    double naa = s.naa;
    double npp = spp.npp;
    double dnaa = s.dnaa_stat;
    double dnpp = spp.dnpp_stat;
    //cout << "correl b-NJpsi = " << spp.bfracpp << " ; " << s.correlaa << endl;
    if (!ctauCut) {
      if (doprompt) {
	naa = s.naa*(1.-s.bfracaa);
	npp = spp.npp*(1.-spp.bfracpp);
	dnaa = naa*sqrt(pow(s.dnaa_stat/s.naa,2)
			- 2.*s.correlaa*s.dnaa_stat*s.dbfracaa/(s.naa*s.bfracaa)
			+ pow(s.dbfracaa/s.bfracaa,2));
	dnpp = npp*sqrt(pow(spp.dnpp_stat/spp.npp,2)
			- 2.*spp.correlpp*spp.dnpp_stat*spp.dbfracpp/(spp.npp*spp.bfracpp)
			+ pow(spp.dbfracpp/spp.bfracpp,2));
      }
      if (dononprompt) {
	naa = s.naa*s.bfracaa;
	npp = spp.npp*spp.bfracpp;
	dnaa = naa*sqrt(pow(s.dnaa_stat/s.naa,2)
			+ 2.*s.correlaa*s.dnaa_stat*s.dbfracaa/(s.naa*s.bfracaa)
			+ pow(s.dbfracaa/s.bfracaa,2));
	dnpp = npp*sqrt(pow(spp.dnpp_stat/spp.npp,2)
			+ 2.*spp.correlpp*spp.dnpp_stat*spp.dbfracpp/(spp.npp*spp.bfracpp)
			+ pow(spp.dbfracpp/spp.bfracpp,2));
      }
    }

    if (dopp && doPbPb) {
      naa *= normfactoraa;
      npp *= normfactorpp;
      dnaa *= normfactoraa;
      dnpp *= normfactorpp;
    }
    
    double njj;
    double dnjj;
    double syst_low;
    double syst_high; 
    
    if (doPbPb && dopp) {
      njj = npp>0 ? naa / npp : 0;
      dnjj = njj>0 ? njj*sqrt(pow(dnaa/naa,2) + pow(dnpp/npp,2)) : 0;
      syst_low = njj*sqrt(pow(spp.systpp,2)+pow(s.systaa,2));
      syst_high = syst_low;
    }
    else if (doPbPb) {
      njj = naa;
      dnjj = dnaa;
      syst_low = spp.systaa;
      syst_high = syst_low;
    }
    else {
      njj = npp;
      dnjj = dnpp;
      syst_low = spp.systpp;
      syst_high = syst_low;
    }

    thebinOrig.print();

    theVarsBinned[thebin].push_back(njj);
    theVarsBinned_stat[thebin].push_back(dnjj);
    theVarsBinned_syst_low[thebin].push_back(syst_low);
    theVarsBinned_syst_high[thebin].push_back(syst_high);
  }
  
  // systematics
  vector< map<anabin, syst> > all_glb_low;
  //all_glb_low.push_back(syst_PP);
  //all_glb_low.push_back(stat_PP);
  //all_glb_low.push_back(syst_lumipp);
  //all_glb.push_back(syst_Nmb);
  //syst_glb_low = combineSyst(all_glb_low,"global_low");
  //syst_glb_high = syst_glb_low;
  // vector< map<anabin, syst> > all_glb_low;
    //all_glb_low.push_back(syst_taa_low);
    //all_glb_low.push_back(syst_lumipp);
    //all_glb_low.push_back(syst_Nmb);
    //syst_glb_low = combineSyst(all_glb_low,"global");
    
    //vector< map<anabin, syst> > all_glb_high;
    //all_glb_high.push_back(syst_taa_high);
    //all_glb_high.push_back(syst_lumipp);
    //all_glb_high.push_back(syst_Nmb);
    //syst_glb_high = combineSyst(all_glb_high,"global");
  
  // make TGraphAsymmErrors
  int cnt=0;
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    int n = theBins[*it].size();
    if(n==0) {
    cout << "Error, nothing found for category" << endl;
    theGraphs[*it] = NULL;
    continue;
    }
    
    theGraphs[*it] = new TGraphAsymmErrors(n);
    theGraphs[*it]->SetName(Form("bin_%i",cnt));
    theGraphs_syst[*it] = new TGraphAsymmErrors(n);
    theGraphs_syst[*it]->SetName(Form("bin_%i_syst",cnt));
    
    double integ = 0;
    double intErr = 0;
    double intSyst = 0;

    if (plotUnfolded) {
      for (int i=0; i<n; i++) {
	anabin thebin = theBins[*it][i];
	double x = thebin.zbin().low() + thebin.zbin().high();
	integ = integ + unfHist->GetBinContent(unfHist->FindBin(x/2));
	intErr = intErr + pow(unfHist->GetBinError(unfHist->FindBin(x/2)),2);
	intSyst = intSyst + pow(unfSyst->GetBinError(unfSyst->FindBin(x/2)),2);
      }
      intErr = sqrt(intErr);
      intSyst = sqrt(intSyst);
      if (integ == 0) cout<<"[ERROR] the integral of the graph is 0."<<endl;
    }
    for (int i=0; i<n; i++) {
      double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
      double exsyst=0, eysyst_low=0,eysyst_high=0;
      double low=0, high=0;
      anabin thebin = theBins[*it][i];
      if (xaxis=="z") {
	low = thebin.zbin().low();
	high = thebin.zbin().high();
      }
      y = theVarsBinned[*it][i];
      x = (low+high)/2.;
      if (plotUnfolded) {y = unfHist->GetBinContent(unfHist->FindBin(x)); y = y*(1.0/(high-low))*(1.0/integ);}
      exh = (high-low)/2.;
      exl = (high-low)/2.;
      exsyst = (xaxis=="pt") ? 0.5 : 0.05;
      eysyst_low = theVarsBinned_syst_low[*it][i];
      eysyst_high = theVarsBinned_syst_high[*it][i];
      eyl = fabs(theVarsBinned_stat[*it][i]);
      if (plotUnfolded) {eyl = unfHist->GetBinError(unfHist->FindBin(x)); eyl = eyl*1.0/integ;}
      eyh = eyl;      
      eysyst_low = y*eysyst_low;
      eysyst_high = y*eysyst_high;
      if (plotUnfolded) {
	eysyst_low = unfSyst->GetBinError(unfSyst->FindBin(x))/unfSyst->GetBinContent(unfSyst->FindBin(x));
	eysyst_low = y*eysyst_low;//*1.0/integ;
	eysyst_high = eysyst_low;
      }
      theGraphs[*it]->SetPoint(i,x,y);
      theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
      theGraphs_syst[*it]->SetPoint(i,x,y);
      theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst_low,eysyst_high);
      cout << "final : x = " << x << " ,y = " << y << " ,eyl = " << eyl << " ,eyh = " << eyh << " ,syst = " << eysyst_low << endl;
      if (y > histMax) histMax = y;
      // theGraphs[*it]->Sort();
      // theGraphs_syst[*it]->Sort();
    }
    cnt++;
  }
  
  // plot
  plotGraphNJJ(theGraphs, theGraphs_syst, xaxis, outputDir, syst_glb_low, syst_glb_high);
}

void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs, map<anabin, TGraphAsymmErrors*> theGraphs_syst, string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high) {
  setTDRStyle();
  
  const char* ylabel = "N(J/#psi)";
  if (plotUnfolded) ylabel = "1/N dN/dz";  
  int intervals2Plot = theGraphs.size();
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  haxes = new TH1F("haxes","haxes", 5, 0, 1);
  haxes->GetXaxis()->SetNdivisions(505);
  //line = TLine(0,1,1,1);
  if (plotUnfolded) haxes->GetYaxis()->SetRangeUser(0, 8);
  else
    haxes->GetYaxis()->SetRangeUser(0, histMax*2);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = "z";//(xaxis=="pt") ? "p_{T} (GeV)" : ((xaxis=="rap") ? "|y|" : "N_{part}");
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->GetYaxis()->CenterTitle(true);
  if (plotUnfolded)
    haxes->GetYaxis()->SetTitleOffset(0.8);
  haxes->Draw();
//  line.Draw();

  TFile* mcFile = TFile::Open(Form("Output/MCResults/mcResult_%s_underflowOff.root", doprompt?"prompt":"nonprompt"));
  TFile* bFile = TFile::Open("Output/MCResults/mcResult_nonprompt_underflowOff_bHadronPt.root");
  if (!mcFile) {
    cout<<"[INFO] MC file not found"<<endl;
  }
  TH1F* mcHist(0x0);
  if (mcON)
    mcHist =  (TH1F*) mcFile->Get(Form("zDist_%s",plotMid?"mid":"fwd"));
  TH1F* bHist(0x0); 
  if (bffON)
    bHist =  (TH1F*) bFile->Get(Form("zDist_%s",plotMid?"mid":"fwd"));

  //double integ=0;
  double intErr=0;
  if (mcON)
    {
      //integ = mcHist->Integral();
      //for (int i=0; i<mcHist->GetSize(); i++)
      //intErr = intErr + pow(mcHist->GetBinError(i),2);
      //intErr=sqrt(intErr);
      
      mcHist->Scale(1.0/mcHist->Integral("width")); 
      bHist->Scale(1.0/bHist->Integral("width"));
    }

  if (plotUnfolded && mcON)
    {
      //gStyle->SetEndErrorSize(5);
      mcHist->SetMarkerColor(color(3));
      mcHist->SetLineColor(color(3));
      mcHist->SetLineWidth(2);
      mcHist->SetMarkerStyle(kOpenSquare);
      mcHist->SetMarkerColor(color(3));
      mcHist->SetMarkerSize(1.5);
      //mcHist->Draw("same");
    }
  if (plotUnfolded && bffON)
    {
      //mcHist->SetMarkerColor(color(8));
      //mcHist->SetLineColor(color(8));
      bHist->SetMarkerColor(color(5));
      bHist->SetLineColor(color(5));
      bHist->SetLineWidth(3);
      bHist->SetLineStyle(2);
    }
  double xshift=0.025;
  TLegend *tleg(0x0);
  tleg = new TLegend(0.70,0.62,0.98,0.72);
  if (bffON) {
    tleg = new TLegend(0.70,0.57,0.98,0.76);
    if (plotMid) tleg = new TLegend(0.20,0.35,0.5,0.55);
  }
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  tleg->SetTextFont(42);
  tleg->SetTextSize(0.04);
  if (bffON) tleg->SetTextSize(0.04);
  bool drawLegend(true);
  
  // prepare for the printing of the result tables
  const char* xname = (xaxis=="cent") ? "Centrality" : (xaxis=="pt" ? "\\pt" : (xaxis=="z" ? "\\z" : "$|y|$"));
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_JPSI_NJJ_%s%s.tex",outputDir.c_str() ,xaxis.c_str(), nameTag.c_str());
  string yname("$\\N (\\Jpsi-Jet)$");
  inittex(texname, xname, yname);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg || !tg_syst) continue;
    
    theCats.push_back(thebin);
    
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;
    //if (plotUnfolded) colorF = 17;
    if (intervals2Plot==2)
    {
      style = 4;
      cnt==0 ? colorI = 9 : colorI =4;
      colorF = color(colorI)-11;
    }
    if (intervals2Plot==3)
    {
      style = 0;
      colorI = cnt+6;
      colorF = color(colorI)-10;
    }
    if (intervals2Plot>3)
    {
      style = cnt+1;
      colorI = cnt+1;
      colorF = color(colorI)-11;
    }
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    tg_syst->SetLineColor(color(colorI));
    tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    bool plot = true;
    
    // in the case where the centrality dependence is plotted: treat the PP uncertainties as global systematics
    // if (xaxis == "cent") {
    if (plot)
    {
      double x(0.), dx(0.), y(0.), dy_low(0.), dy_high(0.);
      double rightA = 0.;
      int centminGlob(0),centmaxGlob(200);
      if (xaxis=="cent") {
        dx = 10;
        rightA = 420.;
      } else if (xaxis=="pt") {
        dx = 0.625;
        if (intervals2Plot == 1)
        {
          rightA = 50.;
          dx = 1.25;
        }
        else
        {
          rightA = 30.;
          dx = 0.65;
        }
        if (intervals2Plot == 3)
        {
          centminGlob = it->first.centbin().low();
          centmaxGlob = it->first.centbin().high();
        }
      } else if (xaxis=="rap") {
        dx = 0.06;
        rightA = 2.4;
      }
      else if (xaxis=="z") {
	dx = 0.06;
	rightA = 1;
      }
      x = rightA - (2*dx*cnt + dx);
      //if (plotFwdMid && cnt>0) x = rightA - (2*dx + dx);
      y = 1;
      
      anabin thebinglb(it->first.zbin().low(),
		       it->first.zbin().high(),
		       it->first.rapbin().low(),
                       it->first.rapbin().high(),
                       it->first.ptbin().low(),
                       it->first.ptbin().high(),
                       centminGlob,centmaxGlob);
      thebinglb.print();
      dy_low = gsyst_low[thebinglb].value;
      dy_high = gsyst_high[thebinglb].value;
      cout << "global syst: " << "+" << dy_high << " -" << dy_low << endl;
      TBox *tbox = new TBox(x-dx,y-dy_low,x+dx,y+dy_high);
      tbox->SetFillColorAlpha(colorF, 1);
      tbox->SetLineColor(color(colorI));
      //tbox->Draw("l");
      TBox *tboxl = (TBox*) tbox->Clone("tboxl");
      tboxl->SetFillStyle(0);
      //tboxl->Draw("l");
    }
    // }
    
    // Plot graphs after uncertainties to avoid overlap

    if (plot)
    {
      tg_syst->Draw("5");
      gStyle->SetEndErrorSize(5);
      tg->Draw("P");
      // tg->Draw("[]");
    }
  
    if (plotUnfolded && mcON){
      mcHist->Draw("E1 same");
    }
    if (plotUnfolded && bffON)
      bHist->Draw("same hist");
    TLatex tl;
    double tlx = 0.25; //0.92;
    double tly = 0.80; //0.69;
    tl.SetTextFont(42); // consistent font for symbol and plain text

    if (plot)
    {
      TString zlabel = Form("%.1f < z < %.1f",it->first.zbin().low(),it->first.zbin().high());
      TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
      if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
      TString ptlabel = Form("%g < p_{T} < %g GeV",it->first.ptbin().low(), it->first.ptbin().high());
      TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
      
      if (mcON)
	{
	  tleg->AddEntry(tg, "Data", "lp");
	  if (bffON){
	    tleg->AddEntry(mcHist, "PYTHIA 8", "lp");
	    tleg->AddEntry(bHist, "#splitline{PYTHIA 8}{b hadron}", "lp");
	  }
	  else 
	    tleg->AddEntry(mcHist, "PYTHIA 8", "lp");
	}
      if (xaxis == "pt")
        {
	  if (intervals2Plot > 3) tleg->AddEntry(tg, raplabel, "p");
	  else if (intervals2Plot == 3) tleg->AddEntry(tg, Form("Cent. %s",centlabel.Data()), "p");
	  else if (intervals2Plot == 1) drawLegend = false;
	  else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",raplabel.Data(),centlabel.Data()));
	  //else tleg->AddEntry(tg, Form("%s, Cent. %s",raplabel.Data(),centlabel.Data()), "p");
        }
      if (xaxis == "cent")
        {
	  if (intervals2Plot > 3 ) tleg->AddEntry(tg, raplabel, "p");
	  else if (intervals2Plot == 2) tleg->AddEntry(tg, ptlabel, "p");
	  else if (intervals2Plot == 1) drawLegend = false;
	  else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()));
          }
      if (xaxis == "rap")
        {
          if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
          else if (intervals2Plot == 1) drawLegend = false;
          else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
        }
      if (xaxis == "z" && !mcON) drawLegend = false; 
    }
    
    // print tex
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
    //if (xaxis == "cent" || xaxis == "rap") 
    oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";
    
    addline(texname,oss.str());
    printGraph(tg, tg_syst, texname);
    
    // for the centrality dependence: we want Npart plotted, not the centrality
    if (xaxis == "cent") {
      centrality2npart(tg, false, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
      centrality2npart(tg_syst, true, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
    }
    
    cnt++;
    it_syst++;
  }
  
  
  if (drawLegend) tleg->Draw();
  if(xaxis=="pt" && intervals2Plot == 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="rap" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }

  if (xaxis=="z" && plotMid)
    {
      if (plotMid && plotFwd)
	{
	  TLatex *tex = new TLatex(0.2,0.78,"|y_{J/#psi}| < 2.4");
	  tex->SetNDC();
	  tex->SetTextSize(0.044);
	  tex->SetTextFont(42);
	  tex->SetLineWidth(2);
	  tex->Draw();
	  
	  TLatex *tex0 = new TLatex(0.2,0.66,"|y_{jet}| < 2.4");
	  tex0->SetNDC();
	  tex0->SetTextSize(0.044);
	  tex0->SetTextFont(42);
	  tex0->SetLineWidth(2);
	  tex0->Draw();
	  
	  TLatex *tex1 = new TLatex(0.2,0.72,"6.5 < p_{T,J/#psi} < 100 GeV");
	  tex1->SetNDC();
	  tex1->SetTextSize(0.044);
	  tex1->SetTextFont(42);
	  tex1->SetLineWidth(2);
	  tex1->Draw();
	  
	  TLatex *tex2 = new TLatex(0.2,0.60,Form("%d < p_{T,jet} < %d GeV",jtPtRange == 0? 30:(jtPtRange == -1? 20:40), jtPtRange == 0? 40:(jtPtRange == -1? 30:50)));
	  tex2->SetNDC();
	  tex2->SetTextSize(0.044);
	  tex2->SetTextFont(42);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	}
      else if (plotMid && !plotFwd) 
	{
	  TLatex *tex = new TLatex(0.2,0.78,"|y_{J/#psi}| < 1.6");
	  tex->SetNDC();
	  tex->SetTextSize(0.044);
	  tex->SetTextFont(42);
	  tex->SetLineWidth(2);
	  tex->Draw();
	  
	  TLatex *tex0 = new TLatex(0.2,0.66,"|y_{jet}| < 2.4");
	  tex0->SetNDC();
	  tex0->SetTextSize(0.044);
	  tex0->SetTextFont(42);
	  tex0->SetLineWidth(2);
	  tex0->Draw();
	  
	  TLatex *tex1 = new TLatex(0.2,0.72,"6.5 < p_{T,J/#psi} < 100 GeV");
	  tex1->SetNDC();
	  tex1->SetTextSize(0.044);
	  tex1->SetTextFont(42);
	  tex1->SetLineWidth(2);
	  tex1->Draw();
	  
	  TLatex *tex2 = new TLatex(0.2,0.60,Form("%d < p_{T,jet} < %d GeV",jtPtRange == 0? 30:(jtPtRange == -1? 20:40), jtPtRange == 0? 40:(jtPtRange == -1? 30:50)));
	  tex2->SetNDC();
	  tex2->SetTextSize(0.044);
	  tex2->SetTextFont(42);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	}
      else if (plotFwd)
	{
	  TLatex *tex = new TLatex(0.2,0.78,"1.6 < |y_{J/#psi}| < 2.4");
	  tex->SetNDC();
	  tex->SetTextSize(0.044);
	  tex->SetTextFont(42);
	  tex->SetLineWidth(2);
	  tex->Draw();
	  
	  TLatex *tex0 = new TLatex(0.2,0.66,"|y_{jet}| < 2.4");
	  tex0->SetNDC();
	  tex0->SetTextSize(0.044);
	  tex0->SetTextFont(42);
	  tex0->SetLineWidth(2);
	  tex0->Draw();
	  
	  TLatex *tex1 = new TLatex(0.2,0.72,"3 < p_{T,J/#psi} < 100 GeV");
	  tex1->SetNDC();
	  tex1->SetTextSize(0.044);
	  tex1->SetTextFont(42);
	  tex1->SetLineWidth(2);
	  tex1->Draw();
	  
	  TLatex *tex2 = new TLatex(0.2,0.60,Form("%d < p_{T,jet} < %d GeV",jtPtRange == 0? 30:(jtPtRange == -1? 20:40), jtPtRange == 0? 40:(jtPtRange == -1? 30:50)));
	  tex2->SetNDC();
	  tex2->SetTextSize(0.044);
	  tex2->SetTextFont(42);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	}
    }
  
  TLatex tl;
  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;
  tl.SetTextFont(62); // consistent font for symbol and plain text
  tl.SetTextSize(0.044); 
  if (doprompt) tl.DrawLatexNDC(tlx,tly,"Prompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  tl.SetTextSize(0.046);
  
  int iPos = 33;
  //if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  //else CMS_lumi( (TPad*) gPad, 107, iPos, "", isPreliminary );
  if ( dopp && doPbPb ) CMS_lumi( (TPad*) gPad, 111, iPos, "", isPreliminary );
  else if ( dopp ) CMS_lumi( (TPad*) gPad, 109, iPos, "", isPreliminary );
  else if ( doPbPb ) CMS_lumi( (TPad*) gPad, 110, iPos, "", isPreliminary );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/RESULT/plot/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/RESULT/plot/root/result_JPSI_NJJ_%s%s.root",outputDir.c_str(), xaxis.c_str(), nameTag.c_str()));
  gSystem->mkdir(Form("Output/%s/RESULT/plot/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/RESULT/plot/png/result_JPSI_NJJ_%s%s.png",outputDir.c_str(), xaxis.c_str(), nameTag.c_str()));
  gSystem->mkdir(Form("Output/%s/RESULT/plot/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/RESULT/plot/pdf/result_JPSI_NJJ_%s%s.pdf",outputDir.c_str(), xaxis.c_str(), nameTag.c_str()));
  gSystem->mkdir(Form("Output/%s/RESULT/plot/cMacro/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/RESULT/plot/cMacro/result_JPSI_NJJ_%s%s.C",outputDir.c_str(), xaxis.c_str(), nameTag.c_str()));
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}


void centrality2npart(TGraphAsymmErrors* tg, bool issyst, double xshift) {
  int n = tg->GetN();
  for (int i=0; i<n; i++) {
    double x, y, exl, exh, exlC, exhC, eyl, eyh;
    x = tg->GetX()[i];
    if (x>0)
    {
      y = tg->GetY()[i];
      exlC = tg->GetErrorXlow(i);
      exhC = tg->GetErrorXhigh(i);
      eyl = tg->GetErrorYlow(i);
      eyh = tg->GetErrorYhigh(i);
      if (!issyst) {
        exl = HI::findNpartSyst_low(2.*(x-exlC),2.*(x+exhC));//0.;
        exh = HI::findNpartSyst_high(2.*(x-exlC),2.*(x+exhC));//0.;
      } else {
        exl = 5;
        exh = exl;
      }
      x = HI::findNpartAverage(2.*(x-exlC),2.*(x+exhC));
      tg->SetPoint(i,x+xshift,y);
      tg->SetPointError(i,exl,exh,eyl,eyh);
    }
    else
    {
      tg->SetPoint(i,-1000,0);
      tg->SetPointError(i,0,0,0,0);
    }
  }
}

int color(int i) {
  if (i==0) return kMagenta+2;
  else if (i==1) return kBlue+2;
  else if (i==2) return kRed+2;
  else if (i==3) return kCyan+2;
  else if (i==4) return kGreen+2;
  else if (i==5) return kOrange+2;
  else if (i==6) return kRed+1;
  else if (i==7) return kYellow+1;
  else if (i==8) return kAzure+1;
  else if (i==9) return kBlack;
  else return kBlack;
}

int markerstyle(int i) {
  if (i==0) return kFullSquare;
  else if (i==1) return kFullCircle;
  else if (i==2) return kFullStar;
  else if (i==3) return kFullCross;
  else if (i==4) return kFullDiamond;
  else if (i==5) return kOpenSquare;
  else if (i==6) return kOpenCircle;
  else if (i==7) return kOpenStar;
  else if (i==8) return kFullTriangleDown;
  else return kOpenCross;
}

void setOptions(bool adoPbPb, bool adopp, bool adoprompt, bool adononprompt, bool aplotMid, bool aplotFwd, bool aexcludeNonFitSyst, string anameTag_base, int ajtPtRange, bool aplotUnfolded, bool aunderflowOff, bool amcON, bool abffON, bool actauCut) {
  doPbPb = adoPbPb;
  dopp = adopp;
  doprompt = adoprompt;
  dononprompt = adononprompt;
  plotMid = aplotMid;
  plotFwd = aplotFwd;
  excludeNonFitSyst = aexcludeNonFitSyst;
  jtPtRange = ajtPtRange;
  plotUnfolded = aplotUnfolded;
  underflowOff = aunderflowOff;
  nameTag = anameTag_base;
  mcON = amcON;
  bffON = abffON;  
  ctauCut = actauCut;

  if (doPbPb && dopp) nameTag += "_Raa";
  else if (doPbPb) nameTag += "_PbPb";
  else nameTag += "_pp";
  if (doprompt) nameTag += "_prompt";
  if (dononprompt) nameTag += "_nonprompt";
  nameTag+= (plotMid)?(plotFwd?"_024":"016"):"_1624";
  if (jtPtRange == 0) nameTag += "_midJtPt";
  else if (jtPtRange == 1) nameTag += "_highJtPt";
  else if (jtPtRange == -1) nameTag += "_lowJtPt";
  nameTag+= (plotUnfolded)?"_unfolded":"";
  nameTag+= (underflowOff)?"_underflowOff":"";
  nameTag+= (mcON)?"_mcON":"";
  nameTag+= (bffON)?"_bffON":"";
  nameTag+= (ctauCut)?"_ctauCut":"";
  histMax = 0;
}

void printOptions() {
  cout <<
    "doPbPb = "<< doPbPb << ", " <<
    "dopp = "<< dopp << ", " <<
    "doprompt = " << doprompt << ", " <<
    "dononprompt = " << dononprompt << ", " <<
    "plotMid = " << plotMid << ", " <<
    "plotFwd = " << plotFwd << ", " <<
    "jtPtRange = " << jtPtRange << ", " <<
    "unfolded = " << plotUnfolded << ", " <<
    "underflowOff = " << underflowOff << ", " <<
    "mcON = " << mcON << ", " <<
    "bffON = " << bffON << ", " <<
    "ctauCut = " << ctauCut << ", " <<
    "excludeNonFitSyst = " << excludeNonFitSyst << ", " <<
    "nameTag_base = \"" << nameTag_base << "\"" <<
    endl;
}

map<anabin, njj_input > readResults(const char* resultsFile)
{
  cout <<"reading results from file "<<resultsFile<<endl;
  map<anabin, njj_input> ans;
  njj_input theresult;
  
  ifstream file(resultsFile);
  if (!(file.good())) return ans;
  
  string resultname; getline(file,resultname);
  
  cout << "[INFO] Reading results : " << resultname.c_str() << endl;
  
  string line;
  double zmin=0, zmax=0, rapmin=0, rapmax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;
  
  while (file.good()) {
    getline(file,line);
    if (line.size()==0) break;
    TString tline(line.c_str());
    TString t; Int_t from = 0, cnt=0;
    while (tline.Tokenize(t, from , ",")) {
      t.Strip(TString::kBoth,' ');
      value = atof(t.Data());
      if (cnt==0) rapmin = atof(t.Data());
      else if (cnt==1) rapmax = value;
      else if (cnt==2) ptmin = value;
      else if (cnt==3) ptmax = value;
      else if (cnt==4) zmin = value;
      else if (cnt==5) zmax = value;
      else if (cnt==6) centmin = value;
      else if (cnt==7) centmax = value;
      else if (cnt==8) theresult.npp = value;
      else if (cnt==9) theresult.dnpp_stat = value;
      else if (cnt==10) theresult.systpp = value;
      else if (cnt > 10) {
        cout << "Warning, too many fields, I'll take the last one." << endl;
        continue;
      }
      cnt++;
    }
    anabin thebin(zmin, zmax, rapmin, rapmax, ptmin, ptmax, centmin, centmax);
    ans[thebin] = theresult;
  }
  
  file.close();
  
  return ans;
}

void drawArrow(double x, double ylow, double yhigh, double dx, Color_t color) {
  TArrow *arrow = new TArrow(x,yhigh,x,ylow<=0. ? 0.01 : ylow,0.03,ylow<=0. ? ">" : "<>");
  arrow->SetLineColor(color);
  arrow->Draw();
  if (ylow<=0.) {
    TLine *line = new TLine(x-dx,yhigh,x+dx,yhigh);
    line->SetLineColor(color);
    line->Draw();
  }
}
