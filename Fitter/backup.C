#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/monster.h"
#include "Systematics/syst.h"

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
//const char* poiname       = "N_Jpsi"; // for NJJ (will correct automatically for efficiency)
//#endif
//const char* ylabel        = "N_{JJ}";

bool  doprompt      = false;  // prompt Jpsi
bool  dononprompt   = true;  // nonprompt Jpsi
bool  is18XXX       = true; //plot results in JPsiJet analysis
bool  plot14005     = false;
bool  plotPsi2S     = false; // plot Psi2S
bool  applyEff      = false;
bool  applyAcc      = false;
bool  doLogPt       = false;
bool  plotFwdMid    = false;
bool  isPreliminary = false;
string nameTag_base = "_prompt";    // can put here e.g. "_prompt", "_nonprompt", ...

const bool useNcoll = false; // false -> use TAA / NMB, true -> use Ncoll / lumiPbPb

//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void setOptions(bool adoprompt, bool adononprompt, bool aplotFwdMid, bool ais18XXX, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, string anameTag_base="");
void plotNJJ(vector<anabin> thecats, string xaxis, string outputDir);
void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs,/* map<anabin, TGraphAsymmErrors*> theGraphs_syst,*/ string xaxis, string outputDir);
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
  double effpp;
  double effaa;
  double accpp;
  double accaa;
  double effppP; // Only used for the b fraction correction
  double effaaP; // Only used for the b fraction correction
  double accppP; // Only used for the b fraction correction
  double accaaP; // Only used for the b fraction correction
  double effppNP; // Only used for the b fraction correction
  double effaaNP; // Only used for the b fraction correction
  double accppNP; // Only used for the b fraction correction
  double accaaNP; // Only used for the b fraction correction
  double systeffppP; // Only used for the b fraction correction
  double systeffaaP; // Only used for the b fraction correction
  double systeffppNP; // Only used for the b fraction correction
  double systeffaaNP; // Only used for the b fraction correction
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


void plotZed(string workDirName, string poiname, int iplot) {
  // poiname is the observable to be plotted = "NJJ", "BF" and "XS" (For NJJ, b-fraction and cross-section)
  
  string xaxis = "zed";
  vector<anabin> theCats;  

  // 10 z bins in 3 rap intervals 
  if (iplot==0) {
    theCats.push_back(anabin(0.0,1.0,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.0,0.1,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.1,0.2,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.2,0.3,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.3,0.4,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.4,0.5,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.5,0.6,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.6,0.7,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.7,0.8,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.8,0.9,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.9,1.0,0,2.4,6.5,50,0,200));
  }

  if (iplot==1) {
    theCats.push_back(anabin(0.0,1.0,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.0,0.1,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.1,0.2,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.2,0.3,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.3,0.4,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.4,0.5,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.5,0.6,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.6,0.7,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.7,0.8,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.8,0.9,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0.9,1.0,0,1.6,6.5,50,0,200));
  }

  if (iplot==2) {
    theCats.push_back(anabin(0.0,1.0,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.0,0.1,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.1,0.2,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.2,0.3,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.3,0.4,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.4,0.5,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.5,0.6,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.6,0.7,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.7,0.8,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.8,0.9,1.6,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0.9,1.0,1.6,2.4,6.5,50,0,200));
  }
  
  
  nameTag = nameTag_base + Form("_%i",iplot);
  
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
  
};

void plotPt(string workDirName, string poiname, int iplot) {
  // poiname is the observable to be plotted = "NJJ", "BF" and "XS" (For NJJ, b-fraction and cross-section)
  
  string xaxis = "pt";
  vector<anabin> theCats;  

  // 3 rap intervals integrated in z
  if (iplot==0) {
    theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0,1,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0,1,1.6,2.4,3,50,0,200));
  }
  
  
  nameTag = nameTag_base + Form("_%i",iplot);
  
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
  
};

void plotCent(string workDirName, string poiname, int iplot) {
  // poiname is the observable to be plotted = "NJJ", "BF" and "XS" (For NJJ, b-fraction and cross-section)
  
  string xaxis = "cent";
  vector<anabin> theCats;
  
  // 3 rapidity intervals intwgrated in z
  if (iplot==0) {
    theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0,1,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0,1,1.6,2.4,3,50,0,200));
  }
   
  nameTag = nameTag_base + Form("_%i",iplot);
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};

void plotRap(string workDirName, string poiname) {
  // poiname is the observable to be plotted = "NJJ", "BF" and "XS" (For NJJ, b-fraction and cross-section)
  
  string xaxis = "rap";
  vector<anabin> theCats;
  
  theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
  
  nameTag = nameTag_base;
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};

void plotAll(string workDirName, string poiname) {
  // poiname is the observable to be plotted = "NJJ", "BF" and "XS" (For NJJ, b-fraction and cross-section)
//  if (dononprompt) nameTag_base = "_nonprompt";
  if (!doprompt && !dononprompt) nameTag_base = "";
  
  if (is18XXX)
  {
    plotZed(workDirName,poiname,0);
    plotPt(workDirName,poiname,0);
  }
};

void doAllplots(bool is18XXX=false) {
  
  if (is18XXX)
  {
    setOptions(true,false,true,true, false, false, false, false, false);
    printOptions();
    plotAll("DataFits_18XXX_2D_2CB_polBkg_nominal","NJJ");
  }

};


/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plotNJJ(vector<anabin> thecats, string xaxis, string outputDir){
  // thecats contains the categories. eg 0<y<1.6 and 1.6<y<2.4
  // xaxis is the variable to be plotted. "pt", "rap" or "cent"
  // outputDir is the directory to save the plots
  outputDir= outputDir +"/mass/DATA";  
  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  
  TString poi("NJpsi");
  if (doprompt) poi = "NJpsi_prompt";
  if (dononprompt) poi = "NJpsi_nonprompt";
  
  TFile *f = new TFile(treeFileName(outputDir.c_str(),""));
  if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),"");
    f = new TFile(treeFileName(outputDir.c_str(),""));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  
  TString sTag("16025");
  if (is18XXX) sTag = "18XXX";
  
  map<anabin, njj_input> theVars_inputs;

  //map<anabin, syst> syst_PP = readSyst_all("PP",poi.Data(),sTag.Data(),false,false);
  //map<anabin, syst> syst_PbPb = readSyst_all("PbPb",poi.Data(),sTag.Data(),false,false);
  //map<anabin, syst> syst_taa_low = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_low.csv",sTag.Data()),true);
  //map<anabin, syst> syst_taa_high = readSyst(Form("Systematics/csv/syst_%s_PbPb_taa_high.csv",sTag.Data()),true);
  //map<anabin, syst> syst_Nmb = readSyst(Form("Systematics/csv/syst_%s_PbPb_Nmb.csv",sTag.Data()),true);
  //map<anabin, syst> syst_lumipp = readSyst(Form("Systematics/csv/syst_%s_PP_lumi.csv",sTag.Data()),true);
  //map<anabin, syst> stat_PP; // for PP statistics
  //map<anabin, syst> syst_glb_low; // for the boxes at 1
  //map<anabin, syst> syst_glb_high; // for the boxes at 1
    
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
  tr->SetBranchAddress("N_Jpsi_val",&val);
  tr->SetBranchAddress("N_Jpsi_errL",&errL);
  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  if (!dononprompt)
  {
    tr->SetBranchAddress("eff_val",&eff);
    tr->SetBranchAddress("acc_val",&acc);
  }
  else
  {
    tr->SetBranchAddress("effnp_val",&eff);
    tr->SetBranchAddress("accnp_val",&acc);
  }
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
    
    if (xaxis=="rap" && ((ymin==0 && ymax<=0.61 && ymax>=0.59 ) || (ymin>=0.59 && ymin<=0.61 && ymax>=1.19 && ymax <=1.21) || (ymin>=1.19 && ymin<=1.21 && ymax>=1.79 && ymax <=1.81) || (ymin>=1.79 && ymin<=1.81 && ymax>=2.39 && ymax <=2.41))) continue;
    
    anabin thebin(zmin, zmax, ymin, ymax, ptmin, ptmax, centmin, centmax);
    
    bool ispp = (TString(collSystem)=="PP");
    
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp_stat = errL;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_errL;
      //theVars_inputs[thebin].systpp = syst_PP[thebin].value;
      //syst thestat_PP; thestat_PP.name = "stat_PP"; thestat_PP.value = errL/val;
      //stat_PP[thebin] = thestat_PP;
      //theVars_inputs[thebin].statpp = thestat_PP;
      theVars_inputs[thebin].lumipp = lumi;
      theVars_inputs[thebin].effpp = eff;
      theVars_inputs[thebin].accpp = acc;
      theVars_inputs[thebin].correlpp = correl;
    
  }
  
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned;
  map<anabin, vector<double> > theVarsBinned_stat;
  map<anabin, vector<double> > theVarsBinned_syst_low;
  map<anabin, vector<double> > theVarsBinned_syst_high;
  map<anabin, TGraphAsymmErrors* > theGraphs;
  //map<anabin, TGraphAsymmErrors* > theGraphs_syst;
  
  map<anabin, njj_input > theResults18XXX; //In case they are needed to compute the Psi2S NJJ
  if (is18XXX && plotPsi2S) theResults18XXX = readResults("njj_18XXX.csv");
  
  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned[*it] = vector<double>();
  }
  
  for (map<anabin, njj_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    njj_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    njj_input spp = theVars_inputs[thebinPP];
           
    double naa = s.naa;
    double npp = spp.npp;
    double dnaa = s.dnaa_stat;
    double dnpp = spp.dnpp_stat;
    
    
    double njj = npp;
    double dnjj = njj>0 ? dnpp:0;//njj*sqrt(pow(dnaa/naa,2) + pow(dnpp/npp,2)) : 0;
    //double syst_low = njj*sqrt(pow(spp.systpp,2)+pow(s.systaa,2));
    //double syst_high = syst_low;
    
    theVarsBinned[thebin].push_back(njj);
    theVarsBinned_stat[thebin].push_back(dnjj);
  }
    
  // make TGraphAsymmErrors
  //int cnt=0;
  //for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
  //int n = theBins[*it].size();
  //if(n==0) {
  //  cout << "Error, nothing found for category" << endl;
  //  theGraphs[*it] = NULL;
  //  continue;
  //}
  //
  // theGraphs[*it] = new TGraphAsymmErrors(n);
  //theGraphs[*it]->SetName(Form("bin_%i",cnt));
  //theGraphs_syst[*it] = new TGraphAsymmErrors(n);
  //theGraphs_syst[*it]->SetName(Form("bin_%i_syst",cnt));
  //
  //for (int i=0; i<n; i++) {
  //  double x=0, exl=0, exh=0, y=0, eyl=0, eyh=0;
  //  double exsyst=0, eysyst_low=0,eysyst_high=0;
  //  double low=0, high=0;
  //  anabin thebin = theBins[*it][i];
  //  y = theVarsBinned[*it][i];
  //  if (xaxis=="pt" || xaxis=="rap" || xaxis=="zed") {
  //    if (xaxis=="pt") {
  //      low= thebin.ptbin().low();
  //      high = thebin.ptbin().high();
  //    } else if (xaxis=="rap"){
  //      low= thebin.rapbin().low();
  //      high = thebin.rapbin().high();
  //    } else {
  //      low= thebin.zbin().low();
  //      high = thebin.zbin().high();
  //    }
  //    x = (low+high)/2.;
  //    exh = (high-low)/2.;
  //    exl = (high-low)/2.;
  //    exsyst = (xaxis=="pt") ? 0.5 : 0.05;
  //    eysyst_low = theVarsBinned_syst_low[*it][i];
  //    eysyst_high = theVarsBinned_syst_high[*it][i];
  //  }
  //  if (xaxis=="cent") {
  //    low= thebin.centbin().low();
  //    high = thebin.centbin().high();
  //    // exl = 0.;
    //    // exh = 0.;
  //  x = (low+high)/2./2.;
  //    exh = (high-low)/2./2.;
  //    exl = (high-low)/2./2.;
  //    exsyst = exl;
  //    eysyst_low = theVarsBinned_syst_low[*it][i];
  //    eysyst_high = theVarsBinned_syst_high[*it][i];
  //  }
  //  eyl = fabs(theVarsBinned_stat[*it][i]);
  //  eyh = eyl;
      
      // eysyst = y*eysyst;
      
  //  theGraphs[*it]->SetPoint(i,x,y);
  //  theGraphs[*it]->SetPointError(i,exl,exh,eyl,eyh);
  //  theGraphs_syst[*it]->SetPoint(i,x,y);
  //  theGraphs_syst[*it]->SetPointError(i,exsyst,exsyst,eysyst_low,eysyst_high);
//       cout << "final = " << x << " " << y << " " << eyl << " " << eyh << " " << eysyst << endl;
      
      // theGraphs[*it]->Sort();
      // theGraphs_syst[*it]->Sort();
  //}
  //cnt++;
  //}
  
  // plot
  plotGraphNJJ(theGraphs, xaxis, outputDir);
}

void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs,/* map<anabin, TGraphAsymmErrors*> theGraphs_syst,*/ string xaxis, string outputDir, map<anabin, syst> gsyst_low, map<anabin, syst> gsyst_high) {
  setTDRStyle();
  
  const char* ylabel = "N_{J/#psi-Jet}";
  
  int intervals2Plot = theGraphs.size();
  
  //if (intervals2Plot>1 && plot14005) return;
  
  if (plotFwdMid && intervals2Plot<=3) return;
  
  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  
  // in the case of the centrality dependence, we need the minimum bias panel on the right
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="zed"){
    haxes = new TH1F("haxes","haxes",1,0,1);
    //haxes->GetXaxis()->SetNdivisions(306,false);
    line = TLine(0,1,1,1);
  } 
  if (xaxis=="pt") {
    if (intervals2Plot != 1) {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0, 30);
      line = TLine(doLogPt ? 3 : 0, 1, 30, 1);
    } else {
      haxes = new TH1F("haxes","haxes",1,doLogPt ? 3 : 0,is18XXX ? 30 : 50);
      line = TLine(doLogPt ? 3 : 0,1,is18XXX ? 30 : 50,1);
    }
    if (doLogPt) c1->SetLogx();
  }
  if (xaxis=="rap") {
    haxes = new TH1F("haxes","haxes",1,0,2.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
    line = TLine(0,1,2.4,1);
  }
  if (xaxis=="cent") {
    haxes = new TH1F("haxes","haxes",1,0,420);
    haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X"));
    line = TLine(0,1,420,1);
  }
  haxes->GetYaxis()->SetRangeUser(0,1.5);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : ((xaxis=="rap") ? "|y|" : "N_{part}");
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();
//  line.Draw();
  
  double xshift=0.025;
//  if (xaxis=="cent") xshift=0.05;
//  TLegend *tleg = new TLegend(0.16+xshift,0.67,0.50, ((xaxis == "cent") && is18XXX) ? 0.89 : 0.83);
  TLegend *tleg(0x0);
  if (xaxis!="cent" && intervals2Plot == 2) tleg = new TLegend(0.44,0.50,0.76,0.62);
  else if (xaxis=="cent" && intervals2Plot == 2) tleg = new TLegend(0.19,0.16,0.51,0.28);
  else if ((xaxis=="cent" || xaxis=="rap" || xaxis=="zed")  && intervals2Plot == 1) tleg = new TLegend(0.51,0.47,0.83,0.62);
  else if (xaxis=="pt" && intervals2Plot == 1) tleg = new TLegend(0.19,0.49,0.51,0.64);
  else if (dononprompt && intervals2Plot == 3) tleg = new TLegend(0.56,0.47,0.88,0.62);
  else if (doprompt && intervals2Plot == 3) tleg = new TLegend(0.19,0.49,0.51,0.64);
  else if (plotFwdMid && intervals2Plot == 4) tleg = new TLegend(0.56,0.47,0.88,0.62);
  else tleg = new TLegend(0.56,0.42,0.88,0.62);
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  tleg->SetTextFont(42);
  tleg->SetTextSize(0.04);
  bool drawLegend(true);
  
  // prepare for the printing of the result tables
  const char* xname = (xaxis=="cent") ? "Centrality" : (xaxis=="pt" ? "\\pt" : (xaxis=="zed"? "z":"$|y|$"));
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s_NJJ_%s%s%s_%s.tex",outputDir.c_str(),plotPsi2S ? "Psi2S" : "JPsi",xaxis.c_str(),nameTag.c_str(), (xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr");
  string yname(plotPsi2S ? "$\\njj (\\psiP)$" : "$\\njj (\\Jpsi)$");
  inittex(texname, xname, yname);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  //map<anabin, TGraphAsymmErrors*>::const_iterator it_syst=theGraphs_syst.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    //TGraphAsymmErrors* tg_syst = it_syst->second;
    if (!tg) continue;
    
    theCats.push_back(thebin);
    
    TGraphErrors *tg_18XXX(0x0);
    
    int style =  cnt;
    int colorI =  cnt;
    int colorF = color(colorI)-11;
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
    //tg_syst->SetLineColor(color(colorI));
    //tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    if (xaxis=="cent") {
      // do not plot wide centrality bins
      prune(tg);
    }
    if (is18XXX && plotPsi2S)
    {
      Int_t nPoints = tg->GetN();
      double x, y, exl, exh;
//      double x_sys, y_sys;
      double liml, limh;
      for (Int_t i = 0 ; i < nPoints ; i++)
      {
        tg->GetPoint(i,x,y);
        if (y == 0)
        {
          exl = tg->GetErrorXlow(i);
          exh = tg->GetErrorXhigh(i);
          
          liml = tg->GetErrorYhigh(i);
//          tg_syst->GetPoint(i,x_sys,y_sys);
          limh = tg->GetErrorYhigh(i);;
          
          double dx(0.);
          if (xaxis=="cent")
          {
            x = HI::findNpartAverage(2.*(x-exl),2.*(x+exh));
            dx =10;
          }
          if (xaxis=="pt") dx = 0.5;
          
          drawArrow(x, liml, limh, dx, color(cnt));
          tg->SetPoint(i,-1000.,0);
          tg->SetPointError(i,0,0,0,0);
          //tg_syst->SetPoint(i,-1000.,0);
          //tg_syst->SetPointError(i,0,0,0,0);
            
          //delete point in res and syst
        }
      }
    }
    
    bool plot = true;
    if (plotFwdMid && (cnt!=0 && cnt!=3)) plot = false;
    
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
      x = rightA - (2*dx*cnt + dx);
      if (plotFwdMid && cnt>0) x = rightA - (2*dx + dx);
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
      tbox->Draw("l");
      TBox *tboxl = (TBox*) tbox->Clone("tboxl");
      tboxl->SetFillStyle(0);
      tboxl->Draw("l");
    }
    // }
    
    // Plot graphs after uncertainties to avoid overlap

    if (plot)
    {
      //tg_syst->Draw("5");
      gStyle->SetEndErrorSize(5);
      tg->Draw("P");
      // tg->Draw("[]");
    }
  
    TLatex tl;
    double tlx = 0.25; //0.92;
    double tly = 0.80; //0.69;
    tl.SetTextFont(42); // consistent font for symbol and plain text

    if (plot)
    {
        TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
        if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
        TString ptlabel = Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
        TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
        
        if (xaxis == "pt")
        {
          if (is18XXX) tleg->AddEntry(tg, raplabel, "p");
          else
          {
            if (intervals2Plot > 3) tleg->AddEntry(tg, raplabel, "p");
            else if (intervals2Plot == 3) tleg->AddEntry(tg, Form("Cent. %s",centlabel.Data()), "p");
            else if (intervals2Plot == 1) drawLegend = false;
            else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",raplabel.Data(),centlabel.Data()));
            //else tleg->AddEntry(tg, Form("%s, Cent. %s",raplabel.Data(),centlabel.Data()), "p");
          }
        }
        if (xaxis == "cent")
        {
          if (is18XXX) tleg->AddEntry(tg, Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()), "p");
          else
          {
            if (intervals2Plot > 3 ) tleg->AddEntry(tg, raplabel, "p");
            else if (intervals2Plot == 2) tleg->AddEntry(tg, ptlabel, "p");
            else if (intervals2Plot == 1) drawLegend = false;
            else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()));
            //else tleg->AddEntry(tg, Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()), "p");
          }
        }
        if (xaxis == "rap")
        {
          if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
          else if (intervals2Plot == 1) drawLegend = false;
          else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
          //else tleg->AddEntry(tg, Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()), "p");
        }
    }
    
    // print tex
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.centbin().low()/2.) << "\\% - " << (int) (it->first.centbin().high()/2.) << "\\%";
    if (xaxis == "cent" || xaxis == "rap") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";
    
    addline(texname,oss.str());
    printGraph(tg, texname);
    
    // for the centrality dependence: we want Npart plotted, not the centrality
    //if (xaxis == "cent") {
      //centrality2npart(tg, false, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
      //centrality2npart(tg_syst, true, (intervals2Plot > 3 && !plotFwdMid) ? ((30./1.8)*it->first.rapbin().low()) : 0.);
      //}
    
    cnt++;
    //it_syst++;
  }
  
  
  if (drawLegend) tleg->Draw();
  line.Draw();
  
  if(xaxis!="cent" && intervals2Plot > 3)
  {
    //      TLatex *tex = new TLatex(0.21,0.86,"Cent. 0-100%");
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 2 && !is18XXX)
  {
    TLatex *tex = new TLatex(0.2,0.78,"1.8 < |y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot > 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"6.5 < p_{T} < 50 GeV/c");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
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
  if(xaxis=="pt" && is18XXX)
  {
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="rap" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
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
  
  TLatex tl;
//  double tlx = 0.54; //0.92;
  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;
//  tl.SetTextAlign(32); // right adjusted
  tl.SetTextFont(42); // consistent font for symbol and plain text
  tl.SetTextSize(0.057); 
  if (doprompt) tl.DrawLatexNDC(tlx,tly,plotPsi2S ? "Prompt #psi(2S)" : "Prompt J/#psi");
  //if (doprompt) tl.DrawLatexNDC(tlx-0.08,tly,plotPsi2S ? "Prompt #psi(2S)" : "Prompt J/#psi");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"Nonprompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi from b hadrons");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi (b hadrons)");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi #leftarrow B");
//  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"B #rightarrow J/#psi");
  tl.SetTextSize(0.046);
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  // CMS_lumi( (TPad*) gPad, 103, iPos, "" );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_NJJ_%s%s%s_%s.root",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_NJJ_%s%s%s_%s.png",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_NJJ_%s%s%s_%s.pdf",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  
  if (plotFwdMid)
  {
    gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_NJJ_%s%s%s_2Rapranges_%s.root",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_NJJ_%s%s%s_2Rapranges_%s.png",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
    gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
    c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_NJJ_%s%s%s_2Rapranges_%s.pdf",outputDir.c_str(), plotPsi2S ? "Psi2S" : "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "",applyEff ? (applyAcc ? "accEffCorr" : "effCorr") : "noCorr"));
  }
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}

void centrality2npart(TGraphAsymmErrors* tg, bool issyst, double xshift){
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
  if (i==0) return kRed+2;
  else if (i==1) return kBlue+2;
  else if (i==2) return kMagenta+2;
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

void setOptions(bool adoprompt, bool adononprompt, bool aplotFwdMid, bool ais18XXX, bool aplotPsi2S, bool aplot14005, bool aapplyEff, bool aapplyAcc, bool adoLogPt, string anameTag_base){
  doprompt = adoprompt;
  dononprompt = adononprompt;
  plotFwdMid = aplotFwdMid;
  is18XXX = ais18XXX;
  plotPsi2S = aplotPsi2S;
  plot14005 = aplot14005;
  applyEff = aapplyEff;
  applyAcc = aapplyAcc;
  doLogPt = adoLogPt;
  
  nameTag_base = anameTag_base;
  
  if (plotPsi2S && !is18XXX)
  {
    cout << "[ERROR] Options not set: you want to plot Psi2 but bins are not 18XXX ones" << endl;
    printOptions();
    return;
  }
  if (doprompt) nameTag_base += "_prompt";
  if (dononprompt) nameTag_base += "_nonprompt";
  
//  if (aplotPsi2S) nameTag_base += "_Psi2S";
//  else nameTag_base += "_JPsi";
  
  if (plotFwdMid) nameTag_base += "_2RapRanges";
    
  if (is18XXX)  nameTag_base += "_18XXX";
  if (plot14005) nameTag_base += "_14005";
}

void printOptions() {
  cout <<
  "doprompt = " << doprompt << ", " <<
  "dononprompt = " << dononprompt << ", " <<
  "is18XXX = " << is18XXX << ", " <<
  "plotPsi2S = " << plotPsi2S << ", " <<
  "plot14005 = " << plot14005 << ", " <<
  "applyEff = " << applyEff << ", " <<
  "applyAcc = " << applyAcc << ", " <<
  "doLogPt = " << doLogPt << ", " <<
  "nameTag_base = \"" << nameTag_base << "\"" <<
  endl;
}

map<anabin, njj_input > readResults(const char* resultsFile)
{
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
      else if (cnt==4) centmin = value;
      else if (cnt==5) centmax = value;
      else if (cnt==6) theresult.naa = value;
      else if (cnt==7) theresult.dnaa_stat = value;
      else if (cnt==8) theresult.systaa = value;
      else if (cnt>8) {
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
