#include "../plotVars.C"
#include "syst.C"

#include "TString.h"
#include <sstream>

using namespace std;

TString spoiname(""); // It can be NJpsi, BJpsi, NJpsi_prompt or NJpsi_nonprompt
//bool bins16004=false; // If false use 16-025 bins, if true use 16-004 bins
TString nameTag(""); // It can be 16025 or 16004 or...
bool midBins=true;
// DECLARATIONS
void plotSysts(anabin thebin, string xaxis, string collTag, bool plotEffSyst=false, bool plotSigSyst=true, bool plotGlobalSysts=false);
void plotSystsAll(const char* apoiname="NJpsi", bool plotEffSyst=false, bool plotSigSyst=true, bool plotGlobalSysts=false);


void plotSystsAll(const char* apoiname, bool plotEffSyst, bool plotSigSyst, bool plotGlobalSysts) {
  spoiname = apoiname;
  if (!spoiname.CompareTo("NJpsi") && !spoiname.CompareTo("BJpsi") && !spoiname.CompareTo("NJpsi_prompt") && !spoiname.CompareTo("NJpsi_nonprompt"))
  {
    cout << "[ERROR] : unknown systematic" << endl;
    return;
  }

  //mid JtPt
  nameTag = "midJtPt";
  //plotSysts(anabin(0.064,1.0,0,2.4,6.5,100,0,200),"z","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
  plotSysts(anabin(0.22,1.0,0,2.4,6.5,100,0,200),"z","PP",plotEffSyst,plotSigSyst,plotGlobalSysts);
  //plotSysts(anabin(0.064,1.0,0,2.4,6.5,100,0,180),"z","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
  plotSysts(anabin(0.22,1.0,0,2.4,6.5,100,0,180),"z","PbPb",plotEffSyst,plotSigSyst,plotGlobalSysts);
}

void plotSysts(anabin thebin, string xaxis, string collTag, bool plotEffSyst, bool plotSigSyst, bool plotGlobalSysts) {
  //cout << "[INFO] at the biginning of plotSyst" << endl;
  float zmin = thebin.zbin().low();
  float zmax = thebin.zbin().high();
  float rapmin = thebin.rapbin().low();
  float rapmax = thebin.rapbin().high();
  float ptmin = thebin.ptbin().low();
  float ptmax = thebin.ptbin().high();
  int centmin = thebin.centbin().low();
  int centmax = thebin.centbin().high();
  
  vector<string> tags;
  vector<TGraphErrors*> graphs;
  vector<map<anabin, syst> > systs;
  
  //cout<<"[INFO] read the systematics"<<endl;
  // total
  systs.push_back(readSyst_all(collTag.c_str(),spoiname.Data(),nameTag.Data(),plotEffSyst,plotEffSyst,plotSigSyst,"../"));
  tags.push_back(systs.back().begin()->second.name);
  //tags.push_back("Total");
  // detail
  //if (collTag=="PbPb") {
    
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_massBkg.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_massSig.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauErr.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauTrue.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauRes.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauBkg.csv",nameTag.Data(),spoiname.Data())));
  //tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PbPb_fulltnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PbPb_muidtnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PbPb_statnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
  //} else {
  if (plotSigSyst) {
    if (plotEffSyst) {
      if (collTag=="PbPb") {
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_fullAccEff.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
      }
      else {
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_fullAccEff.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      }
    }
    if (collTag=="PbPb") {
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_massBkg.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_massSig.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauErr.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauTrue.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauRes.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_ctauBkg.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
    }
    else{
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_massBkg.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_massSig.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_ctauErr.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_ctauTrue.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_ctauRes.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
      systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_ctauBkg.csv",nameTag.Data(),spoiname.Data())));
      tags.push_back(systs.back().begin()->second.name);
    }
    //systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_jetEnergyScale.csv",nameTag.Data(),spoiname.Data())));
    //tags.push_back(systs.back().begin()->second.name);
    //systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_jetEnergyRes.csv",nameTag.Data(),spoiname.Data())));
    //tags.push_back(systs.back().begin()->second.name);

//    systs.push_back(readSyst("csv/syst_PP_fulltnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PP_muidtnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PP_statnp.csv"));
//    tags.push_back(systs.back().begin()->second.name);
  }
  else if (plotEffSyst)
    {
      if (collTag=="PbPb") {
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnptagSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnpmuidSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnpmuidStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnptrgSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnptrgStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnptrkSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_tnptrkStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_EffStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);

	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_AccStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PbPb_EffMisMod.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
      }
      else {
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_tnpglbSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_tnpglbStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_tnpmuidtrgSyst.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_tnpmuidtrgStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_EffStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name)
;
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_AccStat.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name)
;
	systs.push_back(readSyst(Form("csv/syst_%s_%s_PP_EffMisMod.csv",nameTag.Data(),spoiname.Data())));
	tags.push_back(systs.back().begin()->second.name);
      }
      //else cout << "[WARNING] you have to specify the type of uncertainties to plot";
    }
  // global systs
  //if (plotGlobalSysts) {
  //systs.push_back(readSyst("csv/syst_PbPb_taa.csv"));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst("csv/syst_PbPb_Nmb.csv"));
  //tags.push_back(systs.back().begin()->second.name);
  //systs.push_back(readSyst("csv/syst_PP_lumi.csv"));
  //tags.push_back(systs.back().begin()->second.name);
  //    systs.push_back(readSyst("csv/syst_PP_trk.csv"));
//    tags.push_back(systs.back().begin()->second.name);
//    systs.push_back(readSyst("csv/syst_PbPb_trk.csv"));
//    tags.push_back(systs.back().begin()->second.name);
  //}
  
  set<anabin> sb;
  //if (bins16004) sb = allbins16004();
  //else sb = allbins();
  cout<<"[INFO] setting the bins"<<endl;
  //if (midBins) sb = midbins18012();
  //else sb = forbins18012();
  sb = allbins19xxx();
  for (unsigned int i=0; i<systs.size(); i++) {
    map<anabin, syst> thesyst = systs[i];
    vector<double> x, y, dx, dy;
    
    double valmax=-1e99;
    
    for (set<anabin>::const_iterator it=sb.begin(); it!=sb.end(); it++) {
      anabin it2 = *it;
      
      if (!binok(thebin,xaxis,it2,false)) continue;
      if (thesyst.find(*it)==thesyst.end()) continue;
      
      double low, high;
      if (xaxis=="z") {
	low= it->zbin().low();
	high = it->zbin().high();
        if (high<=0.22) continue;
      }
      else if (xaxis=="pt") {
        low= it->ptbin().low();
        high = it->ptbin().high();
      } else if (xaxis=="rap") {
        low= it->rapbin().low();
        high = it->rapbin().high();
        if ((low==0 && high<=0.61 && high>=0.59 ) || (low>=0.59 && low<=0.61 && high>=1.19 && high <=1.21) || (low>=1.19 && low<=1.21 && high>=1.79 && high <=1.81) || (low>=1.79 && low<=1.81 && high>=2.39 && high <=2.41)) continue;
      } else {
        low= it->centbin().low()/2;
        high = it->centbin().high()/2;
      }
      x.push_back((low+high)/2.);
      dx.push_back((high-low)/2.);
      y.push_back(0);
      dy.push_back(thesyst[*it].value);
      valmax = max(valmax,dy.back());
    }

    TGraphErrors *thegraph = new TGraphErrors(x.size(),x.data(),y.data(),dx.data(),dy.data());
    TH1F *haxes=NULL;
    if (xaxis=="z") {
      haxes = new TH1F(tags[i].c_str(),Form(";z(J/#psi);Syst_%s",spoiname.Data()),5,0,1);
    } else if (xaxis=="pt") {
      haxes = new TH1F(tags[i].c_str(),Form(";p_{T} (GeV/c);Syst_%s",spoiname.Data()),1,0,50);
    } else if (xaxis=="cent") {
      haxes = new TH1F(tags[i].c_str(),Form(";Centrality percentile;Syst_%s",spoiname.Data()),1,0,100);
    } else { // if (xaxis=="rap")
      haxes = new TH1F(tags[i].c_str(),Form(";|y|;Syst_%s",spoiname.Data()),1,0,2.4);
    }
    haxes->GetYaxis()->SetLimits(-1.5*valmax, 2.5*valmax);
    haxes->GetYaxis()->SetRangeUser(-1.5*valmax, 2.5*valmax);
    //haxes->GetYaxis()->SetLimits(-1,1);
    //haxes->GetYaxis()->SetRangeUser(-1,1);
    haxes->GetYaxis()->SetTitleOffset(1.4);
    thegraph->SetHistogram(haxes);

    thegraph->Sort();
    graphs.push_back(thegraph);
  }
  
  cout << "[INFO] plotting the graphs" <<endl;  
  plotGraphs(graphs, tags, "systematics", collTag,
             //Form("z%i%i_pt%i%i_rap%i%i_cent%i%i_%s",(int)zmin*10,(int)zmax*10,(int)ptmin*10,(int)ptmax*10,(int)rapmin*10,(int)rapmax*10,centmin,centmax,nameTag.Data()),
	     Form("z%i%i_pt%i%i_rap%i%i%s_%s",(int)(zmin*10),(int)(zmax*10),(int)(ptmin*10),(int)(ptmax*10),(int)(rapmin*10),(int)(rapmax*10), plotSigSyst?(plotEffSyst?"_SigEff":"_Sig"):(plotEffSyst?"_Eff":""), nameTag.Data()),
             //Form("%.1f<z<%.1f, %.1f<|y|<%.1f, %.1f<p_{T}<%.1f, %i-%i%s",zmin,zmax,rapmin,rapmax,ptmin,ptmax,centmin/2,centmax/2,"%"));
	     Form("%.1f<|y|<%.1f, %.1f<p_{T}<%.1f, %i-%i%s",rapmin,rapmax,ptmin,ptmax,centmin/2,centmax/2,"%"));
}
