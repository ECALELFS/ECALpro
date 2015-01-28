#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooChebychev.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>
//Root Stuff
#include "TLatex.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TFormula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::pair;
using std::stringstream;

static const int MAX_IETA = 85;
static const int MAX_IPHI = 360;
static const int MIN_IETA = 1;
static const int MIN_IPHI = 1;

using namespace RooFit;

//5_3_6:  gROOT->ProcessLine(".include /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/include/")
//7_1_0:  gROOT->ProcessLine(".include /afs/cern.ch/cms/slc6_amd64_gcc481/lcg/roofit/5.34.22-cms2/include/")
//Usage: .x Total_fit.C+("2012C_epsilonPlots.root", "allEpsilon_EB", true, true )
void Total_fit( TString File, TString Hname, bool isEB, bool Are_pi0_=true ){

    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    //Files
    cout<<"Opening: "<<File.Data()<<endl;
    TFile* fin = TFile::Open(File.Data());
    //Histos
    TH1F *h     = (TH1F*) fin->Get( Hname.Data() );
    TString outName = "plots/FIT" + Hname + ".png";

    //Fit Method
    int ngaus=1; //1: simple Gaussian, 2: two Gaussian
    float xlo = Are_pi0_? 0.08:0.4, xhi = Are_pi0_? 0.22:0.65;
    RooRealVar x("x","#gamma#gamma invariant mass",xlo, xhi, "GeV/c^2");
    RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h);

    RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.116:0.57,  Are_pi0_? 0.105:0.5, Are_pi0_? 0.150:0.62,"GeV/c^{2}");
    RooRealVar sigma("sigma","#pi^{0} core #sigma",0.013, 0.005,0.020,"GeV/c^{2}");


    if(!isEB)  {
	  mean.setRange( Are_pi0_? 0.10:0.45, Are_pi0_? 0.140:0.62); // 0.200
	  mean.setVal(Are_pi0_? 0.120:0.55);
	  sigma.setRange(0.005, 0.060);
    }
    if(isEB){
	  mean.setRange(Are_pi0_? 0.105:0.47, Are_pi0_? 0.155:0.62);
	  sigma.setRange(0.003, 0.030);
    }

    RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
    Nsig.setVal( h->GetSum()*0.1);

    RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

    RooRealVar sigmaTail("sigmaTail","#pi^{0} tail #sigma",0.040, 0.020,0.065,"GeV/c^{2}");

    RooGaussian gaus2("gaus2","Tail Gaussian",x, mean,sigmaTail);

    RooRealVar fcore("fcore","f_{core}",0.9,0.,1.);
    RooAddPdf  signal("signal","signal model",RooArgList(gaus,gaus2),fcore);

    RooRealVar p0("p0","p0", 1000.,-1.e5,1.e5);
    RooRealVar p1("p1","p1", -3000.,-1.e5,1.e5);
    RooRealVar p2("p2","p2", 10000.,-1.e5,1.e5);
    RooRealVar p3("p3","p3", -10000.,-1.e5,1.e5);
    RooRealVar p4("p4","p4",-4000.,-1.e5,1.e5);
    RooRealVar p5("p5","p5", 5.,-1.e5,1.e5);
    RooRealVar p6("p6","p6", 6.,-1.e5,1.e5);

    RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
    RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
    RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
    RooRealVar cb3("cb3","cb3",-0.1, -0.5,0.5);
    RooRealVar cb4("cb4","cb4", 0.1, -1.,1.);
    RooRealVar cb5("cb5","cb5", 0.1, -1.,1.);
    RooRealVar cb6("cb6","cb6", 0.3, -1.,1.);
    RooArgList cbpars(cb0,cb1,cb2);
    if(Are_pi0_) cbpars.add( cb3);
//    if(isEB){
//	  cb3.setRange(-1,1.);
//	  cb4.setRange(-0.3,0.3);
//	  cbpars.add( cb4 );
//    }
    RooChebychev bkg("bkg","bkg model", x, cbpars );
    RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,1.e8);
    Nbkg.setVal( h->GetSum()*0.8 );

    RooAbsPdf* model=0;

    RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
    RooAddPdf model2("model","sig+bkg",RooArgList(signal,bkg),RooArgList(Nsig,Nbkg));

    if(ngaus==1)      model = &model1;
    else if(ngaus==2) model = &model2;

    RooNLLVar nll("nll","log likelihood var",*model,dh);
    RooMinuit m(nll);
    m.setVerbose(kFALSE);
    m.migrad();
    //RooFitResult* res = m.save() ;

    x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
    RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));
    RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

    RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
    //int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();
    x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());

    RooPlot*  xframe = x.frame(h->GetNbinsX());
    xframe->SetTitle(h->GetTitle());
    dh.plotOn(xframe);
    model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
    model->plotOn(xframe);

    xframe->Draw();
    TLatex lat;
    char line[300];
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.SetTextColor(1);

    float xmin(0.55), yhi(0.80), ypass(0.05);
    if(!Are_pi0_) yhi=0.30;
    sprintf(line,"Yield: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
    lat.DrawLatex(xmin,yhi, line);
    sprintf(line,"m_{#gamma#gamma}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
    lat.DrawLatex(xmin,yhi-ypass, line);
    sprintf(line,"#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
    lat.DrawLatex(xmin,yhi-2.*ypass, line);
    sprintf(line,"S/B (3#sigma): %.2f", (integralSig->getVal()*Nsig.getVal())/(integralBkg->getVal()*Nbkg.getVal()) );
    lat.DrawLatex(xmin,yhi-3.*ypass, line);
    sprintf(line,"S: %.2f", (integralSig->getVal()*Nsig.getVal()) );
    lat.DrawLatex(xmin,yhi-4.*ypass, line);
    sprintf(line,"B: %.2f", (integralBkg->getVal()*Nbkg.getVal()) );
    lat.DrawLatex(xmin,yhi-5.*ypass, line);
    myc1->SaveAs(outName.Data());
    cout<<"Histo Entries: "<<h->GetEntries()<<" and integral: "<<h->Integral()<<endl;
    cout<<"integralSig->getVal()*Nsig.getVal(): "<<integralSig->getVal()*Nsig.getVal()<<" : "<<integralSig->getVal()<<" * "<<Nsig.getVal()<<endl;;
    cout<<"integralBkg->getVal()*Nbkg.getVal(): "<<integralBkg->getVal()*Nbkg.getVal()<<" : "<<integralBkg->getVal()<<" * "<<Nbkg.getVal()<<endl;;
}
