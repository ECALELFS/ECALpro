#ifndef utility_h
#define utility_h

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <map>
#include <iomanip> //for input/output manipulators

#include <algorithm>  // to use the "reverse" function to reverse the order in the array
#include <Rtypes.h> // to use kColor

//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TStyle.h>
#include <TString.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"

using namespace RooFit;
using namespace std;

static string PhpToCopy = "/afs/cern.ch/user/m/mciprian/public/index.php";
// some parameters used for the fit

static bool Are_pi0_ = false;  // fixme: to be passed as an option
static double upper_bound_pi0mass_EB = 0.15;
static double upper_bound_pi0mass_EE = 0.16;
static double upper_bound_etamass_EB = 0.62;
static double upper_bound_etamass_EE = 0.62;
struct Pi0FitResult {
  RooFitResult* res;
  float Nsig;     // signal total            
  float Nsigerr;
  float Nbkg;     // bkg total              
  float Nbkgerr;
  float S;     // signal in 3 sigma region            
  float Serr;
  float B;     // bkg in 3 sigma region              
  float Berr;
  float SoB;  // S/B   
  float SoBerr;
  float mean;     // fit mean (peak) 
  float meanerr;
  float sigma;     // fit sigma              
  float sigmaerr;
  float chi2;
  float chi2red;
  int nFitParam;
  int dof;  // after subtracting fit parameters
  float probchi2;
};

//======================================================

void createPlotDirAndCopyPhp(const string& outputDIR) {

  if (outputDIR != "./") {
    system(("mkdir -p " + outputDIR).c_str());
    system(("cp "+ PhpToCopy + " " + outputDIR).c_str());
  }

}

//======================================================                                                                                                 

string getStringFromDouble(const Double_t& num = 1.0, const Double_t epsilon = 0.00001) {

  Int_t i = (Int_t) num;
  Int_t int_decim = (Int_t) (1000 * (num - (Double_t) i + epsilon));
  if      (int_decim%1000 == 0) return string(Form("%dp%d",i,int_decim/1000));
  else if (int_decim%100 == 0)  return string(Form("%dp%d",i,int_decim/100));
  else if (int_decim%10 == 0)   return string(Form("%dp%d",i,int_decim/10));
  else                          return string(Form("%dp%d",i,int_decim));

}

//======================================================                                                                                                                     

void myAddOverflowInLastBin(TH1 *h) {

  Int_t lastBinNumber = h->GetNbinsX();
  Int_t overflowBinNumber = 1 + lastBinNumber;
  Double_t lastBinContent = h->GetBinContent(lastBinNumber);
  Double_t overflowBinContent = h->GetBinContent(overflowBinNumber);
  Double_t lastBinError = h->GetBinError(lastBinNumber);
  Double_t overflowBinError = h->GetBinError(overflowBinNumber);

  // add content of overflow bin in last bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)                           
  h->SetBinContent(lastBinNumber, lastBinContent + overflowBinContent);
  h->SetBinError(lastBinNumber, sqrt(lastBinError * lastBinError + overflowBinError * overflowBinError));
  // deleting content of last bin (safer, since I might be using that bin to add it again somewhere and I want it to be empty)                                               
  h->SetBinContent(overflowBinNumber,0.0);
  h->SetBinError(overflowBinNumber,0.0);

}


//======================================================                                               

void myRebinHisto(TH1 *h, const Int_t rebinFactor = 1) {

  if (rebinFactor != 1) {
    h->Rebin(rebinFactor);
    if ( (h->GetNbinsX() % rebinFactor) != 0) myAddOverflowInLastBin(h);
  }

}

//=================================================

Pi0FitResult drawHisto(TH1* hist = NULL, 
	 	       const bool isEB = true, 
	 	       const string& outDir = "./", 
	 	       const string& hName = "", 
	 	       const double lumi = 8.6
		       )
{

  createPlotDirAndCopyPhp(outDir);

  TGaxis::SetMaxDigits(3); 

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetRightMargin(0.06);

  hist->SetStats(0);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);

  hist->SetTitle(0);
  
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  hist->GetXaxis()->SetTitleSize(0.05);  
  hist->GetXaxis()->SetTitleOffset(0.9);
  hist->GetXaxis()->SetRangeUser(Are_pi0_ ? 0.05 : 0.35, Are_pi0_ ? 0.25 : 0.7);

  double maxY = hist->GetBinContent(hist->GetMaximumBin());
  hist->GetYaxis()->SetRangeUser(0.0, 1.2*maxY);
  hist->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->Draw("EP");

  /////////////////////////
  /////////////////////////
  // fit and draw result

  // better not to use extreme values of histogram as the fit range, because the actual range is shorter
  // see example here --> http://mciprian.web.cern.ch/mciprian/test_plot/pi0Mass_EB_h_xtal_iter0.png
  // I suggest using 0.080 and 0.21 for pi0
  RooRealVar x("x","#gamma#gamma invariant mass", Are_pi0_? 0.07:0.4, Are_pi0_? 0.21:0.65, "GeV/c^2");
  if (Are_pi0_ && not isEB) x.setRange(0.075, 0.24);

  RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),hist);

  RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.5, Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB,"GeV/c^{2}");
  RooRealVar sigma("sigma","#pi^{0} core #sigma",0.011, 0.005,Are_pi0_ ? 0.015 : 0.03,"GeV/c^{2}");
  if(not isEB)  {
    mean.setRange( Are_pi0_? 0.1:0.45, Are_pi0_? upper_bound_pi0mass_EE:upper_bound_etamass_EE);
    mean.setVal(Are_pi0_? 0.13:0.55);
    sigma.setRange(0.005, Are_pi0_ ? 0.020: 0.035);
  }

  RooRealVar Nsig("Nsig","#pi^{0} yield", hist->Integral()*0.15,0.,hist->Integral()*10.0);
  //Nsig.setVal( hist->Integral()*0.1);

  //sig model
  RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

  // bkg model
  RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
  RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
  RooRealVar cb2("cb2","cb2", 0.1,  -1.,1.);
  RooRealVar cb3("cb3","cb3", -0.1, -0.5,0.5);
  // try to use a second order polynomial, if the fit is bad add other terms                                                                                               
  // if you start with many terms, the fit creates strange curvy shapes trying to fit the statistical fluctuations                                                         
  // 2nd order means a curve with no change of concavity 
  //RooArgList cbpars(cb0,cb1,cb2,cb3);
  // FIXME: should try to repeat fit with more free parameters in B if fit not satisfactory
  RooArgList cbpars(cb0,cb1,cb2);
  RooChebychev bkg("bkg","bkg model", x, cbpars );

  RooRealVar Nbkg("Nbkg","background yield",hist->Integral()*0.85,0.,hist->Integral()*10.0);
  //Nbkg.setVal( hist->Integral()*0.8 );

  RooAbsPdf* model = 0;
  // can use many models
  RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
  // modelXXX ...
  model = &model1;

  RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
  //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same                                                                                      

  // FIT 1
  // copied from ECALpro
  // RooMinuit m(nll);
  // m.setVerbose(kFALSE);
  // //m.setVerbose(kTRUE);                                                                                                       
  // m.migrad();
  // //m.hesse();                                                                                                                          
  // RooFitResult* res = m.save() ;

  // FIT2
  // copied from Raffaele Gerosa
  RooMinimizer mfit(nll);
  mfit.setVerbose(kFALSE);
  mfit.setPrintLevel(-1);
  cout << "######### Minimize" << endl;
  mfit.minimize("Minuit2","minimize");
  cout << "######### Minimize hesse " << endl;
  mfit.minimize("Minuit2","hesse");
  cout<<"######### Estimate minos errors for all parameters"<<endl;
  mfit.minos(RooArgSet(Nsig,Nbkg));
  RooFitResult* res = mfit.save("res") ;

  // FIT 1 and FIT 2 yields practically the same result, using the second

  cout << "print fit result" << endl;
  res->Print();

  RooChi2Var chi2("chi2","chi2 var",*model,dh, true);

  int nFitParam = res->floatParsFinal().getSize();
  // FIXME: must use bins in the fit range
  int ndof = hist->GetNbinsX() - nFitParam; // nBins - floating parameters in model after fit

  //compute S/B and chi2                 
  // use 3 sigma range around mean to get S/B                                                                                       
  x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
  RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));

  RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

  float normSig = integralSig->getVal();
  float normBkg = integralBkg->getVal();

  RooPlot*  xframe = x.frame(hist->GetNbinsX());
  xframe->SetTitle(0);
  dh.plotOn(xframe,Name("data")); 
  model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
  model->plotOn(xframe,Name("model"));

  Pi0FitResult pi0res; // this is the output value of this method                                                                                                          
  pi0res.res = res;
  pi0res.Nsig = Nsig.getVal();
  pi0res.Nsigerr = Nsig.getError();
  pi0res.Nbkg = Nbkg.getVal();
  pi0res.Nbkgerr = Nbkg.getError();
  pi0res.S = normSig*Nsig.getVal();
  pi0res.Serr = normSig*Nsig.getError();
  pi0res.B = normBkg*Nbkg.getVal();
  pi0res.Berr = normBkg*Nbkg.getError();
  pi0res.mean = mean.getVal();
  pi0res.meanerr = mean.getError();
  pi0res.sigma = sigma.getVal();
  pi0res.sigmaerr = sigma.getError();
  pi0res.SoB =  pi0res.S/pi0res.B;
  pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) +
				    pow(pi0res.Berr/pi0res.B,2) ) ;
  pi0res.dof = ndof;
  pi0res.nFitParam = nFitParam;
  pi0res.chi2red = xframe->chiSquare("model","data",nFitParam);
  pi0res.chi2 = pi0res.chi2red * pi0res.dof;
  pi0res.probchi2 = TMath::Prob(pi0res.chi2, ndof);

  xframe->Draw("same");


  TLatex lat;
  char line[300];
  lat.SetNDC();
  lat.SetTextSize(0.035);
  lat.SetTextColor(1);

  float xmin(0.6), yhi(0.85), ypass(0.05);
  if(isEB and not Are_pi0_) yhi=0.30;
  sprintf(line,"Nsig: %.0f #pm %.0f", pi0res.Nsig, pi0res.Nsigerr );
  lat.DrawLatex(xmin,yhi, line);

  sprintf(line,"m_{#gamma#gamma}: %.2f #pm %.2f", pi0res.mean*1000., pi0res.meanerr*1000. );
  lat.DrawLatex(xmin,yhi-ypass, line);

  sprintf(line,"#sigma: %.2f #pm %.2f (%.2f%s)", pi0res.sigma*1000., pi0res.sigmaerr*1000., pi0res.sigma*100./pi0res.mean, "%" );
  lat.DrawLatex(xmin,yhi-2.*ypass, line);

  sprintf(line,"S/B(3#sigma): %.2f #pm %.2f", pi0res.SoB, pi0res.SoBerr );                                                                                               
  lat.DrawLatex(xmin,yhi-3.*ypass, line);

  sprintf(line,"#Chi^{2}: %.2f (%d dof)", pi0res.chi2, pi0res.dof );
  lat.DrawLatex(xmin,yhi-4.*ypass, line);

  // end of fit part
  /////////////////////////
  /////////////////////////

  canvas->RedrawAxis("sameaxis");
  if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),false,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),false,false);
  setTDRStyle();

  string title = Are_pi0_ ? "pi0mass" : "eta0mass";
  title += isEB ? "_EB_" : "_EE_";
  string canvasTitle = outDir + title + hName;
  canvas->SaveAs((canvasTitle + ".png").c_str());

  // save fit parameters in file named as the canvas but with txt extension
  
  string fitParameterFileName = canvasTitle + ".txt";
  ofstream fitParameterFile(fitParameterFileName.c_str(),ios::out);
  if ( !fitParameterFile.is_open() ) {  
    cout<<"Error: unable to open file " << fitParameterFileName <<" !"<<endl;
    exit(EXIT_FAILURE);
  } else {
    fitParameterFile << "FIT PARAMETERS:" << endl;
    fitParameterFile << "-------------------" << endl;
    fitParameterFile << setw(15) <<  " Nsig: " << pi0res.Nsig << " +/- " << pi0res.Nsigerr << endl;
    fitParameterFile << setw(15) <<  " Nbkg: " << pi0res.Nbkg << " +/- " << pi0res.Nbkgerr << endl;
    fitParameterFile << setw(15) <<  " Nsig(3sigma): " << pi0res.S << " +/- " << pi0res.Serr << endl;
    fitParameterFile << setw(15) <<  " Nbkg(3sigma): " << pi0res.B << " +/- " << pi0res.Berr << endl;
    fitParameterFile << setw(15) <<  " S/B(3sigma): " << pi0res.SoB << " +/- " << pi0res.SoBerr << endl;
    fitParameterFile << setw(15) <<  " mean: " << pi0res.mean*1000. << " +/- " << pi0res.meanerr*1000. << endl;
    fitParameterFile << setw(15) <<  " sigma: " << pi0res.sigma*1000. << " +/- " << pi0res.sigmaerr*1000. << endl;
    fitParameterFile << setw(15) <<  " chi2: " << pi0res.chi2 << endl;
    fitParameterFile << setw(15) <<  " DOF (fit.par.sub): " << pi0res.dof << endl;
    fitParameterFile << setw(15) <<  " nFitParam: " << pi0res.nFitParam << endl;
    fitParameterFile << setw(15) <<  " chi2(reduced): " << pi0res.chi2red << endl;
    fitParameterFile << setw(15) <<  " prob(chi2): " << pi0res.probchi2 << endl;
    fitParameterFile.close();
  }

  return pi0res;

}


//=============================================================


Pi0FitResult fitMassSingleHisto(TH1* hist) {

  /////////////////////////
  /////////////////////////
  // fit and draw result

  // better not to use extreme values of histogram as the fit range, because the actual range is shorter
  // see example here --> http://mciprian.web.cern.ch/mciprian/test_plot/pi0Mass_EB_h_xtal_iter0.png
  // I suggest using 0.080 and 0.21 for pi0
  RooRealVar x("x","#gamma#gamma invariant mass", Are_pi0_? 0.07:0.4, Are_pi0_? 0.22:0.65, "GeV/c^2");

  RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),hist);

  RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.5, Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB,"GeV/c^{2}");
  RooRealVar sigma("sigma","#pi^{0} core #sigma",0.011, 0.005,Are_pi0_? 0.015: 0.03,"GeV/c^{2}");

  string hname = hist->GetName();
  bool isEB = (hname.find("EB") != string::npos) ? true : false;
 
  if(not isEB)  {
    mean.setRange( Are_pi0_? 0.1:0.45, Are_pi0_? upper_bound_pi0mass_EE:upper_bound_etamass_EE);
    mean.setVal(Are_pi0_? 0.13:0.55);
    sigma.setRange(0.005, Are_pi0_ ? 0.020 : 0.030);
  }

  RooRealVar Nsig("Nsig","#pi^{0} yield", hist->Integral()*0.15,0.,hist->Integral()*100.0);
  //Nsig.setVal( hist->Integral()*0.1);

  //sig model
  RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

  // bkg model
  RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
  RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
  RooRealVar cb2("cb2","cb2", 0.1,  0.,1.);
  // try to use a second order polynomial, if the fit is bad add other terms                                                                                               
  // if you start with many terms, the fit creates strange curvy shapes trying to fit the statistical fluctuations                                                         
  // 2nd order means a curve with no change of concavity 
  RooArgList cbpars(cb0,cb1,cb2);
  RooChebychev bkg("bkg","bkg model", x, cbpars );

  RooRealVar Nbkg("Nbkg","background yield",hist->Integral()*0.85,0.,hist->Integral()*100.0);
  //Nbkg.setVal( hist->Integral()*0.8 );

  RooAbsPdf* model = 0;
  // can use many models
  RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
  // modelXXX ...
  model = &model1;

  RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
  //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same                                                                                     


  // FIT 1
  // copied from ECALpro
  // RooMinuit m(nll);
  // m.setVerbose(kFALSE);
  // //m.setVerbose(kTRUE);                                                    
  // m.migrad();
  // //m.hesse();                                                                 
  // RooFitResult* res = m.save() ;
  // FIT2
  // copied from Raffaele Gerosa
  RooMinimizer mfit(nll);
  mfit.setVerbose(kFALSE);
  mfit.setPrintLevel(-1);
  cout << "######### Minimize" << endl;
  mfit.minimize("Minuit2","minimize");
  cout << "######### Minimize hesse " << endl;
  mfit.minimize("Minuit2","hesse");
  cout<<"######### Estimate minos errors for all parameters"<<endl;
  mfit.minos(RooArgSet(Nsig,Nbkg));
  RooFitResult* res = mfit.save("res") ;

  // FIT 1 and FIT 2 yields practically the same result, using the second

  cout << "print fit result" << endl;
  res->Print();

  RooChi2Var chi2("chi2","chi2 var",*model,dh, true);

  int nFitParam = res->floatParsFinal().getSize();
  int ndof = hist->GetNbinsX() - nFitParam; // nBins - floating parameters in model after fit

  //compute S/B and chi2                 
  // use 3 sigma range around mean to get S/B                                        
                                               
  x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
  RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));

  RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

  float normSig = integralSig->getVal();
  float normBkg = integralBkg->getVal();
  RooPlot*  xframe = x.frame(hist->GetNbinsX());
  xframe->SetTitle(0);
  dh.plotOn(xframe,Name("data"));  // already drawn
  model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
  model->plotOn(xframe,Name("model"));

  Pi0FitResult pi0res; // this is the output value of this method   
                                                                                    
  pi0res.res = res;
  pi0res.Nsig = Nsig.getVal();
  pi0res.Nsigerr = Nsig.getError();
  pi0res.Nbkg = Nbkg.getVal();
  pi0res.Nbkgerr = Nbkg.getError();
  pi0res.S = normSig*Nsig.getVal();
  pi0res.Serr = normSig*Nsig.getError();
  pi0res.B = normBkg*Nbkg.getVal();
  pi0res.Berr = normBkg*Nbkg.getError();
  pi0res.mean = mean.getVal();
  pi0res.meanerr = mean.getError();
  pi0res.sigma = sigma.getVal();
  pi0res.sigmaerr = sigma.getError();
  pi0res.SoB =  pi0res.S/pi0res.B;
  pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) +
				    pow(pi0res.Berr/pi0res.B,2) ) ;
  pi0res.dof = ndof;
  pi0res.nFitParam = nFitParam;
  pi0res.chi2red = xframe->chiSquare("model","data",nFitParam);
  pi0res.chi2 = pi0res.chi2red * (double) pi0res.dof;
  pi0res.probchi2 = TMath::Prob(pi0res.chi2, ndof);

  return pi0res;


}

//=============================================================

void fitHistoSavePar(const vector<TH1*>& vecHist1d_orig = {},
		     const string& canvasName = "default", 
		     const string& outputDIR = "./",
		     const vector<string>& vecLegEntry = {""}		     
		     ) 
{

  // save fit parameters in file named as the canvas but with txt extension
  vector<Pi0FitResult> vfit;
  for (UInt_t i = 0; i < vecHist1d_orig.size(); i++) {
    vfit.push_back( fitMassSingleHisto(vecHist1d_orig[i]) );
  }
  
  string fitParameterFileName = outputDIR + canvasName + ".txt";
  ofstream fitParameterFile(fitParameterFileName.c_str(),ios::out);
  if ( !fitParameterFile.is_open() ) {  
    cout<<"Error: unable to open file " << fitParameterFileName <<" !"<<endl;
    exit(EXIT_FAILURE);
  } else {
    fitParameterFile << "DUMPING DETAILED LIST OF FIT PARAMETERS FOR EACH PLOT" << endl;
    fitParameterFile << endl;
    for (UInt_t i = 0; i < vecHist1d_orig.size(); i++) {
      fitParameterFile << "-------------------" << endl;
      fitParameterFile << vecHist1d_orig[i]->GetName() << "    " << vecLegEntry[i] << endl;    
      fitParameterFile << "-------------------" << endl;
      fitParameterFile << setw(15) <<  " Nsig: " << vfit[i].Nsig << " +/- " << vfit[i].Nsigerr << endl;
      fitParameterFile << setw(15) <<  " Nbkg: " << vfit[i].Nbkg << " +/- " << vfit[i].Nbkgerr << endl;
      fitParameterFile << setw(15) <<  " Nsig(3sigma): " << vfit[i].S << " +/- " << vfit[i].Serr << endl;
      fitParameterFile << setw(15) <<  " Nbkg(3sigma): " << vfit[i].B << " +/- " << vfit[i].Berr << endl;
      fitParameterFile << setw(15) <<  " S/B(3sigma): " << vfit[i].SoB << " +/- " << vfit[i].SoBerr << endl;
      fitParameterFile << setw(15) <<  " mean: " << vfit[i].mean*1000. << " +/- " << vfit[i].meanerr*1000. << endl;
      fitParameterFile << setw(15) <<  " sigma: " << vfit[i].sigma*1000. << " +/- " << vfit[i].sigmaerr*1000. << endl;
      fitParameterFile << setw(15) <<  " chi2: " << vfit[i].chi2 << endl;
      fitParameterFile << setw(15) <<  " DOF (fit.par.sub): " << vfit[i].dof << endl;
      fitParameterFile << setw(15) <<  " nFitParam: " << vfit[i].nFitParam << endl;
      fitParameterFile << setw(15) <<  " chi2(reduced): " << vfit[i].chi2red << endl;
      fitParameterFile << setw(15) <<  " prob(chi2): " << vfit[i].probchi2 << endl;
      fitParameterFile << endl;

    }
    fitParameterFile << endl;
    fitParameterFile << "SUMMARY TABLE OF MAIN PARAMETERS" << endl;
    fitParameterFile << "-------------------" << endl;
    /* fitParameterFile << setw(25) << " " << "FIT PARAMETERS:" << endl; */
    /* fitParameterFile << setw(25) << " " << "S and B evaluated in 3 sigma window around mean" << endl; */
    fitParameterFile << "FIT PARAMETERS:" << endl;
    fitParameterFile << "S and B evaluated in 3 sigma window around mean" << endl;
    fitParameterFile << "-------------------" << endl;
    fitParameterFile << setw(25) << "label" << setw(16) << "S" << setw(16) << "B" << setw(16) << "S/B" << setw(16) << "mean" << setw(16) << "sigma" << endl;
    for (UInt_t i = 0; i < vecHist1d_orig.size(); i++) {
      fitParameterFile << setw(25) << vecLegEntry[i];
      fitParameterFile << fixed << setprecision(0) << setw(10) << vfit[i].S << "+/-" << setw(4) << vfit[i].Serr;
      fitParameterFile << fixed << setprecision(0) << setw(10) << vfit[i].B << "+/-" << setw(4) << vfit[i].Berr;
      fitParameterFile << fixed << setprecision(2) << setw(10) << vfit[i].SoB << "+/-" << setw(4) << vfit[i].SoBerr;
      fitParameterFile << fixed << setprecision(2) << setw(10) << 1000.*vfit[i].mean << "+/-" << setw(4) << 1000.*vfit[i].meanerr;
      fitParameterFile << fixed << setprecision(2) << setw(10) << 1000.*vfit[i].sigma << "+/-" << setw(4) << 1000.*vfit[i].sigmaerr;
      fitParameterFile << endl;
    }

    fitParameterFile.close();

    

  }



}


//=============================================================


Bool_t getAxisRangeFromUser(string& xAxisName, Double_t& xmin, Double_t& xmax, 
			    const string& xAxisNameTmp = "", 
			    const string& separator = "::", 
			    const string& rangeSeparator = ","
			    ) {
  
  Bool_t setXAxisRangeFromUser = false;
  size_t pos = xAxisNameTmp.find(separator);
    
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    xAxisName.assign(xAxisNameTmp, 0, pos);
    xrange.assign(xAxisNameTmp, pos + separator.size(), string::npos);
    pos = xrange.find(rangeSeparator);
    string numString = "";
    numString.assign(xrange,0,pos);
    xmin = std::stod(numString);
    numString.assign(xrange,pos + rangeSeparator.size(), string::npos);
    xmax = std::stod(numString);
  } else {
    xAxisName = xAxisNameTmp;
  }

  return setXAxisRangeFromUser;

}

//=============================================================

void draw_nTH1(const vector<TH1*>& vecHist1d_orig = {}, 
	       const string& xAxisNameTmp = "", 
	       const string& yAxisName = "Events", 
	       const string& canvasName = "default", 
	       const string& outputDIR = "./", 
	       const vector<string>& vecLegEntry = {""},
	       const string& ratioPadYaxisName = "var/nominal",
	       const Double_t lumi = -1.0, 
	       const Int_t rebinFactorConst = 1, 
	       const Bool_t drawPlotLogY = true,
	       const Bool_t drawRatioWithNominal = false,
	       const vector<Int_t>& rebinFactorPerHisto = {1}
	       ) 
{

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  // assume the "nominal histogram is the first one

  // this is needed if we rebin and use again the histogram outside the function (we would be modifying the histogram passed to the function)
  vector<TH1*> vecHist1d;
  vector<Double_t> integral_orig;
  for (UInt_t i = 0; i < vecHist1d_orig.size(); i++) {
    vecHist1d.push_back(new TH1F(*((TH1F*) vecHist1d_orig[i]->Clone()) ) );
    integral_orig.push_back(vecHist1d_orig[i]->Integral());
  }

  double legLowY = 0.75;
  if (vecHist1d.size() > 5) legLowY = max( 0.5, legLowY - 0.03 * (vecHist1d.size() - 5) );


  // if using the rebinning for each histogram, the global one is ignored
  Int_t rebinFactor = rebinFactorConst;

  string xAxisName = "";
  Double_t xmin = 0;
  Double_t xmax = 0;
  Bool_t setXAxisRangeFromUser = getAxisRangeFromUser(xAxisName, xmin, xmax, xAxisNameTmp);

  // FIXME
  /* Double_t intNum, intDen, errNum, errDen; */
  /* intNum = h1->IntegralAndError(1,h1->GetNbinsX(),errNum); */
  /* intDen = h2->IntegralAndError(1,h2->GetNbinsX(),errDen); */
  /* Double_t IntegralRatio = intNum/intDen; */
  /* Double_t ratioError = IntegralRatio * sqrt(errNum*errNum/(intNum*intNum) + errDen*errDen/(intDen*intDen)); */


  // cout << "xAxisName = " << xAxisName << "   xmin = " << xmin << "  xmax = " << xmax << endl;

  // rebin
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    if (rebinFactorPerHisto.size() == vecHist1d.size()) {
      myRebinHisto(vecHist1d[i],rebinFactorPerHisto[i]);
      if (yAxisName == "a.u.") vecHist1d[i]->Scale(1./(rebinFactorPerHisto[i]*vecHist1d[i]->Integral()));
    } else {
      myRebinHisto(vecHist1d[i],rebinFactor);
      if (yAxisName == "a.u.") vecHist1d[i]->Scale(1./vecHist1d[i]->Integral());
    }
    vecHist1d[i]->SetStats(0);
  }

  // first rebin and then fit
  if (xAxisNameTmp.find("#gamma#gamma invariant mass (GeV/c^{2})") != string::npos) fitHistoSavePar(vecHist1d_orig, canvasName, outputDIR, vecLegEntry );

  // rescale if required (only after fit)
  if (yAxisName == "a.u.")  {
    for (UInt_t i = 0; i < vecHist1d.size(); i++) {
      if (rebinFactorPerHisto.size() == vecHist1d.size()) vecHist1d[i]->Scale(1./(rebinFactorPerHisto[i]*vecHist1d[i]->Integral()));
      else                                                vecHist1d[i]->Scale(1./vecHist1d[i]->Integral());
    }
  }


  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  if (drawRatioWithNominal) canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) vecHist1d[0]->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  Int_t colorList[] = {kBlack, kBlue, kRed+1, kGreen+2, kOrange+1, kCyan+1, kViolet, kGreen, kCyan, kGray+1, kYellow+2};
  vector<Int_t> histColor;
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)         
    vecHist1d[i]->SetLineColor(colorList[i]);
    vecHist1d[i]->SetLineWidth(2);
    vecHist1d[i]->SetFillColor(0);
  }

  if (drawRatioWithNominal) {
    vecHist1d[0]->GetXaxis()->SetLabelSize(0);
    vecHist1d[0]->GetXaxis()->SetTitle(0);
  } else {
    vecHist1d[0]->GetXaxis()->SetTitle(xAxisName.c_str());
    // vecHist1d[0]->GetXaxis()->SetTitleOffset(0.8);
    vecHist1d[0]->GetXaxis()->SetLabelSize(0.04);
    vecHist1d[0]->GetXaxis()->SetTitleSize(0.05);    
  }
  vecHist1d[0]->GetYaxis()->SetTitle(yAxisName.c_str());
  vecHist1d[0]->GetYaxis()->SetTitleOffset(1.1);
  // vecHist1d[0]->GetYaxis()->SetTitleOffset(0.8);  // was 1.03 without setting also the size
  vecHist1d[0]->GetYaxis()->SetTitleSize(0.05);
  //vecHist1d[0]->GetYaxis()->SetRangeUser(0.0, max(vecHist1d[0]->GetMaximum(),h2->GetMaximum()) * 1.2);

  //////////////////////////////
  // set X and Y axis range

  // search for maximum Y and for minimum > 0 (latter only if using log scale for Y axis
  Double_t maxY = -999.0;
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    if ( vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMaximumBin()) > maxY ) maxY = vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMaximumBin());
  }

  Double_t minY = 1e34;

  if (drawPlotLogY) {

    // quick check if there are no empty bins
    for (UInt_t i = 0; i < vecHist1d.size(); i++) {
      if ( vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMinimumBin()) < minY ) minY = vecHist1d[i]->GetBinContent(vecHist1d[i]->GetMinimumBin());
    }

    if (fabs(minY) < 0.00000001) {

      minY = 1e34;
      
      for (UInt_t ihist = 0; ihist < vecHist1d.size(); ihist++) {

	for (Int_t ibin = 0; ibin <= vecHist1d[ihist]->GetNbinsX(); ibin++ ) {
	  if (vecHist1d[ihist]->GetBinContent(ibin) > 0.0000001 && minY > vecHist1d[ihist]->GetBinContent(ibin)) minY = vecHist1d[ihist]->GetBinContent(ibin);
	}      
      
      }

    }

  }

  vecHist1d[0]->GetYaxis()->SetRangeUser(0.0, maxY * 1.2);

  if (setXAxisRangeFromUser) vecHist1d[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  //////////////////////

  vecHist1d[0]->Draw("Hist");
  vecHist1d[0]->SetFillColor(0);
  vecHist1d[0]->SetMarkerStyle(0);
  for (UInt_t i = 1; i < vecHist1d.size(); i++) {
    vecHist1d[i]->Draw("hist same");
  }

  TLegend leg (0.58,legLowY,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    leg.AddEntry(vecHist1d[i],vecLegEntry[i].c_str(),"L");
  }
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  TLegend legNevents(0.15,legLowY-0.03,0.55,0.93);
  if (yAxisName == "a.u.") {
    legNevents.SetFillColor(0);
    legNevents.SetFillStyle(0);
    legNevents.SetBorderSize(0);
    legNevents.SetHeader("# events");
    for (UInt_t i = 0; i < vecHist1d.size(); i++) {
      legNevents.AddEntry(vecHist1d[i],Form("%1.0f",integral_orig[i]),"L");
    }
    legNevents.Draw("same");
    canvas->RedrawAxis("sameaxis");
  }

  //  CMS_lumi(canvas,Form("%.1f",lumi));
  bool cmsPreliminaryIsUp = false;
  if (yAxisName == "a.u.") cmsPreliminaryIsUp = true;

  // FIXME
  /* TPaveText *pvtxt = NULL; */
  /* if (yAxisName == "a.u.") { */
  /*   pvtxt = new TPaveText(0.5,0.6,0.90,0.7, "BR NDC"); */
  /*   pvtxt->SetFillColor(0); */
  /*   pvtxt->SetFillStyle(0); */
  /*   pvtxt->SetBorderSize(0); */
  /*   pvtxt->AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError)); */
  /*   pvtxt->Draw(); */
  /* } */

  if (lumi < 0) CMS_lumi(canvas,"",cmsPreliminaryIsUp,false);
  else {
    if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),cmsPreliminaryIsUp,false);
    else CMS_lumi(canvas,Form("%.1f",lumi),cmsPreliminaryIsUp,false);
  }

  setTDRStyle();

  if (drawRatioWithNominal) {
    pad2->Draw();
    pad2->cd();
  
    frame->Reset("ICES");
    if (canvasName.find("comparisonMassVariation") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99, 1.01);
      /* if      (outputDIR.find("/eta_0/") != string::npos) frame->GetYaxis()->SetRangeUser(0.99, 1.01); */
      /* else if (outputDIR.find("/eta_1/") != string::npos) frame->GetYaxis()->SetRangeUser(0.98, 1.02); */
      /* else if (outputDIR.find("/eta_2/") != string::npos) frame->GetYaxis()->SetRangeUser(0.98, 1.02); */
    } else if (canvasName.find("elescale") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99,1.01);
    } else if (canvasName.find("elescale") != string::npos) {
      frame->GetYaxis()->SetRangeUser(0.99,1.01);
    } else frame->GetYaxis()->SetRangeUser(0.9,1.1);
    frame->GetYaxis()->SetNdivisions(5);
    frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
    frame->GetYaxis()->SetTitleOffset(1.2);
    // frame->GetYaxis()->SetTitleSize(0.15);
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitle(xAxisName.c_str());
    if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
    // frame->GetXaxis()->SetTitleOffset(0.8);
    frame->GetXaxis()->SetTitleSize(0.05);

    vector<TH1D*> ratio;
    for (UInt_t ivar = 1; ivar < vecHist1d.size(); ivar++) 
      ratio.push_back( (TH1D*) vecHist1d[ivar]->Clone(Form("ratio_%d",ivar)) );

    TH1D* den_noerr = (TH1D*) vecHist1d[0]->Clone("den_noerr");
    TH1D* den = (TH1D*) vecHist1d[0]->Clone("den");
    for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++)
      den_noerr->SetBinError(iBin,0.);

    den->Divide(den_noerr);
    den->SetFillColor(kGray);
    frame->Draw();
    den->Draw("E2same");
    for (UInt_t ir = 0; ir < ratio.size(); ir++) {
      ratio[ir]->Divide(den_noerr);
      // ratio[ir]->SetMarkerSize(0.65);
      // ratio[ir]->Draw("EPsame");
      ratio[ir]->SetMarkerStyle(0);
      ratio[ir]->SetLineWidth(2);
      ratio[ir]->Draw("Hist same");
    }
 

    TF1* line = new TF1("horiz_line","1",den->GetXaxis()->GetBinLowEdge(1),den->GetXaxis()->GetBinLowEdge(den->GetNbinsX()+1));
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("Lsame");
    // for (UInt_t ir = 0; ir < ratio.size(); ir++)
    //   ratio[ir]->Draw("EPsame");
    pad2->RedrawAxis("sameaxis");

  }  // end of ratio plot settings

  if (canvasName.find("tmpToBeRemoved") == string::npos) { 
    canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
    canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());
  }

  if (drawPlotLogY) {

    if (yAxisName == "a.u.") vecHist1d[0]->GetYaxis()->SetRangeUser(minY*0.05, maxY*100);
    else vecHist1d[0]->GetYaxis()->SetRangeUser(minY*0.05, maxY*100);
    canvas->SetLogy();
    /* if (lumi < 0) CMS_lumi(canvas,"",true,false); */
    /* else CMS_lumi(canvas,Form("%.1f",lumi),true,false); */
    if (canvasName.find("tmpToBeRemoved") == string::npos) { 
      canvas->SaveAs((outputDIR + canvasName + "_logY.png").c_str());
      canvas->SaveAs((outputDIR + canvasName + "_logY.pdf").c_str());
    }
    canvas->SetLogy(0);

  }    

  delete canvas;
  for (UInt_t i = 0; i < vecHist1d.size(); i++) {
    delete vecHist1d[i];
  }

}



//=============================================================


void drawCorrelationPlot(TH2* h2D, 
			 const string & labelX_input = "xaxis", const string & labelY_input = "yaxis", const string& labelZ_input = "zaxis", 
			 const string& canvasName = "default", const string& plotLabel = "", const string & outputDIR = "./", 
			 const Int_t rebinFactorY = 1,
			 const Bool_t smoothPlot = true,
			 const Bool_t plotProfileX = true)
{

  createPlotDirAndCopyPhp(outputDIR);

  if (rebinFactorY > 1) h2D->RebinY(rebinFactorY);
  string labelX = labelX_input;
  string labelY = labelY_input;
  string labelZ = labelZ_input;

  if (labelX_input.find("log::") != string::npos) labelX = labelX_input.substr(5,string::npos); //remove log: from title
  if (labelY_input.find("log::") != string::npos) labelY = labelY_input.substr(5,string::npos); //remove log: from title
  if (labelZ_input.find("log::") != string::npos) labelZ = labelZ_input.substr(5,string::npos); //remove log: from title
  
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  h2D->SetTitle(0);
  h2D->SetStats(0);

  system(("mkdir -p "+outputDIR).c_str());
  // normalize to 1
  canvas->SetRightMargin(0.18);

  if (labelZ == "a.u.") h2D->Scale(1./h2D->Integral());

  TGraph2D* h2DGraph = NULL;

  TH2* h2DPlot = NULL;
  if (not smoothPlot) h2DPlot = h2D;
  else {
    h2DGraph = new TGraph2D();
    h2DGraph->SetNpx(300);
    h2DGraph->SetNpy(300);
    int nPoint = 0;
    for(int iBinX = 0; iBinX < h2D->GetNbinsX() ; iBinX++){
      for(int iBinY = 0; iBinY < h2D->GetNbinsY() ; iBinY++){
	h2DGraph->SetPoint(nPoint,h2D->GetXaxis()->GetBinCenter(iBinX+1),h2D->GetYaxis()->GetBinCenter(iBinY+1),h2D->GetBinContent(iBinX+1,iBinY+1));      
	nPoint++;
      }
    }
    h2DPlot = h2DGraph->GetHistogram();
  }

  h2DPlot->GetXaxis()->SetTitle(labelX.c_str());
  h2DPlot->GetYaxis()->SetTitle(labelY.c_str());
  h2DPlot->GetZaxis()->SetTitle(labelZ.c_str());
  h2DPlot->Draw("colz");

  TProfile* h2DProfile = h2D->ProfileX(Form("%s_pfx",h2D->GetName()));
  h2DProfile->SetMarkerColor(kBlack);
  h2DProfile->SetMarkerStyle(20);
  h2DProfile->SetMarkerSize(1);
  if (plotProfileX) h2DProfile->Draw("EPsame");

  CMS_lumi(canvas,"",true,false);
  setTDRStyle();

  TLegend leg(0.4,0.6,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,plotLabel.c_str(),"");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",h2DPlot->GetCorrelationFactor()),"");
  leg.Draw("same");

  if (labelX_input.find("log::") != string::npos) canvas->SetLogx();
  if (labelY_input.find("log::") != string::npos) canvas->SetLogy();
  /* canvas->SaveAs((outputDIR + canvasName + ".png").c_str(),"png"); */
  /* canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str(),"pdf"); */
  if (labelZ_input.find("log::") != string::npos) {
    canvas->SetLogz();
    canvas->SaveAs((outputDIR + canvasName + "_logZ.png").c_str(),"png");
    canvas->SaveAs((outputDIR + canvasName + "_logZ.pdf").c_str(),"pdf");
    canvas->SetLogz(0);
  } else {
    canvas->SaveAs((outputDIR + canvasName + ".png").c_str(),"png");
    canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str(),"pdf");
  }

  delete canvas;

}



//======================================================

#endif
