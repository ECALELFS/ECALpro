#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h    
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++    

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

#include "./CMS_lumi.h"

using namespace std;
using namespace RooFit;


void drawRooPlotFromFile(const string& inputDir = "", const bool isEB = true, const string &inputFileName = "") {

  double lumi = 0.18; //in fb-1

  TFile* f = TFile::Open((inputDir+inputFileName).c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< inputDir+inputFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  string plotInFile = isEB ? "Fit_n_30003" : "Fit_n_6397";

  RooPlot * xframe = (RooPlot*) f->Get(plotInFile.c_str());
  if (!xframe) {
    cout << "Warning: RooPlot object not found in file";
    cout << inputFileName;
    cout << "Abort" << endl;
    exit(EXIT_FAILURE);
  }

  string canvasname = isEB ? "pi0massEBxtal" : "pi0massEExtal";
  TCanvas *canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  // canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
  xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");

  xframe->GetXaxis()->SetLabelSize(0.04);
  xframe->GetXaxis()->SetTitleSize(0.05);
  xframe->GetYaxis()->SetTitleOffset(1.1);
  xframe->GetYaxis()->SetTitleSize(0.05);
  xframe->Draw();

  CMS_lumi(canvas,Form("%.2f",lumi),true,false);
  setTDRStyle();
  
  canvas->SaveAs((inputDir + canvasname + ".pdf").c_str());
  canvas->SaveAs((inputDir + canvasname + ".png").c_str());

  delete canvas;

  f->Close();
  delete f;

}



void makeRooPlotFromFile() { 

  string dirName = "AlCaP0_Run2017A_runs296966to296980_v2_EB7xtal_EE2xtal";
  string inputDirEB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_0/fitResPlots/Barrel/";
  string inputDirEE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_0/fitResPlots/Endcap/";

  string inputFileName = "pi0Mass_singleXtal_Rooplot.root";

  // it seems that the first time CMS_lumi is used the settings are screwed up
  // produce a dummy plot (either do not save it or remove it)
  double lumi = 0.18; //in fb-1
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,Form("%.2f",lumi),false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;

  // here we go with the real part
  drawRooPlotFromFile(inputDirEB, true, inputFileName);
  drawRooPlotFromFile(inputDirEE, false, inputFileName);

}
