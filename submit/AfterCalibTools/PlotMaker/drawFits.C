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

using namespace std;
using namespace RooFit;

void createFits(const string& Barrel_or_Endcap = "Barrel", 
		const Int_t& min_fitIndex = 0, 
		const Int_t& max_fitIndex = 100, 
		const string& wwwPath = "", 
		const string& eosPath = "", 
		const string& dirName = "", 
		const string& iterNumber = "", 
		const string& tagName = "") {

  string storePath = wwwPath + dirName + "/" + iterNumber + "/fitResPlots/";
  string filename = "root://eoscms//eos/cms" + eosPath + dirName + "/" + iterNumber + "/" + tagName + Barrel_or_Endcap + "_fitRes.root";
  storePath += (Barrel_or_Endcap + "/");

  string rooplotbasename = "Fit_n_";

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }


  TCanvas *c = NULL;

  for (Int_t i = min_fitIndex; i < max_fitIndex; i++) {

    ostringstream convert;   // stream used for the conversion                                    
    convert << i;      // insert the textual representation of 'i' in the characters in the stream                                          
    string Result = convert.str();  

    string rooplotname = rooplotbasename + Result;
    c = new TCanvas("c",rooplotname.c_str());
    string canvasname = rooplotname + ".png";
    
    RooPlot * xframe = (RooPlot*) f->Get(rooplotname.c_str());
    if (!xframe) {
      cout << "Warning: RooPlot object named \"" << rooplotname << "\" not found in file." <<endl;
    } else {
      xframe->Draw();
      c->SaveAs((storePath + canvasname).c_str());
    }
    delete c;
    
  }


}


void drawFits(string wwwPath, string eosPath, string dirName, string iterNumber, string tagName) {

  // examples of arguments
  // string wwwPath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/";
  // string eosPath = "/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/";
  // string dirName = "AlcaP0_fromRun273158_2016_v2";
  // string iterNumber = "iter_6";
  // string tagName = "AlcaP0_fromRun273158_2016_v2_";

  cout << "Producing fits for Barrel." << endl;
  createFits("Barrel",48300,48350,wwwPath,eosPath,dirName,iterNumber,tagName);
  cout << endl;
  cout << endl;
  cout << "Producing fits for Endcap." << endl;
  createFits("Endcap",400,450,wwwPath,eosPath,dirName,iterNumber,tagName);
  

}
