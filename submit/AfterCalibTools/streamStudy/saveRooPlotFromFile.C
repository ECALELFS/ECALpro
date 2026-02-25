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

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

using namespace std;
using namespace RooFit;

void getRooplotWithIndex(const string& fitResFileOnEos = "", 
			 const bool isEB = true, 
			 const string& outputDIR = "./", 
			 const Int_t rooplotIndex = 1, 
			 const string& outputFileName = ""
			 ) {

  TFile* f = TFile::Open(fitResFileOnEos.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fitResFileOnEos<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  // stringstream ss_index;
  // ss_index << rooplotIndex;
  // string str_index = ss_index.str();

  RooPlot * xframe = (RooPlot*) f->Get(Form("Fit_n_%d",rooplotIndex));
  if (!xframe) {
    cout << "Warning: RooPlot object not found in file. Skipping and going on with next object" <<endl;
    exit(EXIT_FAILURE);
  }

  if (xframe) {
    xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
    xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  }

  TFile* outputFile = new TFile((outputDIR + outputFileName).c_str(),"UPDATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  outputFile->cd();
  
  xframe->Write();
  outputFile->Close();
  delete outputFile;

  f->Close();
  delete f;

}


void saveRooPlotFromFile() {

  string dirName = "AlCaP0_Run2017A_runs296966to296980_v2_EB7xtal_EE2xtal";

  string outputDIR_EB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_0/fitResPlots/Barrel/";
  string outputDIR_EE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_0/fitResPlots/Endcap/";
  string inputDIR_EB = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/" + dirName + "/iter_0/" + dirName + "_Barrel_15_fitRes.root";
  string inputDIR_EE = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/" + dirName + "/iter_0/" + dirName + "_Endcap_3_fitRes.root";

  string outputFileName = "pi0Mass_singleXtal_Rooplot.root";

  getRooplotWithIndex(inputDIR_EB, true, outputDIR_EB, 30003, outputFileName);
  getRooplotWithIndex(inputDIR_EE, false, outputDIR_EE, 6397, outputFileName);

}
