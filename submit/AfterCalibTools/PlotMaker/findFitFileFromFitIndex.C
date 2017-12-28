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

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

using namespace std;
using namespace RooFit;

// example: 
//
//root -l -b -q 'findFitFileFromFitIndex.C+(5000,"/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/","AlcaP0_Run2016G_sel17optim_reg12",0,true)'
//
// returns number of fit file that contains fit with index 5000. This can be used with drawFitsSingleFile.sh
// It looks inside /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlcaP0_Run2016G_sel17optim_reg12/iter_0/ and search for fits in EB

void findFitFileFromFitIndex(const Int_t& fitIndex = 5000,
			     const string& eosPath = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/",
			     const string& dirName = "AlcaP0_Run2016G_sel17optim_reg12", 
			     const Int_t&  iterNumber = 0,
			     const Bool_t& isEB = true 
			     ) 
{

  Int_t nMaxFile = isEB ? 30 : 7; // fitRes files goes from 0 to 30 (7) in EB (EE)

  string BarrelOrEndcap = isEB ? "Barrel" : "Endcap";

  for (Int_t iFile = 0; iFile <= nMaxFile; iFile++ ) {

    string fitResFileOnEos = eosPath + dirName + Form("/iter_%d/",iterNumber) + dirName + Form("_%s_%d_fitRes.root",BarrelOrEndcap.c_str(),iFile);
    TFile* f = TFile::Open(fitResFileOnEos.c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<fitResFileOnEos<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }

    string fit_rooplot = Form("Fit_n_%d_attempt0_rp",fitIndex);
    Bool_t fitFound = f->GetListOfKeys()->Contains(fit_rooplot.c_str());

    f->Close();
    delete f;

    if (fitFound) {
      cout << "Fit " << fit_rooplot << " found in file " << fitResFileOnEos << endl;
      break;
    }

  }

  cout << "Sorry, I didn't find the fit with index " << fitIndex << " anywhere." << endl;

}
