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
#include <TH2D.h>
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

int Xtal_ID[170][360]={0};
int Xtal_Ieta[61200]={0};
int Xtal_Iphi[61200]={0};

int Xtal_Ix[14648]={0};
int Xtal_Iy[14648]={0};
int Xtal_Iz[14648]={0};

// this macro creates a root file with TH2D histograms, that can be used to get iphi(iX) and ieta(iy) in EB(EE) given the fit index used by ECALpro
// it can also be used to get a given fit index given the coordinates

void convert_fitIndex_iphiix_ietaiy() {

  string filename= "convert_fitIndex_iphiix_ietaiy.root";

  TFile* f = TFile::Open(filename.c_str(),"RECREATE");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH2D* fitIndex_vs_ietaiphi_EB = new TH2D("fitIndex_vs_ietaiphi_EB","fit index in ECALpro given i#eta and i#phi",171,-85.5,85.5,360,0.5,360.5);
  TH2D* fitIndex_vs_ixiy_EEp = new TH2D("fitIndex_vs_ixiy_EEp","fit index in ECALpro given ix and iy in EE+",100,0.5,100.5,100,0.5,100.5);
  TH2D* fitIndex_vs_ixiy_EEm = new TH2D("fitIndex_vs_ixiy_EEm","fit index in ECALpro given ix and iy in EE-",100,0.5,100.5,100,0.5,100.5);


  for(int i = 0; i < 61200; i++) {

    int det_ID = EBDetId::detIdFromDenseIndex(i);
    EBDetId ebseed(det_ID);
    int ieta = ebseed.ieta();
    int iphi = ebseed.iphi();		
    fitIndex_vs_ietaiphi_EB->Fill(ieta,iphi,(Double_t)i);

  }
  
  // EE
  for(int i = 1; i <= fitIndex_vs_ixiy_EEp->GetNbinsX(); i++) {
    for(int j = 1; j <= fitIndex_vs_ixiy_EEp->GetNbinsX(); j++) {
      fitIndex_vs_ixiy_EEp->SetBinContent(i,j,-1.0);
      fitIndex_vs_ixiy_EEm->SetBinContent(i,j,-1.0);
    }
  }

  for(int i = 0; i < 14648; i++) {

    int det_ID = EEDetId::detIdFromDenseIndex(i);    
    EEDetId eeseed(det_ID);
    int ix = eeseed.ix();
    int iy = eeseed.iy();		       
    int iz = eeseed.zside();		
    if (iz > 0) fitIndex_vs_ixiy_EEp->Fill(ix,iy,(Double_t)i);
    else        fitIndex_vs_ixiy_EEm->Fill(ix,iy,(Double_t)i);

  }
  
  f->Write();
  f->Close();
  delete f;


}
