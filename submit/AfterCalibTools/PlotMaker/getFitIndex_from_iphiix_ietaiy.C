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


using namespace std;

// this macro opens a root file with TH2D histograms, that can be used to get iphi(iX) and ieta(iy) in EB(EE) given the fit index used by ECALpro
// it can also be used to get a given fit index given the coordinates
// example -->  root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,12,0)'
// Returns the fit index associated to iphi=10, ieta=12 in EB

void getFitIndex_from_iphiix_ietaiy(const Int_t& iphiiy = 50, const Int_t& ietaix = 10, const Int_t iz = 0) {

  // iz = 0 for EB (default), -1,+1 for EE-,EE+

  string filename= "convert_fitIndex_iphiix_ietaiy.root";

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH2D* fitIndex_vs_ietaiphi_EB = (TH2D*) f->Get("fitIndex_vs_ietaiphi_EB");
  TH2D* fitIndex_vs_ixiy_EEp    = (TH2D*) f->Get("fitIndex_vs_ixiy_EEp");
  TH2D* fitIndex_vs_ixiy_EEm    = (TH2D*) f->Get("fitIndex_vs_ixiy_EEm");

  Int_t index = -1;
  if (iz == 0)     index = fitIndex_vs_ietaiphi_EB->GetBinContent(fitIndex_vs_ietaiphi_EB->FindFixBin(ietaix,iphiiy));
  else if (iz > 0) index = fitIndex_vs_ixiy_EEp->GetBinContent(fitIndex_vs_ixiy_EEp->FindFixBin(ietaix,iphiiy));
  else if (iz < 0) index = fitIndex_vs_ixiy_EEm->GetBinContent(fitIndex_vs_ixiy_EEm->FindFixBin(ietaix,iphiiy));

  cout << "Fit index is " << index << endl;

  f->Close();
  delete f;


}
