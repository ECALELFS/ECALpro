#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2.h>
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

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                                                                                
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++                                       

using namespace std;

// use as:
// root -b -q 'convert_ix_iy_eta.C+("<filename>.dat","<filename>.root")'

void convert_ix_iy_eta_to_TH2(const string datfileName = "../../../common/geometry_ietaix_iphiiy_0iz_eta.dat", const string rootfileName = "../../../common/geometry_encap_ix_iy_iz_eta.root") {
  
  ifstream inputFile(datfileName.c_str());

  Int_t NbinsX_2Dmap = 100;
  Double_t lowerX_2Dmap = 0.5;
  Double_t upperX_2Dmap = 100.5;
  Int_t NbinsY_2Dmap = 100;
  Double_t lowerY_2Dmap = 0.5;
  Double_t upperY_2Dmap = 100.5;

  // open file to store histogram
  TFile *rootFile = new TFile((rootfileName).c_str(),"RECREATE");
  if (!rootFile || !rootFile->IsOpen()) {
    cout << "Error: file \"" << rootfileName << "\" was not opened." << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *hEEp = new TH2F("hEEp_eta","cystals #eta map in EE+",NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  TH2F *hEEm = new TH2F("hEEm_eta","cystals #eta map in EE-",NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);

  // fill histograms with some value that will not be used for eta
  for (Int_t ix = 1; ix <= 100; ix++) {
    for (Int_t iy = 1; iy <= 100; iy++) {
      hEEp->SetBinContent(ix, iy, -1);
      hEEm->SetBinContent(ix, iy, -1);
    }
  }
  
  // file.dat format is --> a b c d, that is --> iX, iY, Z side, eta                 
  Int_t a,b, c;
  Double_t d;  

  if (inputFile.is_open()) {

    while (inputFile >> a >> b >> c >> d) {

      if (c > 0)      hEEp->SetBinContent(a,b,d);  
      else if (c < 0) hEEm->SetBinContent(a,b,d);

    }

  } else {
    std::cout << "Error: could not open file " << datfileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();

  rootFile->Write();
  rootFile->Close();
  delete rootFile;

}
