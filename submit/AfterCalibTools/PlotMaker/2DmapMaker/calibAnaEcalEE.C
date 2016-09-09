#define calibAnaEcalEE_cxx

#include "calibAnaEcalEE.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

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

#ifdef calibAnaEcalEE_cxx

calibAnaEcalEE::calibAnaEcalEE(TTree *tree) : calibAnaEcal(tree) {

  ////////////////////////////                                                                                                                                         
  //initializing data members                                                                                                                                          
  ///////////////////////////                                                                                                                                          
  
  //////////////////////////                                                                                                                                           
  // protected data members                                                                                                                                            

  EBorEE = "EB";

  ////////////////////////////////                                                                                                                                     
  // public data members                    

  // iX and iY are integer, but we create bins for histograms such that, e.g., the bin with (iX = 2 ) goes from 1.5 to 2.5       
  // note that 0 < iX <= 100, same for iY. 
  // However, we make maps with 102 bins on both x and y axis so that there is a 1 bin white margin in the plots (this is just a style choice
  NbinsX_2Dmap = 102;
  lowerX_2Dmap = -0.5;
  upperX_2Dmap = 101.5;
  NbinsY_2Dmap = 102;
  lowerY_2Dmap = -0.5;
  upperY_2Dmap = 101.5;

  // using eta ring for EE
  NbinsX_etaProfile = 40;
  lowerX_etaProfile = 0.0;
  upperX_etaProfile = 40.0;

  etaRing = -1;

  Init(tree);

}

//===============================================                                                                                                                      

void calibAnaEcalEE::setHistograms() {

  calibAnaEcal::setHistograms();

  // I was using these values starting for SoverB, SoverSqrtSplusB, SigmaMeanOverMean, mean,                                                                          
  // set only some of them with h->GetBinContent(h->GetMinimumBin()) and h->GetBinContent(h->GetMaximumBin())                                                          
  th2dMinZaxisVector.push_back(0.0);  // was set to hSignal->GetBinContent(hSignal->GetMinimumBin())                                                                  
  th2dMinZaxisVector.push_back(0.0);  // hBackground->GetBinContent(hBackground->GetMinimumBin())                                                                      
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0);//0.0                                                                                                                            
  th2dMinZaxisVector.push_back(0.13);
  th2dMinZaxisVector.push_back(0.005);

}

//===============================================                                                                                                                      

void calibAnaEcalEE::set2DmapMaxZaxisVector() {

  th2dMaxZaxisVector.push_back(hSignal->GetBinContent(hSignal->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(hBackground->GetBinContent(hBackground->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(10e9);
  th2dMaxZaxisVector.push_back(10e9); // when this value is very large (bigger than the default) use the default to plot axis                  
  th2dMaxZaxisVector.push_back(0.1);                                                
  th2dMaxZaxisVector.push_back(0.16);
  th2dMaxZaxisVector.push_back(0.020);

}

//===============================================                                                                                                                      

void calibAnaEcalEE::draw2Dmap(TH2D* hist2d) {

  calibAnaEcal::draw2Dmap(hist2d);

}

//===============================================                                                                                                                      

void calibAnaEcalEE::drawProfile(TProfile *profile, const string& yAxisName) {
 
  calibAnaEcal::drawProfile(profile, yAxisName);
   
}

//===============================================                                                                                                                      

Int_t calibAnaEcalEE::getEtaRingFromIxIyZside(const Int_t &ix, const Int_t &iy, const Int_t &zside) {

  Int_t thisEtaRing = -1;

  string fileName = "eerings_modified.dat";
  ifstream inputFile(fileName.c_str());

  // file format is --> a b c d                                                                                                                         
  Int_t a,b,c,d; 

  if (inputFile.is_open()) {

    while ((thisEtaRing == -1) && (inputFile >> a >> b >> c >> d)) {

      // first check zside, because it will remove half of the line to check                                                   
      if (c == zside && a == ix && b == iy) {  
        thisEtaRing = d;
      }

    }

  } else {
    std::cout << "Error in calibAnaEcalEE::getEtaRing(): could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();

  return thisEtaRing;

}


//===============================================                                                                                                                      

void calibAnaEcalEE::Init(TTree* tree) {

  calibAnaEcal::Init(tree);

  // EE                                                                                                                                                            
  fChain->SetBranchAddress("ix", &ix, &b_ix);                                                                                                                      
  fChain->SetBranchAddress("iy", &iy, &b_iy);                                                                                                                      
  fChain->SetBranchAddress("zside", &zside, &b_zside);                                                                                                             
  fChain->SetBranchAddress("sc", &sc, &b_sc);                                                                                                                      
  fChain->SetBranchAddress("isc", &isc, &b_isc);                                                                                                                   
  fChain->SetBranchAddress("ic", &ic, &b_ic);                                                                                                                      
  fChain->SetBranchAddress("iquadrant", &iquadrant, &b_iquadrant);  
  Notify();    

}

//===============================================                                                                                                                      

void calibAnaEcalEE::Loop()
{  

  if (fChain == 0) return;

  setHistograms();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;                                                                                                                              
      
    if (jentry % 50000 == 0) cout << jentry << endl;

    if (EBorEE == "EEp") {
      if (zside < 0) continue;
    } else {
      if (zside > 0) continue;
    }

    if ((abs(Backgr) > 0.00001) && (abs(Signal) > 0.00001)) { // avoid empty crystals due to masked towers or whatever                                                 

      normalizedS = Signal * fit_Snorm;
      normalizedB = Backgr * fit_Bnorm;

      etaRing = getEtaRingFromIxIyZside(ix,iy,zside);  // will return -1 if eta ring is not found for any reason 

      // to avoid that in 2D maps points below lower threshold in z axis are drawn white (as if they are empty), fill with the maximum between threshold and value     
      hSignal->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[0],(Double_t)normalizedS));
      hBackground->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[1],(Double_t)normalizedB));
      SoverB->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[2],(Double_t)normalizedS/normalizedB));
      SoverSqrtSplusB->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[3],(Double_t)normalizedS/sqrt(normalizedS + normalizedB)));
      SigmaMeanOverMean->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[4],(Double_t)fit_mean_err/fit_mean));
      mean->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[5],(Double_t)fit_mean));
      sigma->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[6],(Double_t)fit_sigma));

      hSignal_etaProfile->Fill((Double_t)etaRing,normalizedS);
      hBackground_etaProfile->Fill((Double_t)etaRing,normalizedB);
      SoverB_etaProfile->Fill((Double_t)etaRing,normalizedS/normalizedB);
      SoverSqrtSplusB_etaProfile->Fill((Double_t)etaRing, normalizedS/sqrt(normalizedS + normalizedB));
      SigmaMeanOverMean_etaProfile->Fill((Double_t)etaRing, fit_mean_err/fit_mean);
      mean_etaProfile->Fill((Double_t)etaRing, fit_mean);
      sigma_etaProfile->Fill((Double_t)etaRing, fit_sigma);

    }

  }

  // set preference for max value in the vertical scale
  set2DmapMaxZaxisVector();  
  // now set the vertical axis maximum value based on user input (will use the least between the default and the user input). 
  setVerticalRangeInHisto();

  for ( UInt_t i = 0; i < th2dVector.size(); i++ ) {

    draw2Dmap(th2dVector[i]);
    drawProfile(profileEtaVector[i], profileYaxisTitle[i]);

  }

}

#endif
