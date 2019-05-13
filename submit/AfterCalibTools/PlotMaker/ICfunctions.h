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
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
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
#include <fstream>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#define EPSILON 0.000001  // used to assess whether a float or double is equal to or close enough to another value. The tolerance depends on the numbers you manage

using namespace std;

//====================================================

int getArrayIndexOfFoldedSMfromIetaIphi(const int ieta = 1, const int iphi = 1) {

  // note that the index in SM returned by this function in not the same as the index returned by EBDetId::ic()                 
  // the difference is mainly in the folding of EB+ on EB-                                                                                      
  // In our case, we overlay crystals such that ieta,iphi=20,40 goes on ieta,iphi=-20,40,                                                          
  // i.e. the iphi coordinate is preserved when we consider two facing SM in EB+ and EB-                                                        
  // The usual CMSSW numbering scheme for crystals in SM is such that, looking at the center of the barrel, the crystal number 1 is always on the left                      
  // which means that the folding would overlay ieta,iphi=20,40 on ieta,iphi=-20,21                                                                                         
  // first 85 crystals correspond to iphi = 1 (in a SM)                                                                                          
  return (fabs(ieta) - 1) + EBDetId::kCrystalsInEta * ((iphi - 1) % EBDetId::kCrystalsInPhi);

}


//====================================================

Bool_t isGoodIC(const Double_t val) {

  if (val > EPSILON && fabs(val -1.0) > EPSILON) return true;
  else return false;

}

//====================================================


void makeICdistributionFromMap(TH1* h, const TH2* mapEB = NULL, const Bool_t noBadXtals = true) {

  Double_t binContent = 0.0;

  for (Int_t x = 1; x <= mapEB->GetNbinsX(); ++x) {

    for (Int_t y = 1; y <= mapEB->GetNbinsY(); ++y) {      

      binContent = mapEB->GetBinContent(x,y);
      if (noBadXtals) {
	if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) h->Fill(binContent);
      } else 
	h->Fill(binContent);

    }

  }

}

//================================================

void makeICprofileIphiFromMap(TProfile* h, const TH2* mapEB = NULL, const Bool_t iphiOnXaxis = true, const Bool_t noBadXtals = true, 
			      const Bool_t justOneSM = true, Int_t imodule = 0) {

  Double_t binContent = 0.0;
  Int_t maxIetaBin = justOneSM ? 85 : 171; 
  Int_t maxIphiBin = justOneSM ? 20 : 360;

  Int_t absieta = -99;

  for (Int_t x = 1; x <= maxIphiBin; ++x) {

    for (Int_t y = 1; y <= maxIetaBin; ++y) {      

      absieta = justOneSM ? y : fabs(y-86);
      if (absieta == 0) continue;

      if      (imodule == 1 &&  absieta > 25 ) continue;	
      else if (imodule == 2 && (absieta < 26 || absieta > 45)) continue;	
      else if (imodule == 3 && (absieta < 46 || absieta > 65)) continue;	
      else if (imodule == 4 &&  absieta < 66) continue;	

      if (iphiOnXaxis) binContent = mapEB->GetBinContent(x,y);
      else             binContent = mapEB->GetBinContent(y,x);

      if (noBadXtals) {
	if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) {
	  // get abscissa corresponding to bin x (likely the two corresponds if using iphi)
	  h->Fill(h->GetBinCenter(x),binContent);
	}
      } else {
	h->Fill(h->GetBinCenter(x),binContent);
      }
    }

  }

}

//================================================

 
void copyMapAllEB(TH2* hnew = NULL, const TH2* hold = NULL, const Bool_t iphiOnXaxis = true) {

  // hnew and hold must have been defined already, and have the same binning on x and y (we check it explicitly)

  if (!hnew) {
    cout << "Warning in copyMapAllEB(): hnew is NULL. Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (!hold) {
    cout << "Warning in copyMapAllEB(): hold is NULL. Exit." << endl;
    exit(EXIT_FAILURE);
  }

  Double_t binContent = 0.0;
  Int_t bin = 0;

  Int_t hnew_nBinsX = hnew->GetNbinsX();
  Int_t hnew_nBinsY = hnew->GetNbinsY();
  Int_t hold_nBinsX = hnew->GetNbinsX();
  Int_t hold_nBinsY = hnew->GetNbinsY();

  if (hnew_nBinsX != hold_nBinsX) {
    cout << "Warning in copyMapAllEB(): hnew and hold have different number of bins in X ("<< hnew_nBinsX << "," << hold_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (hnew_nBinsY != hold_nBinsY) {
    cout << "Warning in copyMapAllEB(): hnew and hold have different number of bins in Y ("<< hnew_nBinsY << "," << hold_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;

      if (iphiOnXaxis) bin = hold->FindFixBin(iphi,ieta);
      else             bin = hold->FindFixBin(ieta,iphi);
      binContent = hold->GetBinContent( bin );
      hnew->SetBinContent(bin, binContent );

    }
    
  }


}

//================================================

 
void copyMapEE(TH2* hnew = NULL, const TH2* hold = NULL) {

  // hnew and hold must have been defined already, and have the same binning on x and y (we check it explicitly)

  if (!hnew) {
    cout << "Warning in copyMapEE(): hnew is NULL. Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (!hold) {
    cout << "Warning in copyMapEE(): hold is NULL. Exit." << endl;
    exit(EXIT_FAILURE);
  }

  Double_t binContent = 0.0;
  Int_t bin = 0;

  Int_t hnew_nBinsX = hnew->GetNbinsX();
  Int_t hnew_nBinsY = hnew->GetNbinsY();
  Int_t hold_nBinsX = hnew->GetNbinsX();
  Int_t hold_nBinsY = hnew->GetNbinsY();

  if (hnew_nBinsX != hold_nBinsX) {
    cout << "Warning in copyMapEE(): hnew and hold have different number of bins in X ("<< hnew_nBinsX << "," << hold_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (hnew_nBinsY != hold_nBinsY) {
    cout << "Warning in copyMapEE(): hnew and hold have different number of bins in Y ("<< hnew_nBinsY << "," << hold_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }

  for (Int_t ix = 1; ix <= 100; ix++) {

    for (Int_t iy = 1; iy <= 100; iy++) {      

      bin = hold->FindFixBin(ix,iy);
      binContent = hold->GetBinContent( bin );
      hnew->SetBinContent(bin, binContent );

    }
    
  }


}

//==================================================

void checkICnormalizedTo1_inEtaRing(TH2* h = NULL,  const Bool_t isEB = true, const Bool_t iphiOnXaxis = true, const Bool_t printOnlyBadMean = true) {

  cout << "####################################" << endl;
  cout << "### Checking that IC are normalized to 1 in etaring" << endl;
  cout << "### Histogram: " << h->GetName() << endl;
  cout << "####################################" << endl;

  Int_t bin = 0;
  Double_t binContent = 0.0;
  Double_t mean = 0.0;
  Double_t nGoodIC = 0;

  if (isEB) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      mean = 0.0;
      nGoodIC = 0;
      
      for (Int_t iphi = 1; iphi <= 360; iphi++) {

	if (ieta == 0) continue;
	
	if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
	else             bin = h->FindFixBin(ieta,iphi);
	binContent = h->GetBinContent(bin);
	if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) {
	  mean += binContent;
	  nGoodIC += 1.0;
	}

      }

      mean /= nGoodIC;
      if (printOnlyBadMean && fabs(mean-1.0) > EPSILON) cout << "ieta: " << ieta << " --> mean = " << setprecision(4) << mean << endl;
      else if (not printOnlyBadMean)                    cout << "ieta: " << ieta << " --> mean = " << setprecision(4) << mean << endl;
    }

  }


}

//==================================================

void normalizeEBMapTo1_inEtaRing(TH2* h = NULL, const Bool_t iphiOnXaxis = true, const Bool_t excludeMod2EBm16 = false) {

  vector<Double_t> meanInEtaRing(171, 0.0); // 171 eta rings (we keep ieta = 0 for simplicity)
  vector<Double_t> nGoodXtalsInEtaRing(171, 0.0);

  Double_t binContent = 0.0;
  Int_t bin = 0;
  Int_t arrayIndex = -1;

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;

      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);
      
      arrayIndex = ieta + 85;
      binContent = h->GetBinContent(bin);
      if (excludeMod2EBm16) {
	if (iphi > 300 && iphi <= 320 && ieta < -25 && ieta >= -45) binContent = 1;
      }

      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) {
	nGoodXtalsInEtaRing[arrayIndex]++;
	meanInEtaRing[arrayIndex] += binContent;
      }

    }

  }

  // now loop again and fill histogram with normalized value 
  // we divide by the mean (which was not divided by number of entries yet)

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;

      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);
      
      arrayIndex = ieta + 85;
      binContent = h->GetBinContent(bin);

      // do not change IC at 1, because they belong to crystals that should not be updated (bad fit or whatever)
      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) 
	h->SetBinContent(bin, binContent * nGoodXtalsInEtaRing[arrayIndex] / meanInEtaRing[arrayIndex]);

    }

  }  

}


//==================================================

void loadMapFromFile(TH2* h, const string& inputFile = "", const string& mapNameInFile = "") {


  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error in loadMapFromFile() when opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2 *map = NULL;  // I expect it to have iphi on X and ieta on Y (IC map from DrawIC.py, not from TH2F in *_calibMap.root)

  map = (TH2*) f->Get(mapNameInFile.c_str());
  if (!map) {
    cout << "Error in loadMapFromFile(): could not get '" << mapNameInFile << "' histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  Int_t bin = 0;

  for (Int_t i = 1; i <= map->GetNbinsX(); ++i) {

    for (Int_t j = 1; j <= map->GetNbinsY(); ++j) {

      bin = h->FindFixBin(map->GetXaxis()->GetBinCenter(i),map->GetYaxis()->GetBinCenter(j));
      h->SetBinContent(bin,map->GetBinContent(i,j));

    }

  }

  f->Close();

}

//==================================================


void normalizeEEMapTo1_inEtaRing(TH2* h = NULL) {

  // eta-ring number goes from 0 to 38 included (0 at lower eta)
  // some people uses convention from 1 to 39, but it doesn't matter for us
  // the histogram content is -1 for bins not associated to any crystal (if you draw it, you would see them white if the Z scale goes from 0 to 38)
  // the etaRing is not dependent on EE+ or EE-, just on iX and iY (etaring is like the radius of EE)

  TH2* map_EE_ixiyToIetaRing = new TH2F("map_EE_ixiyToIetaRing","Conversion map (iX,iY)-->#eta-ring",100,0.5,100.5,100,0.5,100.5);
  loadMapFromFile(map_EE_ixiyToIetaRing,"/afs/cern.ch/user/m/mciprian/public/ECALproTools/EE_xyzToEtaRing/eerings_modified.root","hEEp");

  vector<Double_t> meanInEtaRing(39, 0.0); // 39 eta rings
  vector<Double_t> nGoodXtalsInEtaRing(39, 0.0);

  Double_t binContent = 0.0;
  Int_t bin = 0;
  Int_t arrayIndex = -1;

  for (Int_t ix = 1; ix <= 100; ix++) {

    for (Int_t iy = 1; iy <= 100; iy++) {

      bin = h->FindFixBin(ix,iy);
      arrayIndex = map_EE_ixiyToIetaRing->GetBinContent(ix,iy);  // already goes from 0 to 38
      binContent = h->GetBinContent(bin);

      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) {
	nGoodXtalsInEtaRing[arrayIndex]++;
	meanInEtaRing[arrayIndex] += binContent;
      }

    }

  }

  // now loop again and fill histogram with normalized value
  // we divide by the mean (which was not divided by number of entries yet)

  for (Int_t ix = 1; ix <= 100; ix++) {

    for (Int_t iy = 1; iy <= 100; iy++) {
      
      bin = h->FindFixBin(ix,iy);
      arrayIndex = map_EE_ixiyToIetaRing->GetBinContent(ix,iy);;
      binContent = h->GetBinContent(bin);
      // do not change IC at 1, because they belong to crystals that should not be updated (bad fit or whatever)
      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) h->SetBinContent(bin, binContent * nGoodXtalsInEtaRing[arrayIndex] / meanInEtaRing[arrayIndex]);

    }

  }

}



//==================================================

void normalizeEBMapTo1_inEachModule(TH2* h = NULL, const Bool_t iphiOnXaxis = true, const Bool_t excludeMod2EBm16 = false) {

  vector<Double_t> meanInModule(36*4,0.0); // 36 SM * 4 modules each, initialized to 0
  vector<Double_t> nGoodXtalsInModule(36*4,0.0); // 36 SM * 4 modules each, initialized to 0

  Double_t binContent = 0.0;
  Int_t bin = 0;
  Int_t arrayIndex = -1;

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;
      EBDetId id(ieta,iphi);
      Int_t ism = id.ism();
      Int_t im  = id.im();

      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);
      
      arrayIndex = 4 * (ism-1) + im-1; // ind = 0 for ism = im = 0 and so on (first 4 elements are 4 modules of SM_1, and so on)
      binContent = h->GetBinContent(bin);
      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) {
	nGoodXtalsInModule[arrayIndex]++;
	meanInModule[arrayIndex] += binContent;
      }

    }

  }

  // now loop again and fill histogram with normalized value 
  // we divide by the mean (which was not divided by number of entries yet)

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;
      EBDetId id(ieta,iphi);
      Int_t ism = id.ism();
      Int_t im  = id.im();

      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);
      
      arrayIndex = 4 * (ism-1) + im-1; // ind = 0 for ism = im = 0 and so on (first 4 elements are 4 modules of SM_1, and so on)
      binContent = h->GetBinContent(bin);

      // do not change IC at 1, because they belong to crystals that should not be updated (bad fit or whatever)
      if (binContent > EPSILON && fabs(binContent -1.0) > EPSILON) h->SetBinContent(bin, binContent * nGoodXtalsInModule[arrayIndex] / meanInModule[arrayIndex]);

      // in this case, if we are in the 2nd module of EB-16, just put 1 everywhere
      if (excludeMod2EBm16) {
	if (iphi > 300 && iphi <= 320 && ieta < -25 && ieta >= -45) {
	  h->SetBinContent(bin,1.);
	  cout << ">>> WARNING in normalizeEBMapTo1_inEachModule(): Setting IC of mod2 in EB-16 to 1" << endl;
	}
      }


    }

  }  

}

//==================================================

//take h, fold in SM averaging content of each crystal in the 36 SM and repeat the structure in all EB
// keep also the map in one SM to better visualize it
// all the maps must have been created outside the function

void foldEBMapInSM(const TH2* h = NULL, TH2* hSM = NULL, TH2* hSM_allEB = NULL, const Bool_t iphiOnXaxis = true, const Bool_t iphiOnXaxisSM = false, const Bool_t useEBDetId_ic_scheme = true) {


  vector<Int_t> nGoodXtalsInSM(1700,0); // 1 cell for each crystal in SM, to keep track of how many time we are filling that cell, index is given by EBDetId::ic()
  vector<Double_t> meanInSM(1700,0.0);

  Double_t binContent = 0.0;
  Int_t bin = 0;
  Int_t arrayIndex = -1;

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;
      EBDetId id(ieta,iphi);	
      /* Int_t iphi_foldSM = id.iphiSM(); */
      /* Int_t ieta_foldSM = id.ietaSM(); */
      Int_t ic_foldSM = id.ic();  // goes from 1 to 1700, but the array index goes from 0 to 1699, remember to subtract 1

      Int_t foldIndex = useEBDetId_ic_scheme ? (ic_foldSM -1) : getArrayIndexOfFoldedSMfromIetaIphi(ieta,iphi);

      // avoid empty cells (dead channels)
      // can keep crystals with IC = 1, they do not change the average
      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);

      binContent = h->GetBinContent( bin );

      if ( binContent > EPSILON) {
	nGoodXtalsInSM[foldIndex]++;
	meanInSM[foldIndex] += binContent;
	//cout << Form("ieta,iphi,ic,IC = %d %d %d %d %.3g",ieta,iphi,ic_foldSM,binContent) << endl; 
      } // else {
      // 	badICxtalFile << Form("ieta,iphi,ic,IC = %d %d %d %.3g",ieta,iphi,ic_foldSM,binContent) << endl; 
      // }

    }
      
  }

  // now fill each cell's content with the mean of the ~36 values after folding the SM
  // if there are no dead cells, this number is 36 (36 SM)

  for (Int_t iphi = 1; iphi <= 20; iphi++) {

    for (Int_t ieta = 1; ieta <= 85; ieta++) {      

      EBDetId id(ieta,iphi);	
      Int_t ic_foldSM = id.ic();  // goes from 1 to 1700, but the array index goes from 0 to 1699, remember to subtract 1
      Int_t foldIndex = useEBDetId_ic_scheme ? (ic_foldSM -1) : getArrayIndexOfFoldedSMfromIetaIphi(ieta,iphi);
      if (iphiOnXaxisSM) hSM->SetBinContent(iphi, ieta, meanInSM[foldIndex] / nGoodXtalsInSM[foldIndex] );
      else               hSM->SetBinContent(ieta, iphi, meanInSM[foldIndex] / nGoodXtalsInSM[foldIndex]);
      //cout << Form("ieta,iphi,ic,nEntries,IC = %d %d %d %d %.3g",ieta,iphi,ic_foldSM,nGoodXtalsInSM[foldIndex], mapEB_binContent) << endl; 

    }

  }

  // now we fill the full EB map repeating the SM
  // note that the single SM was chosen as EB+1, so iphi,ieta = 1,1 for the folded map should go into the same crystal for full EB
  // this means I should not use iphiSM() when using useEBDetId_ic_scheme, but rather (21 - id.iphiSM()) or equivalently (((iphi -1) % EBDetId::kCrystalsInPhi) +1)
  if (hSM_allEB != NULL) {

    for (Int_t iphi = 1; iphi <= 360; iphi++) {

      for (Int_t ieta = -85; ieta <= 85; ieta++) {      

	if (ieta == 0) continue;
	EBDetId id(ieta,iphi);	
	// Int_t iphi_foldSM = useEBDetId_ic_scheme ? id.iphiSM() : (((iphi -1) % EBDetId::kCrystalsInPhi) +1);
	Int_t iphi_foldSM = ((iphi -1) % EBDetId::kCrystalsInPhi) +1;
	Int_t ieta_foldSM = id.ietaSM();

	if (iphiOnXaxisSM) binContent = hSM->GetBinContent(iphi_foldSM, ieta_foldSM);
	else               binContent = hSM->GetBinContent(ieta_foldSM, iphi_foldSM);

	if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
	else             bin = h->FindFixBin(ieta,iphi);

	hSM_allEB->SetBinContent(bin, binContent);

      }

    }

  }

}


//==================================================


void foldEBMapInSM_plusMinusSeparate(const TH2* h = NULL, TH2* hSM_allEB = NULL, const Bool_t iphiOnXaxis = true, const Bool_t useEBDetId_ic_scheme = true) {


  vector<Int_t> nGoodXtalsInSM_plus(1700,0); // 1 cell for each crystal in SM, to keep track of how many time we are filling that cell, index is given by EBDetId::ic()
  vector<Double_t> meanInSM_plus(1700,0.0);

  vector<Int_t> nGoodXtalsInSM_minus(1700,0); // 1 cell for each crystal in SM, to keep track of how many time we are filling that cell, index is given by EBDetId::ic()
  vector<Double_t> meanInSM_minus(1700,0.0);

  Double_t binContent = 0.0;
  Int_t bin = 0;
  Int_t arrayIndex = -1;

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;
      EBDetId id(ieta,iphi);	
      /* Int_t iphi_foldSM = id.iphiSM(); */
      /* Int_t ieta_foldSM = id.ietaSM(); */
      Int_t ic_foldSM = id.ic();  // goes from 1 to 1700, but the array index goes from 0 to 1699, remember to subtract 1

      Int_t foldIndex = useEBDetId_ic_scheme ? (ic_foldSM -1) : getArrayIndexOfFoldedSMfromIetaIphi(ieta,iphi);

      // avoid empty cells (dead channels)
      // can keep crystals with IC = 1, they do not change the average
      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);

      binContent = h->GetBinContent( bin );

      if ( binContent > EPSILON) {
	if (ieta > 0) {
	  nGoodXtalsInSM_plus[foldIndex]++;
	  meanInSM_plus[foldIndex] += binContent;
	} else {
	  nGoodXtalsInSM_minus[foldIndex]++;
	  meanInSM_minus[foldIndex] += binContent;
	}
	//cout << Form("ieta,iphi,ic,IC = %d %d %d %d %.3g",ieta,iphi,ic_foldSM,binContent) << endl; 
      } // else {
      // 	badICxtalFile << Form("ieta,iphi,ic,IC = %d %d %d %.3g",ieta,iphi,ic_foldSM,binContent) << endl; 
      // }

    }
      
  }

  for (UInt_t i = 0; i < nGoodXtalsInSM_plus.size(); ++i) {
    meanInSM_plus[i] /= nGoodXtalsInSM_plus[i];
    meanInSM_minus[i] /= nGoodXtalsInSM_minus[i];
  }


  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {      

      if (ieta == 0) continue;
      EBDetId id(ieta,iphi);	
      Int_t ic_foldSM = id.ic();  // goes from 1 to 1700, but the array index goes from 0 to 1699, remember to subtract 1
      Int_t foldIndex = useEBDetId_ic_scheme ? (ic_foldSM -1) : getArrayIndexOfFoldedSMfromIetaIphi(ieta,iphi);

      if (ieta < 0) binContent = meanInSM_minus[foldIndex];
      else          binContent = meanInSM_plus[foldIndex];

      if (iphiOnXaxis) bin = h->FindFixBin(iphi,ieta);
      else             bin = h->FindFixBin(ieta,iphi);

      hSM_allEB->SetBinContent(bin, binContent);

    }

  }

}


//==================================================


void drawMap(TH2* map2D = NULL, 
	     const string& xaxisName = "i#phi", const string& yaxisName = "i#eta", 
	     const string& canvasName = "",
	     const string& outDir = "",
	     const Double_t& mapMin = 0.95, const Double_t& mapMax = 1.05,
	     const Int_t& canvasSizeX = 600, const Int_t& canvasSizeY = 700,
	     const Int_t& nPaletteContours = 0
	     ) 
{
  
  // try to have color steps of the order of the uncertainty you expect (stat. is ~ 0.1% for pi0 in EB)
  gStyle->SetNumberContours((nPaletteContours > 0) ? nPaletteContours : (mapMax - mapMin)/0.001); 

  // draw map folded in SM
  TCanvas* cEB = new TCanvas("cEB","",canvasSizeX,canvasSizeY);
  // cEB->SetLeftMargin(0.16);
  cEB->SetRightMargin(0.14);
  cEB->cd();
  map2D->Draw("COLZ");
  map2D->GetXaxis()->SetTitle(xaxisName.c_str());
  map2D->GetXaxis()->SetTitleSize(0.06);
  map2D->GetXaxis()->SetTitleOffset(0.7);
  map2D->GetYaxis()->SetTitle(yaxisName.c_str());
  map2D->GetYaxis()->SetTitleSize(0.06);
  map2D->GetYaxis()->SetTitleOffset(0.8);
  map2D->GetZaxis()->SetRangeUser(mapMin,mapMax);
  map2D->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/%s.pdf",outDir.c_str(),canvasName.c_str()));
  cEB->SaveAs(Form("%s/%s.png",outDir.c_str(),canvasName.c_str()));
  delete cEB;

}

//====================================================

void drawDistribution(TH1* h = NULL, 
		      const string& xaxisName = "IC", const string& yaxisName = "events", 
		      const string& canvasName = "",
		      const string& outDir = "",
		      const Double_t& xmin = 0.95, const Double_t& xmax = 1.05,
		      const Int_t& canvasSizeX = 700, const Int_t& canvasSizeY = 600,
		      const Bool_t noStatBox = false
		      ) 
{

  // draw map folded in SM
  TCanvas* cEB = new TCanvas("cEB","",canvasSizeX,canvasSizeY);
  cEB->SetTickx(1);
  cEB->SetTicky(1);
  cEB->SetGrid();
  cEB->cd();
  h->Draw("HE");
  if (noStatBox) h->SetStats(0);
  else {
    h->SetStats(1);
    gStyle->SetOptStat(111110);
  }          
  h->GetXaxis()->SetTitle(xaxisName.c_str());
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(0.7);
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->GetYaxis()->SetTitle(yaxisName.c_str());
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.8);
  //h->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/%s.pdf",outDir.c_str(),canvasName.c_str()));
  cEB->SaveAs(Form("%s/%s.png",outDir.c_str(),canvasName.c_str()));
  delete cEB;

}


//====================================================


void runOnTree(TH2* map_IC = NULL,
	       TH2* map_IC_err = NULL,
	       const string& treeRootFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_6/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_calibMap.root", 
	       const Int_t izside = 0,
	       const Bool_t iphiOnXaxis = true)
{

  const Bool_t isEB = (izside == 0) ? true : false; 

  string treeName = isEB ? "calibEB" : "calibEE";
  TChain* chain = new TChain(treeName.c_str());
  chain->Add(TString(treeRootFile.c_str()));
  TTreeReader reader (chain);

  TTreeReaderValue<Int_t>   ieta_ix     (reader, isEB ? "ieta_" : "ix_");
  TTreeReaderValue<Int_t>   iphi_iy     (reader, isEB ? "iphi_" : "iy_");
  TTreeReaderValue<Int_t>*  iz  = isEB ? NULL : new TTreeReaderValue<Int_t>(reader,"zside_");
  TTreeReaderValue<Float_t> fit_mean    (reader,"fit_mean_");     // == 0 for dead channels
  TTreeReaderValue<Float_t> fit_mean_err(reader,"fit_mean_err_"); // == 0 for dead channels
  TTreeReaderValue<Float_t> IC          (reader,"coeff_");        // == 1 for dead channels

  Int_t bin = -1;

  while(reader.Next()) {
 
    if (isEB) {
      if (iphiOnXaxis) bin = map_IC->FindFixBin(*iphi_iy,*ieta_ix);
      else             bin = map_IC->FindFixBin(*ieta_ix,*iphi_iy);
    } else {
      if (izside != **iz) continue;
      bin = map_IC->FindFixBin(*ieta_ix, *iphi_iy);
    }      

    if (*fit_mean > EPSILON) {
      map_IC->    SetBinContent(bin, (Double_t) *IC);
      map_IC_err->SetBinContent(bin, (Double_t) *IC * *fit_mean_err / *fit_mean);
    } else {
      map_IC->    SetBinContent(bin, (Double_t) -1.0);
      map_IC_err->SetBinContent(bin, (Double_t)  0.0);
    }

  }

  delete chain;  

}



//==============================================

void divideEBmap(TH2* h, const TH2* num, const TH2* den, const Bool_t noBadXtals = true, const Double_t ICforBadXtals = -1.0) {

  // set h equal to num/den
  // if noBadXtals = true, set any bad xtal to ICforBadXtals
  // in general dead xtals should have IC 0 or 1 or maybe -1
  // in all other functions, tals with IC <= 0 or equal to 1 are not used for normalization

  Int_t num_nBinsX = num->GetNbinsX();
  Int_t num_nBinsY = num->GetNbinsY();
  Int_t den_nBinsX = den->GetNbinsX();
  Int_t den_nBinsY = den->GetNbinsY();
  Int_t h_nBinsX   = h->GetNbinsX();
  Int_t h_nBinsY   = h->GetNbinsY();

  if (num_nBinsX != den_nBinsX) {
    cout << "Warning in divideEBmap(): num and den have different number of bins in X ("<< num_nBinsX << "," << den_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (num_nBinsY != den_nBinsY) {
    cout << "Warning in divideEBmap(): num and den have different number of bins in Y ("<< num_nBinsY << "," << den_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (h_nBinsX != den_nBinsX) {
    cout << "Warning in divideEBmap(): h and den have different number of bins in X ("<< h_nBinsX << "," << den_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (h_nBinsY != den_nBinsY) {
    cout << "Warning in divideEBmap(): h and den have different number of bins in Y ("<< h_nBinsY << "," << den_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }


  Int_t bin = -1;
  Double_t ratio = 0.0;

  for (Int_t ix = 1; ix <= num_nBinsX; ++ix) {

    for (Int_t iy = 1; iy <= num_nBinsY; ++iy) {

      bin = num->GetBin(ix,iy);
      if (noBadXtals) {
	if (num->GetBinContent(bin) > EPSILON && fabs(num->GetBinContent(bin) -1.0) > EPSILON)
	  ratio = (den->GetBinContent(bin) != 0.0) ? (num->GetBinContent(bin) / den->GetBinContent(bin)) : ICforBadXtals;
	else
	  ratio = ICforBadXtals;
      } else {
	ratio = (den->GetBinContent(bin) != 0.0) ? (num->GetBinContent(bin) / den->GetBinContent(bin)) : ICforBadXtals;
      }

      h->SetBinContent(bin, ratio);

    }

  }


}  


//=============================================


void divideEEmap(TH2* h, const TH2* num, const TH2* den, const Bool_t noBadXtals = true, const Double_t ICforBadXtals = -1.0) {

  // set h equal to num/den
  // if noBadXtals = true, set any bad xtal to ICforBadXtals
  // in general dead xtals should have IC 0 or 1 or maybe -1
  // in all other functions, tals with IC <= 0 or equal to 1 are not used for normalization

  Int_t num_nBinsX = num->GetNbinsX();
  Int_t num_nBinsY = num->GetNbinsY();
  Int_t den_nBinsX = den->GetNbinsX();
  Int_t den_nBinsY = den->GetNbinsY();
  Int_t h_nBinsX   = h->GetNbinsX();
  Int_t h_nBinsY   = h->GetNbinsY();

  if (num_nBinsX != den_nBinsX) {
    cout << "Warning in divideEBmap(): num and den have different number of bins in X ("<< num_nBinsX << "," << den_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (num_nBinsY != den_nBinsY) {
    cout << "Warning in divideEBmap(): num and den have different number of bins in Y ("<< num_nBinsY << "," << den_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (h_nBinsX != den_nBinsX) {
    cout << "Warning in divideEBmap(): h and den have different number of bins in X ("<< h_nBinsX << "," << den_nBinsX << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }
  if (h_nBinsY != den_nBinsY) {
    cout << "Warning in divideEBmap(): h and den have different number of bins in Y ("<< h_nBinsY << "," << den_nBinsY << "). Exit." << endl;
    exit(EXIT_FAILURE);
  }


  Int_t bin = -1;
  Double_t ratio = 0.0;

  for (Int_t ix = 1; ix <= num_nBinsX; ++ix) {

    for (Int_t iy = 1; iy <= num_nBinsY; ++iy) {

      bin = num->GetBin(ix,iy);
      if (noBadXtals) {
	if (num->GetBinContent(bin) > EPSILON && fabs(num->GetBinContent(bin) -1.0) > EPSILON)
	  ratio = (den->GetBinContent(bin) != 0.0) ? (num->GetBinContent(bin) / den->GetBinContent(bin)) : ICforBadXtals;
	else
	  ratio = ICforBadXtals;
      } else {
	ratio = (den->GetBinContent(bin) != 0.0) ? (num->GetBinContent(bin) / den->GetBinContent(bin)) : ICforBadXtals;
      }

      h->SetBinContent(bin, ratio);

    }

  }


}  
