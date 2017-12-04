#ifndef calibAnaEcal_h
#define calibAnaEcal_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TTree.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include "calibAnaEcal_base.h"

class calibAnaEcal : public calibAnaEcal_base {
 public:

  calibAnaEcal(TTree *tree); 
  virtual ~calibAnaEcal() { std::cout<<"~calibAnaEcal() called"<<std::endl; }

  ///////////////////////////////////////
  // public member functions
  virtual void setHistograms();
  virtual void draw2Dmap(TH2D*);
  virtual void drawProfile(TProfile*, const std::string& );
  virtual void drawChisquare(TH2D*, const Bool_t drawProfileX);
  virtual void setVerticalRangeInHisto();
  virtual void Init(TTree *);
  virtual void Loop();

  TH2D *hSignal = NULL;
  TH2D *hBackground = NULL;
  TH2D *SoverB = NULL;                      
  TH2D *SoverSqrtSplusB = NULL;
  TH2D *SigmaMeanOverMean = NULL; 
  TH2D *mean = NULL;
  TH2D *sigma = NULL;
  TH2D *chisquare = NULL;

  TH2D *chisquare_vs_etaring = NULL; // 2D plots chi^2 vs etaring

  // # of bins, lower and upper edges for TH2D. Will be set differently in the derived class for EB and EE
  // X is iphi or iX, Y is ieta or iY
  Int_t NbinsX_2Dmap;
  Double_t lowerX_2Dmap;
  Double_t upperX_2Dmap;
  Int_t NbinsY_2Dmap;
  Double_t lowerY_2Dmap;
  Double_t upperY_2Dmap;

  // Profiles of 2D maps along eta (for EB) or eta ring (for EE)
  TProfile *hSignal_etaProfile;
  TProfile *hBackground_etaProfile;
  TProfile *SoverB_etaProfile;
  TProfile *SoverSqrtSplusB_etaProfile;
  TProfile *SigmaMeanOverMean_etaProfile;
  TProfile *mean_etaProfile;
  TProfile *sigma_etaProfile;
  TProfile *chisquare_etaProfile;

  // # of bins, lower and upper edges for profile. Will be set differently in the derived class for EB and EE
  Int_t NbinsX_etaProfile;
  Double_t lowerX_etaProfile;
  Double_t upperX_etaProfile;

  std::vector<TH2D*> th2dVector;  // to keep the list TH2D above
  std::vector<Double_t> th2dMinZaxisVector;
  std::vector<Double_t> th2dMaxZaxisVector;  // code will use default when this is lower than the maximum   

  std::vector<TProfile*> profileEtaVector;  // to keep the list of TProfiles
  std::vector<std::string> profileYaxisTitle;  // names for y axis in profile plots

  ////////////////////////////////////////////////////
  // variables used in member functions, such as Loop()
  Double_t normalizedS;
  Double_t normalizedB;


  ///////////////////////////////////////////////////
  // member functions to access protected data member

  // getter
  std::string getEBorEE() const { return EBorEE; }
  std::string getPi0orEta() const { return Pi0orEta; }
  //std::string getWhichEE() const { return whichEE; }
  std::string getDirName() const { return dirName; }
  std::string getIterNumber() const { return iterNumber; }
  std::string getWwwPath() const { return wwwPath; }
  // setter
  void setEBorEE(const std::string input) { EBorEE = input; }
  void setPi0orEta(const std::string input) { Pi0orEta = input; }
  //  void setWhichEE(const std::string input) { whichEE = input; }
  void setDirName(const std::string input) { dirName = input; }
  void setIterNumber(const std::string input) { iterNumber = input; }
  void setWwwPath(const std::string input) { wwwPath = input; }
  
 protected:
  std::string EBorEE; // can be : EB, EEp, EEm
  std::string Pi0orEta; // can be : Pi0, Eta
  //  std::string whichEE;  // if EB, whichEE is "", else if EE then whichEE can be EEp or EEm
  std::string dirName;
  std::string iterNumber;
  std::string wwwPath; // to store the path to directory associated to a website (to store and display plots)

};


#endif
