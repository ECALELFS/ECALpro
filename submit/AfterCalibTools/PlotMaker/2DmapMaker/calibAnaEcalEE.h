#ifndef calibAnaEcalEE_h
#define calibAnaEcalEE_h

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

#include "calibAnaEcal.h"

class calibAnaEcalEE : public calibAnaEcal {
 public:

  calibAnaEcalEE(TTree *tree); 
  virtual ~calibAnaEcalEE() { std::cout<<"~calibAnaEcalEE() called"<<std::endl; }

  ///////////////////////////////////////
  // public member functions
  virtual void setHistograms();
  virtual void draw2Dmap(TH2D*, const Bool_t);
  virtual void drawProfile(TProfile*, const std::string&, const Bool_t );
  virtual void Init(TTree *);
  virtual void Loop();
  virtual void set2DmapMaxZaxisVector();
  virtual Int_t getEtaRingFromIxIyZside(const Int_t &, const Int_t &, const Int_t &);

  ////////////////////////////////////////////////////
  // variables used in member functions, such as Loop()


  ///////////////////////////////////////////////////
  // member functions to access protected data member

  ///////////////////////////////////////////////////
  // private or protected data members
 protected:
  Int_t etaRing;

};


#endif
