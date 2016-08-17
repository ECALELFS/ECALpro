#ifndef calibAnaEcalEB_h
#define calibAnaEcalEB_h

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

class calibAnaEcalEB : public calibAnaEcal {
 public:

  calibAnaEcalEB(TTree *tree); 
  virtual ~calibAnaEcalEB() { std::cout<<"~calibAnaEcalEB() called"<<std::endl; }

  ///////////////////////////////////////
  // public member functions
  virtual void setHistograms();
  virtual void draw2Dmap(TH2D*);
  virtual void drawProfile(TProfile*, const std::string& );
  virtual void Init(TTree *);
  virtual void Loop();
  virtual void set2DmapMaxZaxisVector();

  TProfile *mean_iphiProfile = NULL;

  ////////////////////////////////////////////////////
  // variables used in member functions, such as Loop()


  ///////////////////////////////////////////////////
  // member functions to access protected data member

  ///////////////////////////////////////////////////
  // private or protected data members


};


#endif
