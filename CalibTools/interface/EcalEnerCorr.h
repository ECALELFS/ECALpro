#ifndef EcalEnerCorr_H
#define EcalEnerCorr_H

//
// $Id: EcalEnerCorr.h,v 1.0 2011/06/30 18:45:02 montanin

#include <vector>
#include<iostream>
#include <stdexcept>
using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

#define ENBINSEB  12
#define ETABINSEB 85

#define ENBINSEE  15
#define ETABINSEE 30

//----------------------------------------
class EcalEnerCorr {
//----------------------------------------
  public:
    EcalEnerCorr();

    ~EcalEnerCorr() {
    }
    
    double getContainmentCorrectionsEB( double energy, int ieta );
    bool loadContainmentCorrectionsEB( const char* cfile);
    double getContainmentCorrectionsEE( double energy, double ieta );
    bool loadContainmentCorrectionsEE( const char* cfile);
    double getContainmentPointCorrectionsEE( double energy, double ieta );
    bool loadContainmentPointCorrectionsEE( const char* cfile);
    bool loadContainmentMinvCorrections(const char* cfile);
    double getContainmentMinvCorrections( double ieta );
    bool etaBorderS(int i);
    bool etaBorderM(int i);
    bool uniqueFunc(int i);

// modify default varibles
    void setVerbosity(bool n);

// definition of energy bins 
    int getEnBinsEB();
    double enBinBoundEB[ENBINSEB+1]; 
    void fillEnergyEdgeEB();
    int getEnBinsEE();
    double enBinBoundEE[ENBINSEE+1];
    void fillEnergyEdgeEE();
 
// definition of eta bins
    int getEtaBinsEB();
    int getEtaBinsEE();
    double etaBinBoundEE[ETABINSEE+1];
    void fillEtaEdgeEE();
    
    static const int MaxEnergyBins = 50;

  private:

   // energy corrections EB
   TF1   *extCorrFunc1SupModBorders[ENBINSEB];
   TF1   *extCorrFunc1ModBorders[ENBINSEB];
   TF1   *extCorrFunc1Bulk[ENBINSEB];
   //  TF1   *extCorrFunc2SupModBorders[ENBINSEB];
   TF1   *extCorrFunc2ModBorders[ENBINSEB];
   TF1   *extCorrFunc2Bulk[ENBINSEB]; 
   TH1F  *extCorrPointsSupModBorders[ENBINSEB];
   TH1F  *extCorrMinvFunctionEB;

   // energy corrections EE
   TF1   *extCorrFunctionEE[ENBINSEE];
   //TH1F  *extCorrPointEE[ENBINSEE];  
   TH1F  **extCorrPointEE;  
   vector<TH1F*> extVCorrPointEE;
   TH1F  *extCorrMinvFunctionEE; 

   bool noCorrectionsEB_;            // cache to decide whether we have corrections to apply
   bool noCorrectionsEE_;            // cache to decide whether we have corrections to apply

// default variables

   bool  vdebug;
   float EBedge;
   float EEedge;

};

#endif
