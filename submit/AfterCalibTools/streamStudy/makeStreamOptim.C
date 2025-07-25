#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
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

// #include "DataFormats/DetId/interface/DetId.h"
// #include "DataFormats/EcalDetId/interface/EBDetId.h"
// #include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "./CMS_lumi.h"
//#include "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_8_0_25/src/WmassAnalysis/macros/utility.h"
#include "./utility.h"

#define CHECK_EVERY_N 10000

using namespace std;
using namespace RooFit;

static string endcap_ix_iy_zside_ietaRing = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/EE_xyzToEtaRing/eerings_modified.root";
static string deadXtalFileName = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/test_DeadXtal_AlCaP0_Run2017B_3July_upToRun297723/h_DeadXtal.root";
static bool drawAllMassPlot = false;

static double EBetaRegionBoundary = 1.0; // warning, 1.0 not exactly equal to gap between module 3 and 4
static double EEetaRegionBoundary = 1.8; // 1.8 is boundary between region 1 and 2 in EE, region 3 basically equal to region 2

//static vector<Float_t> ptGamCut = {0.8, 1.0, 1.2, 1.4};
static vector<Float_t> ptGamCut = {0.5, 0.65, 0.8, 1.4};
//static vector<Float_t> ptPairCut = {1.0, 1.5, 2.0, 2.5};
static vector<Float_t> ptPairCut = {2.0, 2.5, 2.7, 3.0};
//static vector<Float_t> ptPairOverMCut = {10, 12, 15, 18};
static vector<Float_t> ptPairOverMCut = {15, 30, 40, 50};
// static vector<Float_t> s4s9Cut = {0.8, 0.85, 0.9, 0.95};
static vector<Int_t> nXtalCut = {4, 5, 6, 7};
static vector<Float_t> clusIsoCut = {0.1, 0.2, 0.3, 0.5};
static vector<Float_t> s4s9Cut = {0.8, 0.85, 0.9, 0.95};


// following objects are used to do mass in bins of nXtal
static vector<Int_t> nXtal1Bins = {4, 5, 6, 7, 8, 9};
static vector<Int_t> nXtal2Bins = {4, 5, 6, 7, 8, 9};
static Int_t nBinsXtal1 = (Int_t) nXtal1Bins.size();
static Int_t nBinsXtal2 = (Int_t) nXtal2Bins.size();


//static bool useHLTisoCalibForComparison = false;
//static bool useStream2017ForComparison = false;

// note
// in EB I would use pt(pi0) > 2.0-2.5 and pT/M > 25-30, must see other cuts, and must also reduce cluster iso to 0.1 or 0.2 instead of 0.5
// I would then use s4s9 > 0.88 in EB, as in the stream
// also for EE it should be as tight as possible, like 

//=================================================

class vectorManager {

public:
  vectorManager() { };

  vectorManager(const vector<TH1*> & histPtrs,
		const vector<string> & histNames,
		const vector<string> & histLegends,
		const vector<Int_t> & histRebins
		)
  {
    histPtrs_    = vector<TH1*  >(histPtrs);
    histNames_   = vector<string>(histNames);
    histLegends_ = vector<string>(histLegends);
    histRebins_ = vector<Int_t>(histRebins);
  };

  ~vectorManager() {};

  vector<TH1*>   getHistPtrs()    const { return histPtrs_;    };
  vector<string> getHistNames()   const { return histNames_;   };
  vector<string> getHistLegends() const { return histLegends_; };
  vector<Int_t> getHistRebins() const { return histRebins_; };

  void addComponent(TH1* histPtr = NULL, const string& histName = "name", const string& histLegend = "leg", const Int_t& histRebin = 1)  { 
    histPtrs_.push_back(histPtr);
    histNames_.push_back(histName);
    histLegends_.push_back(histLegend);
    histRebins_.push_back(histRebin);
  };

private:

  vector<TH1*> histPtrs_;
  vector<string> histNames_;
  vector<string> histLegends_;
  vector<Int_t> histRebins_;

};

//=================================================

bool noDeadXtalIn3x3matrixSeededByThisXtal(const TH2F* hDeadXtals = NULL, const int x = 1, const int y = 1) {

  // WARNING: it is assumed that the seed is already not adjacent to a gap, therefore we won't have abs(eta)=0 or abs(ieta)=85 or iphi=1 or iphi=360 for the seed 

  int nDeadXtals = 0;

  for (int xspan = x-1; xspan <= x+1 && nDeadXtals == 0; xspan++) {
    for (int yspan = y-1; yspan <= y+1 && nDeadXtals == 0; yspan++) {
      nDeadXtals += (int) (0.5 + hDeadXtals->GetBinContent(xspan,yspan)); // histogram returns float, to avoid bad truncation sum 0.5 and then round to integer 
    }
  }

  return (nDeadXtals == 0) ? true : false;

}

//======================================================                                                                                                                     

Int_t getHistIndexByXY_int(const Int_t& varX = -1, const Int_t& varY = -1, const vector<Int_t>& xBins = {0}, const vector<Int_t>& yBins = {0}) {
  
  Bool_t xFound = false;
  Bool_t yFound = false;

  Int_t ix = 0;
  Int_t iy = 0;

  if (xBins.size() < 2) {
    cout << "Error in getHistIndexByXY_int(): x variable has " << xBins.size() << " bins (at least 2 required). Please check. Exit." << endl;
    exit(EXIT_FAILURE);    
  } else if (yBins.size() < 2) {
    cout << "Error in getHistIndexByXY_int(): y variable has " << yBins.size() << " bins (at least 2 required). Please check. Exit." << endl;
    exit(EXIT_FAILURE);    
  }
  
  if (varX < xBins.front() || varX > xBins.back()) return -1;
  if (varY < yBins.front() || varY > yBins.back()) return -2;
  
  while ( !xFound && (ix < (Int_t) xBins.size()) ) {

    if (varX == xBins[ix]) xFound = true;
    ix++;

  }

  while ( !yFound && (iy < (Int_t) yBins.size()) ) {

    if (varY == yBins[iy]) yFound = true;
    iy++;

  }

  // subtract 1 from varX and varY so that bin 0 of vector (first vector bin) is returned when the value is in the first varX and varY bin        
  return (ix -1) * ((Int_t) yBins.size()) + (iy -1);

}
  
//======================================================    

void makeStreamOptim(const bool isEB = true,
		     const bool isPi0 = true,
		     const bool useHLTisoCalibForComparison = false,
		     const bool useStream2017ForComparison = false,
		     const double lumi = 7.5,
		     const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/Optimization/streamSelection/",
		     const string& pathToNtuples = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/",
		     const TString& dirList = "AlCaP0_fromHLTPhysics_2017AB_TreeOptim,AlCaP0_from_ParkingHLTPhysics_TreeOptim"
		     ) 
{

  createPlotDirAndCopyPhp(outDir);

  TObjArray* array = dirList.Tokenize(",");
  vector<string> eosDir;
  for (Int_t j = 0; j < array->GetEntries(); j++) {
    TString str = ((TObjString *) array->At(j))->String();
    eosDir.push_back(str.Data());
    cout << j << " --> " << eosDir[j] << endl;
  }

  // might want to pass it as an option (comma separated list to be parsed, for example)
  // if (isPi0) {
  //   eosDir.push_back("AlCaP0_fromHLTPhysics_2017AB_TreeOptim");
  //   eosDir.push_back("AlCaP0_from_ParkingHLTPhysics_TreeOptim");
  // } else {
  //   eosDir.push_back("AlCaEta_fromHLTPhysics_2017AB_TreeOptim");
  //   eosDir.push_back("AlCaEta_fromParkingHLTPhysics_2017AB_TreeOptim");
  // }

  //  double lumi = 7.5; // to be updated

  TChain* chain = new TChain("Tree_Optim");

  for (UInt_t i = 0; i < eosDir.size(); i++) {

    string pathToNtuplesComplete = pathToNtuples;
    if (pathToNtuplesComplete.back() != '/') pathToNtuplesComplete += "/";
    pathToNtuplesComplete += eosDir[i];
    if (pathToNtuplesComplete.back() != '/') pathToNtuplesComplete += "/";
    pathToNtuplesComplete += "/iter_0/";

    // add specific files to the chain                        
    // only good ntuples were selected using 
    // ls -l /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_fromHLTPhysics_2017AB_TreeOptim/iter_0/ | awk '$5 > 10000 {print $9}' > AlCaP0_fromHLTPhysics_2017AB_TreeOptim_goodNtuples.txt
    string goodNtuplesFileName = eosDir[i] + "_goodNtuples.txt";
    ifstream infile(goodNtuplesFileName.c_str());
    string line;
    if(infile.is_open()){
      while(getline(infile,line)){
	if(line == "") continue;
	if(not TString(line).Contains("root")) continue;
	chain->Add((pathToNtuplesComplete+line).c_str());      
      }
    } else {
      cout << "Error opening file " << goodNtuplesFileName << endl;
      exit(EXIT_FAILURE);
    }
    infile.close();

  }

  long int nTotal = chain->GetEntries();
  long int nEvents = 0;

  cout << "Chain has: " << endl;
  cout << "##  " << chain->GetNtrees() << " trees" << endl;
  cout << "##  " << nTotal << " events" << endl;

  TTreeReader reader (chain);

  // trigger, depends on pi0 or eta
  TTreeReaderValue<Bool_t> * HLT_EB = NULL; 
  TTreeReaderValue<Bool_t> * HLT_EE = NULL; 
  if (isPi0) {
    HLT_EB = new TTreeReaderValue<Bool_t>(reader,"AlCa_EcalPi0EBonly");
    HLT_EE = new TTreeReaderValue<Bool_t> (reader,"AlCa_EcalPi0EEonly");
  } else {
    HLT_EB = new TTreeReaderValue<Bool_t>(reader,"AlCa_EcalEtaEBonly");
    HLT_EE = new TTreeReaderValue<Bool_t> (reader,"AlCa_EcalEtaEEonly");
  }

  // scalar variables
  TTreeReaderValue<ULong64_t> Event (reader,"Event");
  TTreeReaderValue<Int_t> LumiBlock (reader,"LumiBlock");
  TTreeReaderValue<Int_t> Run (reader,"Run");
  TTreeReaderValue<Int_t> STr2_NPi0_rec (reader,"STr2_NPi0_rec");
  TTreeReaderValue<Short_t> L1_SingleEG5 (reader,"L1_SingleEG5");
  TTreeReaderValue<Short_t> L1_SingleEG10 (reader,"L1_SingleEG10");
  TTreeReaderValue<Short_t> L1_SingleEG15 (reader,"L1_SingleEG15");
  TTreeReaderValue<Short_t> L1_SingleEG18 (reader,"L1_SingleEG18");
  TTreeReaderValue<Short_t> L1_SingleEG24 (reader,"L1_SingleEG24");
  TTreeReaderValue<Short_t> L1_SingleEG26 (reader,"L1_SingleEG26");
  TTreeReaderValue<Short_t> L1_SingleEG28 (reader,"L1_SingleEG28");
  TTreeReaderValue<Short_t> L1_SingleEG30 (reader,"L1_SingleEG30");
  TTreeReaderValue<Short_t> L1_SingleEG32 (reader,"L1_SingleEG32");
  TTreeReaderValue<Short_t> L1_SingleEG34 (reader,"L1_SingleEG34");
  TTreeReaderValue<Short_t> L1_SingleEG36 (reader,"L1_SingleEG36");
  TTreeReaderValue<Short_t> L1_SingleEG38 (reader,"L1_SingleEG38");
  TTreeReaderValue<Short_t> L1_SingleEG40 (reader,"L1_SingleEG40");
  TTreeReaderValue<Short_t> L1_SingleEG45 (reader,"L1_SingleEG45");
  TTreeReaderValue<Short_t> L1_SingleIsoEG18 (reader,"L1_SingleIsoEG18");
  TTreeReaderValue<Short_t> L1_SingleIsoEG20 (reader,"L1_SingleIsoEG20");
  TTreeReaderValue<Short_t> L1_SingleIsoEG22 (reader,"L1_SingleIsoEG22");
  TTreeReaderValue<Short_t> L1_SingleIsoEG24 (reader,"L1_SingleIsoEG24");
  TTreeReaderValue<Short_t> L1_SingleIsoEG26 (reader,"L1_SingleIsoEG26");
  TTreeReaderValue<Short_t> L1_SingleIsoEG28 (reader,"L1_SingleIsoEG28");
  TTreeReaderValue<Short_t> L1_SingleIsoEG30 (reader,"L1_SingleIsoEG30");
  TTreeReaderValue<Short_t> L1_SingleIsoEG32 (reader,"L1_SingleIsoEG32");
  TTreeReaderValue<Short_t> L1_SingleIsoEG34 (reader,"L1_SingleIsoEG34");
  TTreeReaderValue<Short_t> L1_SingleIsoEG36 (reader,"L1_SingleIsoEG36");
  TTreeReaderValue<Short_t> L1_SingleIsoEG18er2p1 (reader,"L1_SingleIsoEG18er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG20er2p1 (reader,"L1_SingleIsoEG20er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG22er2p1 (reader,"L1_SingleIsoEG22er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG24er2p1 (reader,"L1_SingleIsoEG24er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG26er2p1 (reader,"L1_SingleIsoEG26er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG28er2p1 (reader,"L1_SingleIsoEG28er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG30er2p1 (reader,"L1_SingleIsoEG30er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG32er2p1 (reader,"L1_SingleIsoEG32er2p1");
  TTreeReaderValue<Short_t> L1_SingleIsoEG34er2p1 (reader,"L1_SingleIsoEG34er2p1");
  TTreeReaderValue<Short_t> L1_DoubleEG_15_10 (reader,"L1_DoubleEG_15_10");
  TTreeReaderValue<Short_t> L1_DoubleEG_18_17 (reader,"L1_DoubleEG_18_17");
  TTreeReaderValue<Short_t> L1_DoubleEG_20_18 (reader,"L1_DoubleEG_20_18");
  TTreeReaderValue<Short_t> L1_DoubleEG_22_10 (reader,"L1_DoubleEG_22_10");
  TTreeReaderValue<Short_t> L1_DoubleEG_22_12 (reader,"L1_DoubleEG_22_12");
  TTreeReaderValue<Short_t> L1_DoubleEG_22_15 (reader,"L1_DoubleEG_22_15");
  TTreeReaderValue<Short_t> L1_DoubleEG_23_10 (reader,"L1_DoubleEG_23_10");
  TTreeReaderValue<Short_t> L1_DoubleEG_24_17 (reader,"L1_DoubleEG_24_17");
  TTreeReaderValue<Short_t> L1_DoubleEG_25_12 (reader,"L1_DoubleEG_25_12");
  TTreeReaderValue<Short_t> L1_SingleJet16 (reader,"L1_SingleJet16");
  TTreeReaderValue<Short_t> L1_SingleJet20 (reader,"L1_SingleJet20");
  TTreeReaderValue<Short_t> L1_SingleJet35 (reader,"L1_SingleJet35");
  TTreeReaderValue<Short_t> L1_SingleJet60 (reader,"L1_SingleJet60");
  TTreeReaderValue<Short_t> L1_SingleJet90 (reader,"L1_SingleJet90");
  TTreeReaderValue<Short_t> L1_SingleJet120 (reader,"L1_SingleJet120");
  TTreeReaderValue<Short_t> L1_SingleJet140 (reader,"L1_SingleJet140");
  TTreeReaderValue<Short_t> L1_SingleJet150 (reader,"L1_SingleJet150");
  TTreeReaderValue<Short_t> L1_SingleJet160 (reader,"L1_SingleJet160");
  TTreeReaderValue<Short_t> L1_SingleJet170 (reader,"L1_SingleJet170");
  TTreeReaderValue<Short_t> L1_SingleJet180 (reader,"L1_SingleJet180");
  TTreeReaderValue<Short_t> L1_SingleJet200 (reader,"L1_SingleJet200");
  TTreeReaderValue<Short_t> L1_DoubleJet40er3p0 (reader,"L1_DoubleJet40er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet50er3p0 (reader,"L1_DoubleJet50er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet60er3p0 (reader,"L1_DoubleJet60er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet80er3p0 (reader,"L1_DoubleJet80er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet100er3p0 (reader,"L1_DoubleJet100er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet112er3p0 (reader,"L1_DoubleJet112er3p0");
  TTreeReaderValue<Short_t> L1_DoubleJet120er3p0 (reader,"L1_DoubleJet120er3p0");
  TTreeReaderValue<Short_t> L1_TripleJet_84_68_48_VBF (reader,"L1_TripleJet_84_68_48_VBF");
  TTreeReaderValue<Short_t> L1_TripleJet_88_72_56_VBF (reader,"L1_TripleJet_88_72_56_VBF");
  TTreeReaderValue<Short_t> L1_TripleJet_92_76_64_VBF (reader,"L1_TripleJet_92_76_64_VBF");
  TTreeReaderValue<Short_t> L1_QuadJet40er3p0 (reader,"L1_QuadJet40er3p0");
  TTreeReaderValue<Short_t> L1_QuadJet50er3p0 (reader,"L1_QuadJet50er3p0");
  TTreeReaderValue<Short_t> L1_QuadJet60er3p0 (reader,"L1_QuadJet60er3p0");
  TTreeReaderValue<Short_t> L1_HTT120er (reader,"L1_HTT120er");
  TTreeReaderValue<Short_t> L1_HTT160er (reader,"L1_HTT160er");
  TTreeReaderValue<Short_t> L1_HTT200er (reader,"L1_HTT200er");
  TTreeReaderValue<Short_t> L1_HTT220er (reader,"L1_HTT220er");
  TTreeReaderValue<Short_t> L1_HTT240er (reader,"L1_HTT240er");
  TTreeReaderValue<Short_t> L1_HTT255er (reader,"L1_HTT255er");
  TTreeReaderValue<Short_t> L1_HTT270er (reader,"L1_HTT270er");
  TTreeReaderValue<Short_t> L1_HTT280er (reader,"L1_HTT280er");
  TTreeReaderValue<Short_t> L1_HTT300er (reader,"L1_HTT300er");
  TTreeReaderValue<Short_t> L1_HTT320er (reader,"L1_HTT320er");
  TTreeReaderValue<Short_t> L1_IsolatedBunch (reader,"L1_IsolatedBunch");
  TTreeReaderValue<Short_t> L1_AlwaysTrue (reader,"L1_AlwaysTrue");

  // vector variables
  TTreeReaderArray<Int_t> STr2_Pi0recIsEB (reader,"STr2_Pi0recIsEB");
  TTreeReaderArray<Float_t> STr2_IsoPi0_rec (reader,"STr2_IsoPi0_rec");
  TTreeReaderArray<Float_t> STr2_HLTIsoPi0_rec (reader,"STr2_HLTIsoPi0_rec");
  TTreeReaderArray<Int_t> STr2_n1CrisPi0_rec (reader,"STr2_n1CrisPi0_rec");
  TTreeReaderArray<Int_t> STr2_n2CrisPi0_rec (reader,"STr2_n2CrisPi0_rec");
  TTreeReaderArray<Float_t> STr2_mPi0_rec (reader,"STr2_mPi0_rec");
  TTreeReaderArray<Float_t> STr2_enG1_rec (reader,"STr2_enG1_rec");
  TTreeReaderArray<Float_t> STr2_enG2_rec (reader,"STr2_enG2_rec");
  TTreeReaderArray<Float_t> STr2_etaPi0_rec (reader,"STr2_etaPi0_rec");
  TTreeReaderArray<Float_t> STr2_ptPi0_rec (reader,"STr2_ptPi0_rec");
  TTreeReaderArray<Float_t> STr2_ptPi0_nocor (reader,"STr2_ptPi0_nocor");
  TTreeReaderArray<Float_t> STr2_enG1_nocor (reader,"STr2_enG1_nocor");
  TTreeReaderArray<Float_t> STr2_enG2_nocor (reader,"STr2_enG2_nocor");
  TTreeReaderArray<Float_t> STr2_mPi0_nocor (reader,"STr2_mPi0_nocor");
  TTreeReaderArray<Float_t> STr2_DeltaRG1G2 (reader,"STr2_DeltaRG1G2");
  TTreeReaderArray<Float_t> STr2_Es_e1_1 (reader,"STr2_Es_e1_1");
  TTreeReaderArray<Float_t> STr2_Es_e1_2 (reader,"STr2_Es_e1_2");
  TTreeReaderArray<Float_t> STr2_Es_e2_1 (reader,"STr2_Es_e2_1");
  TTreeReaderArray<Float_t> STr2_Es_e2_2 (reader,"STr2_Es_e2_2");
  TTreeReaderArray<Float_t> STr2_S4S9_1 (reader,"STr2_S4S9_1");
  TTreeReaderArray<Float_t> STr2_S4S9_2 (reader,"STr2_S4S9_2");
  TTreeReaderArray<Float_t> STr2_S2S9_1 (reader,"STr2_S2S9_1");
  TTreeReaderArray<Float_t> STr2_S2S9_2 (reader,"STr2_S2S9_2");
  TTreeReaderArray<Float_t> STr2_S1S9_1 (reader,"STr2_S1S9_1");
  TTreeReaderArray<Float_t> STr2_S1S9_2 (reader,"STr2_S1S9_2");
  TTreeReaderArray<Float_t> STr2_Eta_1 (reader,"STr2_Eta_1");
  TTreeReaderArray<Float_t> STr2_Eta_2 (reader,"STr2_Eta_2");
  TTreeReaderArray<Float_t> STr2_Phi_1 (reader,"STr2_Phi_1");
  TTreeReaderArray<Float_t> STr2_Phi_2 (reader,"STr2_Phi_2");
  TTreeReaderArray<Float_t> STr2_Time_1 (reader,"STr2_Time_1");
  TTreeReaderArray<Float_t> STr2_Time_2 (reader,"STr2_Time_2");
  TTreeReaderArray<Int_t> STr2_iEtaiX_1 (reader,"STr2_iEtaiX_1");
  TTreeReaderArray<Int_t> STr2_iEtaiX_2 (reader,"STr2_iEtaiX_2");
  TTreeReaderArray<Int_t> STr2_iPhiiY_1 (reader,"STr2_iPhiiY_1");
  TTreeReaderArray<Int_t> STr2_iPhiiY_2 (reader,"STr2_iPhiiY_2");
  TTreeReaderArray<Int_t> STr2_iEta_1on5 (reader,"STr2_iEta_1on5");
  TTreeReaderArray<Int_t> STr2_iEta_2on5 (reader,"STr2_iEta_2on5");
  TTreeReaderArray<Int_t> STr2_iPhi_1on2 (reader,"STr2_iPhi_1on2");
  TTreeReaderArray<Int_t> STr2_iPhi_2on2 (reader,"STr2_iPhi_2on2");
  TTreeReaderArray<Int_t> STr2_iEta_1on2520 (reader,"STr2_iEta_1on2520");
  TTreeReaderArray<Int_t> STr2_iEta_2on2520 (reader,"STr2_iEta_2on2520");
  TTreeReaderArray<Int_t> STr2_iPhi_1on20 (reader,"STr2_iPhi_1on20");
  TTreeReaderArray<Int_t> STr2_iPhi_2on20 (reader,"STr2_iPhi_2on20");


  float massMinEB = isPi0 ? 0.06 : 0.2;
  float massMaxEB = isPi0 ? 0.22 : 0.8;
  int nBinsMassEB = isPi0 ? 160 : 300;
  float massMinEE = isPi0 ? 0.05 : 0.2;
  float massMaxEE = isPi0 ? 0.25 : 0.8;
  int nBinsMassEE = isPi0 ? 200 : 300;

  ///////////////////////
  // histograms

  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////
  //// Work in progress
  vector<string> regionTagName;
  regionTagName.push_back("region1EB");
  regionTagName.push_back("region2EB");
  regionTagName.push_back("region1EE");
  regionTagName.push_back("region2EE");
  regionTagName.push_back("region3EE");

  // used in the plot to write some details                                      
  vector<string> plotTag;
  plotTag.push_back("EB |#eta|<1.0");
  plotTag.push_back("EB |#eta|>1.0");
  plotTag.push_back("EE |#eta|< 1.8");
  plotTag.push_back("EE 1.8<|#eta|<2.0");
  plotTag.push_back("EE |#eta|>2.0");

  vector <TH1F*> hmass_ntpSel;
  vector <TH1F*> hmass_2016Sel;
  vector <TH1F*> hmass_2017Sel;
  vector <TH1F*> hmass_2016SelCalib;
  vector <TH1F*> hmass_2017SelCalib;
  vector <TH1F*> hmass_2017SelCalibLowPtGam;
  vector <TH1F*> hmass_2017SelCalibHLTiso;

  vector <TH1F*> hmass_2017SelCalibOptim5nXtal;
  vector <TH1F*> hmass_2017SelCalibOptim6nXtal;
  vector <TH1F*> hmass_2017SelCalibOptim7nXtal;

  vector <TH1F*> hmass_2017SelCalibOptim65nXtal;
  vector <TH1F*> hmass_2017SelCalibOptim75nXtal;
  vector <TH1F*> hmass_2017SelCalibOptim76nXtal;
  //  vector <TH1F*> hmass_2017SelCalibOptim0p3clusIso;

  vector<vectorManager*> vm;
  vector<vectorManager*> vm_bis;
  vector<vectorManager*> vmOptim;
  vector<vectorManager*> vmOptimNxtal;

  vector<Int_t> rebinFactorEBin;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorEBin.push_back(1);
  rebinFactorEBin.push_back(1);
  rebinFactorEBin.push_back(1);
  rebinFactorEBin.push_back(1);
  rebinFactorEBin.push_back(1);

  // used only for region 2
  vector<Int_t> rebinFactorEB;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(2);
  rebinFactorEB.push_back(2);

  vector<Int_t> rebinFactorEE;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(8);
  rebinFactorEE.push_back(5);

  vector<Int_t> rebinFactorBisEBin;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorBisEBin.push_back(1);
  rebinFactorBisEBin.push_back(1);
  rebinFactorBisEBin.push_back(1);
  rebinFactorBisEBin.push_back(1);
  rebinFactorBisEBin.push_back(1);

  vector<Int_t> rebinFactorBisEB;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorBisEB.push_back(1);
  rebinFactorBisEB.push_back(2);
  rebinFactorBisEB.push_back(2);
  rebinFactorBisEB.push_back(2);
  rebinFactorBisEB.push_back(2);

  vector<Int_t> rebinFactorBisEE;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorBisEE.push_back(2);
  rebinFactorBisEE.push_back(5);
  rebinFactorBisEE.push_back(8);
  rebinFactorBisEE.push_back(8);
  rebinFactorBisEE.push_back(5);

  vector<Int_t> rebinFactorOptimEBin;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorOptimEBin.push_back(1);
  rebinFactorOptimEBin.push_back(1);
  rebinFactorOptimEBin.push_back(1);
  rebinFactorOptimEBin.push_back(1);
  rebinFactorOptimEBin.push_back(1);

  // used only for region 2
  vector<Int_t> rebinFactorOptimEB;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorOptimEB.push_back(1);
  rebinFactorOptimEB.push_back(2);
  rebinFactorOptimEB.push_back(1);
  rebinFactorOptimEB.push_back(2);
  rebinFactorOptimEB.push_back(2);

  vector<Int_t> rebinFactorOptimEE;  // rebins associated to 5 selections above, not the 5 regions of ECAL  
  rebinFactorOptimEE.push_back(2);
  rebinFactorOptimEE.push_back(8);
  rebinFactorOptimEE.push_back(5);
  rebinFactorOptimEE.push_back(5);
  rebinFactorOptimEE.push_back(8);

  vector <TH1F*> hmass_2017Sel_SingleEGX;
  vector <TH1F*> hmass_2017Sel_SingleIsoEGX;
  vector <TH1F*> hmass_2017Sel_SingleIsoEGXer2p1;
  vector <TH1F*> hmass_2017Sel_DoubleEGX;                                                                                                              
  vector <TH1F*> hmass_2017Sel_SingleJetX;                                                                                                             
  vector <TH1F*> hmass_2017Sel_DoubleJetXer3p0;                                                                                                        
  vector <TH1F*> hmass_2017Sel_TripleJetX;                                                                                                             
  vector <TH1F*> hmass_2017Sel_QuadJetXer3p0;                                                                                                          
  vector <TH1F*> hmass_2017Sel_HTTXer;

  vector <TH1F*> hmass_2016Sel_1nearGap;  // both photons with seed close to gap or dead region (seed facing gap)
  vector <TH1F*> hnXtal1_2016Sel_1nearGap; 
  vector <TH1F*> hnXtal2_2016Sel_1nearGap; 
  vector <TH1F*> hmass_2016Sel_1far1nearGap;
  vector <TH1F*> hnXtal1_2016Sel_1far1nearGap; 
  vector <TH1F*> hnXtal2_2016Sel_1far1nearGap; 
  vector <TH1F*> hmass_2016Sel_1farGap;  // both photons with seed far from gap or dead region (at least 1 crystal) 

  vector<vectorManager*> vm_nearfar;
  vector<vectorManager*> vm_nearfar_nxtal;

  vector <TH1F*> hmass_2016Sel_notDoubleEGX;
  vector <TH1F*> hmass_2016Sel_notSingleEGX;
  vector <TH1F*> hmass_2017Sel_notDoubleEGX;
  vector <TH1F*> hmass_2017Sel_notSingleEGX;
    
  vector<vectorManager*> vm_L1seedEG;
  vector<vectorManager*> vm_L1seedJet;

  vector<vectorManager*> vm_L1seed_SingleOrDoubleEG;

  vector< vector<TH1F*> > hmass_2017SelCalib_ptGam(regionTagName.size()); // 5 object of type vector<TH1F*>
  vector< vector<TH1F*> > hmass_2017SelCalib_ptPair(regionTagName.size());
  vector< vector<TH1F*> > hmass_2017SelCalib_ptPairOverM(regionTagName.size());
  vector< vector<TH1F*> > hmass_2017SelCalib_nXtal(regionTagName.size());
  vector< vector<TH1F*> > hmass_2017SelCalib_clusIso(regionTagName.size());
  vector< vector<TH1F*> > hmass_2017SelCalib_s4s9(regionTagName.size());

  vector<vectorManager*> vm_ptGam;
  vector<vectorManager*> vm_ptPair;
  vector<vectorManager*> vm_ptPairOverM;
  vector<vectorManager*> vm_nXtal;
  vector<vectorManager*> vm_clusIso;
  vector<vectorManager*> vm_s4s9;

  vector <TH1F*> hnXtal1_ntpSel;
  vector <TH1F*> hnXtal1_2016Sel;
  vector <TH1F*> hnXtal1_2016SelCalib;
  vector <TH1F*> hnXtal2_ntpSel;
  vector <TH1F*> hnXtal2_2016Sel;
  vector <TH1F*> hnXtal2_2016SelCalib;

  vector<vectorManager*> vm_nXtalDistr;

  vector< vector <TH1F*> > hmass_nXtal1_nXtal2(regionTagName.size()); // use a 1D array (one for each region), then I will build a 2D map 
  vector<TH2F*> h2_massPeak_nXtal1_nXtal2;
  vector<TH2F*> h2_massSigma_nXtal1_nXtal2;
  vector<TH2F*> h2_massReso_nXtal1_nXtal2;
  string selForXtalStudy = "calibSel2016"; // "calibSel2016" or "streamSel2016"

  // pT distributions of photons and pi0
  // using base selection (unbiased wrt pT) or stream in 2016 and 2017 (difference is only in nXtal) but removing pT cuts
  vector <TH1F*> hptGam1_ntpSel;
  vector <TH1F*> hptGam1_2016Sel_noPtCuts;
  vector <TH1F*> hptGam1_2017Sel_noPtCuts;

  vector <TH1F*> hptGam2_ntpSel;
  vector <TH1F*> hptGam2_2016Sel_noPtCuts;
  vector <TH1F*> hptGam2_2017Sel_noPtCuts;

  vector <TH1F*> hptPair_ntpSel;
  vector <TH1F*> hptPair_2016Sel_noPtCuts;
  vector <TH1F*> hptPair_2017Sel_noPtCuts;

  vector <TH1F*> hptPair_2017SelCalib_optim2018;

  vector<vectorManager*> vm_ptGam1Distr;
  vector<vectorManager*> vm_ptGam2Distr;
  vector<vectorManager*> vm_ptPairDistr;

  // energy
  TH1F* henergyGam1_EBreg1_calibOptim = new TH1F("henergyGam1_EBreg1_calibOptim","",250,0,25.0);
  TH1F* henergyGam1_EBreg2_calibOptim = new TH1F("henergyGam1_EBreg2_calibOptim","",250,0,25.0);;

  vectorManager* vm_energyGam1Distr = new vectorManager();
  vm_energyGam1Distr->addComponent(henergyGam1_EBreg1_calibOptim, henergyGam1_EBreg1_calibOptim->GetName(), "EB |#eta|<1.0");
  vm_energyGam1Distr->addComponent(henergyGam1_EBreg2_calibOptim, henergyGam1_EBreg2_calibOptim->GetName(), "EB |#eta|>1.0");

  string histName = "";
  string doubleToStr = "";
  
  for (UInt_t i = 0; i < regionTagName.size(); i++) {

    string name = regionTagName[i];
    double massMin = (i < 2) ? massMinEB : massMinEE;
    double massMax = (i < 2) ? massMaxEB : massMaxEE;
    int nBinsMass = (i < 2) ? nBinsMassEB : nBinsMassEE;

    hmass_ntpSel.push_back( new TH1F(Form("hmass_%s_ntpSel",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2016Sel.push_back( new TH1F(Form("hmass_%s_2016Sel",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel.push_back( new TH1F(Form("hmass_%s_2017Sel",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2016SelCalib.push_back( new TH1F(Form("hmass_%s_2016SelCalib",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalib.push_back( new TH1F(Form("hmass_%s_2017SelCalib",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibLowPtGam.push_back( new TH1F(Form("hmass_%s_2017SelCalibLowPtGam",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibHLTiso.push_back( new TH1F(Form("hmass_%s_2017SelCalibHLTiso",name.c_str()),"",nBinsMass,massMin,massMax) );   

    hmass_2017SelCalibOptim5nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim5nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibOptim6nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim6nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibOptim65nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim65nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibOptim75nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim75nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibOptim76nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim76nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017SelCalibOptim7nXtal.push_back( new TH1F(Form("hmass_%s_2017SelCalibOptim7nXtal",name.c_str()),"",nBinsMass,massMin,massMax) );

    vm.push_back( new vectorManager() );
    vm[i]->addComponent(hmass_ntpSel[i], hmass_ntpSel[i]->GetName(), Form("%s, loose",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2016Sel[i], hmass_2016Sel[i]->GetName(), Form("%s, 2016",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso"); // shorter form			

    vm_bis.push_back( new vectorManager() );
    //vm_bis[i]->addComponent(hmass_ntpSel[i], hmass_ntpSel[i]->GetName(), Form("%s, loose",plotTag[i].c_str()));
    vm_bis[i]->addComponent(hmass_2016Sel[i], hmass_2016Sel[i]->GetName(), Form("%s, 2016",plotTag[i].c_str()));
    vm_bis[i]->addComponent(hmass_2016SelCalib[i], hmass_2016SelCalib[i]->GetName(), Form("%s, 2016 calib",plotTag[i].c_str()));
    //vm_bis[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso"); // shorter form			
    vm_bis[i]->addComponent(hmass_2017SelCalibLowPtGam[i], hmass_2017SelCalibLowPtGam[i]->GetName(),  "2017 calib, p_{T}(#gamma)>0.8");
    vm_bis[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));
    vm_bis[i]->addComponent(hmass_2017SelCalibOptim7nXtal[i], hmass_2017SelCalibOptim7nXtal[i]->GetName(), "2017 calib Optim nXtal>=7");

    vmOptim.push_back( new vectorManager() );
    vmOptim[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));
    vmOptim[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));
    vmOptim[i]->addComponent(hmass_2017SelCalibOptim5nXtal[i], hmass_2017SelCalibOptim5nXtal[i]->GetName(), "2017 calib Optim nXtal>=5");
    vmOptim[i]->addComponent(hmass_2017SelCalibOptim6nXtal[i], hmass_2017SelCalibOptim6nXtal[i]->GetName(), "2017 calib Optim nXtal>=6");
    vmOptim[i]->addComponent(hmass_2017SelCalibOptim7nXtal[i], hmass_2017SelCalibOptim7nXtal[i]->GetName(), "2017 calib Optim nXtal>=7");

    vmOptimNxtal.push_back( new vectorManager() );
    vmOptimNxtal[i]->addComponent(hmass_2017SelCalibOptim65nXtal[i], hmass_2017SelCalibOptim65nXtal[i]->GetName(), Form("%s, nXtal>=6,5",plotTag[i].c_str()));
    vmOptimNxtal[i]->addComponent(hmass_2017SelCalibOptim6nXtal[i], hmass_2017SelCalibOptim6nXtal[i]->GetName(), "2017 calib Optim nXtal>=6");
    vmOptimNxtal[i]->addComponent(hmass_2017SelCalibOptim75nXtal[i], hmass_2017SelCalibOptim75nXtal[i]->GetName(), "2017 calib Optim nXtal>=7,5");
    vmOptimNxtal[i]->addComponent(hmass_2017SelCalibOptim76nXtal[i], hmass_2017SelCalibOptim76nXtal[i]->GetName(), "2017 calib Optim nXtal>=7,6");
    vmOptimNxtal[i]->addComponent(hmass_2017SelCalibOptim7nXtal[i], hmass_2017SelCalibOptim7nXtal[i]->GetName(), "2017 calib Optim nXtal>=7");

    hmass_2017Sel_SingleEGX.push_back( new TH1F(Form("hmass_%s_2017Sel_SingleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel_SingleIsoEGX.push_back( new TH1F(Form("hmass_%s_2017Sel_SingleIsoEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel_SingleIsoEGXer2p1.push_back( new TH1F(Form("hmass_%s_2017Sel_SingleIsoEGXer2p1",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel_DoubleEGX.push_back( new TH1F(Form("hmass_%s_2017Sel_DoubleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );       
    hmass_2017Sel_SingleJetX.push_back( new TH1F(Form("hmass_%s_2017Sel_SingleJetX",name.c_str()),"",nBinsMass,massMin,massMax) );                          
    hmass_2017Sel_DoubleJetXer3p0.push_back( new TH1F(Form("hmass_%s_2017Sel_DoubleJetXer3p0",name.c_str()),"",nBinsMass,massMin,massMax) );      
    hmass_2017Sel_TripleJetX.push_back( new TH1F(Form("hmass_%s_2017Sel_TripleJetX",name.c_str()),"",nBinsMass,massMin,massMax) );     
    hmass_2017Sel_QuadJetXer3p0.push_back( new TH1F(Form("hmass_%s_2017Sel_QuadJetXer3p0",name.c_str()),"",nBinsMass,massMin,massMax) );                                 
    hmass_2017Sel_HTTXer.push_back( new TH1F(Form("hmass_%s_2017Sel_HTTXer",name.c_str()),"",nBinsMass,massMin,massMax) );

    hmass_2016Sel_notDoubleEGX.push_back( new TH1F(Form("hmass_%s_2016Sel_notDoubleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2016Sel_notSingleEGX.push_back( new TH1F(Form("hmass_%s_2016Sel_notSingleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel_notDoubleEGX.push_back( new TH1F(Form("hmass_%s_2017Sel_notDoubleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    hmass_2017Sel_notSingleEGX.push_back( new TH1F(Form("hmass_%s_2017Sel_notSingleEGX",name.c_str()),"",nBinsMass,massMin,massMax) );
    
    vm_L1seedEG.push_back( new vectorManager() );
    vm_L1seedEG[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));
    vm_L1seedEG[i]->addComponent(hmass_2017Sel_SingleEGX[i], hmass_2017Sel_SingleEGX[i]->GetName(), "SingleEGX");
    vm_L1seedEG[i]->addComponent(hmass_2017Sel_SingleIsoEGX[i], hmass_2017Sel_SingleIsoEGX[i]->GetName(), "SingleIsoEGX");
    vm_L1seedEG[i]->addComponent(hmass_2017Sel_SingleIsoEGXer2p1[i], hmass_2017Sel_SingleIsoEGXer2p1[i]->GetName(), "SingleIsoEGXer2p1");
    vm_L1seedEG[i]->addComponent(hmass_2017Sel_DoubleEGX[i], hmass_2017Sel_DoubleEGX[i]->GetName(), "DoubleEGX");

    vm_L1seedJet.push_back( new vectorManager() );
    vm_L1seedJet[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));
    vm_L1seedJet[i]->addComponent(hmass_2017Sel_SingleJetX[i], hmass_2017Sel_SingleJetX[i]->GetName(), "SingleJetX");
    vm_L1seedJet[i]->addComponent(hmass_2017Sel_DoubleJetXer3p0[i], hmass_2017Sel_DoubleJetXer3p0[i]->GetName(), "DoubleJetXer3p0");
    vm_L1seedJet[i]->addComponent(hmass_2017Sel_TripleJetX[i], hmass_2017Sel_TripleJetX[i]->GetName(), "TripleJetX");
    vm_L1seedJet[i]->addComponent(hmass_2017Sel_QuadJetXer3p0[i], hmass_2017Sel_QuadJetXer3p0[i]->GetName(), "QuadJetXer3p0");
    vm_L1seedJet[i]->addComponent(hmass_2017Sel_HTTXer[i], hmass_2017Sel_HTTXer[i]->GetName(), "HTTXer");

    // barrel only
    if (i < 2) {
      hmass_2016Sel_1nearGap.push_back( new TH1F(Form("hmass_%s_2016Sel_1nearGap",name.c_str()),"",nBinsMass,massMin,massMax) );    
      hnXtal1_2016Sel_1nearGap.push_back( new TH1F(Form("hnXtal1_%s_2016Sel_1nearGap",name.c_str()),"",6,3.5,9.5) );    
      hnXtal2_2016Sel_1nearGap.push_back( new TH1F(Form("hnXtal2_%s_2016Sel_1nearGap",name.c_str()),"",6,3.5,9.5) );    
      hmass_2016Sel_1far1nearGap.push_back( new TH1F(Form("hmass_%s_2016Sel_1far1nearGap",name.c_str()),"",nBinsMass,massMin,massMax) );    
      hnXtal1_2016Sel_1far1nearGap.push_back( new TH1F(Form("hnXtal1_%s_2016Sel_1far1nearGap",name.c_str()),"",6,3.5,9.5) );    
      hnXtal2_2016Sel_1far1nearGap.push_back( new TH1F(Form("hnXtal2_%s_2016Sel_1far1nearGap",name.c_str()),"",6,3.5,9.5) );    
      hmass_2016Sel_1farGap.push_back( new TH1F(Form("hmass_%s_2016Sel_1farGap",name.c_str()),"",nBinsMass,massMin,massMax) );
      
      vm_nearfar.push_back( new vectorManager() ); 
      vm_nearfar[i]->addComponent(hmass_2016Sel[i], hmass_2016Sel[i]->GetName(), Form("%s, 2016",plotTag[i].c_str()));
      vm_nearfar[i]->addComponent(hmass_2016Sel_1nearGap[i], hmass_2016Sel_1nearGap[i]->GetName(), "both #gamma far from gap");
      vm_nearfar[i]->addComponent(hmass_2016Sel_1far1nearGap[i], hmass_2016Sel_1far1nearGap[i]->GetName(), ">0 #gamma close to gap");
      //vm_nearfar[i]->addComponent(hmass_2016Sel_1farGap[i], hmass_2016Sel_1farGap[i]->GetName(), "both #gamma close to gap");      

      vm_nearfar_nxtal.push_back( new vectorManager() ); 
      vm_nearfar_nxtal[i]->addComponent(hnXtal1_2016Sel_1nearGap[i], hnXtal1_2016Sel_1nearGap[i]->GetName(), "#gamma1, both far from gap");
      vm_nearfar_nxtal[i]->addComponent(hnXtal2_2016Sel_1nearGap[i], hnXtal2_2016Sel_1nearGap[i]->GetName(), "#gamma2, both far from gap");
      vm_nearfar_nxtal[i]->addComponent(hnXtal1_2016Sel_1far1nearGap[i], hnXtal1_2016Sel_1far1nearGap[i]->GetName(), "#gamma1, >0#gamma close to gap");
      vm_nearfar_nxtal[i]->addComponent(hnXtal2_2016Sel_1far1nearGap[i], hnXtal2_2016Sel_1far1nearGap[i]->GetName(), "#gamma2, >0#gamma close to gap");

    }

    vm_L1seed_SingleOrDoubleEG.push_back( new vectorManager() );
    vm_L1seed_SingleOrDoubleEG[i]->addComponent(hmass_2016Sel[i], hmass_2016Sel[i]->GetName(), Form("%s, 2016",plotTag[i].c_str()));
    vm_L1seed_SingleOrDoubleEG[i]->addComponent(hmass_2016Sel_notDoubleEGX[i], hmass_2016Sel_notDoubleEGX[i]->GetName(), "! DoubleEGX, nXtal>=4");
    vm_L1seed_SingleOrDoubleEG[i]->addComponent(hmass_2016Sel_notSingleEGX[i], hmass_2016Sel_notSingleEGX[i]->GetName(), "! SingleEGX, nXtal>=4");
    vm_L1seed_SingleOrDoubleEG[i]->addComponent(hmass_2017Sel_notDoubleEGX[i], hmass_2016Sel_notDoubleEGX[i]->GetName(), "! DoubleEGX, nXtal>=6");
    vm_L1seed_SingleOrDoubleEG[i]->addComponent(hmass_2017Sel_notSingleEGX[i], hmass_2016Sel_notSingleEGX[i]->GetName(), "! SingleEGX, nXtal>=6");

    vm_ptGam.push_back( new vectorManager() );
    vm_ptPair.push_back( new vectorManager() );
    vm_ptPairOverM.push_back( new vectorManager() );
    vm_nXtal.push_back( new vectorManager() );
    vm_clusIso.push_back( new vectorManager() );
    vm_s4s9.push_back( new vectorManager() );

    if (useStream2017ForComparison) {
      vm_ptGam[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
      vm_ptPair[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
      vm_ptPairOverM[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
      vm_nXtal[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
      vm_clusIso[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
      vm_s4s9[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s, 2017",plotTag[i].c_str()));  
    } else if (useHLTisoCalibForComparison)  {
      vm_ptGam[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
      vm_ptPair[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
      vm_ptPairOverM[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
      vm_nXtal[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
      vm_clusIso[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
      vm_s4s9[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso");  
    } else {
      vm_ptGam[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
      vm_ptPair[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
      vm_ptPairOverM[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
      vm_nXtal[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
      vm_clusIso[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
      vm_s4s9[i]->addComponent(hmass_2017SelCalib[i], hmass_2017SelCalib[i]->GetName(), Form("%s, 2017 calib",plotTag[i].c_str()));  
    }

    for (UInt_t jcut = 0; jcut < ptGamCut.size(); jcut++) {
      doubleToStr = getStringFromDouble(ptGamCut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_ptGam%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_ptGam[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_ptGam[i]->addComponent(hmass_2017SelCalib_ptGam[i][jcut], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[jcut]));  
    }
    for (UInt_t jcut = 0; jcut < ptPairCut.size(); jcut++) {
      doubleToStr = getStringFromDouble(ptPairCut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_ptPair%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_ptPair[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_ptPair[i]->addComponent(hmass_2017SelCalib_ptPair[i][jcut], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[jcut]));  
    }
    for (UInt_t jcut = 0; jcut < ptPairOverMCut.size(); jcut++) {
      doubleToStr = getStringFromDouble(ptPairOverMCut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_ptPairOverM%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_ptPairOverM[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_ptPairOverM[i]->addComponent(hmass_2017SelCalib_ptPairOverM[i][jcut], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[jcut]));  
    }
    for (UInt_t jcut = 0; jcut < nXtalCut.size(); jcut++) {
      doubleToStr = getStringFromDouble(nXtalCut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_nXtal%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_nXtal[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_nXtal[i]->addComponent(hmass_2017SelCalib_nXtal[i][jcut], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[jcut]));  
    }
    for (UInt_t jcut = 0; jcut < clusIsoCut.size(); jcut++) {
      doubleToStr = getStringFromDouble(clusIsoCut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_clusIso%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_clusIso[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_clusIso[i]->addComponent(hmass_2017SelCalib_clusIso[i][jcut], histName.c_str(), Form("clusIso > %1.2f",clusIsoCut[jcut]));  
    }
    for (UInt_t jcut = 0; jcut < s4s9Cut.size(); jcut++) {
      doubleToStr = getStringFromDouble(s4s9Cut[jcut]);
      histName = string(Form("hmass_%s_2017SelCalib_s4s9%s",name.c_str(),doubleToStr.c_str()));
      hmass_2017SelCalib_s4s9[i].push_back(new TH1F(histName.c_str(),"",nBinsMass,massMin,massMax));
      vm_s4s9[i]->addComponent(hmass_2017SelCalib_s4s9[i][jcut], histName.c_str(), Form("s4s9 > %1.2f",s4s9Cut[jcut]));  
    }

    hnXtal1_ntpSel.push_back( new TH1F(Form("hnXtal1_%s_ntpSel",name.c_str()),"", 6, 3.5, 9.5) );
    hnXtal1_2016Sel.push_back( new TH1F(Form("hnXtal1_%s_2016Sel",name.c_str()),"", 6, 3.5, 9.5) );
    hnXtal1_2016SelCalib.push_back( new TH1F(Form("hnXtal1_%s_2016SelCalib",name.c_str()),"", 6, 3.5, 9.5) );
    hnXtal2_ntpSel.push_back( new TH1F(Form("hnXtal2_%s_ntpSel",name.c_str()),"", 6, 3.5, 9.5) );
    hnXtal2_2016Sel.push_back( new TH1F(Form("hnXtal2_%s_2016Sel",name.c_str()),"", 6, 3.5, 9.5) );
    hnXtal2_2016SelCalib.push_back( new TH1F(Form("hnXtal2_%s_2016SelCalib",name.c_str()),"", 6, 3.5, 9.5) );

    vm_nXtalDistr.push_back( new vectorManager() );
    vm_nXtalDistr[i]->addComponent(hnXtal1_ntpSel[i], hnXtal1_ntpSel[i]->GetName(), Form("%s, loose #gamma_{1}",plotTag[i].c_str()));
    vm_nXtalDistr[i]->addComponent(hnXtal2_ntpSel[i], hnXtal2_ntpSel[i]->GetName(), Form("%s, loose #gamma_{2}",plotTag[i].c_str()));
    vm_nXtalDistr[i]->addComponent(hnXtal1_2016Sel[i], hnXtal1_2016Sel[i]->GetName(), Form("%s, stream #gamma_{1}",plotTag[i].c_str()));
    vm_nXtalDistr[i]->addComponent(hnXtal2_2016Sel[i], hnXtal2_2016Sel[i]->GetName(), Form("%s, stream #gamma_{2}",plotTag[i].c_str()));
    vm_nXtalDistr[i]->addComponent(hnXtal1_2016SelCalib[i], hnXtal1_2016SelCalib[i]->GetName(), Form("%s, calib #gamma_{1}",plotTag[i].c_str())); 
    vm_nXtalDistr[i]->addComponent(hnXtal2_2016SelCalib[i], hnXtal2_2016SelCalib[i]->GetName(), Form("%s, calib #gamma_{2}",plotTag[i].c_str())); 

    for (Int_t xtal1bin = 0; xtal1bin < nBinsXtal1; xtal1bin++) {
      for (Int_t xtal2bin = 0; xtal2bin < nBinsXtal2; xtal2bin++) {
	hmass_nXtal1_nXtal2[i].push_back(new TH1F(Form("hmass_%s_nXtalLead%d_nXtalTrail%d",name.c_str(),nXtal1Bins[xtal1bin],nXtal2Bins[xtal2bin]),
						  "",nBinsMass,massMin,massMax)
					 );
      }
    }
    h2_massPeak_nXtal1_nXtal2.push_back(new TH2F(Form("h2_massPeak_%s_nXtal1_nXtal2",name.c_str()),"",
						 nBinsXtal1,nXtal1Bins.front()-0.5,nXtal1Bins.back()+0.5,
						 nBinsXtal2,nXtal2Bins.front()-0.5,nXtal2Bins.back()+0.5)
					);    
    h2_massSigma_nXtal1_nXtal2.push_back(new TH2F(Form("h2_massSigma_%s_nXtal1_nXtal2",name.c_str()),"",
						  nBinsXtal1,nXtal1Bins.front()-0.5,nXtal1Bins.back()+0.5,
						  nBinsXtal2,nXtal2Bins.front()-0.5,nXtal2Bins.back()+0.5)
					 );    
    h2_massReso_nXtal1_nXtal2.push_back(new TH2F(Form("h2_massReso_%s_nXtal1_nXtal2",name.c_str()),"",
						  nBinsXtal1,nXtal1Bins.front()-0.5,nXtal1Bins.back()+0.5,
						  nBinsXtal2,nXtal2Bins.front()-0.5,nXtal2Bins.back()+0.5)
					 );    

    hptGam1_ntpSel.push_back( new TH1F(Form("hptGam1_%s_ntpSel",name.c_str()),"", 100, 0.0, 10) );
    hptGam1_2016Sel_noPtCuts.push_back( new TH1F(Form("hptGam1_%s_2016Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );
    hptGam1_2017Sel_noPtCuts.push_back( new TH1F(Form("hptGam1_%s_2017Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );
    hptGam2_ntpSel.push_back( new TH1F(Form("hptGam2_%s_ntpSel",name.c_str()),"", 100, 0.0, 10) );
    hptGam2_2016Sel_noPtCuts.push_back( new TH1F(Form("hptGam2_%s_2016Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );
    hptGam2_2017Sel_noPtCuts.push_back( new TH1F(Form("hptGam2_%s_2017Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );
    hptPair_ntpSel.push_back( new TH1F(Form("hptPair_%s_ntpSel",name.c_str()),"", 100, 0.0, 10) );
    hptPair_2016Sel_noPtCuts.push_back( new TH1F(Form("hptPair_%s_2016Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );
    hptPair_2017Sel_noPtCuts.push_back( new TH1F(Form("hptPair_%s_2017Sel_noPtCuts",name.c_str()),"", 100, 0.0, 10) );

    hptPair_2017SelCalib_optim2018.push_back( new TH1F(Form("hptPair_%s_2017SelCalib_optim2018",name.c_str()),"", 150, 0.0, 15) );

    vm_ptGam1Distr.push_back( new vectorManager() );
    vm_ptGam1Distr[i]->addComponent(hptGam1_ntpSel[i], hptGam1_ntpSel[i]->GetName(), Form("%s, loose",plotTag[i].c_str()));
    vm_ptGam1Distr[i]->addComponent(hptGam1_2016Sel_noPtCuts[i], hptGam1_2016Sel_noPtCuts[i]->GetName(), "stream 2016 no p_{T} cuts");
    vm_ptGam1Distr[i]->addComponent(hptGam1_2017Sel_noPtCuts[i], hptGam1_2017Sel_noPtCuts[i]->GetName(), "stream 2017 no p_{T} cuts");
    vm_ptGam2Distr.push_back( new vectorManager() );
    vm_ptGam2Distr[i]->addComponent(hptGam2_ntpSel[i], hptGam2_ntpSel[i]->GetName(), Form("%s, loose",plotTag[i].c_str()));
    vm_ptGam2Distr[i]->addComponent(hptGam2_2016Sel_noPtCuts[i], hptGam2_2016Sel_noPtCuts[i]->GetName(), "stream 2016 no p_{T} cuts");
    vm_ptGam2Distr[i]->addComponent(hptGam2_2017Sel_noPtCuts[i], hptGam2_2017Sel_noPtCuts[i]->GetName(), "stream 2017 no p_{T} cuts");
    vm_ptPairDistr.push_back( new vectorManager() );
    vm_ptPairDistr[i]->addComponent(hptPair_ntpSel[i], hptPair_ntpSel[i]->GetName(), Form("%s, loose",plotTag[i].c_str()));
    vm_ptPairDistr[i]->addComponent(hptPair_2016Sel_noPtCuts[i], hptPair_2016Sel_noPtCuts[i]->GetName(), "stream 2016 no p_{T} cuts");
    vm_ptPairDistr[i]->addComponent(hptPair_2017Sel_noPtCuts[i], hptPair_2017Sel_noPtCuts[i]->GetName(), "stream 2017 no p_{T} cuts");
  }


  // end of work in progress
  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////

  Int_t nRegion = -1;

  TH2F* hDeadXtalEB = NULL;
  TH2F* hDeadXtalEEm = NULL;
  TH2F* hDeadXtalEEp = NULL;
  
  TFile* deadXtalFile = TFile::Open(deadXtalFileName.c_str(),"READ");
  if (!deadXtalFile || !deadXtalFile->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< deadXtalFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  
  hDeadXtalEB = (TH2F*) deadXtalFile->Get("rms_EB_r"); // #eta on x #phi on y                                                                                                
  hDeadXtalEEm = (TH2F*) deadXtalFile->Get("rms_EEm");
  hDeadXtalEEp = (TH2F*) deadXtalFile->Get("rms_EEp");

  if (!hDeadXtalEB || hDeadXtalEB == NULL) {
    cout << "Error: histogram rms_EB not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEm || hDeadXtalEEm == NULL) {
    cout << "Error: histogram rms_EEm not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEp || hDeadXtalEEp == NULL) {
    cout << "Error: histogram rms_EEp not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }



  while (reader.Next()) {

    //reader.SetLocalEntry(nEvents);
    
    cout.flush();
    if(nEvents % CHECK_EVERY_N == 0) {
      cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      //cout << "entry : " << nEvents << endl;     
    }
    nEvents++;

    Bool_t passL1seedExpression2016 = false;
    Bool_t passL1seed_SingleEGX = false;                                                                                                              
    Bool_t passL1seed_SingleIsoEGX = false;                                                                                                           
    Bool_t passL1seed_SingleIsoEGXer2p1 = false;                                                                                                      
    Bool_t passL1seed_DoubleEGX = false;                                                                                                              
    Bool_t passL1seed_SingleJetX = false;                                                                                                             
    Bool_t passL1seed_DoubleJetXer3p0 = false;                                                                                                        
    Bool_t passL1seed_TripleJetX = false;                                                                                                             
    Bool_t passL1seed_QuadJetXer3p0 = false;                                                                                                          
    Bool_t passL1seed_HTTXer = false;

    if ( *L1_SingleEG5 || *L1_SingleEG10 || *L1_SingleEG15 || *L1_SingleEG18 || *L1_SingleEG24 || *L1_SingleEG26 || *L1_SingleEG28 || *L1_SingleEG30 || *L1_SingleEG32 || *L1_SingleEG34 || *L1_SingleEG36 || *L1_SingleEG38 || *L1_SingleEG40 || *L1_SingleEG45 ) 
      {
  	passL1seed_SingleEGX = true;
  	//nPassL1seed_SingleEGX++;
      }
    
    if ( *L1_SingleIsoEG18 || *L1_SingleIsoEG20 || *L1_SingleIsoEG22 || *L1_SingleIsoEG24 || *L1_SingleIsoEG26 || *L1_SingleIsoEG28 || *L1_SingleIsoEG30 || *L1_SingleIsoEG32 || *L1_SingleIsoEG34 || *L1_SingleIsoEG36 ) 
      {
  	passL1seed_SingleIsoEGX = true;
  	//nPassL1seed_SingleIsoEGX++;
      }

    if ( *L1_SingleIsoEG18er2p1 || *L1_SingleIsoEG20er2p1 || *L1_SingleIsoEG22er2p1 || *L1_SingleIsoEG24er2p1 || *L1_SingleIsoEG26er2p1 || *L1_SingleIsoEG28er2p1 || *L1_SingleIsoEG30er2p1 || *L1_SingleIsoEG32er2p1 || *L1_SingleIsoEG34er2p1 ) 
      {
  	passL1seed_SingleIsoEGXer2p1 = true;
  	//nPassL1seed_SingleIsoEGXer2p1++;
      }

    if ( *L1_DoubleEG_15_10 || *L1_DoubleEG_18_17 || *L1_DoubleEG_20_18 || *L1_DoubleEG_22_10 || *L1_DoubleEG_23_10 || *L1_DoubleEG_22_12 || *L1_DoubleEG_22_15 || *L1_DoubleEG_24_17 || *L1_DoubleEG_25_12 ) 
      {
  	passL1seed_DoubleEGX = true;
  	//nPassL1seed_DoubleEGX++;
      }

    if ( *L1_SingleJet16 || *L1_SingleJet20 || *L1_SingleJet35 || *L1_SingleJet60 || *L1_SingleJet90 || *L1_SingleJet120 || *L1_SingleJet140 || *L1_SingleJet150 || *L1_SingleJet160 || *L1_SingleJet170 || *L1_SingleJet180 || *L1_SingleJet200 ) 
      {
  	passL1seed_SingleJetX = true;
  	//nPassL1seed_SingleJetX++;
      }

    if ( *L1_DoubleJet40er3p0 || *L1_DoubleJet50er3p0 || *L1_DoubleJet60er3p0 || *L1_DoubleJet80er3p0 || *L1_DoubleJet100er3p0 || *L1_DoubleJet112er3p0 || *L1_DoubleJet120er3p0 ) 
      {
  	passL1seed_DoubleJetXer3p0 = true;
  	//nPassL1seed_DoubleJetXer3p0++;
      }

    if ( *L1_TripleJet_88_72_56_VBF || *L1_TripleJet_84_68_48_VBF || *L1_TripleJet_92_76_64_VBF ) 
      {
  	passL1seed_TripleJetX = true;
  	//nPassL1seed_TripleJetX++;
      }

    if ( *L1_QuadJet40er3p0 || *L1_QuadJet50er3p0 || *L1_QuadJet60er3p0 ) 
      {
  	passL1seed_QuadJetXer3p0 = true;
  	//nPassL1seed_QuadJetXer3p0++;
      }

    if ( *L1_HTT120er || *L1_HTT160er || *L1_HTT200er || *L1_HTT240er || *L1_HTT255er || *L1_HTT270er || *L1_HTT280er || *L1_HTT300er || *L1_HTT320er || *L1_HTT220er ) 
      {
  	passL1seed_HTTXer = true;
  	//nPassL1seed_HTTXer++;
      }

    if ( *L1_AlwaysTrue || *L1_IsolatedBunch || passL1seed_SingleEGX || passL1seed_SingleIsoEGX || passL1seed_SingleIsoEGXer2p1 || passL1seed_DoubleEGX || passL1seed_SingleJetX || passL1seed_DoubleJetXer3p0 || passL1seed_TripleJetX || passL1seed_QuadJetXer3p0 || passL1seed_HTTXer ) 
      {
  	passL1seedExpression2016 = true;
  	//nPassL1seedExpression2016++;
      }

    

    if (not passL1seedExpression2016) continue;
    if (not **HLT_EB && not **HLT_EE) continue;

    for (Int_t i = 0; i < *STr2_NPi0_rec; i++) {

      // regional thresholds, to be set below
      double ptPairThrStream = 0.0;
      double ptGammaThrStream = 0.0;
      double s4s9ThrStream = 0.0;
      int nXtalThrStream = 0;

      double ptPairThrCalib = 0.0;
      double ptGammaThrCalib = 0.0;
      double s4s9ThrCalib = 0.0;
      int nXtalThrCalib = 0;

      double ptPairThrCalibOptim = 0.0;
      double ptGammaThrCalibOptim = 0.0;
      double s4s9ThrCalibOptim = 0.0;
      int nXtalThrCalibOptim = 0;
      double ptPairOverMThrCalibOptim = 0.0;
      double clusIsoThrCalibOptim = 0.0;

      double pi0mass = STr2_mPi0_nocor[i];
      double pi0Pt = STr2_ptPi0_nocor[i];
      double ptOverM = pi0Pt/pi0mass;

      // note: quantities in the ntuples are ordered by energy of photons (actually energy of the seed)
      // eta, phi and other variables are associate to photon before correcting the energy
      // the actual order in energy of two photons can be swapped if a correction on one is high (or if the seed ordering is not equal to the cluster energy one)
      // therefore it is possible that ptGam1 < ptGam2, and therefore we must correct the order if ptGam must be the leading pT photon
      double ptGam1 = STr2_enG1_nocor[i]/cosh(STr2_Eta_1[i]);
      double ptGam2 = STr2_enG2_nocor[i]/cosh(STr2_Eta_2[i]);
      double ptGammaMin = min(ptGam1, ptGam2);
      double s4s9min = min(STr2_S4S9_1[i],STr2_S4S9_2[i]);
      int nXtalMin = min(STr2_n1CrisPi0_rec[i],STr2_n2CrisPi0_rec[i]);
      bool passStreamSel2016 = false;
      bool passStreamSel2017 = false;
      bool passStreamSel2016_noPtCuts = false;
      bool passStreamSel2017_noPtCuts = false;
      bool passCalibSel2016 = false;
      bool passCalibSel2017 = false;
      bool passCalibSel2017LowPtGam = false;
      bool passCalibSel2017Optim = false;
      bool passCalibSel2017HLTiso = false;
      
      if (STr2_Pi0recIsEB[i]) {

	if (pi0mass < massMinEB || pi0mass > massMaxEB) continue;

	if (fabs(STr2_etaPi0_rec[i]) < EBetaRegionBoundary) {
	  nRegion = 0;
	  //stream
	  ptPairThrStream = 2.0;
	  ptGammaThrStream = 0.65;	  
	  s4s9ThrStream = 0.88;
	  nXtalThrStream = 6;
	  //calib
	  ptPairThrCalib = 2.6;
	  ptGammaThrCalib = 1.3;	  
	  s4s9ThrCalib = 0.83;
	  nXtalThrCalib = 7;	

	  // under test
	  ptPairThrCalibOptim = 2.0;
	  ptGammaThrCalibOptim = 0.65;	  
	  s4s9ThrCalibOptim = 0.88;
	  nXtalThrCalibOptim = 5;	
	  ptPairOverMThrCalibOptim = 0.0; //25.0;
	  clusIsoThrCalibOptim = 0.2;

	} else {
	  if (nXtalMin < 4) continue; // for region 2 in EB I forgot to set the minimum nXtal to 4 when producing the nuples, so I require it here 
	  nRegion = 1;
	  //stream
	  ptPairThrStream = 1.75;
	  ptGammaThrStream = 0.65;	  
	  s4s9ThrStream = 0.9;
	  nXtalThrStream = 6;
	  //calib
	  ptPairThrCalib = 2.6;
	  ptGammaThrCalib = 1.3;	  
	  s4s9ThrCalib = 0.83;
	  nXtalThrCalib = 7;

	  // under test
	  ptPairThrCalibOptim = 1.5;
	  ptGammaThrCalibOptim = 0.65;	  
	  s4s9ThrCalibOptim = 0.9;
	  nXtalThrCalibOptim = 5;	
	  ptPairOverMThrCalibOptim = 0.0; //20.0;
	  clusIsoThrCalibOptim = 0.2;

	}

      } else {
    
        if (pi0mass < massMinEE || pi0mass > massMaxEE) continue;
	
	// there should be 3 regions, but region3 has same thresholds as region2 
	if (fabs(STr2_etaPi0_rec[i]) < EEetaRegionBoundary) {
	  nRegion = 2;
	  //stream
	  ptPairThrStream = 3.75;
	  ptGammaThrStream = 1.1;	  
	  s4s9ThrStream = 0.85;
	  nXtalThrStream = 6;
	  //calib	  
	  ptPairThrCalib = 3.75;
	  ptGammaThrCalib = 1.1;	  
	  s4s9ThrCalib = 0.95;
	  nXtalThrCalib = 7;         

	  // under test
	  ptPairThrCalibOptim = 2.0;
	  ptGammaThrCalibOptim = 0.6;	  
	  s4s9ThrCalibOptim = 0.92;
	  nXtalThrCalibOptim = 5;	
	  ptPairOverMThrCalibOptim = 0.0; //30.0;
	  clusIsoThrCalibOptim = 0.2;

	} else {
	  nRegion = 3;
	  //stream
	  ptPairThrStream = 2.0;
	  ptGammaThrStream = 0.95;	  
	  s4s9ThrStream = 0.92;
	  nXtalThrStream = 6;
	  //calib
	  ptPairThrCalib = 2.0;
	  ptGammaThrCalib = 0.95;	  
	  s4s9ThrCalib = 0.95;
	  nXtalThrCalib = 7;

	  // under test
	  ptPairThrCalibOptim = 2.0;
	  ptGammaThrCalibOptim = 0.6;	  
	  s4s9ThrCalibOptim = 0.92;
	  nXtalThrCalibOptim = 5;	
	  ptPairOverMThrCalibOptim = 0.0;// 25.0;
	  clusIsoThrCalibOptim = 0.15;

	  if (fabs(STr2_etaPi0_rec[i]) > 2.0) nRegion = 4;
	}

      }
      
      if (STr2_HLTIsoPi0_rec[i] < 0.5 && s4s9min > s4s9ThrStream) passStreamSel2016_noPtCuts = true;
      if (passStreamSel2016_noPtCuts && ptGammaMin > ptGammaThrStream && pi0Pt > ptPairThrStream) passStreamSel2016 = true;
      if (passStreamSel2016_noPtCuts && nXtalMin >= nXtalThrStream) passStreamSel2017_noPtCuts = true;
      if (passStreamSel2016 && nXtalMin >= nXtalThrStream) passStreamSel2017 = true;
      if (STr2_IsoPi0_rec[i] > 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib) passCalibSel2016 = true;
      if (passCalibSel2016 && nXtalMin >= nXtalThrCalib) passCalibSel2017 = true;	
      if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017HLTiso = true;
      if (STr2_IsoPi0_rec[i] > 0.5 && ptGammaMin > 0.8 && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017LowPtGam = true;
 
      if (STr2_IsoPi0_rec[i] > clusIsoThrCalibOptim && ptGammaMin > ptGammaThrCalibOptim && s4s9min > s4s9ThrCalibOptim && pi0Pt > ptPairThrCalibOptim && nXtalMin >= nXtalThrCalibOptim && ptOverM > ptPairOverMThrCalibOptim) passCalibSel2017Optim = true;

      if (nRegion == 0) {
	if (passStreamSel2016 && nXtalMin >= 7 && s4s9min >= 0.88 && STr2_IsoPi0_rec[i] >= 0.2) {
	  hptPair_2017SelCalib_optim2018[nRegion]->Fill(pi0Pt);
	}
      }

      hmass_ntpSel[nRegion]->Fill(pi0mass);
      hnXtal1_ntpSel[nRegion]->Fill(STr2_n1CrisPi0_rec[i]);
      hnXtal2_ntpSel[nRegion]->Fill(STr2_n2CrisPi0_rec[i]);
      hptGam1_ntpSel[nRegion]->Fill(ptGam1);
      hptGam2_ntpSel[nRegion]->Fill(ptGam2);
      hptPair_ntpSel[nRegion]->Fill(pi0Pt);

      if (passStreamSel2016_noPtCuts) {
	hptGam1_2016Sel_noPtCuts[nRegion]->Fill(ptGam1);
	hptGam2_2016Sel_noPtCuts[nRegion]->Fill(ptGam2);
	hptPair_2016Sel_noPtCuts[nRegion]->Fill(pi0Pt);
	if (passStreamSel2017_noPtCuts) {
	  hptGam1_2017Sel_noPtCuts[nRegion]->Fill(ptGam1);
	  hptGam2_2017Sel_noPtCuts[nRegion]->Fill(ptGam2);
	  hptPair_2017Sel_noPtCuts[nRegion]->Fill(pi0Pt);
	}
      }

      if (passStreamSel2016) {
	hmass_2016Sel[nRegion]->Fill(pi0mass);
	if (not passL1seed_DoubleEGX) hmass_2016Sel_notDoubleEGX[nRegion]->Fill(pi0mass);
	if (not passL1seed_SingleEGX) hmass_2016Sel_notSingleEGX[nRegion]->Fill(pi0mass);

	// barrel only
	if (nRegion < 2) {
	  if (STr2_iPhi_1on20[i] > 2 && STr2_iEta_1on2520[i] > 2 && STr2_iPhi_2on20[i] > 2 && STr2_iEta_2on2520[i] > 2 && 
	      noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEB,STr2_iEtaiX_1[i]+86,STr2_iPhiiY_1[i]) && 
	      noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEB,STr2_iEtaiX_2[i]+86,STr2_iPhiiY_2[i])
	      ) { 
	    hmass_2016Sel_1nearGap[nRegion]->Fill(pi0mass);	    
	    hnXtal1_2016Sel_1nearGap[nRegion]->Fill(STr2_n1CrisPi0_rec[i]);	    
	    hnXtal2_2016Sel_1nearGap[nRegion]->Fill(STr2_n2CrisPi0_rec[i]);	    
	  } else {
	    hmass_2016Sel_1far1nearGap[nRegion]->Fill(pi0mass);	    	    
	    hnXtal1_2016Sel_1far1nearGap[nRegion]->Fill(STr2_n1CrisPi0_rec[i]);	    	    
	    hnXtal2_2016Sel_1far1nearGap[nRegion]->Fill(STr2_n2CrisPi0_rec[i]);	    	    
	  }
	}

	hnXtal1_2016Sel[nRegion]->Fill(STr2_n1CrisPi0_rec[i]);
	hnXtal2_2016Sel[nRegion]->Fill(STr2_n2CrisPi0_rec[i]);

      }

      //if (passStreamSel2016) {
      if ((selForXtalStudy == "streamSel2016" && passStreamSel2016) || (selForXtalStudy == "calibSel2016" && passCalibSel2016) ) {
	int bin = getHistIndexByXY_int(STr2_n1CrisPi0_rec[i], STr2_n2CrisPi0_rec[i], nXtal1Bins, nXtal2Bins);
	if (bin >= 0) hmass_nXtal1_nXtal2[nRegion][bin]->Fill(pi0mass);
      }

      if (passStreamSel2017) {
	hmass_2017Sel[nRegion]->Fill(pi0mass);
	if (passL1seed_SingleEGX) hmass_2017Sel_SingleEGX[nRegion]->Fill(pi0mass);
	if (passL1seed_SingleIsoEGX) hmass_2017Sel_SingleIsoEGX[nRegion]->Fill(pi0mass);
	if (passL1seed_SingleIsoEGXer2p1) hmass_2017Sel_SingleIsoEGXer2p1[nRegion]->Fill(pi0mass);
	if (passL1seed_DoubleEGX) hmass_2017Sel_DoubleEGX[nRegion]->Fill(pi0mass);
	if (passL1seed_SingleJetX) hmass_2017Sel_SingleJetX[nRegion]->Fill(pi0mass);
	if (passL1seed_DoubleJetXer3p0) hmass_2017Sel_DoubleJetXer3p0[nRegion]->Fill(pi0mass);
	if (passL1seed_TripleJetX) hmass_2017Sel_TripleJetX[nRegion]->Fill(pi0mass);
	if (passL1seed_QuadJetXer3p0) hmass_2017Sel_QuadJetXer3p0[nRegion]->Fill(pi0mass);
	if (passL1seed_HTTXer) hmass_2017Sel_HTTXer[nRegion]->Fill(pi0mass);
	if (not passL1seed_DoubleEGX) hmass_2017Sel_notDoubleEGX[nRegion]->Fill(pi0mass);
	if (not passL1seed_SingleEGX && not passL1seed_SingleIsoEGX && not passL1seed_SingleIsoEGXer2p1) hmass_2017Sel_notSingleEGX[nRegion]->Fill(pi0mass);

      }
      if (passCalibSel2016) {
	hmass_2016SelCalib[nRegion]->Fill(pi0mass);
	hnXtal1_2016SelCalib[nRegion]->Fill(STr2_n1CrisPi0_rec[i]);
	hnXtal2_2016SelCalib[nRegion]->Fill(STr2_n2CrisPi0_rec[i]);

      }
      if (passCalibSel2017) hmass_2017SelCalib[nRegion]->Fill(pi0mass);
      if (passCalibSel2017LowPtGam) hmass_2017SelCalibLowPtGam[nRegion]->Fill(pi0mass);
      if (passCalibSel2017HLTiso) hmass_2017SelCalibHLTiso[nRegion]->Fill(pi0mass);
      if (passCalibSel2017Optim) {
	hmass_2017SelCalibOptim5nXtal[nRegion]->Fill(pi0mass);
	if (STr2_n1CrisPi0_rec[i] >= 6) hmass_2017SelCalibOptim65nXtal[nRegion]->Fill(pi0mass);
	if (STr2_n1CrisPi0_rec[i] >= 7) hmass_2017SelCalibOptim75nXtal[nRegion]->Fill(pi0mass);
	if (nXtalMin >= 6) {
	  hmass_2017SelCalibOptim6nXtal[nRegion]->Fill(pi0mass);	
	  if (STr2_n1CrisPi0_rec[i] >= 7) hmass_2017SelCalibOptim76nXtal[nRegion]->Fill(pi0mass);
	}
	if (nXtalMin >= 7) {
	  hmass_2017SelCalibOptim7nXtal[nRegion]->Fill(pi0mass);
	  if (nRegion == 0) 
	    henergyGam1_EBreg1_calibOptim->Fill(STr2_enG1_nocor[i]);
	  else if (nRegion == 1) 
	    henergyGam1_EBreg2_calibOptim->Fill(STr2_enG1_nocor[i]);
	}
      }

      // base calib sel, no Pt cuts
      // since HLTiso < 0.5 already in the ntuples, just take or or STr2_IsoPi0_rec[i] > 0.5 and useHLTisoCalibForComparison
      // if latter is true, than it is like not cutting on cluster isolation
      // if it is false, then the cluster isolation selection decides whether to keep or reject the event

      double s4s9CutToUse = 0.0;
      int nXtalCutToUse = 0;

      if (useStream2017ForComparison) {
	s4s9CutToUse = s4s9ThrStream;
	nXtalCutToUse = nXtalThrStream;
      } else {
	s4s9CutToUse = s4s9ThrCalib;
	nXtalCutToUse = nXtalThrCalib;	
      }

      if ( (useStream2017ForComparison || (useHLTisoCalibForComparison || STr2_IsoPi0_rec[i] > 0.5))) {

	if (nXtalMin >= nXtalCutToUse && s4s9min > s4s9CutToUse) {

	    for (UInt_t ipt = 0; ipt < ptGamCut.size(); ipt++) {
	      if (ptGammaMin > ptGamCut[ipt]) {
		hmass_2017SelCalib_ptGam[nRegion][ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairCut.size(); ipt++) {
	      if (pi0Pt > ptPairCut[ipt]) {
		hmass_2017SelCalib_ptPair[nRegion][ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairOverMCut.size(); ipt++) {
	      if ( ptOverM > ptPairOverMCut[ipt]) {
		hmass_2017SelCalib_ptPairOverM[nRegion][ipt]->Fill(pi0mass);
	      }
	    }
	    	    
	}

	if (s4s9min > s4s9CutToUse) {

	  for (UInt_t ipt = 0; ipt < nXtalCut.size(); ipt++) {
	    if ( nXtalMin >= nXtalCut[ipt]) {
	      hmass_2017SelCalib_nXtal[nRegion][ipt]->Fill(pi0mass);
	    }
	  }
	  
	}
	
	if (nXtalMin >= nXtalCutToUse) {

	  for (UInt_t ipt = 0; ipt < s4s9Cut.size(); ipt++) {
	    if ( s4s9min >= s4s9Cut[ipt]) {
	      hmass_2017SelCalib_s4s9[nRegion][ipt]->Fill(pi0mass);
	    }
	  }

	}

      }

      // for clus iso, if comparing to stream2017, apply that selection, otherwise, just cut on Nxtal and S4S9 (that will be those of the calibration)
      if ( (useStream2017ForComparison && passStreamSel2017) || (nXtalMin >= nXtalCutToUse && s4s9min > s4s9CutToUse)) {

	for (UInt_t ipt = 0; ipt < clusIsoCut.size(); ipt++) {
	  if (STr2_IsoPi0_rec[i] > clusIsoCut[ipt]) {
	    hmass_2017SelCalib_clusIso[nRegion][ipt]->Fill(pi0mass);
	  }
	}
	  
      }


    }

  }

  cout << endl;

  // it seems that the first time CMS_lumi is used the settings are screwed up
  // produce a dummy plot (either do not save it or remove it)                
  //double lumi = 0.18; //in fb-1              
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,Form("%.2f",lumi),false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;


  vector<Int_t> *rebinPtr;
  vector<Int_t> *rebinOptimPtr;
  vector<Int_t> *rebinBisPtr;

  //  draw_nTH1(vm_energyGam1Distr->getHistPtrs(),"energy(#gamma1) (leading seed) [GeV]","a.u.", "energyGam1_"+name+"_optimSel7xtal", outDir, vm_energyGam1Distr->getHistLegends(), "", lumi, 1, true, false);

  ///////////////
  /// hardcoded distribution
  vector<TH1*> hist; hist.push_back(hptPair_2017SelCalib_optim2018[0]);
  vector<string> histLeg; histLeg.push_back("EB |#eta|<1.0 sel 2018");
  vector<Int_t> histRebin; histRebin.push_back(1);
  draw_nTH1(hist,"pT(#gamma#gamma) [GeV]","Events","ptPair_region1EB_2017SelCalib_optim2018",outDir,histLeg,"",lumi,1,false,false,histRebin);

  TFile*f = TFile::Open("hptPair_2017SelCalib_optim2018.root","RECREATE");
  if (f->IsOpen()) {
    f->cd();
    hptPair_2017SelCalib_optim2018[0]->Write();
  } 
  f->Close();

  
  for (UInt_t i = 0; i < regionTagName.size(); i++) {

    if (i == 0) rebinPtr = &rebinFactorEBin; 
    else if (i == 1) rebinPtr = &rebinFactorEB;
    else if (i > 1) rebinPtr = &rebinFactorEE;

    if (i == 0) rebinOptimPtr = &rebinFactorOptimEBin; 
    else if (i == 1) rebinOptimPtr = &rebinFactorOptimEB;
    else if (i > 1) rebinOptimPtr = &rebinFactorOptimEE;

    if (i == 0) rebinBisPtr = &rebinFactorBisEBin; 
    else if (i == 1) rebinBisPtr = &rebinFactorBisEB;
    else if (i > 1) rebinBisPtr = &rebinFactorBisEE;

    int globalRebin = -1; 
    if (i == 0) globalRebin = 2;
    else if (i == 1) globalRebin = 4;
    else globalRebin = 8;

    string name = regionTagName[i];

    draw_nTH1(vm[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name, outDir, vm[i]->getHistLegends(), "", lumi, 1, false, false, *rebinPtr);
    draw_nTH1(vm_bis[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_bis", outDir, vm_bis[i]->getHistLegends(), "", lumi, 1, false, false, *rebinBisPtr);
    draw_nTH1(vmOptim[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_optim", outDir, vmOptim[i]->getHistLegends(), "", lumi, 1, false, false, *rebinOptimPtr);
    draw_nTH1(vmOptimNxtal[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_optimNxtal", outDir, vmOptimNxtal[i]->getHistLegends(), "", lumi, 1, false, false, *rebinOptimPtr);

    draw_nTH1(vm_L1seedEG[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_L1seedEG", outDir, vm_L1seedEG[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_L1seedJet[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_L1seedJet", outDir, vm_L1seedJet[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_L1seed_SingleOrDoubleEG[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_L1seed_SingleEGXorDoubleEGX", outDir, vm_L1seed_SingleOrDoubleEG[i]->getHistLegends(), "", lumi, globalRebin, false, false);

    draw_nTH1(vm_ptGam[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_ptGam", outDir, vm_ptGam[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_ptPair[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_ptPair", outDir, vm_ptPair[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_ptPairOverM[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_ptPairOverM", outDir, vm_ptPairOverM[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_nXtal[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_nXtal", outDir, vm_nXtal[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_clusIso[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_clusIso", outDir, vm_clusIso[i]->getHistLegends(), "", lumi, globalRebin, false, false);
    draw_nTH1(vm_s4s9[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_s4s9", outDir, vm_s4s9[i]->getHistLegends(), "", lumi, globalRebin, false, false);

    draw_nTH1(vm_nXtalDistr[i]->getHistPtrs(),"number of crystals","Events", "nXtal_"+name, outDir, vm_nXtalDistr[i]->getHistLegends(), "", lumi, 1, false, false);

    draw_nTH1(vm_ptGam1Distr[i]->getHistPtrs(),"p_{T}(#gamma1) (leading seed) [GeV]","Events", "ptGam1_leadingSeed_"+name, outDir, vm_ptGam1Distr[i]->getHistLegends(), "", lumi, 1, true, false);
    draw_nTH1(vm_ptGam2Distr[i]->getHistPtrs(),"p_{T}(#gamma2) (trailing seed) [GeV]","Events", "ptGam2_trailingSeed_"+name, outDir, vm_ptGam2Distr[i]->getHistLegends(), "", lumi, 1, true, false);
    draw_nTH1(vm_ptPairDistr[i]->getHistPtrs(),"p_{T}(#gamma#gamma) [GeV]","Events", "ptPair_"+name, outDir, vm_ptPairDistr[i]->getHistLegends(), "", lumi, 1, true, false);

    if (i < 2) {
      draw_nTH1(vm_nearfar[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+name+"_nearFarGaps", outDir, vm_nearfar[i]->getHistLegends(), "", lumi, 1, false, false);     
      draw_nTH1(vm_nearfar_nxtal[i]->getHistPtrs(),"number of crystals","a.u.", "nXtal_"+name+"_nearFarGaps", outDir, vm_nearfar_nxtal[i]->getHistLegends(), "", lumi, 1, false, false);     
    }

    // plots for nXtal only up to outer EE (hard to see a peak in high eta EE)
    if (i < 3) {

      for (Int_t xtal1bin = 0; xtal1bin < nBinsXtal1; xtal1bin++) {
	for (Int_t xtal2bin = 0; xtal2bin < nBinsXtal2; xtal2bin++) {
	  
	  int bin = xtal1bin + nBinsXtal2 * xtal2bin;
	  Pi0FitResult res = drawHisto(hmass_nXtal1_nXtal2[i][bin],(i<2),outDir+"nXtalStudy/"+selForXtalStudy+"/",
				       Form("pi0Mass_%s_nXtalLead%d_nXtalTrail%d",name.c_str(),nXtal1Bins[xtal1bin],nXtal2Bins[xtal2bin]),lumi);
	  // double mean = res.mean;
	  // cout << "res.mean = " << res.mean << endl;
	  // cout << "### CHECK 1" << endl;
	  h2_massPeak_nXtal1_nXtal2[i]->SetBinContent(xtal1bin+1, xtal2bin+1, res.mean);
	  h2_massSigma_nXtal1_nXtal2[i]->SetBinContent(xtal1bin+1, xtal2bin+1, res.sigma);
	  h2_massReso_nXtal1_nXtal2[i]->SetBinContent(xtal1bin+1, xtal2bin+1, 100.0 * res.sigma/res.mean);
	  //h2_massPeak_nXtal1_nXtal2[i]->Fill(nXtal1Bins[xtal1bin]+0.5, nXtal2Bins[xtal2bin]+0.5, mean);
	  //cout << "### CHECK 2" << endl;

	}	
      }
      drawCorrelationPlot(h2_massPeak_nXtal1_nXtal2[i],
      			  "nXtal leading photon (seed)", "nXtal trailing photon (seed)","#gamma#gamma mass peak (GeV/c^{2})",
			  "pi0Mass_peak_"+name+"_nXtal1_nXtal2",
			  plotTag[i],outDir+"nXtalStudy/"+selForXtalStudy+"/" ,1,false,false);
      drawCorrelationPlot(h2_massSigma_nXtal1_nXtal2[i],
      			  "nXtal leading photon (seed)", "nXtal trailing photon (seed)","#gamma#gamma mass sigma (MeV/c^{2})",
			  "pi0Mass_sigma_"+name+"_nXtal1_nXtal2",
			  plotTag[i],outDir+"nXtalStudy/"+selForXtalStudy+"/",1,false,false);
      drawCorrelationPlot(h2_massReso_nXtal1_nXtal2[i],
      			  "nXtal leading photon (seed)", "nXtal trailing photon (seed)","#gamma#gamma mass resolution (#sigma/#mu) [%]",
			  "pi0Mass_reso_"+name+"_nXtal1_nXtal2",
			  plotTag[i],outDir+"nXtalStudy/"+selForXtalStudy+"/",1,false,false);

    }

  }

  cout << "THE END!" << endl;

}
