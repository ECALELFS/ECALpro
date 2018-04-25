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
#include "RooMinuit.h"

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

static bool useHLTisoCalibForComparison = false;


// note
// in EB I would use pt(pi0) > 2.0-2.5 and pT/M > 25-30, must see other cuts, and must also reduce cluster iso to 0.1 or 0.2 instead of 0.5

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

  void addComponent(TH1* histPtr = NULL, const string& histName = "name", const string& histLegend = "leg", const Int_t& histRebins = 1)  { 
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

//=================================================


void makeStreamOptim(const bool isEB = true,
		     const bool isPi0 = true,
		     const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/Optimization/streamSelection/",
		     const string& eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
		     ) 
{

  createPlotDirAndCopyPhp(outDir);

  // might want to pass it as an option (comma separated list to be parsed, for example)
  vector<string> eosDir;
  if (isPi0) {
    eosDir.push_back("AlCaP0_fromHLTPhysics_2017AB_TreeOptim");
    eosDir.push_back("AlCaP0_from_ParkingHLTPhysics_TreeOptim");
  } else {
    eosDir.push_back("AlCaEta_fromHLTPhysics_2017AB_TreeOptim");
    eosDir.push_back("AlCaEta_fromParkingHLTPhysics_2017AB_TreeOptim");
  }

  double lumi = 7.5; // to be updated

  TChain* chain = new TChain("Tree_Optim");

  for (UInt_t i = 0; i < eosDir.size(); i++) {

    string eosPathComplete = eosPath;
    if (eosPathComplete.back() != '/') eosPathComplete += "/";
    eosPathComplete += eosDir[i];
    if (eosPathComplete.back() != '/') eosPathComplete += "/";
    eosPathComplete += "/iter_0/";

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
	chain->Add((eosPathComplete+line).c_str());      
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
  float massMinEE = isPi0 ? 0.05 : 0.2;
  float massMaxEE = isPi0 ? 0.25 : 0.8;


  ///////////////////////
  // histograms

  // EB
  TH1F* hmass_allEB_ntpSel = new TH1F("hmass_allEB_ntpSel","",160,0.060,0.220);
  TH1F* hmass_allEB_2016Sel = new TH1F("hmass_allEB_2016Sel","",160,0.060,0.220); // selection used in 2016
  TH1F* hmass_allEB_2017Sel = new TH1F("hmass_allEB_2017Sel","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_allEB_2017SelCalib = new TH1F("hmass_allEB_2017SelCalib","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)
  TH1F* hmass_allEB_2017SelCalibHLTiso = new TH1F("hmass_allEB_2017SelCalibHLTiso","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)

  vectorManager* vm_allEB = new vectorManager();
  vm_allEB->addComponent(hmass_allEB_ntpSel, "hmass_allEB_ntpSel", "all EB loose");
  vm_allEB->addComponent(hmass_allEB_2016Sel, "hmass_allEB_2016Sel", "all EB 2016");
  vm_allEB->addComponent(hmass_allEB_2017Sel, "hmass_allEB_2017Sel", "all EB 2017");
  vm_allEB->addComponent(hmass_allEB_2017SelCalib, "hmass_allEB_2017SelCalib", "all EB 2017 calib");
  vm_allEB->addComponent(hmass_allEB_2017SelCalibHLTiso, "hmass_allEB_2017SelCalibHLTiso", "EB 2017 calib HLTiso");

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
  plotTag.push_back("EB |#eta| < 1.0");
  plotTag.push_back("EB |#eta| > 1.0");
  plotTag.push_back("EE |#eta| < 1.8");
  plotTag.push_back("EE 1.8 < |#eta| < 2.0");
  plotTag.push_back("EE |#eta| > 2.0");

  vector <TH1F*> hmass_ntpSel;
  vector <TH1F*> hmass_2016Sel;
  vector <TH1F*> hmass_2017Sel;
  vector <TH1F*> hmass_2017SelCalib;
  vector <TH1F*> hmass_2017SelCalibHLTiso;

  vector<vectorManager*> vm;
  
  for (UInt_t i = 0; i < regionTagName.size(); i++) {
    string name = regionTagName[i];
    double massMin = (i < 2) ? massMinEB : massMinEE;
    double massMax = (i < 2) ? massMaxEB : massMaxEE;

    hmass_ntpSel.push_back( new TH1F(Form("hmass_%s_ntpSel",name.c_str()),"",160,0.060,0.220) );
    hmass_2016Sel.push_back( new TH1F(Form("hmass_%s_2016Sel",name.c_str()),"",160,0.060,0.220) );
    hmass_2017Sel.push_back( new TH1F(Form("hmass_%s_2017Sel",name.c_str()),"",160,0.060,0.220) );
    hmass_2017SelCalib.push_back( new TH1F(Form("hmass_%s_2017SelCalib",name.c_str()),"",160,0.060,0.220) );
    hmass_2017SelCalibHLTiso.push_back( new TH1F(Form("hmass_%s_2017SelCalibHLTiso",name.c_str()),"",160,0.060,0.220) );   

    vm.push_back( new vectorManager() );
    vm[i]->addComponent(hmass_ntpSel[i], hmass_ntpSel[i]->GetName(), Form("%s loose",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2016Sel[i], hmass_2016Sel[i]->GetName(), Form("%s 2016",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017Sel[i], hmass_2017Sel[i]->GetName(), Form("%s 2017",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017SelCalib[i], 2017SelCalib[1]->GetName(), Form("%s 2017 calib",plotTag[i].c_str()));
    vm[i]->addComponent(hmass_2017SelCalibHLTiso[i], hmass_2017SelCalibHLTiso[i]->GetName(), "2017 calib HLTiso"); // shorter form			

  }

  // end of work in progress
  /////////////////////////////////////
  /////////////////////////////////////
  /////////////////////////////////////

  TH1F* hmass_region1EB_ntpSel = new TH1F("hmass_region1EB_ntpSel","",160,0.060,0.220);
  TH1F* hmass_region1EB_2016Sel = new TH1F("hmass_region1EB_2016Sel","",160,0.060,0.220); // selection used in 2016
  TH1F* hmass_region1EB_2017Sel = new TH1F("hmass_region1EB_2017Sel","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region1EB_2017SelCalib = new TH1F("hmass_region1EB_2017SelCalib","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)
  TH1F* hmass_region1EB_2017SelCalibHLTiso = new TH1F("hmass_region1EB_2017SelCalibHLTiso","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)

  vectorManager* vm_region1EB = new vectorManager();
  vm_region1EB->addComponent(hmass_region1EB_ntpSel, "hmass_region1EB_ntpSel", "EB, |#eta|<1.0, loose");
  vm_region1EB->addComponent(hmass_region1EB_2016Sel, "hmass_region1EB_2016Sel", "EB, |#eta|<1.0, 2016");
  vm_region1EB->addComponent(hmass_region1EB_2017Sel, "hmass_region1EB_2017Sel", "EB, |#eta|<1.0, 2017");
  vm_region1EB->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");
  vm_region1EB->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");

  TH1F* hmass_region2EB_ntpSel = new TH1F("hmass_region2EB_ntpSel","",160,0.060,0.220);
  TH1F* hmass_region2EB_2016Sel = new TH1F("hmass_region2EB_2016Sel","",160,0.060,0.220); // selection used in 2016
  TH1F* hmass_region2EB_2017Sel = new TH1F("hmass_region2EB_2017Sel","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region2EB_2017SelCalib = new TH1F("hmass_region2EB_2017SelCalib","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)
  TH1F* hmass_region2EB_2017SelCalibHLTiso = new TH1F("hmass_region2EB_2017SelCalibHLTiso","",160,0.060,0.220); // selection used in 2017 (as 2016 but with Nxtal cut >= 7)
  // used only for region 2
  vector<Int_t> rebinFactorEB;
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(1);
  rebinFactorEB.push_back(2);
  rebinFactorEB.push_back(2);

  vectorManager* vm_region2EB = new vectorManager();
  vm_region2EB->addComponent(hmass_region2EB_ntpSel, "hmass_region2EB_ntpSel", "EB, |#eta|>1.0, loose");
  vm_region2EB->addComponent(hmass_region2EB_2016Sel, "hmass_region2EB_2016Sel", "EB, |#eta|>1.0, 2016");
  vm_region2EB->addComponent(hmass_region2EB_2017Sel, "hmass_region2EB_2017Sel", "EB, |#eta|>1.0, 2017");
  vm_region2EB->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  vm_region2EB->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");


  string histName = "";
  string doubleToStr = "";

  vector<TH1F*> hmass_region1EB_2017SelCalib_ptGam;
  vector<TH1F*> hmass_region2EB_2017SelCalib_ptGam;

  vectorManager* vm_region1EB_ptGam = new vectorManager();
  vectorManager* vm_region2EB_ptGam = new vectorManager();
  if (useHLTisoCalibForComparison)  {
    vm_region1EB_ptGam->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
    vm_region2EB_ptGam->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
  } else {
    vm_region1EB_ptGam->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
    vm_region2EB_ptGam->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  }
  for (UInt_t i = 0; i < ptGamCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptGamCut[i]);
    histName = string(Form("hmass_region1EB_2017SelCalib_ptGam%s",doubleToStr.c_str()));
    hmass_region1EB_2017SelCalib_ptGam.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region1EB_ptGam->addComponent(hmass_region1EB_2017SelCalib_ptGam[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[i]));  
    histName = string(Form("hmass_region2EB_2017SelCalib_ptGam%s",doubleToStr.c_str()));
    hmass_region2EB_2017SelCalib_ptGam.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region2EB_ptGam->addComponent(hmass_region2EB_2017SelCalib_ptGam[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[i]));  
  }

  vector<TH1F*> hmass_region1EB_2017SelCalib_ptPair;
  vector<TH1F*> hmass_region2EB_2017SelCalib_ptPair;

  vectorManager* vm_region1EB_ptPair = new vectorManager();
  vectorManager* vm_region2EB_ptPair = new vectorManager();
  if (useHLTisoCalibForComparison)  {
    vm_region1EB_ptPair->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
    vm_region2EB_ptPair->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
  } else {
    vm_region1EB_ptPair->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
    vm_region2EB_ptPair->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  }
  for (UInt_t i = 0; i < ptPairCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptPairCut[i]);
    histName = string(Form("hmass_region1EB_2017SelCalib_ptPair%s",doubleToStr.c_str()));
    hmass_region1EB_2017SelCalib_ptPair.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region1EB_ptPair->addComponent(hmass_region1EB_2017SelCalib_ptPair[i], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[i]));  
    histName = string(Form("hmass_region2EB_2017SelCalib_ptPair%s",doubleToStr.c_str()));
    hmass_region2EB_2017SelCalib_ptPair.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region2EB_ptPair->addComponent(hmass_region2EB_2017SelCalib_ptPair[i], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[i]));  
  }

  vector<TH1F*> hmass_region1EB_2017SelCalib_ptPairOverM;
  vector<TH1F*> hmass_region2EB_2017SelCalib_ptPairOverM;

  vectorManager* vm_region1EB_ptPairOverM = new vectorManager();
  vectorManager* vm_region2EB_ptPairOverM = new vectorManager();
  if (useHLTisoCalibForComparison)  {
    vm_region1EB_ptPairOverM->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
    vm_region2EB_ptPairOverM->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
  } else {
    vm_region1EB_ptPairOverM->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
    vm_region2EB_ptPairOverM->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  }
  for (UInt_t i = 0; i < ptPairOverMCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptPairOverMCut[i]);
    histName = string(Form("hmass_region1EB_2017SelCalib_ptPairOverM%s",doubleToStr.c_str()));
    hmass_region1EB_2017SelCalib_ptPairOverM.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region1EB_ptPairOverM->addComponent(hmass_region1EB_2017SelCalib_ptPairOverM[i], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[i]));  
    histName = string(Form("hmass_region2EB_2017SelCalib_ptPairOverM%s",doubleToStr.c_str()));
    hmass_region2EB_2017SelCalib_ptPairOverM.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region2EB_ptPairOverM->addComponent(hmass_region2EB_2017SelCalib_ptPairOverM[i], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[i]));  
  }

  vector<TH1F*> hmass_region1EB_2017SelCalib_nXtal;
  vector<TH1F*> hmass_region2EB_2017SelCalib_nXtal;

  vectorManager* vm_region1EB_nXtal = new vectorManager();
  vectorManager* vm_region2EB_nXtal = new vectorManager();
  if (useHLTisoCalibForComparison)  {
    vm_region1EB_nXtal->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
    vm_region2EB_nXtal->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
  } else {
    vm_region1EB_nXtal->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
    vm_region2EB_nXtal->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  }
  for (UInt_t i = 0; i < nXtalCut.size(); i++) {
    histName = string(Form("hmass_region1EB_2017SelCalib_nXtal%d",nXtalCut[i]));
    hmass_region1EB_2017SelCalib_nXtal.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region1EB_nXtal->addComponent(hmass_region1EB_2017SelCalib_nXtal[i], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[i]));  
    histName = string(Form("hmass_region2EB_2017SelCalib_nXtal%d",nXtalCut[i]));
    hmass_region2EB_2017SelCalib_nXtal.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region2EB_nXtal->addComponent(hmass_region2EB_2017SelCalib_nXtal[i], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[i]));  
  }

  vector<TH1F*> hmass_region1EB_2017SelCalib_clusIso;
  vector<TH1F*> hmass_region2EB_2017SelCalib_clusIso;

  vectorManager* vm_region1EB_clusIso = new vectorManager();
  vectorManager* vm_region2EB_clusIso = new vectorManager();
  if (useHLTisoCalibForComparison)  {
    vm_region1EB_clusIso->addComponent(hmass_region1EB_2017SelCalibHLTiso, "hmass_region1EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
    vm_region2EB_clusIso->addComponent(hmass_region2EB_2017SelCalibHLTiso, "hmass_region2EB_2017SelCalibHLTiso", "EB, 2017 calib HLTiso");  
  } else {
    vm_region1EB_clusIso->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
    vm_region2EB_clusIso->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  }
  for (UInt_t i = 0; i < clusIsoCut.size(); i++) {
    doubleToStr = getStringFromDouble(clusIsoCut[i]);
    histName = string(Form("hmass_region1EB_2017SelCalib_clusIso%s",doubleToStr.c_str()));
    hmass_region1EB_2017SelCalib_clusIso.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region1EB_clusIso->addComponent(hmass_region1EB_2017SelCalib_clusIso[i], histName.c_str(), Form("clus iso > %1.1f",clusIsoCut[i]));  
    histName = string(Form("hmass_region2EB_2017SelCalib_clusIso%s",doubleToStr.c_str()));
    hmass_region2EB_2017SelCalib_clusIso.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
    vm_region2EB_clusIso->addComponent(hmass_region2EB_2017SelCalib_clusIso[i], histName.c_str(), Form("clus iso > %1.1f",clusIsoCut[i]));  
  }



  // vector<TH1F*> hmass_region1EB_2017SelCalib_s4s9;
  // vector<TH1F*> hmass_region2EB_2017SelCalib_s4s9;

  // vectorManager* vm_region1EB_s4s9 = new vectorManager();
  // vectorManager* vm_region2EB_s4s9 = new vectorManager();
  // vm_region1EB_s4s9->addComponent(hmass_region1EB_2017SelCalib, "hmass_region1EB_2017SelCalib", "EB, |#eta|<1.0, 2017 calib");  
  // vm_region2EB_s4s9->addComponent(hmass_region2EB_2017SelCalib, "hmass_region2EB_2017SelCalib", "EB, |#eta|>1.0, 2017 calib");  
  // for (UInt_t i = 0; i < s4s9Cut.size(); i++) {
  //   doubleToStr = getStringFromDouble(s4s9Cut[i]);
  //   histName = string(Form("hmass_region1EB_2017SelCalib_s4s9%s",doubleToStr.c_str()));
  //   hmass_region1EB_2017SelCalib_s4s9.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
  //   vm_region1EB_s4s9->addComponent(hmass_region1EB_2017SelCalib_s4s9[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",s4s9Cut[i]));  
  //   histName = string(Form("hmass_region2EB_2017SelCalib_s4s9%s",doubleToStr.c_str()));
  //   hmass_region2EB_2017SelCalib_s4s9.push_back(new TH1F(histName.c_str(),"",160,0.060,0.220));
  //   vm_region2EB_s4s9->addComponent(hmass_region2EB_2017SelCalib_s4s9[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",s4s9Cut[i]));  
  // }



  ///////////////////////////////////
  // EE
  TH1F* hmass_allEE_ntpSel = new TH1F("hmass_allEE_ntpSel","",200,0.050,0.250);
  TH1F* hmass_allEE_2016Sel = new TH1F("hmass_allEE_2016Sel","",200,0.050,0.250); // selection used in 2016
  TH1F* hmass_allEE_2017Sel = new TH1F("hmass_allEE_2017Sel","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_allEE_2017SelCalib = new TH1F("hmass_allEE_2017SelCalib","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_allEE_2017SelCalibHLTiso = new TH1F("hmass_allEE_2017SelCalibHLTiso","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)

  TH1F* hmass_region1EE_ntpSel = new TH1F("hmass_region1EE_ntpSel","",200,0.050,0.250);
  TH1F* hmass_region1EE_2016Sel = new TH1F("hmass_region1EE_2016Sel","",200,0.050,0.250); // selection used in 2016
  TH1F* hmass_region1EE_2017Sel = new TH1F("hmass_region1EE_2017Sel","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region1EE_2017SelCalib = new TH1F("hmass_region1EE_2017SelCalib","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region1EE_2017SelCalibHLTiso = new TH1F("hmass_region1EE_2017SelCalibHLTiso","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)

  TH1F* hmass_region2EE_ntpSel = new TH1F("hmass_region2EE_ntpSel","",200,0.050,0.250);
  TH1F* hmass_region2EE_2016Sel = new TH1F("hmass_region2EE_2016Sel","",200,0.050,0.250); // selection used in 2016
  TH1F* hmass_region2EE_2017Sel = new TH1F("hmass_region2EE_2017Sel","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region2EE_2017SelCalib = new TH1F("hmass_region2EE_2017SelCalib","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region2EE_2017SelCalibHLTiso = new TH1F("hmass_region2EE_2017SelCalibHLTiso","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)

  TH1F* hmass_region3EE_ntpSel = new TH1F("hmass_region3EE_ntpSel","",200,0.050,0.250);
  TH1F* hmass_region3EE_2016Sel = new TH1F("hmass_region3EE_2016Sel","",200,0.050,0.250); // selection used in 2016
  TH1F* hmass_region3EE_2017Sel = new TH1F("hmass_region3EE_2017Sel","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region3EE_2017SelCalib = new TH1F("hmass_region3EE_2017SelCalib","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)
  TH1F* hmass_region3EE_2017SelCalibHLTiso = new TH1F("hmass_region3EE_2017SelCalibHLTiso","",200,0.050,0.250); // selection used in 2017 (as 2016 but with Nxtal cut >= 6)


  vector<Int_t> rebinFactorEE;
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(2);
  rebinFactorEE.push_back(8);
  rebinFactorEE.push_back(5);

  vectorManager* vm_allEE = new vectorManager();
  vm_allEE->addComponent(hmass_allEE_ntpSel, "hmass_allEE_ntpSel", "all EE loose");
  vm_allEE->addComponent(hmass_allEE_2016Sel, "hmass_allEE_2016Sel", "all EE 2016");
  vm_allEE->addComponent(hmass_allEE_2017Sel, "hmass_allEE_2017Sel", "all EE 2017");
  vm_allEE->addComponent(hmass_allEE_2017SelCalib, "hmass_allEE_2017SelCalib", "all EE 2017 calib");
  vm_allEE->addComponent(hmass_allEE_2017SelCalibHLTiso, "hmass_allEE_2017SelCalibHLTiso", "EE 2017 calib HLTiso");

  vectorManager* vm_region1EE = new vectorManager();
  vm_region1EE->addComponent(hmass_region1EE_ntpSel, "hmass_region1EE_ntpSel", "EE, |#eta|<1.8, loose");
  vm_region1EE->addComponent(hmass_region1EE_2016Sel, "hmass_region1EE_2016Sel", "EE, |#eta|<1.8, 2016");
  vm_region1EE->addComponent(hmass_region1EE_2017Sel, "hmass_region1EE_2017Sel", "EE, |#eta|<1.8, 2017");
  vm_region1EE->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");
  vm_region1EE->addComponent(hmass_region1EE_2017SelCalibHLTiso, "hmass_region1EE_2017SelCalibHLTiso", "EE, 2017 calib HLTiso");

  vectorManager* vm_region2EE = new vectorManager();
  vm_region2EE->addComponent(hmass_region2EE_ntpSel, "hmass_region2EE_ntpSel", "EE, |#eta|<2.0, loose");
  vm_region2EE->addComponent(hmass_region2EE_2016Sel, "hmass_region2EE_2016Sel", "EE, |#eta|<2.0, 2016");
  vm_region2EE->addComponent(hmass_region2EE_2017Sel, "hmass_region2EE_2017Sel", "EE, |#eta|<2.0, 2017");
  vm_region2EE->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib"); 
  vm_region2EE->addComponent(hmass_region2EE_2017SelCalibHLTiso, "hmass_region2EE_2017SelCalibHLTiso", "EE, 2017 calib HLTiso"); 

  vectorManager* vm_region3EE = new vectorManager();
  vm_region3EE->addComponent(hmass_region3EE_ntpSel, "hmass_region3EE_ntpSel", "EE, |#eta|>2.0, loose");
  vm_region3EE->addComponent(hmass_region3EE_2016Sel, "hmass_region3EE_2016Sel", "EE, |#eta|>2.0, 2016");
  vm_region3EE->addComponent(hmass_region3EE_2017Sel, "hmass_region3EE_2017Sel", "EE, |#eta|>2.0, 2017");
  vm_region3EE->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");
  vm_region3EE->addComponent(hmass_region3EE_2017SelCalibHLTiso, "hmass_region3EE_2017SelCalibHLTiso", "EE, 2017 calib HLTiso");

  vector<TH1F*> hmass_region1EE_2017SelCalib_ptGam;
  vector<TH1F*> hmass_region2EE_2017SelCalib_ptGam;
  vector<TH1F*> hmass_region3EE_2017SelCalib_ptGam;

  vectorManager* vm_region1EE_ptGam = new vectorManager();
  vectorManager* vm_region2EE_ptGam = new vectorManager();
  vectorManager* vm_region3EE_ptGam = new vectorManager();
  vm_region1EE_ptGam->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  vm_region2EE_ptGam->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  vm_region3EE_ptGam->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  for (UInt_t i = 0; i < ptGamCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptGamCut[i]);
    histName = string(Form("hmass_region1EE_2017SelCalib_ptGam%s",doubleToStr.c_str()));
    hmass_region1EE_2017SelCalib_ptGam.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region1EE_ptGam->addComponent(hmass_region1EE_2017SelCalib_ptGam[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[i]));  
    histName = string(Form("hmass_region2EE_2017SelCalib_ptGam%s",doubleToStr.c_str()));
    hmass_region2EE_2017SelCalib_ptGam.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region2EE_ptGam->addComponent(hmass_region2EE_2017SelCalib_ptGam[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[i]));  
    histName = string(Form("hmass_region3EE_2017SelCalib_ptGam%s",doubleToStr.c_str()));
    hmass_region3EE_2017SelCalib_ptGam.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region3EE_ptGam->addComponent(hmass_region3EE_2017SelCalib_ptGam[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",ptGamCut[i]));  
  }

  vector<TH1F*> hmass_region1EE_2017SelCalib_ptPair;
  vector<TH1F*> hmass_region2EE_2017SelCalib_ptPair;
  vector<TH1F*> hmass_region3EE_2017SelCalib_ptPair;

  vectorManager* vm_region1EE_ptPair = new vectorManager();
  vectorManager* vm_region2EE_ptPair = new vectorManager();
  vectorManager* vm_region3EE_ptPair = new vectorManager();
  vm_region1EE_ptPair->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  vm_region2EE_ptPair->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  vm_region3EE_ptPair->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  for (UInt_t i = 0; i < ptPairCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptPairCut[i]);
    histName = string(Form("hmass_region1EE_2017SelCalib_ptPair%s",doubleToStr.c_str()));
    hmass_region1EE_2017SelCalib_ptPair.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region1EE_ptPair->addComponent(hmass_region1EE_2017SelCalib_ptPair[i], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[i]));  
    histName = string(Form("hmass_region2EE_2017SelCalib_ptPair%s",doubleToStr.c_str()));
    hmass_region2EE_2017SelCalib_ptPair.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region2EE_ptPair->addComponent(hmass_region2EE_2017SelCalib_ptPair[i], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[i]));  
    histName = string(Form("hmass_region3EE_2017SelCalib_ptPair%s",doubleToStr.c_str()));
    hmass_region3EE_2017SelCalib_ptPair.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region3EE_ptPair->addComponent(hmass_region3EE_2017SelCalib_ptPair[i], histName.c_str(), Form("p_{T}(#pi^{0}) > %1.1f",ptPairCut[i]));  
  }

  vector<TH1F*> hmass_region1EE_2017SelCalib_ptPairOverM;
  vector<TH1F*> hmass_region2EE_2017SelCalib_ptPairOverM;
  vector<TH1F*> hmass_region3EE_2017SelCalib_ptPairOverM;

  vectorManager* vm_region1EE_ptPairOverM = new vectorManager();
  vectorManager* vm_region2EE_ptPairOverM = new vectorManager();
  vectorManager* vm_region3EE_ptPairOverM = new vectorManager();
  vm_region1EE_ptPairOverM->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  vm_region2EE_ptPairOverM->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  vm_region3EE_ptPairOverM->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  for (UInt_t i = 0; i < ptPairOverMCut.size(); i++) {
    doubleToStr = getStringFromDouble(ptPairOverMCut[i]);
    histName = string(Form("hmass_region1EE_2017SelCalib_ptPairOverM%s",doubleToStr.c_str()));
    hmass_region1EE_2017SelCalib_ptPairOverM.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region1EE_ptPairOverM->addComponent(hmass_region1EE_2017SelCalib_ptPairOverM[i], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[i]));  
    histName = string(Form("hmass_region2EE_2017SelCalib_ptPairOverM%s",doubleToStr.c_str()));
    hmass_region2EE_2017SelCalib_ptPairOverM.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region2EE_ptPairOverM->addComponent(hmass_region2EE_2017SelCalib_ptPairOverM[i], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[i]));  
    histName = string(Form("hmass_region3EE_2017SelCalib_ptPairOverM%s",doubleToStr.c_str()));
    hmass_region3EE_2017SelCalib_ptPairOverM.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region3EE_ptPairOverM->addComponent(hmass_region3EE_2017SelCalib_ptPairOverM[i], histName.c_str(), Form("p_{T}/M(#pi^{0}) > %1.1f",ptPairOverMCut[i]));  
  }

  vector<TH1F*> hmass_region1EE_2017SelCalib_nXtal;
  vector<TH1F*> hmass_region2EE_2017SelCalib_nXtal;
  vector<TH1F*> hmass_region3EE_2017SelCalib_nXtal;

  vectorManager* vm_region1EE_nXtal = new vectorManager();
  vectorManager* vm_region2EE_nXtal = new vectorManager();
  vectorManager* vm_region3EE_nXtal = new vectorManager();
  vm_region1EE_nXtal->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  vm_region2EE_nXtal->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  vm_region3EE_nXtal->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  for (UInt_t i = 0; i < nXtalCut.size(); i++) {
    doubleToStr = getStringFromDouble(nXtalCut[i]);
    histName = string(Form("hmass_region1EE_2017SelCalib_nXtal%d",nXtalCut[i]));
    hmass_region1EE_2017SelCalib_nXtal.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region1EE_nXtal->addComponent(hmass_region1EE_2017SelCalib_nXtal[i], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[i]));  
    histName = string(Form("hmass_region2EE_2017SelCalib_nXtal%d",nXtalCut[i]));
    hmass_region2EE_2017SelCalib_nXtal.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region2EE_nXtal->addComponent(hmass_region2EE_2017SelCalib_nXtal[i], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[i]));  
    histName = string(Form("hmass_region3EE_2017SelCalib_nXtal%d",nXtalCut[i]));
    hmass_region3EE_2017SelCalib_nXtal.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region3EE_nXtal->addComponent(hmass_region3EE_2017SelCalib_nXtal[i], histName.c_str(), Form("n(Xtal) >= %d",nXtalCut[i]));  
  }


  vector<TH1F*> hmass_region1EE_2017SelCalib_clusIso;
  vector<TH1F*> hmass_region2EE_2017SelCalib_clusIso;
  vector<TH1F*> hmass_region3EE_2017SelCalib_clusIso;

  vectorManager* vm_region1EE_clusIso = new vectorManager();
  vectorManager* vm_region2EE_clusIso = new vectorManager();
  vectorManager* vm_region3EE_clusIso = new vectorManager();
  vm_region1EE_clusIso->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  vm_region2EE_clusIso->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  vm_region3EE_clusIso->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  for (UInt_t i = 0; i < clusIsoCut.size(); i++) {
    doubleToStr = getStringFromDouble(clusIsoCut[i]);
    histName = string(Form("hmass_region1EE_2017SelCalib_clusIso%s",doubleToStr.c_str()));
    hmass_region1EE_2017SelCalib_clusIso.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region1EE_clusIso->addComponent(hmass_region1EE_2017SelCalib_clusIso[i], histName.c_str(), Form("clus iso > %1.1f",clusIsoCut[i]));  
    histName = string(Form("hmass_region2EE_2017SelCalib_clusIso%s",doubleToStr.c_str()));
    hmass_region2EE_2017SelCalib_clusIso.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region2EE_clusIso->addComponent(hmass_region2EE_2017SelCalib_clusIso[i], histName.c_str(), Form("clus iso > %1.1f",clusIsoCut[i]));  
    histName = string(Form("hmass_region3EE_2017SelCalib_clusIso%s",doubleToStr.c_str()));
    hmass_region3EE_2017SelCalib_clusIso.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
    vm_region3EE_clusIso->addComponent(hmass_region3EE_2017SelCalib_clusIso[i], histName.c_str(), Form("clus iso > %1.1f",clusIsoCut[i]));  
  }

  // vector<TH1F*> hmass_region1EE_2017SelCalib_s4s9;
  // vector<TH1F*> hmass_region2EE_2017SelCalib_s4s9;
  // vector<TH1F*> hmass_region3EE_2017SelCalib_s4s9;

  // vectorManager* vm_region1EE_s4s9 = new vectorManager();
  // vectorManager* vm_region2EE_s4s9 = new vectorManager();
  // vectorManager* vm_region3EE_s4s9 = new vectorManager();
  // vm_region1EE_s4s9->addComponent(hmass_region1EE_2017SelCalib, "hmass_region1EE_2017SelCalib", "EE, |#eta|<1.8, 2017 calib");  
  // vm_region2EE_s4s9->addComponent(hmass_region2EE_2017SelCalib, "hmass_region2EE_2017SelCalib", "EE, |#eta|<2.0, 2017 calib");  
  // vm_region3EE_s4s9->addComponent(hmass_region3EE_2017SelCalib, "hmass_region3EE_2017SelCalib", "EE, |#eta|>2.0, 2017 calib");  
  // for (UInt_t i = 0; i < s4s9Cut.size(); i++) {
  //   doubleToStr = getStringFromDouble(s4s9Cut[i]);
  //   histName = string(Form("hmass_region1EE_2017SelCalib_s4s9%s",doubleToStr.c_str()));
  //   hmass_region1EE_2017SelCalib_s4s9.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
  //   vm_region1EE_s4s9->addComponent(hmass_region1EE_2017SelCalib_s4s9[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",s4s9Cut[i]));  
  //   histName = string(Form("hmass_region2EE_2017SelCalib_s4s9%s",doubleToStr.c_str()));
  //   hmass_region2EE_2017SelCalib_s4s9.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
  //   vm_region2EE_s4s9->addComponent(hmass_region2EE_2017SelCalib_s4s9[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",s4s9Cut[i]));  
  //   histName = string(Form("hmass_region3EE_2017SelCalib_s4s9%s",doubleToStr.c_str()));
  //   hmass_region3EE_2017SelCalib_s4s9.push_back(new TH1F(histName.c_str(),"",200,0.050,0.250));
  //   vm_region3EE_s4s9->addComponent(hmass_region3EE_2017SelCalib_s4s9[i], histName.c_str(), Form("p_{T}(#gamma) > %1.1f",s4s9Cut[i]));  
  // }

  Int_t nRegion = -1;

  while (reader.Next()) {

    //reader.SetLocalEntry(nEvents);
    
    cout.flush();
    if(nEvents % CHECK_EVERY_N == 0) {
      cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      //cout << "entry : " << nEvents << endl;     
    }
    nEvents++;

    Bool_t passL1seedExpression2016 = false;
    if ( *L1_AlwaysTrue || *L1_IsolatedBunch || *L1_SingleEG5 || *L1_SingleEG10 || *L1_SingleEG15 || *L1_SingleEG18 || *L1_SingleEG24 || *L1_SingleEG26 || *L1_SingleEG28 || *L1_SingleEG30 || *L1_SingleEG32 || *L1_SingleEG34 || *L1_SingleEG36 || *L1_SingleEG38 || *L1_SingleEG40 || *L1_SingleEG45 || *L1_SingleIsoEG18 || *L1_SingleIsoEG20 || *L1_SingleIsoEG22 || *L1_SingleIsoEG24 || *L1_SingleIsoEG26 || *L1_SingleIsoEG28 || *L1_SingleIsoEG30 || *L1_SingleIsoEG32 || *L1_SingleIsoEG34 || *L1_SingleIsoEG36 || *L1_SingleIsoEG18er2p1 || *L1_SingleIsoEG20er2p1 || *L1_SingleIsoEG22er2p1 || *L1_SingleIsoEG24er2p1 || *L1_SingleIsoEG26er2p1 || *L1_SingleIsoEG28er2p1 || *L1_SingleIsoEG30er2p1 || *L1_SingleIsoEG32er2p1 || *L1_SingleIsoEG34er2p1 || *L1_DoubleEG_15_10 || *L1_DoubleEG_18_17 || *L1_DoubleEG_20_18 || *L1_DoubleEG_22_10 || *L1_DoubleEG_23_10 || *L1_DoubleEG_22_12 || *L1_DoubleEG_22_15 || *L1_DoubleEG_24_17 || *L1_DoubleEG_25_12 ||  *L1_SingleJet16 || *L1_SingleJet20 || *L1_SingleJet35 || *L1_SingleJet60 || *L1_SingleJet90 || *L1_SingleJet120 || *L1_SingleJet140 || *L1_SingleJet150 || *L1_SingleJet160 || *L1_SingleJet170 || *L1_SingleJet180 || *L1_SingleJet200 || *L1_DoubleJet40er3p0 || *L1_DoubleJet50er3p0 || *L1_DoubleJet60er3p0 || *L1_DoubleJet80er3p0 || *L1_DoubleJet100er3p0 || *L1_DoubleJet112er3p0 || *L1_DoubleJet120er3p0 || *L1_TripleJet_88_72_56_VBF || *L1_TripleJet_84_68_48_VBF || *L1_TripleJet_92_76_64_VBF || *L1_QuadJet40er3p0 || *L1_QuadJet50er3p0 || *L1_QuadJet60er3p0 || *L1_HTT120er || *L1_HTT160er || *L1_HTT200er || *L1_HTT240er || *L1_HTT255er || *L1_HTT270er || *L1_HTT280er || *L1_HTT300er || *L1_HTT320er || *L1_HTT220er ) 
      {
	passL1seedExpression2016 = true;
      }


    if (not passL1seedExpression2016 || not **HLT_EB || not **HLT_EE) continue;

    for (Int_t i = 0; i < *STr2_NPi0_rec; i++) {

      // regional thresholds, to be set below
      double ptPairThrStream = 0.0;
      double ptGammaThrStream = 0.0;
      double s4s9ThrStream = 0.0;
      int nXtalThrStream = 0;
	// for region 2 ask always Nxtal >= 4, it should have been used, but actually I noticed that it wasn't in the ntuples for region2
      int nXtalThrStream2016 = 0;
      bool isRegion1 = false;

      double ptPairThrCalib = 0.0;
      double ptGammaThrCalib = 0.0;
      double s4s9ThrCalib = 0.0;
      int nXtalThrCalib = 0;

      double pi0mass = STr2_mPi0_nocor[i];
      double pi0Pt = STr2_ptPi0_nocor[i];
      double ptOverM = pi0Pt/pi0mass;

      double ptGam1 = STr2_enG1_nocor[i]/cosh(STr2_Eta_1[i]);
      double ptGam2 = STr2_enG2_nocor[i]/cosh(STr2_Eta_2[i]);
      double ptGammaMin = min(ptGam1, ptGam2);
      double s4s9min = min(STr2_S4S9_1[i],STr2_S4S9_2[i]);
      int nXtalMin = min(STr2_n1CrisPi0_rec[i],STr2_n2CrisPi0_rec[i]);
      bool passStreamSel2016 = false;
      bool passStreamSel2017 = false;
      bool passCalibSel2017 = false;
      bool passCalibSel2017HLTiso = false;
      
      if (STr2_Pi0recIsEB[i]) {

	if (pi0mass < massMinEB || pi0mass > massMaxEB) continue;

	if (fabs(STr2_etaPi0_rec[i]) < EBetaRegionBoundary) {
	  ptPairThrStream = 2.0;
	  ptGammaThrStream = 0.65;	  
	  s4s9ThrStream = 0.88;
	  nXtalThrStream = 6;
	  isRegion1 = true;
	  ptPairThrCalib = 2.6;
	  ptGammaThrCalib = 1.3;	  
	  s4s9ThrCalib = 0.83;
	  nXtalThrCalib = 7;
	  nXtalThrStream2016 = 0;
	  nRegion = 0;
	} else {
	  ptPairThrStream = 1.75;
	  ptGammaThrStream = 0.65;	  
	  s4s9ThrStream = 0.9;
	  nXtalThrStream = 6;
	  isRegion1 = false;
	  ptPairThrCalib = 2.6;
	  ptGammaThrCalib = 1.3;	  
	  s4s9ThrCalib = 0.83;
	  nXtalThrCalib = 7;
	  if (nXtalMin < 4) continue; // for region 2 in EB I forgot to set the minimum nXtal to 4 when producing the nuples, so I require it here 
	  nRegion = 1;
	}

	hmass_allEB_ntpSel->Fill(pi0mass);
	if (isRegion1) {
	  hmass_region1EB_ntpSel->Fill(pi0mass);
	} else {	  
	  hmass_region2EB_ntpSel->Fill(pi0mass);
	}

	//if (STr2_IsoPi0_rec[i] < 0.5) continue; // this is actually a DR distance between pi0 and any other cluster (excluding the two of the photons from pi0)
	// if (STr2_HLTIsoPi0_rec[i] > 0.5) continue; // this is the isolation defined also in the HLT

	if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrStream && s4s9min > s4s9ThrStream && pi0Pt > ptPairThrStream) passStreamSel2016 = true;
	if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrStream && s4s9min > s4s9ThrStream && pi0Pt > ptPairThrStream && nXtalMin >= nXtalThrStream) passStreamSel2017 = true;
	if (STr2_IsoPi0_rec[i] > 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017 = true;
	if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017HLTiso = true;
      
        if (passStreamSel2016) { 
	  hmass_allEB_2016Sel->Fill(pi0mass);
	  if (isRegion1) hmass_region1EB_2016Sel->Fill(pi0mass);
	  else hmass_region2EB_2016Sel->Fill(pi0mass);
	}
	
	if (passStreamSel2017) {
	  hmass_allEB_2017Sel->Fill(pi0mass);
	  if (isRegion1) hmass_region1EB_2017Sel->Fill(pi0mass);
	  else hmass_region2EB_2017Sel->Fill(pi0mass);
	}
	
	if (passCalibSel2017) {
	  hmass_allEB_2017SelCalib->Fill(pi0mass);
	  if (isRegion1) hmass_region1EB_2017SelCalib->Fill(pi0mass);
	  else hmass_region2EB_2017SelCalib->Fill(pi0mass);
	}

	if (passCalibSel2017HLTiso) {
	  hmass_allEB_2017SelCalibHLTiso->Fill(pi0mass);
	  if (isRegion1) hmass_region1EB_2017SelCalibHLTiso->Fill(pi0mass);
	  else hmass_region2EB_2017SelCalibHLTiso->Fill(pi0mass);
	}

	// base calib sel, no Pt cuts
	// since HLTiso < 0.5 already in the ntuples, just take or or STr2_IsoPi0_rec[i] > 0.5 and useHLTisoCalibForComparison
	// if latter is true, than it is like not cutting on cluster isolation
	// if it is false, then the cluster isolation selection decides whether to keep or reject the event
	if ( (useHLTisoCalibForComparison || STr2_IsoPi0_rec[i] > 0.5) && s4s9min > s4s9ThrCalib) {

	  if (nXtalMin >= nXtalThrCalib) {

	    for (UInt_t ipt = 0; ipt < ptGamCut.size(); ipt++) {
	      if (ptGammaMin > ptGamCut[ipt]) {
		if (isRegion1) hmass_region1EB_2017SelCalib_ptGam[ipt]->Fill(pi0mass);
		else hmass_region2EB_2017SelCalib_ptGam[ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairCut.size(); ipt++) {
	      if (pi0Pt > ptPairCut[ipt]) {
		if (isRegion1) hmass_region1EB_2017SelCalib_ptPair[ipt]->Fill(pi0mass);
		else hmass_region2EB_2017SelCalib_ptPair[ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairOverMCut.size(); ipt++) {
	      if ( ptOverM > ptPairOverMCut[ipt]) {
		if (isRegion1) hmass_region1EB_2017SelCalib_ptPairOverM[ipt]->Fill(pi0mass);
		else hmass_region2EB_2017SelCalib_ptPairOverM[ipt]->Fill(pi0mass);
	      }
	    }
	    
	  }

	  for (UInt_t ipt = 0; ipt < nXtalCut.size(); ipt++) {
	    if ( nXtalMin >= nXtalCut[ipt]) {
	      if (isRegion1) hmass_region1EB_2017SelCalib_nXtal[ipt]->Fill(pi0mass);
	      else hmass_region2EB_2017SelCalib_nXtal[ipt]->Fill(pi0mass);
	    }
	  }

	}

	if (nXtalMin >= nXtalThrCalib && s4s9min > s4s9ThrCalib) {

	  for (UInt_t ipt = 0; ipt < clusIsoCut.size(); ipt++) {
	    if (STr2_IsoPi0_rec[i] > clusIsoCut[ipt]) {
	      if (isRegion1) hmass_region1EB_2017SelCalib_clusIso[ipt]->Fill(pi0mass);
	      else hmass_region2EB_2017SelCalib_clusIso[ipt]->Fill(pi0mass);
	    }
	  }
	  
	}

      } else {
    
        if (pi0mass < massMinEE || pi0mass > massMaxEE) continue;
	
	// use to distinguish region 3, which is only for EE
	Bool_t isRegion2 = false;

	// there should be 3 regions, but region3 has same thresholds as region2 
	if (fabs(STr2_etaPi0_rec[i]) < EEetaRegionBoundary) {
	  ptPairThrStream = 3.75;
	  ptGammaThrStream = 1.1;	  
	  s4s9ThrStream = 0.85;
	  nXtalThrStream = 6;
	  isRegion1 = true;
	  ptPairThrCalib = 3.75;
	  ptGammaThrCalib = 1.1;	  
	  s4s9ThrCalib = 0.95;
	  nXtalThrCalib = 7;
	  nRegion = 2;
	} else {
	  ptPairThrStream = 2.0;
	  ptGammaThrStream = 0.95;	  
	  s4s9ThrStream = 0.92;
	  nXtalThrStream = 6;
	  isRegion1 = false;
	  ptPairThrCalib = 2.0;
	  ptGammaThrCalib = 0.95;	  
	  s4s9ThrCalib = 0.95;
	  nXtalThrCalib = 7;
	  nRegion = 3;
	  if (fabs(STr2_etaPi0_rec[i]) < 2.0) isRegion2 = true;
	  else nRegion = 4;
	}

	hmass_allEE_ntpSel->Fill(pi0mass);
	if (isRegion1) {
	  hmass_region1EE_ntpSel->Fill(pi0mass);
	} else if (isRegion2) {
	  hmass_region2EE_ntpSel->Fill(pi0mass);
	} else {
	  hmass_region3EE_ntpSel->Fill(pi0mass);
	}

	if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrStream && s4s9min > s4s9ThrStream && pi0Pt > ptPairThrStream) passStreamSel2016 = true;
	if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrStream && s4s9min > s4s9ThrStream && pi0Pt > ptPairThrStream && nXtalMin >= nXtalThrStream) passStreamSel2017 = true;
	if (STr2_IsoPi0_rec[i] > 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017 = true;	
        if (STr2_HLTIsoPi0_rec[i] < 0.5 && ptGammaMin > ptGammaThrCalib && s4s9min > s4s9ThrCalib && pi0Pt > ptPairThrCalib && nXtalMin >= nXtalThrCalib) passCalibSel2017HLTiso = true;
          
        if (passStreamSel2016) { 
	  hmass_allEE_2016Sel->Fill(pi0mass);
	  if (isRegion1) hmass_region1EE_2016Sel->Fill(pi0mass);
	  else if (isRegion2) hmass_region2EE_2016Sel->Fill(pi0mass);
	  else hmass_region3EE_2016Sel->Fill(pi0mass);
	}
	
	if (passStreamSel2017) {
	  hmass_allEE_2017Sel->Fill(pi0mass);
	  if (isRegion1) hmass_region1EE_2017Sel->Fill(pi0mass);
	  else if (isRegion2) hmass_region2EE_2017Sel->Fill(pi0mass);
	  else hmass_region3EE_2017Sel->Fill(pi0mass);
	}
	
	if (passCalibSel2017) {
	  hmass_allEE_2017SelCalib->Fill(pi0mass);
	  if (isRegion1) hmass_region1EE_2017SelCalib->Fill(pi0mass);
	  else if (isRegion2) hmass_region2EE_2017SelCalib->Fill(pi0mass);
	  else hmass_region3EE_2017SelCalib->Fill(pi0mass);
	}

	if (passCalibSel2017HLTiso) {
	  hmass_allEE_2017SelCalibHLTiso->Fill(pi0mass);
	  if (isRegion1) hmass_region1EE_2017SelCalibHLTiso->Fill(pi0mass);
	  else if (isRegion2) hmass_region2EE_2017SelCalibHLTiso->Fill(pi0mass);
	  else hmass_region3EE_2017SelCalibHLTiso->Fill(pi0mass);
	}

	// base calib sel, no Pt cuts
	// see comments in the EB equivalent part
	if ((useHLTisoCalibForComparison || STr2_IsoPi0_rec[i] > 0.5) && s4s9min > s4s9ThrCalib) {

	  if (nXtalMin >= nXtalThrCalib) {

	    for (UInt_t ipt = 0; ipt < ptGamCut.size(); ipt++) {
	      if (ptGammaMin > ptGamCut[ipt]) {
		if (isRegion1) hmass_region1EE_2017SelCalib_ptGam[ipt]->Fill(pi0mass);
		else if (isRegion2) hmass_region2EE_2017SelCalib_ptGam[ipt]->Fill(pi0mass);
		else hmass_region3EE_2017SelCalib_ptGam[ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairCut.size(); ipt++) {
	      if (pi0Pt > ptPairCut[ipt]) {
		if (isRegion1) hmass_region1EE_2017SelCalib_ptPair[ipt]->Fill(pi0mass);
		else if (isRegion2) hmass_region2EE_2017SelCalib_ptPair[ipt]->Fill(pi0mass);
		else hmass_region3EE_2017SelCalib_ptPair[ipt]->Fill(pi0mass);
	      }
	    }
	    for (UInt_t ipt = 0; ipt < ptPairOverMCut.size(); ipt++) {
	      if ( ptOverM > ptPairOverMCut[ipt]) {
		if (isRegion1) hmass_region1EE_2017SelCalib_ptPairOverM[ipt]->Fill(pi0mass);
		else if (isRegion2) hmass_region2EE_2017SelCalib_ptPairOverM[ipt]->Fill(pi0mass);
		else hmass_region3EE_2017SelCalib_ptPairOverM[ipt]->Fill(pi0mass);
	      }
	    }    

	  }
	  
	  for (UInt_t ipt = 0; ipt < nXtalCut.size(); ipt++) {
	    if ( nXtalMin > nXtalCut[ipt]) {
	      if (isRegion1) hmass_region1EE_2017SelCalib_nXtal[ipt]->Fill(pi0mass);
	      else if (isRegion2) hmass_region2EE_2017SelCalib_nXtal[ipt]->Fill(pi0mass);
	      else hmass_region3EE_2017SelCalib_nXtal[ipt]->Fill(pi0mass);
	    }
	  }

	}

	if (nXtalMin >= nXtalThrCalib && s4s9min > s4s9ThrCalib) {

	  for (UInt_t ipt = 0; ipt < clusIsoCut.size(); ipt++) {
	    if (STr2_IsoPi0_rec[i] > clusIsoCut[ipt]) {
	      if (isRegion1) hmass_region1EE_2017SelCalib_clusIso[ipt]->Fill(pi0mass);
	      else if (isRegion2) hmass_region2EE_2017SelCalib_clusIso[ipt]->Fill(pi0mass);
	      else hmass_region3EE_2017SelCalib_clusIso[ipt]->Fill(pi0mass);
	    }
	  }

	}

      }
  
      vm[nRegion]
    
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


  draw_nTH1(vm_allEB->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_allEB", outDir, vm_allEB->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region1EB->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB", outDir, vm_region1EB->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB", outDir, vm_region2EB->getHistLegends(), "", lumi, 1, false, false, rebinFactorEB);

  draw_nTH1(vm_region1EB_ptGam->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB_ptGam", outDir, vm_region1EB_ptGam->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB_ptGam->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB_ptGam", outDir, vm_region2EB_ptGam->getHistLegends(), "", lumi, 2, false, false);
  draw_nTH1(vm_region1EB_ptPair->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB_ptPair", outDir, vm_region1EB_ptPair->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB_ptPair->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB_ptPair", outDir, vm_region2EB_ptPair->getHistLegends(), "", lumi, 2, false, false);
  draw_nTH1(vm_region1EB_ptPairOverM->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB_ptPairOverM", outDir, vm_region1EB_ptPairOverM->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB_ptPairOverM->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB_ptPairOverM", outDir, vm_region2EB_ptPairOverM->getHistLegends(), "", lumi, 2, false, false);
  draw_nTH1(vm_region1EB_nXtal->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB_nXtal", outDir, vm_region1EB_nXtal->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB_nXtal->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB_nXtal", outDir, vm_region2EB_nXtal->getHistLegends(), "", lumi, 2, false, false);
  draw_nTH1(vm_region1EB_clusIso->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region1EB_clusIso", outDir, vm_region1EB_clusIso->getHistLegends(), "", lumi, 1, false, false);
  draw_nTH1(vm_region2EB_clusIso->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", "pi0Mass_region2EB_clusIso", outDir, vm_region2EB_clusIso->getHistLegends(), "", lumi, 1, false, false);



  draw_nTH1(vm_allEE->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_allEE", outDir, vm_allEE->getHistLegends(), "", lumi, 1, false, false, rebinFactorEE);
  draw_nTH1(vm_region1EE->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE", outDir, vm_region1EE->getHistLegends(), "", lumi, 1, false, false, rebinFactorEE);
  draw_nTH1(vm_region2EE->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE", outDir, vm_region2EE->getHistLegends(), "", lumi, 1, false, false, rebinFactorEE);
  draw_nTH1(vm_region3EE->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE", outDir, vm_region3EE->getHistLegends(), "", lumi, 1, false, false, rebinFactorEE);

  draw_nTH1(vm_region1EE_ptGam->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE_ptGam", outDir, vm_region1EE_ptGam->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region2EE_ptGam->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE_ptGam", outDir, vm_region2EE_ptGam->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region3EE_ptGam->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE_ptGam", outDir, vm_region3EE_ptGam->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region1EE_ptPair->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE_ptPair", outDir, vm_region1EE_ptPair->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region2EE_ptPair->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE_ptPair", outDir, vm_region2EE_ptPair->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region3EE_ptPair->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE_ptPair", outDir, vm_region3EE_ptPair->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region1EE_ptPairOverM->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE_ptPairOverM", outDir, vm_region1EE_ptPairOverM->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region2EE_ptPairOverM->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE_ptPairOverM", outDir, vm_region2EE_ptPairOverM->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region3EE_ptPairOverM->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE_ptPairOverM", outDir, vm_region3EE_ptPairOverM->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region1EE_nXtal->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE_nXtal", outDir, vm_region1EE_nXtal->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region2EE_nXtal->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE_nXtal", outDir, vm_region2EE_nXtal->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region3EE_nXtal->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE_nXtal", outDir, vm_region3EE_nXtal->getHistLegends(), "", lumi, 8, false, false);
  draw_nTH1(vm_region1EE_clusIso->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region1EE_clusIso", outDir, vm_region1EE_clusIso->getHistLegends(), "", lumi, 4, false, false);
  draw_nTH1(vm_region2EE_clusIso->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region2EE_clusIso", outDir, vm_region2EE_clusIso->getHistLegends(), "", lumi, 4, false, false);
  draw_nTH1(vm_region3EE_clusIso->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", "pi0Mass_region3EE_clusIso", outDir, vm_region3EE_clusIso->getHistLegends(), "", lumi, 4, false, false);

  for (UInt_t i = 0; i < regionTagName.size(); i++) {

    draw_nTH1(vm[i]->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})","a.u.", "pi0Mass_"+regionTagName[i], outDir, vm[i]->getHistLegends(), "", lumi, 1, false, false, vm[i]->getHistRebins());

  }

  cout << "THE END!" << endl;

}
