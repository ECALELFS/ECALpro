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

//=================================================

class vectorManager {

public:
  vectorManager() { };

  vectorManager(const vector<TH1*> & histPtrs,
		const vector<string> & histNames,
		const vector<string> & histLegends
		)
  {
    histPtrs_    = vector<TH1*  >(histPtrs);
    histNames_   = vector<string>(histNames);
    histLegends_ = vector<string>(histLegends);
  };

  ~vectorManager() {};

  vector<TH1*>   getHistPtrs()    const { return histPtrs_;    };
  vector<string> getHistNames()   const { return histNames_;   };
  vector<string> getHistLegends() const { return histLegends_; };

  void addComponent(TH1* histPtr = NULL, const string& histName = "name", const string& histLegend = "leg")  { 
    histPtrs_.push_back(histPtr);
    histNames_.push_back(histName);
    histLegends_.push_back(histLegend);
  };

private:

  vector<TH1*> histPtrs_;
  vector<string> histNames_;
  vector<string> histLegends_;

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


void makeL1seedStudy(const bool isEB = true,
		     const bool isPi0 = true,
		     const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/Optimization/streamSelection/",
		     const string& eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/",
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


  double lumi = 7.5; // to be updated

  TH1F* L1seedComposition = new TH1F();
  L1seedComposition->SetNameTitle("L1seedComposition","");
  Bool_t isFirstFile = true;

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

	///////////////////
	// open root file and get triggerComposition histogram
	TFile *rootfile = new TFile((eosPathComplete+line).c_str(),"READ");
	if (!rootfile || !rootfile->IsOpen()) {
	  cout << "Error: file \"" << eosPathComplete+line << "\" was not opened." << endl;
	  exit(EXIT_FAILURE);
	}
	TH1F* htmp = (TH1F*) rootfile->Get("triggerComposition");
	if (!htmp || htmp == NULL) {
	  cout << "Error: histogram 'triggerComposition' not found in file ' " << eosPathComplete+line << "'. End of programme." << endl;
	  exit(EXIT_FAILURE);
	}

	// if first histogram, set histogram
	if (isFirstFile) {

	  L1seedComposition->SetBins(htmp->GetNbinsX(),htmp->GetBinLowEdge(1),htmp->GetBinLowEdge(htmp->GetNbinsX()+1));
	  for (Int_t bin = 1; bin <= htmp->GetNbinsX(); bin++) {
	    L1seedComposition->GetXaxis()->SetBinLabel(bin,htmp->GetXaxis()->GetBinLabel(bin));
	    L1seedComposition->SetBinContent(bin,0.0);
	  }

	} 
	// add current histogram
	if (not L1seedComposition->Add((TH1F*) htmp->Clone())) {
	  cout << "Error: L1seedComposition->Add(htmp) failed. End of programme." << endl;
	  exit(EXIT_FAILURE);
	}
	
	isFirstFile = false;
	rootfile->Close();
	delete rootfile;
	///  end of operations on root file
	//////////////////////////////

      }

    } else {
      cout << "Error opening file " << goodNtuplesFileName << endl;
      exit(EXIT_FAILURE);
    }
    infile.close();

  }
  
  if (!L1seedComposition || L1seedComposition == NULL) {
    cout << "Error: histogram 'L1seedComposition' not filled properly. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  long int nTotal = chain->GetEntries();  
  long int nPassL1 = chain->GetEntries();
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

  long int nPassL1seedExpression2016 = 0;
  long int nPassL1seed_SingleEGX = 0;
  long int nPassL1seed_SingleIsoEGX = 0;
  long int nPassL1seed_SingleIsoEGXer2p1 = 0;
  long int nPassL1seed_DoubleEGX = 0;
  long int nPassL1seed_SingleJetX = 0;
  long int nPassL1seed_DoubleJetXer3p0 = 0;
  long int nPassL1seed_TripleJetX = 0;
  long int nPassL1seed_QuadJetXer3p0 = 0;
  long int nPassL1seed_HTTXer = 0;


  // while (reader.Next()) {

  //   //reader.SetLocalEntry(nEvents);
    
  //   cout.flush();
  //   if(nEvents % CHECK_EVERY_N == 0) {
  //     cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
  //     //cout << "entry : " << nEvents << endl;     
  //   }
  //   nEvents++;

  //   Bool_t passL1seedExpression2016 = false;
  //   Bool_t passL1seed_SingleEGX = false;
  //   Bool_t passL1seed_SingleIsoEGX = false;
  //   Bool_t passL1seed_SingleIsoEGXer2p1 = false;
  //   Bool_t passL1seed_DoubleEGX = false;
  //   Bool_t passL1seed_SingleJetX = false;
  //   Bool_t passL1seed_DoubleJetXer3p0 = false;
  //   Bool_t passL1seed_TripleJetX = false;
  //   Bool_t passL1seed_QuadJetXer3p0 = false;
  //   Bool_t passL1seed_HTTXer = false;

  //   vector<Bool_t> L1seedCategories;

  //   if ( *L1_SingleEG5 || *L1_SingleEG10 || *L1_SingleEG15 || *L1_SingleEG18 || *L1_SingleEG24 || *L1_SingleEG26 || *L1_SingleEG28 || *L1_SingleEG30 || *L1_SingleEG32 || *L1_SingleEG34 || *L1_SingleEG36 || *L1_SingleEG38 || *L1_SingleEG40 || *L1_SingleEG45 ) 
  //     {
  // 	passL1seed_SingleEGX = true;
  // 	nPassL1seed_SingleEGX++;
  //     }

  //   if ( *L1_SingleIsoEG18 || *L1_SingleIsoEG20 || *L1_SingleIsoEG22 || *L1_SingleIsoEG24 || *L1_SingleIsoEG26 || *L1_SingleIsoEG28 || *L1_SingleIsoEG30 || *L1_SingleIsoEG32 || *L1_SingleIsoEG34 || *L1_SingleIsoEG36 ) 
  //     {
  // 	passL1seed_SingleIsoEGX = true;
  // 	nPassL1seed_SingleIsoEGX++;
  //     }

  //   if ( *L1_SingleIsoEG18er2p1 || *L1_SingleIsoEG20er2p1 || *L1_SingleIsoEG22er2p1 || *L1_SingleIsoEG24er2p1 || *L1_SingleIsoEG26er2p1 || *L1_SingleIsoEG28er2p1 || *L1_SingleIsoEG30er2p1 || *L1_SingleIsoEG32er2p1 || *L1_SingleIsoEG34er2p1 ) 
  //     {
  // 	passL1seed_SingleIsoEGXer2p1 = true;
  // 	nPassL1seed_SingleIsoEGXer2p1++;
  //     }

  //   if ( *L1_DoubleEG_15_10 || *L1_DoubleEG_18_17 || *L1_DoubleEG_20_18 || *L1_DoubleEG_22_10 || *L1_DoubleEG_23_10 || *L1_DoubleEG_22_12 || *L1_DoubleEG_22_15 || *L1_DoubleEG_24_17 || *L1_DoubleEG_25_12 ) 
  //     {
  // 	passL1seed_DoubleEGX = true;
  // 	nPassL1seed_DoubleEGX++;
  //     }

  //   if ( *L1_SingleJet16 || *L1_SingleJet20 || *L1_SingleJet35 || *L1_SingleJet60 || *L1_SingleJet90 || *L1_SingleJet120 || *L1_SingleJet140 || *L1_SingleJet150 || *L1_SingleJet160 || *L1_SingleJet170 || *L1_SingleJet180 || *L1_SingleJet200 ) 
  //     {
  // 	passL1seed_SingleJetX = true;
  // 	nPassL1seed_SingleJetX++;
  //     }

  //   if ( *L1_DoubleJet40er3p0 || *L1_DoubleJet50er3p0 || *L1_DoubleJet60er3p0 || *L1_DoubleJet80er3p0 || *L1_DoubleJet100er3p0 || *L1_DoubleJet112er3p0 || *L1_DoubleJet120er3p0 ) 
  //     {
  // 	passL1seed_DoubleJetXer3p0 = true;
  // 	nPassL1seed_DoubleJetXer3p0++;
  //     }

  //   if ( *L1_TripleJet_88_72_56_VBF || *L1_TripleJet_84_68_48_VBF || *L1_TripleJet_92_76_64_VBF ) 
  //     {
  // 	passL1seed_TripleJetX = true;
  // 	nPassL1seed_TripleJetX++;
  //     }

  //   if ( *L1_QuadJet40er3p0 || *L1_QuadJet50er3p0 || *L1_QuadJet60er3p0 ) 
  //     {
  // 	passL1seed_QuadJetXer3p0 = true;
  // 	nPassL1seed_QuadJetXer3p0++;
  //     }

  //   if ( *L1_HTT120er || *L1_HTT160er || *L1_HTT200er || *L1_HTT240er || *L1_HTT255er || *L1_HTT270er || *L1_HTT280er || *L1_HTT300er || *L1_HTT320er || *L1_HTT220er ) 
  //     {
  // 	passL1seed_HTTXer = true;
  // 	nPassL1seed_HTTXer++;
  //     }

  //   if ( *L1_AlwaysTrue || *L1_IsolatedBunch || passL1seed_SingleEGX || passL1seed_SingleIsoEGX || passL1seed_SingleIsoEGXer2p1 || passL1seed_DoubleEGX || passL1seed_SingleJetX || passL1seed_DoubleJetXer3p0 || passL1seed_TripleJetX || passL1seed_QuadJetXer3p0 || passL1seed_HTTXer ) 
  //     {
  // 	passL1seedExpression2016 = true;
  // 	nPassL1seedExpression2016++;
  //     }


  //   if (not passL1seedExpression2016) continue;
  //   if (not **HLT_EB && not **HLT_EE) continue;



  // }

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

  TCanvas* canvas = new TCanvas("canvas","",900,600);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetGridx(1);
  canvas->cd();
  canvas->SetRightMargin(0.02);
  canvas->SetLeftMargin(0.05);
  canvas->SetBottomMargin(0.3);
  canvas->SetLogy();
  L1seedComposition->SetStats(0);
  L1seedComposition->SetTitle("L1 seeds composition");

  double minNot0 = 9999.9;
  for (Int_t ibin = 1; ibin <= L1seedComposition->GetNbinsX(); ibin++ ) {
    if (L1seedComposition->GetBinContent(ibin) > 0.0000001 && minNot0 > L1seedComposition->GetBinContent(ibin)) minNot0 = L1seedComposition->GetBinContent(ibin);
  }
  L1seedComposition->GetYaxis()->SetRangeUser( 0.5 * minNot0,
					       5.0 * L1seedComposition->GetBinContent(L1seedComposition->GetMaximumBin())
					       );
  L1seedComposition->Draw("Hist");
  string canvasTag = isPi0 ? "pi0" : "eta";
  canvas->SaveAs((outDir+"L1seedComposition_"+ canvasTag + ".png").c_str());
  canvas->SaveAs((outDir+"L1seedComposition_"+ canvasTag + ".pdf").c_str());

  cout << "THE END!" << endl;

}
