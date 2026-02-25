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

static double EBetaRegionBoundary = 1.0; // warning, 1.0 not exactly equal to gap between module 3 and 4                    
static double EEetaRegionBoundary = 1.8; // 1.8 is boundary between region 1 and 2 in EE, region 3 basically equal to region 2      

void makeStreamCorrelation(const bool isPi0 = true,
			   const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/Optimization/correlation/looseSel/",
			   const string& eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/",
			   const TString& dirList = "AlCaP0_fromHLTPhysics_2017AB_TreeOptim,AlCaP0_from_ParkingHLTPhysics_TreeOptim"
		     ) 
{


  // can also read ntuples from pccmsrm29 (I copied only the good ntuples already)
  // in this case use eosPath = "/u2/mciprian/streamOptimNtuples/"

  createPlotDirAndCopyPhp(outDir);

  TObjArray* array = dirList.Tokenize(",");
  vector<string> eosDir;
  for (Int_t j = 0; j < array->GetEntries(); j++) {
    TString str = ((TObjString *) array->At(j))->String();
    eosDir.push_back(str.Data());
    cout << j << " --> " << eosDir[j] << endl;
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
  // TTreeReaderValue<Short_t> L1_SingleEG5 (reader,"L1_SingleEG5");
  // TTreeReaderValue<Short_t> L1_SingleEG10 (reader,"L1_SingleEG10");
  // TTreeReaderValue<Short_t> L1_SingleEG15 (reader,"L1_SingleEG15");
  // TTreeReaderValue<Short_t> L1_SingleEG18 (reader,"L1_SingleEG18");
  // TTreeReaderValue<Short_t> L1_SingleEG24 (reader,"L1_SingleEG24");
  // TTreeReaderValue<Short_t> L1_SingleEG26 (reader,"L1_SingleEG26");
  // TTreeReaderValue<Short_t> L1_SingleEG28 (reader,"L1_SingleEG28");
  // TTreeReaderValue<Short_t> L1_SingleEG30 (reader,"L1_SingleEG30");
  // TTreeReaderValue<Short_t> L1_SingleEG32 (reader,"L1_SingleEG32");
  // TTreeReaderValue<Short_t> L1_SingleEG34 (reader,"L1_SingleEG34");
  // TTreeReaderValue<Short_t> L1_SingleEG36 (reader,"L1_SingleEG36");
  // TTreeReaderValue<Short_t> L1_SingleEG38 (reader,"L1_SingleEG38");
  // TTreeReaderValue<Short_t> L1_SingleEG40 (reader,"L1_SingleEG40");
  // TTreeReaderValue<Short_t> L1_SingleEG45 (reader,"L1_SingleEG45");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG18 (reader,"L1_SingleIsoEG18");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG20 (reader,"L1_SingleIsoEG20");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG22 (reader,"L1_SingleIsoEG22");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG24 (reader,"L1_SingleIsoEG24");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG26 (reader,"L1_SingleIsoEG26");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG28 (reader,"L1_SingleIsoEG28");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG30 (reader,"L1_SingleIsoEG30");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG32 (reader,"L1_SingleIsoEG32");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG34 (reader,"L1_SingleIsoEG34");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG36 (reader,"L1_SingleIsoEG36");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG18er2p1 (reader,"L1_SingleIsoEG18er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG20er2p1 (reader,"L1_SingleIsoEG20er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG22er2p1 (reader,"L1_SingleIsoEG22er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG24er2p1 (reader,"L1_SingleIsoEG24er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG26er2p1 (reader,"L1_SingleIsoEG26er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG28er2p1 (reader,"L1_SingleIsoEG28er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG30er2p1 (reader,"L1_SingleIsoEG30er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG32er2p1 (reader,"L1_SingleIsoEG32er2p1");
  // TTreeReaderValue<Short_t> L1_SingleIsoEG34er2p1 (reader,"L1_SingleIsoEG34er2p1");
  // TTreeReaderValue<Short_t> L1_DoubleEG_15_10 (reader,"L1_DoubleEG_15_10");
  // TTreeReaderValue<Short_t> L1_DoubleEG_18_17 (reader,"L1_DoubleEG_18_17");
  // TTreeReaderValue<Short_t> L1_DoubleEG_20_18 (reader,"L1_DoubleEG_20_18");
  // TTreeReaderValue<Short_t> L1_DoubleEG_22_10 (reader,"L1_DoubleEG_22_10");
  // TTreeReaderValue<Short_t> L1_DoubleEG_22_12 (reader,"L1_DoubleEG_22_12");
  // TTreeReaderValue<Short_t> L1_DoubleEG_22_15 (reader,"L1_DoubleEG_22_15");
  // TTreeReaderValue<Short_t> L1_DoubleEG_23_10 (reader,"L1_DoubleEG_23_10");
  // TTreeReaderValue<Short_t> L1_DoubleEG_24_17 (reader,"L1_DoubleEG_24_17");
  // TTreeReaderValue<Short_t> L1_DoubleEG_25_12 (reader,"L1_DoubleEG_25_12");
  // TTreeReaderValue<Short_t> L1_SingleJet16 (reader,"L1_SingleJet16");
  // TTreeReaderValue<Short_t> L1_SingleJet20 (reader,"L1_SingleJet20");
  // TTreeReaderValue<Short_t> L1_SingleJet35 (reader,"L1_SingleJet35");
  // TTreeReaderValue<Short_t> L1_SingleJet60 (reader,"L1_SingleJet60");
  // TTreeReaderValue<Short_t> L1_SingleJet90 (reader,"L1_SingleJet90");
  // TTreeReaderValue<Short_t> L1_SingleJet120 (reader,"L1_SingleJet120");
  // TTreeReaderValue<Short_t> L1_SingleJet140 (reader,"L1_SingleJet140");
  // TTreeReaderValue<Short_t> L1_SingleJet150 (reader,"L1_SingleJet150");
  // TTreeReaderValue<Short_t> L1_SingleJet160 (reader,"L1_SingleJet160");
  // TTreeReaderValue<Short_t> L1_SingleJet170 (reader,"L1_SingleJet170");
  // TTreeReaderValue<Short_t> L1_SingleJet180 (reader,"L1_SingleJet180");
  // TTreeReaderValue<Short_t> L1_SingleJet200 (reader,"L1_SingleJet200");
  // TTreeReaderValue<Short_t> L1_DoubleJet40er3p0 (reader,"L1_DoubleJet40er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet50er3p0 (reader,"L1_DoubleJet50er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet60er3p0 (reader,"L1_DoubleJet60er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet80er3p0 (reader,"L1_DoubleJet80er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet100er3p0 (reader,"L1_DoubleJet100er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet112er3p0 (reader,"L1_DoubleJet112er3p0");
  // TTreeReaderValue<Short_t> L1_DoubleJet120er3p0 (reader,"L1_DoubleJet120er3p0");
  // TTreeReaderValue<Short_t> L1_TripleJet_84_68_48_VBF (reader,"L1_TripleJet_84_68_48_VBF");
  // TTreeReaderValue<Short_t> L1_TripleJet_88_72_56_VBF (reader,"L1_TripleJet_88_72_56_VBF");
  // TTreeReaderValue<Short_t> L1_TripleJet_92_76_64_VBF (reader,"L1_TripleJet_92_76_64_VBF");
  // TTreeReaderValue<Short_t> L1_QuadJet40er3p0 (reader,"L1_QuadJet40er3p0");
  // TTreeReaderValue<Short_t> L1_QuadJet50er3p0 (reader,"L1_QuadJet50er3p0");
  // TTreeReaderValue<Short_t> L1_QuadJet60er3p0 (reader,"L1_QuadJet60er3p0");
  // TTreeReaderValue<Short_t> L1_HTT120er (reader,"L1_HTT120er");
  // TTreeReaderValue<Short_t> L1_HTT160er (reader,"L1_HTT160er");
  // TTreeReaderValue<Short_t> L1_HTT200er (reader,"L1_HTT200er");
  // TTreeReaderValue<Short_t> L1_HTT220er (reader,"L1_HTT220er");
  // TTreeReaderValue<Short_t> L1_HTT240er (reader,"L1_HTT240er");
  // TTreeReaderValue<Short_t> L1_HTT255er (reader,"L1_HTT255er");
  // TTreeReaderValue<Short_t> L1_HTT270er (reader,"L1_HTT270er");
  // TTreeReaderValue<Short_t> L1_HTT280er (reader,"L1_HTT280er");
  // TTreeReaderValue<Short_t> L1_HTT300er (reader,"L1_HTT300er");
  // TTreeReaderValue<Short_t> L1_HTT320er (reader,"L1_HTT320er");
  // TTreeReaderValue<Short_t> L1_IsolatedBunch (reader,"L1_IsolatedBunch");
  // TTreeReaderValue<Short_t> L1_AlwaysTrue (reader,"L1_AlwaysTrue");

  // note: quantities in the ntuples are ordered by energy of photons (actually energy of the seed)
  // eta, phi and other variables are associate to photon before correcting the energy
  // the actual order in energy of two photons can be swapped if a correction on one is high (or if the seed ordering is not equal to the cluster energy one)
  // therefore it is possible that ptGam1 < ptGam2, and therefore we must correct the order if ptGam must be the leading pT photon
  
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

  ///////////////////////
  // histograms
  // correlation (we often use the minimum value of the two photons variable, since it is used to cut in the stream)
  vector<TH2F*> h2_ptGam_M;
  vector<TH2F*> h2_ptPair_M;
  vector<TH2F*> h2_ptPairOverM_M;
  vector<TH2F*> h2_s4s9_M;
  vector<TH2F*> h2_nXtal_M;
  vector<TH2F*> h2_nXtal1_nXtal2;
  vector<TH2F*> h2_ptGam1_ptGam2;
  vector<TH2F*> h2_ptGam_ptPair;
  vector<TH2F*> h2_ptGam_ptPairOverM;
  vector<TH2F*> h2_ptPair_ptPairOverM;
  vector<TH2F*> h2_nXtal_ptGam;
  vector<TH2F*> h2_nXtal_ptPair;
  vector<TH2F*> h2_nXtal_DRgam;
  vector<TH2F*> h2_nXtal_s4s9;
  vector<TH2F*> h2_seedDistance_ptPair;
  vector<TH2F*> h2_seedDistance_nXtal;
  vector<TH2F*> h2_seedDistance_ptGam;
  vector<TH2F*> h2_HLTiso_DRiso;
  vector<TH2F*> h2_DRiso_M;
  vector<TH2F*> h2_HLTiso_M;

  float massMinEB = isPi0 ? 0.06 : 0.2;
  float massMaxEB = isPi0 ? 0.22 : 0.8;
  float massMinEE = isPi0 ? 0.05 : 0.2;
  float massMaxEE = isPi0 ? 0.25 : 0.8;

  vector<double> etaBinEdgesEB = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5};
  vector<TH2F*> h2_nXtal_M_EB_eta;
  for (UInt_t i = 0; i < etaBinEdgesEB.size()-1; i++) {
    string str1 = getStringFromDouble(etaBinEdgesEB[i]);
    string str2 = getStringFromDouble(etaBinEdgesEB[i+1]);
    h2_nXtal_M_EB_eta.push_back(new TH2F(Form("h2_nXtal_M_EB_eta%sto%s",str1.c_str(),str2.c_str()),Form("nXtal_M_EB_eta%sto%s",str1.c_str(),str2.c_str()),6,3.5,9.5,80,massMinEB,massMaxEB));
  }

  for (UInt_t i = 0; i < regionTagName.size(); i++) {
    string name = regionTagName[i];
    double massMin = (i < 2) ? massMinEB : massMinEE;
    double massMax = (i < 2) ? massMaxEB : massMaxEE;
    h2_ptGam_M.push_back(new TH2F(Form("h2_ptGam_M_%s",name.c_str()),Form("ptGam_M_%s",name.c_str()),25,0.5,5.5,80,massMin,massMax));
    h2_ptPair_M.push_back(new TH2F(Form("h2_ptPair_M_%s",name.c_str()),Form("ptPair_M_%s",name.c_str()),50,0.5,10.5,80,massMin,massMax));
    h2_ptPairOverM_M.push_back(new TH2F(Form("h2_ptPairOverM_M_%s",name.c_str()),Form("ptPairOverM_M_%s",name.c_str()),24,0.0,48,80,massMin,massMax));
    h2_s4s9_M.push_back(new TH2F(Form("h2_s4s9_M_%s",name.c_str()),Form("s4s9_M_%s",name.c_str()),25,0.75,1.0,80,massMin,massMax));
    h2_nXtal_M.push_back(new TH2F(Form("h2_nXtal_M_%s",name.c_str()),Form("nXtal_M_%s",name.c_str()),6,3.5,9.5,80,massMin,massMax));
    h2_nXtal1_nXtal2.push_back(new TH2F(Form("h2_nXtal1_nXtal2_%s",name.c_str()),Form("nXtal1_nXtal2_%s",name.c_str()),6,3.5,9.5,6,3.5,9.5));
    h2_ptGam1_ptGam2.push_back(new TH2F(Form("h2_ptGam1_ptGam2_%s",name.c_str()),Form("ptGam1_ptGam2_%s",name.c_str()),50,0.5,10.5,25,0.5,5.5));
    h2_ptGam_ptPair.push_back(new TH2F(Form("h2_ptGam_ptPair_%s",name.c_str()),Form("ptGam_ptPair_%s",name.c_str()),25,0.5,5.5,50,0.5,10.5));
    h2_ptGam_ptPairOverM.push_back(new TH2F(Form("h2_ptGam_ptPairOverM_%s",name.c_str()),Form("ptGam_ptPairOverM_%s",name.c_str()),25,0.5,5.5,24,0.0,48));
    h2_ptPair_ptPairOverM.push_back(new TH2F(Form("h2_ptPair_ptPairOverM_%s",name.c_str()),Form("ptPair_ptPairOverM_%s",name.c_str()),50,0.5,10.5,24,0.0,48));
    h2_nXtal_ptGam.push_back(new TH2F(Form("h2_nXtal_ptGam_%s",name.c_str()),Form("nXtal_ptGam_%s",name.c_str()),6,3.5,9.5,25,0.5,5.5));
    h2_nXtal_ptPair.push_back(new TH2F(Form("h2_nXtal_ptPair_%s",name.c_str()),Form("nXtal_ptPair_%s",name.c_str()),6,3.5,9.5,50,0.5,10.5));
    h2_nXtal_DRgam.push_back(new TH2F(Form("h2_nXtal_DRgam_%s",name.c_str()),Form("nXtal_DRgam_%s",name.c_str()),6,3.5,9.5,25,0.0,0.5));
    h2_nXtal_s4s9.push_back(new TH2F(Form("h2_nXtal_s4s9_%s",name.c_str()),Form("nXtal_s4s9_%s",name.c_str()),6,3.5,9.5,25,0.75,1.0));
    h2_seedDistance_ptPair.push_back(new TH2F(Form("h2_seedDistance_ptPair_%s",name.c_str()),Form("seedDistance_ptPair_%s",name.c_str()),20,0.5,20.5,50,0.5,10.5));
    h2_seedDistance_nXtal.push_back(new TH2F(Form("h2_seedDistance_nXtal_%s",name.c_str()),Form("seedDistance_nXtal_%s",name.c_str()),20,0.5,20.5,6,3.5,9.5));
    h2_seedDistance_ptGam.push_back(new TH2F(Form("h2_seedDistance_ptGam_%s",name.c_str()),Form("seedDistance_ptGam_%s",name.c_str()),20,0.5,20.5,25,0.5,5.5));
    h2_HLTiso_DRiso.push_back(new TH2F(Form("h2_HLTiso_DRiso_%s",name.c_str()),Form("HLTiso_DRiso_%s",name.c_str()),25,0.0,0.5,20,0.0,2.0));
    h2_DRiso_M.push_back(new TH2F(Form("h2_DRiso_M_%s",name.c_str()),Form("DRiso_M_%s",name.c_str()),20,0.0,2.0,80,massMin,massMax));
    h2_HLTiso_M.push_back(new TH2F(Form("h2_HLTiso_M_%s",name.c_str()),Form("HLTiso_M_%s",name.c_str()),25,0.0,0.5,80,massMin,massMax));
  }

  Int_t nRegion = -1;

  while (reader.Next()) {

    //reader.SetLocalEntry(nEvents);
    
    cout.flush();
    if(nEvents % CHECK_EVERY_N == 0) {
      cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      //cout << "entry : " << nEvents << endl;     
    }
    nEvents++;

    if (not **HLT_EB && not **HLT_EE) continue;

    for (Int_t i = 0; i < *STr2_NPi0_rec; i++) {

      double pi0mass = STr2_mPi0_nocor[i];
      double pi0Pt = STr2_ptPi0_nocor[i];
      double ptOverM = pi0Pt/pi0mass;

      // note: quantities in the ntuples are ordered by energy of photons (actually energy of the seed)
      // eta, phi and other variables are associate to photon before correcting the energy
      // the actual order in energy of two photons can be swapped if a correction on one is high (or if the seed ordering is not equal to the cluster energy one)
      // therefore it is possible that ptGam1 < ptGam2, and therefore we must correct the order if ptGam must be the leading pT photon
      double ptGam1 = STr2_enG1_nocor[i]/cosh(STr2_Eta_1[i]);
      double ptGam2 = STr2_enG2_nocor[i]/cosh(STr2_Eta_2[i]);
      if (ptGam1 < ptGam2) {
	double tmp = ptGam2;
	ptGam2 = ptGam1;
	ptGam1 = tmp;
      }
      double ptGammaMin = min(ptGam1, ptGam2);
      double s4s9min = min(STr2_S4S9_1[i],STr2_S4S9_2[i]);
      int nXtalMin = min(STr2_n1CrisPi0_rec[i],STr2_n2CrisPi0_rec[i]);
      
      if (STr2_Pi0recIsEB[i]) {

	if (pi0mass < massMinEB || pi0mass > massMaxEB) continue;

	if (fabs(STr2_etaPi0_rec[i]) < EBetaRegionBoundary) {
	  nRegion = 0;
	} else {
	  nRegion = 1;
	  if (nXtalMin < 4) continue; // for region 2 in EB I forgot to set the minimum nXtal to 4 when producing the nuples, so I require it here
	}

	Bool_t binNotFound = true;
	for (UInt_t etaBin = 0; binNotFound && etaBin < etaBinEdgesEB.size()-1; etaBin++) {
	  if (fabs(STr2_etaPi0_rec[i]) < etaBinEdgesEB[etaBin+1]) {
	    h2_nXtal_M_EB_eta[etaBin]->Fill(nXtalMin, pi0mass);
	    binNotFound = false;
	  }
	}


      } else {
    
        if (pi0mass < massMinEE || pi0mass > massMaxEE) continue;
	
	// there should be 3 regions, but region3 has same thresholds as region2 
	if (fabs(STr2_etaPi0_rec[i]) < EEetaRegionBoundary) {
	  nRegion = 2;	 
	} else {
	  nRegion = 3;
	  if (fabs(STr2_etaPi0_rec[i]) > 2.0) nRegion = 4;
	}

      }

      double seedDistance = sqrt(fabs( (STr2_iEtaiX_1[i] - STr2_iEtaiX_2[i]) * (STr2_iPhiiY_1[i] - STr2_iPhiiY_2[i]) ));

      h2_ptGam_M[nRegion]->Fill(ptGammaMin, pi0mass);
      h2_ptPair_M[nRegion]->Fill(pi0Pt, pi0mass);
      h2_ptPairOverM_M[nRegion]->Fill(ptOverM, pi0mass);
      h2_s4s9_M[nRegion]->Fill(s4s9min, pi0mass);
      h2_nXtal_M[nRegion]->Fill(nXtalMin, pi0mass);
      h2_nXtal1_nXtal2[nRegion]->Fill(STr2_n1CrisPi0_rec[i], STr2_n2CrisPi0_rec[i]);
      h2_ptGam1_ptGam2[nRegion]->Fill(ptGam1, ptGam2);
      h2_ptGam_ptPair[nRegion]->Fill(ptGammaMin, pi0Pt);
      h2_ptGam_ptPairOverM[nRegion]->Fill(ptGammaMin, ptOverM);
      h2_ptPair_ptPairOverM[nRegion]->Fill(pi0Pt, ptOverM);
      h2_nXtal_ptGam[nRegion]->Fill(nXtalMin, ptGammaMin);
      h2_nXtal_ptPair[nRegion]->Fill(nXtalMin, pi0Pt);
      h2_nXtal_DRgam[nRegion]->Fill(nXtalMin, STr2_DeltaRG1G2[i]);
      h2_nXtal_s4s9[nRegion]->Fill(nXtalMin, s4s9min);
      h2_seedDistance_ptPair[nRegion]->Fill(seedDistance, pi0Pt);
      h2_seedDistance_nXtal[nRegion]->Fill(seedDistance, nXtalMin);
      h2_seedDistance_ptGam[nRegion]->Fill(seedDistance, ptGammaMin);
      h2_HLTiso_DRiso[nRegion]->Fill(STr2_HLTIsoPi0_rec[i],STr2_IsoPi0_rec[i]);
      h2_DRiso_M[nRegion]->Fill(STr2_IsoPi0_rec[i],pi0mass);
      h2_HLTiso_M[nRegion]->Fill(STr2_HLTIsoPi0_rec[i],pi0mass);

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

  string distanceAxisName = "";
  bool smoothPlot_Nxtal = true;

  for (UInt_t i = 0; i < regionTagName.size(); i++) {

    if (i < 2) distanceAxisName = "#sqrt{#Deltai#eta^{2}+#Deltai#phi^{2}}(seeds)"; // EB
    else       distanceAxisName = "#sqrt{#DeltaiX^{2}+#DeltaiY^{2}}(seeds)";     // EE

    drawCorrelationPlot(h2_ptGam_M[i], "min p_{T}(#gamma) [GeV]", "#gamma#gamma invariant mass [GeV]", h2_ptGam_M[i]->GetTitle(), plotTag[i], outDir);    
    drawCorrelationPlot(h2_ptPair_M[i], "p_{T}(#gamma#gamma) [GeV]", "#gamma#gamma invariant mass [GeV]", h2_ptPair_M[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_ptPairOverM_M[i], "p_{T}/M(#gamma#gamma)", "#gamma#gamma invariant mass [GeV]", h2_ptPairOverM_M[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_s4s9_M[i], "min S4/S9(#gamma)", "#gamma#gamma invariant mass [GeV]", h2_s4s9_M[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_nXtal_M[i], "min nXtal(#gamma)", "#gamma#gamma invariant mass [GeV]", h2_nXtal_M[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_nXtal1_nXtal2[i], "nXtal(leading #gamma)", "nXtal(trailing #gamma)", h2_nXtal1_nXtal2[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_ptGam1_ptGam2[i], "p_{T}(leading #gamma)", "p_{T}(trailing #gamma)", h2_ptGam1_ptGam2[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_ptGam_ptPair[i], "min p_{T}(#gamma) [GeV]", "p_{T}(#gamma#gamma) [GeV]", h2_ptGam_ptPair[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_ptGam_ptPairOverM[i], "min p_{T}(#gamma) [GeV]", "p_{T}/M(#gamma#gamma)", h2_ptGam_ptPairOverM[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_nXtal_ptGam[i], "min nXtal(#gamma)", "min p_{T}(#gamma) [GeV]", h2_nXtal_ptGam[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_nXtal_ptPair[i], "min nXtal(#gamma)", "p_{T}(#gamma#gamma) [GeV]", h2_nXtal_ptPair[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_nXtal_DRgam[i], "min nXtal(#gamma)", "#DeltaR(#gamma1,#gamma2)", h2_nXtal_DRgam[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_nXtal_s4s9[i], "min nXtal(#gamma)", "min S4/S9(#gamma)", h2_nXtal_s4s9[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_seedDistance_ptPair[i], distanceAxisName.c_str(), "p_{T}(#gamma#gamma) [GeV]", h2_seedDistance_ptPair[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_seedDistance_nXtal[i], distanceAxisName.c_str(), "min nXtal(#gamma)", h2_seedDistance_nXtal[i]->GetTitle(), plotTag[i], outDir, 1, smoothPlot_Nxtal);
    drawCorrelationPlot(h2_seedDistance_ptGam[i], distanceAxisName.c_str(), "min p_{T}(#gamma)", h2_seedDistance_ptGam[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_HLTiso_DRiso[i], "HLT isolation", "cluster isolation", h2_HLTiso_DRiso[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_DRiso_M[i], "cluster isolation", "#gamma#gamma invariant mass [GeV]", h2_DRiso_M[i]->GetTitle(), plotTag[i], outDir);
    drawCorrelationPlot(h2_HLTiso_M[i], "HLT isolation", "#gamma#gamma invariant mass [GeV]", h2_HLTiso_M[i]->GetTitle(), plotTag[i], outDir);
  }

  for (UInt_t i = 0; i < etaBinEdgesEB.size()-1; i++) {
    drawCorrelationPlot(h2_nXtal_M_EB_eta[i], "min nXtal(#gamma)", "#gamma#gamma invariant mass [GeV]", h2_nXtal_M_EB_eta[i]->GetTitle(), Form("EB %1.1f < |#eta| < %1.1f",etaBinEdgesEB[i],etaBinEdgesEB[i+1]), outDir, 1 , smoothPlot_Nxtal);
  }


  cout << "THE END!" << endl;

}
