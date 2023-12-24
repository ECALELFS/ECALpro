// -*- C++ -*-
//
// Package:    FitEpsilonPlot
// Class:      FitEpsilonPlot
// 
/**\class FitEpsilonPlot FitEpsilonPlot.cc CalibCode/FitEpsilonPlot/src/FitEpsilonPlot.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Marco Grassi, CMS
//         Created:  Tue Nov  8 17:18:54 CET 2011
// $Id: FitEpsilonPlot.cc,v 1.9 2013/06/17 13:40:42 lpernie Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TStyle.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RooGaussian.h"
#include "RooCBShape.h"
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
#include "RooFitLegacy/RooMinuit.h"
#include "RooMinimizer.h"
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooAbsTestStatistic.h" 

#include "CalibCode/FitEpsilonPlot/interface/FitEpsilonPlot.h"

using std::cout;
using std::endl;
using std::string;
using namespace RooFit;


// in this code the upper boundary of mass for the fit is used to assess the goodness of fit: look for --> if( fitres.chi2 < 5 && fabs(mean-<some_number>)>0.0000001)
// which means "if the Chi2 is good and the mean of the fit is far from the upper boundary ..."
// the upper boundary must be consistent with <some_number>
static double upper_bound_pi0mass_EB = 0.15;
static double upper_bound_pi0mass_EE = 0.16;
static double upper_bound_etamass_EB = 0.6;
static double upper_bound_etamass_EE = 0.6;

static float fitRange_low_pi0 = 0.080; // value used in the fit function to define the fit range
static float fitRange_high_pi0 = 0.212; // value used in the fit function to define the fit range
static float fitRange_high_pi0_ext = 0.222;

static float fitRange_low_eta = 0.380; // value used in the fit function to define the fit range
static float fitRange_low_etaEE = 0.360; // value used in the fit function to define the fit range
static float fitRange_high_eta = 0.680; // value used in the fit function to define the fit range
static float fitRange_high_eta_ext = 0.700;

static float EoverEtrue_integralMin = 25; // require that integral in a given range is > EoverEtrue_integralMin for E/Etrue distribution (used for MC only)

FitEpsilonPlot::FitEpsilonPlot(const edm::ParameterSet& iConfig)

{

    //now do what ever initialization is needed
    currentIteration_ =  iConfig.getUntrackedParameter<int>("CurrentIteration");
    epsilonPlotFileName_ = iConfig.getUntrackedParameter<std::string>("EpsilonPlotFileName");
    outputDir_ = iConfig.getUntrackedParameter<std::string>("OutputDir");
    outfilename_          = iConfig.getUntrackedParameter<std::string>("OutputFile");
    calibMapPath_ = iConfig.getUntrackedParameter<std::string>("calibMapPath");
    inRangeFit_ = iConfig.getUntrackedParameter<int>("NInFit");
    finRangeFit_ = iConfig.getUntrackedParameter<int>("NFinFit");    
    EEoEB_ = iConfig.getUntrackedParameter<std::string>("EEorEB");
    isNot_2010_ = iConfig.getUntrackedParameter<bool>("isNot_2010");
    Are_pi0_ = iConfig.getUntrackedParameter<bool>("Are_pi0");
    StoreForTest_ = iConfig.getUntrackedParameter<bool>("StoreForTest",true);
    Barrel_orEndcap_ = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");
    useMassInsteadOfEpsilon_ = iConfig.getUntrackedParameter<bool>("useMassInsteadOfEpsilon",true);
    isEoverEtrue_ = iConfig.getUntrackedParameter<bool>("isEoverEtrue",false);
    useFit_RooMinuit_ = iConfig.getUntrackedParameter<bool>("useFit_RooMinuit",false);
    foldInSuperModule_ = iConfig.getUntrackedParameter<bool>("foldInSuperModule",false);
    makeFoldedHistograms_ = iConfig.getUntrackedParameter<bool>("makeFoldedHistograms",false);

    // apparently for E/Etrue the fits are much better (I tried RooCMSShape + double-Crystal-Ball)
    // some tuning might be required, though
    // if (isEoverEtrue_) useFit_RooMinuit_ = false;
    // eventually I decided to use FitTo method for E/Etrue

    //foldInSuperModule_ = true;
    fitEoverEtrueWithRooFit_ = true;   // use bare TH1::Fit or RooFit (better, can stay true)

    // read directly folded histograms (the folding is done in this analyzer, so the very first time this option is false)
    readFoldedHistogramFromFile_ = makeFoldedHistograms_ ? false : true; 

    foldEB_all0_onlyPlus1_onlyMinus2_ = 0; // 0 to put all 36 SM in one, 1 for using EB+ only, 2 for using EB- only (but then they are used on all barrel because I only have a single SM map)

    // I should add a code that do the folding before going to the fitting part

    fitFileName_ = outfilename_;
    std::string strToReplace = "calibMap";
    fitFileName_.replace(outfilename_.find(strToReplace.c_str()),strToReplace.size(),"fitRes");
    fitFileName_ = outputDir_ + "/" + fitFileName_;

    inputEpsilonFile_ = nullptr; //  before we open the file, we assign a null value
    outfile_ = nullptr;

    /// setting calibration type
    calibTypeString_ = iConfig.getUntrackedParameter<std::string>("CalibType");
    if(     calibTypeString_.compare("xtal")    == 0 ) {
      calibTypeNumber_ = xtal;
      regionalCalibration_ = &xtalCalib;
      if (isEoverEtrue_) regionalCalibration_g2_ = new EcalRegionalCalibration<EcalCalibType::Xtal>();  //regionalCalibration_g2_ = &xtalCalib_g2;
    } else if(calibTypeString_.compare("tt")      == 0 ) {
      calibTypeNumber_ = tt;
      regionalCalibration_ = &TTCalib;
      if (isEoverEtrue_) regionalCalibration_g2_ = new EcalRegionalCalibration<EcalCalibType::TrigTower>(); // regionalCalibration_g2_ = &TTCalib_g2;
    } else if(calibTypeString_.compare("etaring") == 0 ) {
      calibTypeNumber_ = etaring;
      regionalCalibration_ = &etaCalib;
      if (isEoverEtrue_) regionalCalibration_g2_ = new EcalRegionalCalibration<EcalCalibType::EtaRing>(); // regionalCalibration_g2_ = &etaCalib_g2;
    } else throw cms::Exception("CalibType") << "Calib type not recognized\n";

    cout << "FIT_EPSILON: crosscheck: selected type: " << regionalCalibration_->printType() << endl;

    /// retrieving calibration coefficients of the previous iteration
    // if currentIteration_ = 0, calibMapPath_ contains "iter_-1" unless the current set of ICs was started from another existing set (see parameters.py)
    // therefore, the case with extension is included below
    std::string stringToMatch = "iter_-1";  // used below: this string should not match to trigger true condition 
    if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
    else if(currentIteration_ > 0 || (currentIteration_ == 0 && calibMapPath_.find(stringToMatch)==std::string::npos))
    {
      regionalCalibration_->getCalibMap()->loadCalibMapFromFile(calibMapPath_.c_str(),false);
      if (isEoverEtrue_) regionalCalibration_g2_->getCalibMap()->loadCalibMapFromFile(calibMapPath_.c_str(),true);
    }

    TH1::SetDefaultSumw2(); // all new histograms will automatically activate the storage of the sum of squares of errors (i.e, TH1::Sumw2 is automatically called).

    // load epsilon from current iter
    if (isEoverEtrue_) {

      if ((Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" )) {
	EoverEtrue_g1_EB_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEB()];
	EoverEtrue_g2_EB_h = new TH1F*[regionalCalibration_g2_->getCalibMap()->getNRegionsEB()];
      }
      if ((Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" )) {
	EoverEtrue_g1_EE_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEE()];
	EoverEtrue_g2_EE_h = new TH1F*[regionalCalibration_g2_->getCalibMap()->getNRegionsEE()];
      }
      cout << "FIT_EPSILON: FitEpsilonPlot:: loading EoverEtrue plots from file: " << epsilonPlotFileName_ << endl;
      loadEoverEtruePlot(epsilonPlotFileName_,1);
      loadEoverEtruePlot(epsilonPlotFileName_,2);

      if (foldInSuperModule_ && EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE")) {

	EoverEtrue_g1_EB_SM_hvec.clear();
	EoverEtrue_g2_EB_SM_hvec.clear();

	// 20(phi)*85(ieta) crystals in 1 SM
	// we use a 1D vector ad treat it as a 2D array: the index convention is 
	// index = (fabs(ieta) - 1) + 85 * ((iphi - 1)%20)    ieta = 1, 2, ..., 85     iphi = 1, 2 , ..., 360 
	// ieta,iphi = 1,1 --> index = 0

	// ieta 1 --> 85
	//
	//  * * * * * * * * * * * * * *  iphi 1
	//  * * * * * * * . . . . . . .
	//  . . . . . . . . . . . . . . 
        //
	//  * * * * * * * . . . . . . .  iphi 20

	//cout << "EBDetId::kCrystalsPerSM = " << EBDetId::kCrystalsPerSM << endl;
	for (int iv = 0; iv < EBDetId::kCrystalsPerSM; ++iv) {  // 1700 crystals
	  
	    EoverEtrue_g1_EB_SM_hvec.push_back( new TH1F(Form("EoverEtrue_g1_EB_SM_hvec_%d",iv),
							 "g1 E/E_{true} folded in SM",
							 EoverEtrue_g1_EB_h[inRangeFit_]->GetNbinsX(),
							 EoverEtrue_g1_EB_h[inRangeFit_]->GetBinLowEdge(1),
							 EoverEtrue_g1_EB_h[inRangeFit_]->GetBinLowEdge(1+EoverEtrue_g1_EB_h[inRangeFit_]->GetNbinsX())
							 ) );
	    EoverEtrue_g2_EB_SM_hvec.push_back( new TH1F(Form("EoverEtrue_g2_EB_SM_hvec_%d",iv),
							 "g2 E/E_{true} folded in SM",
							 EoverEtrue_g2_EB_h[inRangeFit_]->GetNbinsX(),
							 EoverEtrue_g2_EB_h[inRangeFit_]->GetBinLowEdge(1),
							 EoverEtrue_g2_EB_h[inRangeFit_]->GetBinLowEdge(1+EoverEtrue_g2_EB_h[inRangeFit_]->GetNbinsX())
							 ) );

	}

	if (readFoldedHistogramFromFile_) {
	  cout << "FIT_EPSILON: reading folded histogram from file" << endl;
	  loadEoverEtruePlotFoldedInSM(1);
	  loadEoverEtruePlotFoldedInSM(2);
	} else {
	  cout << "FIT_EPSILON: folding histograms ..." << endl;
	  addHistogramsToFoldSM(EoverEtrue_g1_EB_SM_hvec,epsilonPlotFileName_,1);
	  addHistogramsToFoldSM(EoverEtrue_g2_EB_SM_hvec,epsilonPlotFileName_,2);
	  cout << "FIT_EPSILON: folding histograms completed..." << endl;
	}

      }

    } else {

      if ((Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" )) {
	epsilon_EB_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEB()];
      }
      if ((Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" )) {
	epsilon_EE_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEE()];
      }
      cout << "FIT_EPSILON: FitEpsilonPlot:: loading epsilon plots from file: " << epsilonPlotFileName_ << endl;
      //loadEpsilonPlot(epsilonPlotFileName_);
      loadEpsilonPlot2D(epsilonPlotFileName_);

      if (foldInSuperModule_ && EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE")) {

	epsilon_EB_SM_hvec.clear();

	// 20(phi)*85(ieta) crystals in 1 SM
	// we use a 1D vector ad treat it as a 2D array: the index convention is 
	// index = (fabs(ieta) - 1) + 85 * ((iphi - 1)%20)    ieta = 1, 2, ..., 85     iphi = 1, 2 , ..., 360 
	// ieta,iphi = 1,1 --> index = 0

	// ieta 1 --> 85
	//
	//  * * * * * * * * * * * * * *  iphi 1
	//  * * * * * * * . . . . . . .
	//  . . . . . . . . . . . . . . 
        //
	//  * * * * * * * . . . . . . .  iphi 20

	//cout << "EBDetId::kCrystalsPerSM = " << EBDetId::kCrystalsPerSM << endl;
	for (int iv = 0; iv < EBDetId::kCrystalsPerSM; ++iv) {  // 1700 crystals
	  
	  epsilon_EB_SM_hvec.push_back( new TH1F(Form("epsilon_EB_SM_hvec_%d",iv),
						 "#pi^{0} mass folded in SM",
						 epsilon_EB_h[inRangeFit_]->GetNbinsX(),
						 epsilon_EB_h[inRangeFit_]->GetBinLowEdge(1),
						 epsilon_EB_h[inRangeFit_]->GetBinLowEdge(1+epsilon_EB_h[inRangeFit_]->GetNbinsX())
						 ) );

	}

	if (readFoldedHistogramFromFile_) {
	  cout << "FIT_EPSILON: reading folded histogram from file" << endl;
	  loadEpsilonPlotFoldedInSM();
	} else {
	  cout << "FIT_EPSILON: folding histograms ..." << endl;
	  addHistogramsToFoldSM(epsilon_EB_SM_hvec,epsilonPlotFileName_,1);
	  cout << "FIT_EPSILON: folding histograms completed..." << endl;
	}

      }

    }

}


FitEpsilonPlot::~FitEpsilonPlot()
{

  cout << "Beginning of destructor" << endl;

  if ((Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" )) {

    if (isEoverEtrue_) {
      deleteEpsilonPlot(EoverEtrue_g1_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
      deleteEpsilonPlot(EoverEtrue_g2_EB_h, regionalCalibration_g2_->getCalibMap()->getNRegionsEB() );
      // delete EoverEtrue_g1_EB_h;
      // delete EoverEtrue_g2_EB_h;
      if (foldInSuperModule_) {
	for (unsigned int i = 0; i < EoverEtrue_g1_EB_SM_hvec.size(); ++i) {
	  delete EoverEtrue_g1_EB_SM_hvec[i];
	  delete EoverEtrue_g2_EB_SM_hvec[i];
	}
	EoverEtrue_g1_EB_SM_hvec.clear();
	EoverEtrue_g2_EB_SM_hvec.clear();
      }
    } else {
      deleteEpsilonPlot(epsilon_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
      if (foldInSuperModule_) {
	for (unsigned int i = 0; i < epsilon_EB_SM_hvec.size(); ++i) {
	  delete epsilon_EB_SM_hvec[i];
	}
	epsilon_EB_SM_hvec.clear();
      }
      // delete epsilon_EB_h;
    }

  }

  cout << "After EB" << endl;


  if ((Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" )) {

    if (isEoverEtrue_) {
      deleteEpsilonPlot(EoverEtrue_g1_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
      deleteEpsilonPlot(EoverEtrue_g2_EE_h, regionalCalibration_g2_->getCalibMap()->getNRegionsEE() );
      // delete EoverEtrue_g1_EE_h;
      // delete EoverEtrue_g2_EE_h;
    } else {
      deleteEpsilonPlot(epsilon_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
      // delete epsilon_EE_h;
    }

  }

  cout << "After EE" << endl;

  // do not close file, it produces a segmentation fault
  if (inputEpsilonFile_ != nullptr) {
    if (inputEpsilonFile_->IsOpen()) inputEpsilonFile_->Close();
  }
  cout << "End of destructor" << endl;

}


//
// member functions
//

int FitEpsilonPlot::getArrayIndexOfFoldedSMfromIetaIphi(const int ieta = 1, const int iphi = 1) {

  // note that the index in SM returned by this function in not the same as the index returned by EBDetId::ic()
  // the difference is mainly in the folding of EB+ on EB-
  // In our case, we overlay crystals such that ieta,iphi=20,40 goes on ieta,iphi=-20,40, 
  // i.e. the iphi coordinate is preserved when we consider two facing SM in EB+ and EB-
  // The usual CMSSW numbering scheme for crystals in SM is such that, looking at the center of the barrel, the crystal number 1 is always on the left
  // which means that the folding would overlay ieta,iphi=20,40 on ieta,iphi=-20,21
  
  // first 85 crystals correspond to iphi = 1 (in a SM)
  return (fabs(ieta) - 1) + EBDetId::kCrystalsInEta * ((iphi - 1) % EBDetId::kCrystalsInPhi);
  
}


int FitEpsilonPlot::getArrayIndexOfFoldedSMfromDenseIndex(const int index = 1, const bool useEBDetId_ic_scheme = true) {

  // when using EBDetId::ic(), the crystal number in SM is given such that ic=1 has iphi=1 in EB+ and iphi=20 in EB- (or the opposite, I don't remember)
  // the idea is that with ic() the crystal number is increased going from left to right

  EBDetId thisEBcrystal(EBDetId::detIdFromDenseIndex( index ));
  if (useEBDetId_ic_scheme) return thisEBcrystal.ic()-1;
  else                      return getArrayIndexOfFoldedSMfromIetaIphi(thisEBcrystal.ietaAbs(),thisEBcrystal.iphi());
  //return getArrayIndexOfFoldedSMfromIetaIphi(thisEBcrystal.ietaAbs(),thisEBcrystal.iphi());
}


void FitEpsilonPlot::addHistogramsToFoldSM(std::vector<TH1F*>& hvec, const std::string& filename, const int whichPhoton = 1) {

  if (hvec.size() == 0) throw cms::Exception("addHistogramsToFoldSM") << "Vector passed to function has size 0\n"; 

  std::string line = "";
  std::string histoNamePattern = isEoverEtrue_ ? Form("%s/EoverEtrue_g%d",EEoEB_.c_str(),whichPhoton) : Form("%s/epsilon",EEoEB_.c_str());

  // open the file if it has not been created so far, otherwise check that it is still open (this would happen on second photon)
  if (inputEpsilonFile_ == nullptr) {
    inputEpsilonFile_ = TFile::Open(filename.c_str(),"READ");
    if(!inputEpsilonFile_) 
      throw cms::Exception("addHistogramsToFoldSM") << "Cannot open file " << filename << "\n"; 
  } else if (not inputEpsilonFile_->IsOpen()) {
    inputEpsilonFile_ = TFile::Open(filename.c_str(),"READ");
  }

  ////////////////////
  // CAVEAT !!
  // If opening the following file, before writing objects we should do TFile::cd() (with the other files where histograms are saved)
  // This is because Root changes the current directory and messes up the filesystem
  //////////////////
  // copy path of current directory
  const char* currFilePath = gDirectory->GetPath();

  // // create file containing folded histograms (could be used later without folding again)
  string foldFileName = filename;
  std::string strToReplace = "epsilonPlots";
  foldFileName.replace(filename.find(strToReplace.c_str()),strToReplace.size(),"histograms_foldedInSM");

  string foldFileOpeningMode = (whichPhoton == 1) ? "RECREATE" : "UPDATE";
  TFile* f = TFile::Open(foldFileName.c_str(),foldFileOpeningMode.c_str());
  if (!f || !f->IsOpen()) {
    throw cms::Exception("FitEpsilonPlot") << "error opening file '" << foldFileName << "' to save folded histogram\n";
  }

  TH1F* htmp = nullptr;

  // if we are here it means we are already in EB, but let's ask again
  if ( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {

    int nRegionsEB = ((whichPhoton == 1) ? regionalCalibration_->getCalibMap()->getNRegionsEB() : regionalCalibration_g2_->getCalibMap()->getNRegionsEB()); 
    
    for (int iR = 0; iR < nRegionsEB; ++iR) {

      EBDetId thisEBcrystal(EBDetId::detIdFromDenseIndex( iR));
      if (foldEB_all0_onlyPlus1_onlyMinus2_ == 1 && thisEBcrystal.ieta() < 0) continue; // if we want to use only EB+
      if (foldEB_all0_onlyPlus1_onlyMinus2_ == 2 && thisEBcrystal.ieta() > 0) continue; // if we want to use only EB-

      line = Form("%s_EB_iR_%d",histoNamePattern.c_str(), iR);
      //if (isTest) line = histoNamePattern;

      htmp = (TH1F*)inputEpsilonFile_->Get(line.c_str());      
      if(!htmp)	throw cms::Exception("addHistogramsToFoldSM") << "FIT_EPSILON: cannot load histogram " << line << "\n";
      
      int crystalIndexInSM = getArrayIndexOfFoldedSMfromDenseIndex(iR);
      if (crystalIndexInSM >= EBDetId::kCrystalsPerSM) {
	std::cout << "FIT_EPSILON: error in SM folding, index = " << crystalIndexInSM << std::endl;
	throw cms::Exception("FitEpsilonPlot") << "crystalIndexInSM >= " << EBDetId::kCrystalsPerSM << "\n";
      }
      if (htmp->GetEntries() > 0) {
	bool AddWasSuccesful = hvec[crystalIndexInSM]->Add(htmp);
	if (not AddWasSuccesful) throw cms::Exception("addHistogramsToFoldSM") << "FIT_EPSILON: failed to add histogram " << line << "\n";
	//if (crystalIndexInSM == 0) std::cout << "EoverEtrue_g1_EB_SM_hvec[0]->Integral = " << hvec[crystalIndexInSM] << std::endl;
      }

    }

  }

  // // save folded histogrmas
  f->cd();
  for (unsigned int i = 0; i < hvec.size(); i++) {
    hvec[i]->Write();
  }
  f->Close();

  // now restore previous path is ROOT filesystem (so to get back to previous file)
  if (currFilePath != nullptr) { 
    gROOT->cd(currFilePath);
  } else {
    throw cms::Exception("addHistogramsToFoldSM") << "Could not restore path in currFilePath: variable was empty\n";
  }

}

void FitEpsilonPlot::loadEoverEtruePlot(const std::string& filename, const int whichPhoton = 1) {

  // here regionalCalibration_ is only used to get the number of regions, which is the same for both photons
  // hence, no need to use regionalCalibration_ or regionalCalibration_g2_

  std::string line = "";
  std::string histoNamePattern = Form("EoverEtrue_g%d",whichPhoton );

  // test the machinary fitting inclusive histogram (otherwise I have no statistics)
  // bool isTest = true;
  // if (isTest) {
  //   if (whichPhoton == 1) {
  //     if (EEoEB_ == "Barrel") histoNamePattern = "allEoverEtrue_g1_EB";
  //     else histoNamePattern = "allEoverEtrue_g1_EE";
  //   } else {
  //     if (EEoEB_ == "Barrel") histoNamePattern = "allEoverEtrue_g2_EB";
  //     else histoNamePattern = "allEoverEtrue_g2_EE";
  //   }
  // }

  // open the file if it has not been created so far, otherwise check that it is still open (this would happen on second photon)
  if (inputEpsilonFile_ == nullptr) {
    inputEpsilonFile_ = TFile::Open(filename.c_str());
    if(!inputEpsilonFile_) 
      throw cms::Exception("loadEpsilonPlot") << "Cannot open file " << filename << "\n"; 
  } else if (not inputEpsilonFile_->IsOpen()) {
    inputEpsilonFile_ = TFile::Open(filename.c_str());
  }

  if ( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    
    for (int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++) {

      line = Form("%s_EB_iR_%d",histoNamePattern.c_str(), iR);
      //if (isTest) line = histoNamePattern;

      if (whichPhoton == 1) {
	EoverEtrue_g1_EB_h[iR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());      
	if(!EoverEtrue_g1_EB_h[iR])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else if(!(iR%1000))
	  cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EB region " << iR << " loaded" << endl;
      } else {
	EoverEtrue_g2_EB_h[iR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());      
	if(!EoverEtrue_g2_EB_h[iR])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else if(!(iR%1000))
	  cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EB region " << iR << " loaded" << endl;
      }
      
    }

  } else if( EEoEB_ == "Endcap" && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {

    for (int jR=inRangeFit_; jR <= finRangeFit_ && jR<EEDetId::kSizeForDenseIndexing; jR++) {
      
      line = Form("%s_EE_iR_%d",histoNamePattern.c_str(), jR);
      //if (isTest) line = histoNamePattern;

      if (whichPhoton == 1) {

	EoverEtrue_g1_EE_h[jR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());
	if(!EoverEtrue_g1_EE_h[jR])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else if(!(jR%1000))
	  cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EE region " << jR << " loaded" << endl;

      } else {

	EoverEtrue_g2_EE_h[jR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());
	if(!EoverEtrue_g2_EE_h[jR])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else if(!(jR%1000))
	  cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EE region " << jR << " loaded" << endl;
	
      }

    }

  }

}

//============================================================

void FitEpsilonPlot::loadEoverEtruePlotFoldedInSM(const int whichPhoton = 1) {

  string foldFileName = epsilonPlotFileName_;
  std::string strToReplace = "epsilonPlots";
  foldFileName.replace(epsilonPlotFileName_.find(strToReplace.c_str()),strToReplace.size(),"histograms_foldedInSM");
  std::string line = "";
  std::string histoNamePattern = Form("EoverEtrue_g%d_EB_SM_hvec",whichPhoton );

  TFile* fileFoldHistogram  = TFile::Open(foldFileName.c_str(),"READ");
  if(!fileFoldHistogram) 
    throw cms::Exception("loadEoverEtruePlotFoldedInSM") << "Cannot open file " << foldFileName << "\n"; 

  if ( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    
    for (int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++) {

      // when folding into SM, we interpret iR as the EBDetId::ic() number (which goes from 1 to 1700, so need to subtract 1)
      //int indexSM = getArrayIndexOfFoldedSMfromDenseIndex(iR);
      int indexSM = iR;
      line = Form("%s_%d",histoNamePattern.c_str(), indexSM);

      if (whichPhoton == 1) {
	EoverEtrue_g1_EB_SM_hvec[indexSM] = (TH1F*)fileFoldHistogram->Get(line.c_str());      
	if(!EoverEtrue_g1_EB_SM_hvec[indexSM])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else {
	  EoverEtrue_g1_EB_SM_hvec[indexSM]->SetDirectory(0);
	  if(!(iR%1000))
	    cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EB region " << iR << " loaded" << endl;
	}
      } else {
	EoverEtrue_g2_EB_SM_hvec[indexSM] = (TH1F*)fileFoldHistogram->Get(line.c_str());      
	if(!EoverEtrue_g2_EB_SM_hvec[indexSM])
	  throw cms::Exception("loadEoverEtruePlot") << "Cannot load histogram " << line << "\n";
	else {
	  EoverEtrue_g2_EB_SM_hvec[indexSM]->SetDirectory(0);
	  if(!(iR%1000))
	    cout << "FIT_EPSILON: EoverEtrue distribution (photon " << whichPhoton << ") for EB region " << iR << " loaded" << endl;
	}
      }

    }

  }
  
  fileFoldHistogram->Close();

}

//============================================================

void FitEpsilonPlot::loadEpsilonPlotFoldedInSM() {

  string foldFileName = epsilonPlotFileName_;
  std::string strToReplace = "epsilonPlots";
  foldFileName.replace(epsilonPlotFileName_.find(strToReplace.c_str()),strToReplace.size(),"histograms_foldedInSM");
  std::string line = "";
  std::string histoNamePattern = "epsilon_EB_SM_hvec";

  TFile* fileFoldHistogram  = TFile::Open(foldFileName.c_str(),"READ");
  if(!fileFoldHistogram) 
    throw cms::Exception("loadEpsilonPlotFoldedInSM") << "Cannot open file " << foldFileName << "\n"; 

  if ( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    
    for (int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++) {

      int indexSM = getArrayIndexOfFoldedSMfromDenseIndex(iR);
      line = Form("%s_%d",histoNamePattern.c_str(), indexSM);

      epsilon_EB_SM_hvec[indexSM] = (TH1F*)fileFoldHistogram->Get(line.c_str());      
      if(!epsilon_EB_SM_hvec[indexSM])
	throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << line << "\n";
      else {
	epsilon_EB_SM_hvec[indexSM]->SetDirectory(0);
	if(!(iR%1000))
	  cout << "FIT_EPSILON: epsilon distribution for EB region " << iR << " loaded" << endl;
      }

    }

  }
  
  fileFoldHistogram->Close();

}

//============================================================

void FitEpsilonPlot::loadEpsilonPlot(const std::string& filename)
{
  std::string line = "";

  inputEpsilonFile_ = TFile::Open(filename.c_str());
  if(!inputEpsilonFile_) 
    throw cms::Exception("loadEpsilonPlot") << "Cannot open file " << filename << "\n"; 
  if( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
    for(int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++)
      {
	line = Form("epsilon_EB_iR_%d",iR);
	epsilon_EB_h[iR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());

	if(!epsilon_EB_h[iR])
	  throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << line << "\n";
	else if(!(iR%1000))
	  cout << "FIT_EPSILON: Epsilon distribution for EB region " << iR << " loaded" << endl;
	epsilon_EB_h[iR]->SetDirectory(0);
      }
  }
  else if( EEoEB_ == "Endcap" && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
    for(int jR=inRangeFit_; jR <= finRangeFit_ && jR<EEDetId::kSizeForDenseIndexing; jR++)
      {
	line = Form("epsilon_EE_iR_%d",jR);
	epsilon_EE_h[jR] = (TH1F*)inputEpsilonFile_->Get(line.c_str());
	if(!epsilon_EE_h[jR])
	  throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << line << "\n";
	else if(!(jR%1000))
	  cout << "FIT_EPSILON: Epsilon distribution for EE region " << jR << " loaded" << endl;
	epsilon_EE_h[jR]->SetDirectory(0);
      }
  }

}

void FitEpsilonPlot::loadEpsilonPlot2D(const std::string& filename)
{
  std::string line = "";

  inputEpsilonFile_ = TFile::Open(filename.c_str());
  if(!inputEpsilonFile_) 
    throw cms::Exception("loadEpsilonPlot2D") << "Cannot open file " << filename << "\n"; 

  if( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ){

    line = "epsilon_EB_iR";
    TH2F* h2_tmp_epsilon = (TH2F*)inputEpsilonFile_->Get(line.c_str());
    if(!h2_tmp_epsilon)
      throw cms::Exception("loadEpsilonPlot2D") << "Cannot load histogram " << line << "\n";    

    for(int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++)
      {
	// take slice of TH2 at index whose bin center is iR. 
	// Since ProjectionX needs the bin number, use iR+1 (iR goes from 0 to N(xtals)-1 ) 
	// range iR+1 to iR+1 will just select bin iR+1
	// do projection to return TH1D*, then clone inot TH1F* (should be the same except for precision loss)
	// but precision loss should be negligible (if not, convert everything to TH1D*)
	epsilon_EB_h[iR] = (TH1F*) (h2_tmp_epsilon->ProjectionX(Form("proj_epsilon_EB_h%d",iR),iR+1,iR+1,"e"))->Clone(Form("epsilon_EB_h%d",iR));

	if(!epsilon_EB_h[iR])
	  throw cms::Exception("loadEpsilonPlot2D") << "Cannot load histogram " << line << "\n";
	else if(!(iR%1000))
	  cout << "FIT_EPSILON: Epsilon distribution for EB region " << iR << " loaded" << endl;
	epsilon_EB_h[iR]->SetDirectory(0);
      }
    

  }
  else if( EEoEB_ == "Endcap" && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){

    line = Form("epsilon_EE_iR");
    TH2F* h2_tmp_epsilon = (TH2F*)inputEpsilonFile_->Get(line.c_str());
    if(!h2_tmp_epsilon)
      throw cms::Exception("loadEpsilonPlot2D") << "Cannot load histogram " << line << "\n";

    for(int jR=inRangeFit_; jR <= finRangeFit_ && jR<EEDetId::kSizeForDenseIndexing; jR++)
      {
	epsilon_EE_h[jR] = (TH1F*) (h2_tmp_epsilon->ProjectionX(Form("proj_epsilon_EE_h%d",jR),jR+1,jR+1,"e"))->Clone(Form("epsilon_EE_h%d",jR));
	if(!epsilon_EE_h[jR])
	  throw cms::Exception("loadEpsilonPlot2D") << "Cannot load histogram " << line << "\n";
	else if(!(jR%1000))
	  cout << "FIT_EPSILON: Epsilon distribution for EE region " << jR << " loaded" << endl;
	epsilon_EE_h[jR]->SetDirectory(0);
      }
  }

}



void  FitEpsilonPlot::deleteEpsilonPlot(TH1F **h, int size)
{
    for(int jR=0; jR<size; jR++)
	  delete h[jR];

    //delete h; // do not delete it, otherwise it makes the code crash and the end of the destructor
}

// void  FitEpsilonPlot::deleteEpsilonPlot2D(TH2F *h)
// {
//   // is it needed? Probably not, and might make the code crash
//   delete h;
// }


void FitEpsilonPlot::saveCoefficients() 
{
  /// output file
  std::string fileName = outputDir_  + "/" + outfilename_;
  outfile_ = new TFile(fileName.c_str(),"RECREATE");
  cout << "FIT_EPSILON: Saving Calibration Coefficients in " << fileName << " ... " << endl;;
  if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << fileName << "\n";
  outfile_->cd();

  // 2D calib map in the barrel
  TH2F* hmap_EB = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y",
			   2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,
			   EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
  hmap_EB->GetXaxis()->SetTitle("i#eta");
  hmap_EB->GetYaxis()->SetTitle("i#phi");
  TH2F* hmap_EEp = new TH2F("calibMap_EEp","EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEp->GetXaxis()->SetTitle("ix");
  hmap_EEp->GetYaxis()->SetTitle("iy");
  TH2F* hmap_EEm = new TH2F("calibMap_EEm","EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEm->GetXaxis()->SetTitle("ix");
  hmap_EEm->GetYaxis()->SetTitle("iy");
  TH1F* hint = new TH1F("hint","Bin1: inRangeFit_ Bin2: finRangeFit_ Bin3: Barrel(0)/Endcap(1)",3,0.,3.);
  hint->SetBinContent(1,inRangeFit_);
  hint->SetBinContent(2,finRangeFit_);
  if( EEoEB_ == "Barrel" ) hint->SetBinContent(3,0);
  else                     hint->SetBinContent(3,1);
  hint->Write();

  /// filling Barrel Map
  for(int j=0; j<regionalCalibration_->getCalibMap()->getNRegionsEB(); ++j)  
    {
      std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
	EBDetId ebid(*iid);
	int ix = ebid.ieta()+EBDetId::MAX_IETA+1;

	float coeffValue = regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ? regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;
	hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
      } // loop over DetId in regions
    }
  hmap_EB->SetMinimum(0.9);
  hmap_EB->SetStats(false);
  hmap_EB->Write();

  for(int jR=0; jR < regionalCalibration_->getCalibMap()->getNRegionsEE(); jR++)
    {
      std::vector<DetId> ids =  regionalCalibration_->allDetIdsInEERegion(jR);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	{ 
	  EEDetId eeid(*iid);
	  float coeffValue =  regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ?  regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;

	  if(eeid.positiveZ())
	    hmap_EEp->Fill(eeid.ix(), eeid.iy(), coeffValue); 
	  else 
	    hmap_EEm->Fill(eeid.ix(), eeid.iy(), coeffValue);
	}
    }

  hmap_EEp->SetMinimum(0.9);
  hmap_EEp->SetStats(false);
  hmap_EEp->Write();

  hmap_EEm->SetMinimum(0.9);
  hmap_EEm->SetStats(false);
  hmap_EEm->Write();

  /*------------- TTREE --------------*/

  uint32_t   rawId;
  int        hashedIndex;
  int        ieta;
  int        iphi;
  int        iSM;
  int        iMod;
  int        iTT;
  int        iTTeta;
  int        iTTphi;
  int        iter = currentIteration_;
  float      regCoeff;
  float      Signal;//#
  float      Backgr; 
  float      Chisqu; 
  float      Ndof; 
  float      fit_mean;
  float      fit_mean_err;
  float      fit_sigma;
  float      fit_Snorm;
  float      fit_b0;
  float      fit_b1;    
  float      fit_b2;    
  float      fit_b3;    
  float      fit_Bnorm; 
  /// endcap variables
  int ix;
  int iy;
  int zside;
  int sc; 
  int isc;
  int ic;
  int iquadrant;

  TTree* treeEB = new TTree("calibEB","Tree of EB Inter-calibration constants");
  TTree* treeEE = new TTree("calibEE","Tree of EE Inter-calibration constants");


  /// barrel
  treeEB->Branch("rawId",&rawId,"rawId/i");
  treeEB->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEB->Branch("ieta",&ieta,"ieta/I");
  treeEB->Branch("iphi",&iphi,"iphi/I");
  treeEB->Branch("iSM",&iSM,"iSM/I");
  treeEB->Branch("iMod",&iMod,"iMod/I");
  treeEB->Branch("iTT",&iTT,"iTT/I");
  treeEB->Branch("iTTeta",&iTTeta,"iTTeta/I");
  treeEB->Branch("iTTphi",&iTTphi,"iTTphi/I");
  treeEB->Branch("iter",&iter,"iter/I");
  treeEB->Branch("coeff",&regCoeff,"coeff/F");
  treeEB->Branch("Signal",&Signal,"Signal/F");//#
  treeEB->Branch("Backgr",&Backgr,"Backgr/F");
  treeEB->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEB->Branch("Ndof",&Ndof,"Ndof/F");
  treeEB->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEB->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEB->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEB->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEB->Branch("fit_b0",&fit_b0,"fit_b0/F");
  treeEB->Branch("fit_b1",&fit_b1,"fit_b1/F");
  treeEB->Branch("fit_b2",&fit_b2,"fit_b2/F");
  treeEB->Branch("fit_b3",&fit_b3,"fit_b3/F");
  treeEB->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

  /// endcap
  treeEE->Branch("ix",&ix,"ix/I");
  treeEE->Branch("iy",&iy,"iy/I");
  treeEE->Branch("zside",&zside,"zside/I");
  treeEE->Branch("sc",&sc,"sc/I");
  treeEE->Branch("isc",&isc,"isc/I");
  treeEE->Branch("ic",&ic,"ic/I");
  treeEE->Branch("iquadrant",&iquadrant,"iquadrant/I");
  treeEE->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEE->Branch("iter",&iter,"iter/I");
  treeEE->Branch("coeff",&regCoeff,"coeff/F");
  treeEE->Branch("Signal",&Signal,"Signal/F");//#
  treeEE->Branch("Backgr",&Backgr,"Backgr/F");
  treeEE->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEE->Branch("Ndof",&Ndof,"Ndof/F");
  treeEE->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEE->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEE->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEE->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEE->Branch("fit_b0",&fit_b0,"fit_b0/F");
  treeEE->Branch("fit_b1",&fit_b1,"fit_b1/F");
  treeEE->Branch("fit_b2",&fit_b2,"fit_b2/F");
  treeEE->Branch("fit_b3",&fit_b3,"fit_b3/F");
  treeEE->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");


  for(int iR=0; iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); ++iR)  {
    std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(iR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
      EBDetId ebid(*iid);
      hashedIndex = ebid.hashedIndex();
      ieta = ebid.ieta();
      iphi = ebid.iphi();
      iSM = ebid.ism();
      iMod = ebid.im();
      iTT  = ebid.tower().hashedIndex();
      iTTeta = ebid.tower_ieta();
      iTTphi = ebid.tower_iphi();
      Signal = EBmap_Signal[ebid.hashedIndex()];//#
      Backgr = EBmap_Backgr[ebid.hashedIndex()];
      Chisqu = EBmap_Chisqu[ebid.hashedIndex()];
      Ndof = EBmap_ndof[ebid.hashedIndex()];
      fit_mean     = EBmap_mean[ebid.hashedIndex()];
      fit_mean_err = EBmap_mean_err[ebid.hashedIndex()];
      fit_sigma  = EBmap_sigma[ebid.hashedIndex()];
      fit_Snorm  = EBmap_Snorm[ebid.hashedIndex()];
      fit_b0     = EBmap_b0[ebid.hashedIndex()];
      fit_b1     = EBmap_b1[ebid.hashedIndex()];
      fit_b2     = EBmap_b2[ebid.hashedIndex()];
      fit_b3     = EBmap_b3[ebid.hashedIndex()];
      fit_Bnorm  = EBmap_Bnorm[ebid.hashedIndex()];

      regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);

      treeEB->Fill();
    } // loop over DetId in regions
  } // loop over regions

  for(int jR=0; jR < regionalCalibration_->getCalibMap()->getNRegionsEE() ; jR++)
    {
      std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	{ 
	  EEDetId eeid(*iid);
	  ix = eeid.ix();
	  iy = eeid.iy();
	  zside = eeid.zside();
	  sc = eeid.sc();
	  isc = eeid.isc();
	  ic = eeid.ic();
	  iquadrant = eeid.iquadrant();
	  hashedIndex = eeid.hashedIndex();
	  regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);
	  Signal = EEmap_Signal[eeid.hashedIndex()];//#
	  Backgr = EEmap_Backgr[eeid.hashedIndex()];
	  Chisqu = EEmap_Chisqu[eeid.hashedIndex()];            
	  Ndof = EEmap_ndof[eeid.hashedIndex()];            
	  fit_mean     = EEmap_mean[eeid.hashedIndex()];
	  fit_mean_err = EEmap_mean_err[eeid.hashedIndex()];
	  fit_sigma  = EEmap_sigma[eeid.hashedIndex()];
	  fit_Snorm  = EEmap_Snorm[eeid.hashedIndex()];
	  fit_b0     = EEmap_b0[eeid.hashedIndex()];
	  fit_b1     = EEmap_b1[eeid.hashedIndex()];
	  fit_b2     = EEmap_b2[eeid.hashedIndex()];
	  fit_b3     = EEmap_b3[eeid.hashedIndex()];
	  fit_Bnorm  = EEmap_Bnorm[eeid.hashedIndex()];

	  treeEE->Fill();
	}
    }

  treeEB->Write();
  treeEE->Write();

  outfile_->Write();
  outfile_->Close();
  cout << "FIT_EPSILON:  done" << endl;

}

//==========================

void FitEpsilonPlot::saveCoefficientsEoverEtrue(const bool isSecondGenPhoton = false) 
{

  // important, if using the second photon the output file is updated, so the call with isSecondGenPhoton = true should be made as the second one
  // otherwise, based on the current implementation, at the time you open the file for the first photon the file would be overwritten due to RECREATE mode

  /// output file
  std::string fileName = outputDir_  + "/" + outfilename_;
  if (isSecondGenPhoton) outfile_ = new TFile(fileName.c_str(),"UPDATE");
  else                   outfile_ = new TFile(fileName.c_str(),"RECREATE");
  cout << "FIT_EPSILON: Saving E/Etrue Coefficients in " << fileName << " ... " << endl;;
  if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << fileName << "\n";
  outfile_->cd();

  // 2D calib map in the barrel
  TH2F* hmap_EB = new TH2F((isSecondGenPhoton ? "calibMap_EB_g2" : "calibMap_EB"),"EB calib coefficients: #eta on x, #phi on y",
			   2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,
			   EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
  hmap_EB->GetXaxis()->SetTitle("i#eta");
  hmap_EB->GetYaxis()->SetTitle("i#phi");
  TH2F* hmap_EEp = new TH2F((isSecondGenPhoton ? "calibMap_EEp_g2" : "calibMap_EEp"),"EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEp->GetXaxis()->SetTitle("ix");
  hmap_EEp->GetYaxis()->SetTitle("iy");
  TH2F* hmap_EEm = new TH2F((isSecondGenPhoton ? "calibMap_EEm_g2" : "calibMap_EEm"),"EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEm->GetXaxis()->SetTitle("ix");
  hmap_EEm->GetYaxis()->SetTitle("iy");
  TH1F* hint = new TH1F("hint","Bin1: inRangeFit_ Bin2: finRangeFit_ Bin3: Barrel(0)/Endcap(1)",3,0.,3.);
  hint->SetBinContent(1,inRangeFit_);
  hint->SetBinContent(2,finRangeFit_);
  if( EEoEB_ == "Barrel" ) hint->SetBinContent(3,0);
  else                     hint->SetBinContent(3,1);
  hint->Write();

  EcalRegionalCalibrationBase* regCalibToUse = (isSecondGenPhoton) ? regionalCalibration_g2_ : regionalCalibration_;
  std::map<int,TFitResultPtr>& EBmap_fitresptrToUse = (isSecondGenPhoton) ? EBmap_fitresptr_g2 : EBmap_fitresptr_g1;
  std::map<int,TFitResultPtr>& EEmap_fitresptrToUse = (isSecondGenPhoton) ? EEmap_fitresptr_g2 : EEmap_fitresptr_g1;

  //filling Barrel Map
  if (foldInSuperModule_) {

    // in this case we only used the index given by ic in a single SM
    // we do not use regCalibToUse->allDetIdsInEBRegion(), we assume we are making xtals

    for(int j = inRangeFit_; j <= finRangeFit_; ++j)  
      {      

	// in this configuration j is interpreted as the xtal number given by EBDetId::ic() (be careful, make sure this is how the number is used in the previous methods)
	// WARNING: j starts from 0, while ic() is always >= 1: need to sum 1 to interpret j as ic() 
	// go to EB+1 from ic(): iphi,ieta = 1,1 has ic = 20 (for EB- iphi = 1 has ic = 1) 
	// we use EBDetId constructor with SM number, xtal number ( i.e. ic() ) and mode = EBDetId::SMCRYSTALMODE == 1
	// SM number is 1 for EB+1 (up to 18), and 19 for EB-1 (up to EB-18 which has iSM = 36)

	EBDetId ebid(1,j+1,1);
	float coeffValue = regCalibToUse->getCalibMap()->coeff(ebid) > 0. ? regCalibToUse->getCalibMap()->coeff(ebid) : 1.;	      
	int ix = ebid.ieta() + EBDetId::MAX_IETA+1;
	hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	// now fill all other SM in the same way
	for (int ism = 2; ism <= 36; ++ism) {
	  EBDetId ebid(ism,j+1,1);
	  int ix = ebid.ieta() + EBDetId::MAX_IETA+1;
	  hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	}
      }

  } else {

    for(int j=0; j<regCalibToUse->getCalibMap()->getNRegionsEB(); ++j)  
      {
	std::vector<DetId> ids = regCalibToUse->allDetIdsInEBRegion(j);
	for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
	  EBDetId ebid(*iid);
	  int ix = ebid.ieta()+EBDetId::MAX_IETA+1;

	  float coeffValue = regCalibToUse->getCalibMap()->coeff(*iid) > 0. ? regCalibToUse->getCalibMap()->coeff(*iid) : 1.;
	  hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	} // loop over DetId in regions
      }

  }

  hmap_EB->SetMinimum(0.9);
  hmap_EB->SetStats(false);
  hmap_EB->Write();

  for(int jR=0; jR < regCalibToUse->getCalibMap()->getNRegionsEE(); jR++)
    {
      std::vector<DetId> ids = regCalibToUse->allDetIdsInEERegion(jR);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
  	{ 
  	  EEDetId eeid(*iid);
  	  float coeffValue =  regCalibToUse->getCalibMap()->coeff(*iid) > 0. ?  regCalibToUse->getCalibMap()->coeff(*iid) : 1.;

  	  if(eeid.positiveZ())
  	    hmap_EEp->Fill(eeid.ix(), eeid.iy(), coeffValue); 
  	  else 
  	    hmap_EEm->Fill(eeid.ix(), eeid.iy(), coeffValue);
  	}
    }

  hmap_EEp->SetMinimum(0.9);
  hmap_EEp->SetStats(false);
  hmap_EEp->Write();

  hmap_EEm->SetMinimum(0.9);
  hmap_EEm->SetStats(false);
  hmap_EEm->Write();

  /*------------- TTREE --------------*/
  uint32_t   rawId;
  int        hashedIndex;
  int        ieta;
  int        iphi;
  int        iSM;
  int        iMod;
  int        iTT;
  int        iTTeta;
  int        iTTphi;
  int        iter = currentIteration_;
  float      regCoeff;
  // float      Signal;//#
  // float      Backgr; 
  float      Chisqu; 
  float      Ndof; 
  float      fit_mean;
  float      fit_mean_err;
  float      fit_sigma;
  // float      fit_Snorm;
  // float      fit_b0;
  // float      fit_b1;    
  // float      fit_b2;    
  // float      fit_b3;    
  // float      fit_Bnorm; 
  /// endcap variables
  int ix;
  int iy;
  int zside;
  int sc; 
  int isc;
  int ic;
  int iquadrant;

  TTree* treeEB = new TTree((isSecondGenPhoton ? "calibEB_g2" : "calibEB"),"Tree of EB Inter-calibration constants");
  TTree* treeEE = new TTree((isSecondGenPhoton ? "calibEE_g2" : "calibEE"),"Tree of EE Inter-calibration constants");

  /// barrel
  treeEB->Branch("rawId",&rawId,"rawId/i");
  treeEB->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEB->Branch("ieta",&ieta,"ieta/I");
  treeEB->Branch("iphi",&iphi,"iphi/I");
  treeEB->Branch("iSM",&iSM,"iSM/I");
  treeEB->Branch("iMod",&iMod,"iMod/I");
  treeEB->Branch("iTT",&iTT,"iTT/I");
  treeEB->Branch("iTTeta",&iTTeta,"iTTeta/I");
  treeEB->Branch("iTTphi",&iTTphi,"iTTphi/I");
  treeEB->Branch("iter",&iter,"iter/I");
  treeEB->Branch("coeff",&regCoeff,"coeff/F");
  // treeEB->Branch("Signal",&Signal,"Signal/F");//#
  // treeEB->Branch("Backgr",&Backgr,"Backgr/F");
  treeEB->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEB->Branch("Ndof",&Ndof,"Ndof/F");
  treeEB->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEB->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEB->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  // treeEB->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  // treeEB->Branch("fit_b0",&fit_b0,"fit_b0/F");
  // treeEB->Branch("fit_b1",&fit_b1,"fit_b1/F");
  // treeEB->Branch("fit_b2",&fit_b2,"fit_b2/F");
  // treeEB->Branch("fit_b3",&fit_b3,"fit_b3/F");
  // treeEB->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

  /// endcap
  treeEE->Branch("ix",&ix,"ix/I");
  treeEE->Branch("iy",&iy,"iy/I");
  treeEE->Branch("zside",&zside,"zside/I");
  treeEE->Branch("sc",&sc,"sc/I");
  treeEE->Branch("isc",&isc,"isc/I");
  treeEE->Branch("ic",&ic,"ic/I");
  treeEE->Branch("iquadrant",&iquadrant,"iquadrant/I");
  treeEE->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEE->Branch("iter",&iter,"iter/I");
  treeEE->Branch("coeff",&regCoeff,"coeff/F");
  // treeEE->Branch("Signal",&Signal,"Signal/F");//#
  // treeEE->Branch("Backgr",&Backgr,"Backgr/F");
  treeEE->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEE->Branch("Ndof",&Ndof,"Ndof/F");
  treeEE->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEE->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEE->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  // treeEE->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  // treeEE->Branch("fit_b0",&fit_b0,"fit_b0/F");
  // treeEE->Branch("fit_b1",&fit_b1,"fit_b1/F");
  // treeEE->Branch("fit_b2",&fit_b2,"fit_b2/F");
  // treeEE->Branch("fit_b3",&fit_b3,"fit_b3/F");
  // treeEE->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

  for(int iR=0; iR < regCalibToUse->getCalibMap()->getNRegionsEB(); ++iR)  {

    std::vector<DetId> ids = regCalibToUse->allDetIdsInEBRegion(iR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {

      EBDetId ebid(*iid);
      hashedIndex = ebid.hashedIndex();
      ieta = ebid.ieta();
      iphi = ebid.iphi();
      iSM = ebid.ism();
      iMod = ebid.im();
      iTT  = ebid.tower().hashedIndex();
      iTTeta = ebid.tower_ieta();
      iTTphi = ebid.tower_iphi();
      // // Signal = EBmap_Signal[ebid.hashedIndex()];//#
      // // Backgr = EBmap_Backgr[ebid.hashedIndex()];
      // Chisqu = EBmap_Chisqu[ebid.hashedIndex()];
      // Ndof = EBmap_ndof[ebid.hashedIndex()];
      // fit_mean     = EBmap_mean[ebid.hashedIndex()];
      // fit_mean_err = EBmap_mean_err[ebid.hashedIndex()];
      // fit_sigma  = EBmap_sigma[ebid.hashedIndex()];
      // // fit_Snorm  = EBmap_Snorm[ebid.hashedIndex()];
      // // fit_b0     = EBmap_b0[ebid.hashedIndex()];
      // // fit_b1     = EBmap_b1[ebid.hashedIndex()];
      // // fit_b2     = EBmap_b2[ebid.hashedIndex()];
      // // fit_b3     = EBmap_b3[ebid.hashedIndex()];
      // // fit_Bnorm  = EBmap_Bnorm[ebid.hashedIndex()];

      if (EBmap_fitresptrToUse[ebid.hashedIndex()] >= 0 && EBmap_fitresptrToUse[ebid.hashedIndex()].Get() != nullptr) {
	Chisqu       = EBmap_fitresptrToUse[ebid.hashedIndex()]->Chi2();
	Ndof         = EBmap_fitresptrToUse[ebid.hashedIndex()]->Ndf();
	fit_mean     = EBmap_fitresptrToUse[ebid.hashedIndex()]->Parameter(1);  // for the double CB the mean is parameter [1] (as for the gaussian)
	fit_mean_err = EBmap_fitresptrToUse[ebid.hashedIndex()]->ParError(1);
	fit_sigma    = EBmap_fitresptrToUse[ebid.hashedIndex()]->Parameter(2);
      } else {
	Chisqu       = -999;
	Ndof         = -999;
	fit_mean     = -999;
	fit_mean_err = -999;
	fit_sigma    = -999;
      }

      regCoeff = regCalibToUse->getCalibMap()->coeff(*iid);

      treeEB->Fill();
    } // loop over DetId in regions
  } // loop over regions

  for(int jR=0; jR < regCalibToUse->getCalibMap()->getNRegionsEE() ; jR++) {

    std::vector<DetId> ids = regCalibToUse->allDetIdsInEERegion(jR);

    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) { 

      EEDetId eeid(*iid);
      ix = eeid.ix();
      iy = eeid.iy();
      zside = eeid.zside();
      sc = eeid.sc();
      isc = eeid.isc();
      ic = eeid.ic();
      iquadrant = eeid.iquadrant();
      hashedIndex = eeid.hashedIndex();
      regCoeff = regCalibToUse->getCalibMap()->coeff(*iid);
      // // Signal = EEmap_Signal[eeid.hashedIndex()];//#
      // // Backgr = EEmap_Backgr[eeid.hashedIndex()];
      // Chisqu = EEmap_Chisqu[eeid.hashedIndex()];            
      // Ndof = EEmap_ndof[eeid.hashedIndex()];            
      // fit_mean     = EEmap_mean[eeid.hashedIndex()];
      // fit_mean_err = EEmap_mean_err[eeid.hashedIndex()];
      // fit_sigma  = EEmap_sigma[eeid.hashedIndex()];
      // // fit_Snorm  = EEmap_Snorm[eeid.hashedIndex()];
      // // fit_b0     = EEmap_b0[eeid.hashedIndex()];
      // // fit_b1     = EEmap_b1[eeid.hashedIndex()];
      // // fit_b2     = EEmap_b2[eeid.hashedIndex()];
      // // fit_b3     = EEmap_b3[eeid.hashedIndex()];
      // // fit_Bnorm  = EEmap_Bnorm[eeid.hashedIndex()];

      if (EEmap_fitresptrToUse[eeid.hashedIndex()] >= 0 && EEmap_fitresptrToUse[eeid.hashedIndex()].Get() != nullptr) {
	Chisqu       = EEmap_fitresptrToUse[eeid.hashedIndex()]->Chi2();
	Ndof         = EEmap_fitresptrToUse[eeid.hashedIndex()]->Ndf();
	fit_mean     = EEmap_fitresptrToUse[eeid.hashedIndex()]->Parameter(1);  // for the double CB the mean is parameter [1]
	fit_mean_err = EEmap_fitresptrToUse[eeid.hashedIndex()]->ParError(1);
	fit_sigma    = EEmap_fitresptrToUse[eeid.hashedIndex()]->Parameter(2);
      } else {
	Chisqu       = -999;
	Ndof         = -999;
	fit_mean     = -999;
	fit_mean_err = -999;
	fit_sigma    = -999;
      }

      treeEE->Fill();

    }

  }

  treeEB->Write();
  treeEE->Write();

  outfile_->Write();
  outfile_->Close();
  cout << "FIT_EPSILON:  done" << endl;

}

//==============================================================

//==========================

void FitEpsilonPlot::saveCoefficientsEoverEtrueRooFit(const bool isSecondGenPhoton = false) 
{

  // when saving the coefficients, if we are folding in SM, we need to fill just one SM and copy on all the other
  // so, from denseIndex we go to ic() (beware, it could be ic() was not used for the folding, so there could be an inconsistency)
  // then we must get iphi and ieta in SM (ieta in 1-85 and iphi in 1-20)

  // important, if using the second photon the output file is updated, so the call with isSecondGenPhoton = true should be made as the second one
  // otherwise, based on the current implementation, at the time you open the file for the first photon the file would be overwritten due to RECREATE mode

  /// output file
  std::string fileName = outputDir_  + "/" + outfilename_;
  if (isSecondGenPhoton) outfile_ = new TFile(fileName.c_str(),"UPDATE");
  else                   outfile_ = new TFile(fileName.c_str(),"RECREATE");
  cout << "FIT_EPSILON: Saving E/Etrue Coefficients in " << fileName << " ... " << endl;;
  if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << fileName << "\n";
  outfile_->cd();

  // 2D calib map in the barrel
  TH2F* hmap_EB = new TH2F((isSecondGenPhoton ? "calibMap_EB_g2" : "calibMap_EB"),"EB calib coefficients: #eta on x, #phi on y",
			   2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,
			   EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
  hmap_EB->GetXaxis()->SetTitle("i#eta");
  hmap_EB->GetYaxis()->SetTitle("i#phi");
  TH2F* hmap_EEp = new TH2F((isSecondGenPhoton ? "calibMap_EEp_g2" : "calibMap_EEp"),"EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEp->GetXaxis()->SetTitle("ix");
  hmap_EEp->GetYaxis()->SetTitle("iy");
  TH2F* hmap_EEm = new TH2F((isSecondGenPhoton ? "calibMap_EEm_g2" : "calibMap_EEm"),"EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEm->GetXaxis()->SetTitle("ix");
  hmap_EEm->GetYaxis()->SetTitle("iy");
  TH1F* hint = new TH1F("hint","Bin1: inRangeFit_ Bin2: finRangeFit_ Bin3: Barrel(0)/Endcap(1)",3,0.,3.);
  hint->SetBinContent(1,inRangeFit_);
  hint->SetBinContent(2,finRangeFit_);
  if( EEoEB_ == "Barrel" ) hint->SetBinContent(3,0);
  else                     hint->SetBinContent(3,1);
  hint->Write();

  EcalRegionalCalibrationBase* regCalibToUse = (isSecondGenPhoton) ? regionalCalibration_g2_ : regionalCalibration_;

  //filling Barrel Map
  if (foldInSuperModule_) {

    // in this case we only used the index given by ic in a single SM
    // we do not use regCalibToUse->allDetIdsInEBRegion(), we assume we are making xtals

    for(int j = inRangeFit_; j <= finRangeFit_; ++j)  
      {      

	// in this configuration j is interpreted as the xtal number given by EBDetId::ic() (be careful, make sure this is how the number is used in the previous methods)
	// WARNING: j starts from 0, while ic() is always >= 1: need to sum 1 to interpret j as ic() 
	// go to EB+1 from ic(): iphi,ieta = 1,1 has ic = 20 (for EB- iphi = 1 has ic = 1) 
	// we use EBDetId constructor with SM number, xtal number ( i.e. ic() ) and mode = EBDetId::SMCRYSTALMODE == 1
	// SM number is 1 for EB+1 (up to 18), and 19 for EB-1 (up to EB-18 which has iSM = 36)

	EBDetId ebid(1,j+1,1);
	float coeffValue = regCalibToUse->getCalibMap()->coeff(ebid) > 0. ? regCalibToUse->getCalibMap()->coeff(ebid) : 1.;	      
	int ix = ebid.ieta() + EBDetId::MAX_IETA+1;
	hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	// now fill all other SM in the same way
	for (int ism = 2; ism <= 36; ++ism) {
	  EBDetId ebid(ism,j+1,1);
	  int ix = ebid.ieta() + EBDetId::MAX_IETA+1;
	  hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	}
      }

  } else {

    //filling Barrel Map
    for(int j=0; j<regCalibToUse->getCalibMap()->getNRegionsEB(); ++j)  
      {
      	std::vector<DetId> ids = regCalibToUse->allDetIdsInEBRegion(j);

	for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {

	  EBDetId ebid(*iid);
	  float coeffValue = regCalibToUse->getCalibMap()->coeff(*iid) > 0. ? regCalibToUse->getCalibMap()->coeff(*iid) : 1.;
	  int ix = ebid.ieta()+EBDetId::MAX_IETA+1;
	  hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );

	} // loop over DetId in regions

      }

  }

  hmap_EB->SetMinimum(0.9);
  hmap_EB->SetStats(false);
  hmap_EB->Write();

  for(int jR=0; jR < regCalibToUse->getCalibMap()->getNRegionsEE(); jR++)
    {
      std::vector<DetId> ids = regCalibToUse->allDetIdsInEERegion(jR);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
  	{ 
  	  EEDetId eeid(*iid);
  	  float coeffValue =  regCalibToUse->getCalibMap()->coeff(*iid) > 0. ?  regCalibToUse->getCalibMap()->coeff(*iid) : 1.;

  	  if(eeid.positiveZ())
  	    hmap_EEp->Fill(eeid.ix(), eeid.iy(), coeffValue); 
  	  else 
  	    hmap_EEm->Fill(eeid.ix(), eeid.iy(), coeffValue);
  	}
    }

  hmap_EEp->SetMinimum(0.9);
  hmap_EEp->SetStats(false);
  hmap_EEp->Write();

  hmap_EEm->SetMinimum(0.9);
  hmap_EEm->SetStats(false);
  hmap_EEm->Write();

  /*------------- TTREE --------------*/
  uint32_t   rawId;
  int        hashedIndex;
  int        ieta;
  int        iphi;
  int        iSM;
  int        iMod;
  int        iTT;
  int        iTTeta;
  int        iTTphi;
  int        iter = currentIteration_;
  float      regCoeff;
  float      Signal;//#
  float      Backgr; 
  float      Chisqu; 
  float      Ndof; 
  float      fit_mean;
  float      fit_mean_err;
  float      fit_sigma;
  float      fit_Snorm;
  float      fit_b0;
  float      fit_b1;    
  float      fit_b2;    
  float      fit_b3;    
  float      fit_Bnorm; 
  /// endcap variables
  int ix;
  int iy;
  int zside;
  int sc; 
  int isc;
  int ic;
  int iquadrant;

  TTree* treeEB = new TTree((isSecondGenPhoton ? "calibEB_g2" : "calibEB"),"Tree of EB Inter-calibration constants");
  TTree* treeEE = new TTree((isSecondGenPhoton ? "calibEE_g2" : "calibEE"),"Tree of EE Inter-calibration constants");

  /// barrel
  treeEB->Branch("rawId",&rawId,"rawId/i");
  treeEB->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEB->Branch("ieta",&ieta,"ieta/I");
  treeEB->Branch("iphi",&iphi,"iphi/I");
  treeEB->Branch("iSM",&iSM,"iSM/I");
  treeEB->Branch("iMod",&iMod,"iMod/I");
  treeEB->Branch("iTT",&iTT,"iTT/I");
  treeEB->Branch("iTTeta",&iTTeta,"iTTeta/I");
  treeEB->Branch("iTTphi",&iTTphi,"iTTphi/I");
  treeEB->Branch("iter",&iter,"iter/I");
  treeEB->Branch("coeff",&regCoeff,"coeff/F");
  treeEB->Branch("Signal",&Signal,"Signal/F");//#
  treeEB->Branch("Backgr",&Backgr,"Backgr/F");
  treeEB->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEB->Branch("Ndof",&Ndof,"Ndof/F");
  treeEB->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEB->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEB->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEB->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEB->Branch("fit_b0",&fit_b0,"fit_b0/F");
  treeEB->Branch("fit_b1",&fit_b1,"fit_b1/F");
  treeEB->Branch("fit_b2",&fit_b2,"fit_b2/F");
  treeEB->Branch("fit_b3",&fit_b3,"fit_b3/F");
  treeEB->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

  /// endcap
  treeEE->Branch("ix",&ix,"ix/I");
  treeEE->Branch("iy",&iy,"iy/I");
  treeEE->Branch("zside",&zside,"zside/I");
  treeEE->Branch("sc",&sc,"sc/I");
  treeEE->Branch("isc",&isc,"isc/I");
  treeEE->Branch("ic",&ic,"ic/I");
  treeEE->Branch("iquadrant",&iquadrant,"iquadrant/I");
  treeEE->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEE->Branch("iter",&iter,"iter/I");
  treeEE->Branch("coeff",&regCoeff,"coeff/F");
  treeEE->Branch("Signal",&Signal,"Signal/F");//#
  treeEE->Branch("Backgr",&Backgr,"Backgr/F");
  treeEE->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEE->Branch("Ndof",&Ndof,"Ndof/F");
  treeEE->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEE->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEE->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEE->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEE->Branch("fit_b0",&fit_b0,"fit_b0/F");
  treeEE->Branch("fit_b1",&fit_b1,"fit_b1/F");
  treeEE->Branch("fit_b2",&fit_b2,"fit_b2/F");
  treeEE->Branch("fit_b3",&fit_b3,"fit_b3/F");
  treeEE->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

  for(int iR=0; iR < regCalibToUse->getCalibMap()->getNRegionsEB(); ++iR)  {

    std::vector<DetId> ids = regCalibToUse->allDetIdsInEBRegion(iR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {

      EBDetId ebid(*iid);
      hashedIndex = ebid.hashedIndex();
      ieta = ebid.ieta();
      iphi = ebid.iphi();
      iSM = ebid.ism();
      iMod = ebid.im();
      iTT  = ebid.tower().hashedIndex();
      iTTeta = ebid.tower_ieta();
      iTTphi = ebid.tower_iphi();
      if (isSecondGenPhoton) {
	Signal = EBmap_Signal_g2[ebid.hashedIndex()];//#
	Backgr = EBmap_Backgr_g2[ebid.hashedIndex()];
	Chisqu = EBmap_Chisqu_g2[ebid.hashedIndex()];
	Ndof = EBmap_ndof_g2[ebid.hashedIndex()];
	fit_mean     = EBmap_mean_g2[ebid.hashedIndex()];
	fit_mean_err = EBmap_mean_err_g2[ebid.hashedIndex()];
	fit_sigma  = EBmap_sigma_g2[ebid.hashedIndex()];
	fit_Snorm  = EBmap_Snorm_g2[ebid.hashedIndex()];
	fit_b0     = EBmap_b0_g2[ebid.hashedIndex()];
	fit_b1     = EBmap_b1_g2[ebid.hashedIndex()];
	fit_b2     = EBmap_b2_g2[ebid.hashedIndex()];
	fit_b3     = EBmap_b3_g2[ebid.hashedIndex()];
	fit_Bnorm  = EBmap_Bnorm_g2[ebid.hashedIndex()];
      } else {
	Signal = EBmap_Signal[ebid.hashedIndex()];//#
	Backgr = EBmap_Backgr[ebid.hashedIndex()];
	Chisqu = EBmap_Chisqu[ebid.hashedIndex()];
	Ndof = EBmap_ndof[ebid.hashedIndex()];
	fit_mean     = EBmap_mean[ebid.hashedIndex()];
	fit_mean_err = EBmap_mean_err[ebid.hashedIndex()];
	fit_sigma  = EBmap_sigma[ebid.hashedIndex()];
	fit_Snorm  = EBmap_Snorm[ebid.hashedIndex()];
	fit_b0     = EBmap_b0[ebid.hashedIndex()];
	fit_b1     = EBmap_b1[ebid.hashedIndex()];
	fit_b2     = EBmap_b2[ebid.hashedIndex()];
	fit_b3     = EBmap_b3[ebid.hashedIndex()];
	fit_Bnorm  = EBmap_Bnorm[ebid.hashedIndex()];
      }

      regCoeff = regCalibToUse->getCalibMap()->coeff(*iid);

      treeEB->Fill();
    } // loop over DetId in regions
  } // loop over regions

  for(int jR=0; jR < regCalibToUse->getCalibMap()->getNRegionsEE() ; jR++) {

    std::vector<DetId> ids = regCalibToUse->allDetIdsInEERegion(jR);

    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) { 

      EEDetId eeid(*iid);
      ix = eeid.ix();
      iy = eeid.iy();
      zside = eeid.zside();
      sc = eeid.sc();
      isc = eeid.isc();
      ic = eeid.ic();
      iquadrant = eeid.iquadrant();
      hashedIndex = eeid.hashedIndex();
      regCoeff = regCalibToUse->getCalibMap()->coeff(*iid);
      if (isSecondGenPhoton) {
	Signal = EEmap_Signal_g2[eeid.hashedIndex()];//#
	Backgr = EEmap_Backgr_g2[eeid.hashedIndex()];
	Chisqu = EEmap_Chisqu_g2[eeid.hashedIndex()];            
	Ndof = EEmap_ndof_g2[eeid.hashedIndex()];            
	fit_mean     = EEmap_mean_g2[eeid.hashedIndex()];
	fit_mean_err = EEmap_mean_err_g2[eeid.hashedIndex()];
	fit_sigma  = EEmap_sigma_g2[eeid.hashedIndex()];
	fit_Snorm  = EEmap_Snorm_g2[eeid.hashedIndex()];
	fit_b0     = EEmap_b0_g2[eeid.hashedIndex()];
	fit_b1     = EEmap_b1_g2[eeid.hashedIndex()];
	fit_b2     = EEmap_b2_g2[eeid.hashedIndex()];
	fit_b3     = EEmap_b3_g2[eeid.hashedIndex()];
	fit_Bnorm  = EEmap_Bnorm_g2[eeid.hashedIndex()];
      }

      treeEE->Fill();

    }

  }

  treeEB->Write();
  treeEE->Write();

  outfile_->Write();
  outfile_->Close();
  cout << "FIT_EPSILON:  done" << endl;

}




// ------------ method called for each event  ------------

void FitEpsilonPlot::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // if we only wanted to fold histograms, like the first time we call this whole code, we don't need to fit anything here
    if (foldInSuperModule_ && makeFoldedHistograms_) {
      cout << "FIT_EPSILON: not doing anything inside analyze(): we only wanted to fold histograms" << endl;
      return;
    }

    TF1 ffit("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
    ffit.SetParameters(100,0,0.1);
    ffit.SetParNames("Constant","Mean_value","Sigma","a","b");

    ffit.SetParLimits(3,-500,500);
    ffit.SetParLimits(2,0.05,0.22);

    cout << "FIT_EPSILON: About to fit epsilon distributions" << endl; 

    /// compute average weight, eps, and update calib constant
    if( (EEoEB_ == "Barrel") && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ){

      for(uint32_t j= (uint32_t)inRangeFit_; j <= (uint32_t)finRangeFit_ && j < (uint32_t)regionalCalibration_->getCalibMap()->getNRegionsEB(); ++j)  
	{
	  cout<<"FIT_EPSILON: Fitting EB Cristal--> "<<j<<endl;

	  if(!(j%1000)) cout << "FIT_EPSILON: fitting EB region " << j << endl;

	  float mean = 0.;
	  float mean_g2 = 0.; // used only for E/Etrue with MC	      

	  if (isEoverEtrue_) {
		  
	    int crystalIndexInSM = foldInSuperModule_ ? j : getArrayIndexOfFoldedSMfromDenseIndex(j);
	    TH1F* histoToFit_g1 = (foldInSuperModule_ ? EoverEtrue_g1_EB_SM_hvec[crystalIndexInSM] : EoverEtrue_g1_EB_h[j]);
	    // first photon 
	    // int iMin = EoverEtrue_g1_EB_h[j]->GetXaxis()->FindFixBin(0.6); 
	    // int iMax = EoverEtrue_g1_EB_h[j]->GetXaxis()->FindFixBin(1.1);
	    // double integral = EoverEtrue_g1_EB_h[j]->Integral(iMin, iMax);  
	    double integral = histoToFit_g1->Integral();  

	    if(integral > EoverEtrue_integralMin) {

	      if (fitEoverEtrueWithRooFit_) {
		Pi0FitResult fitres = FitEoverEtruePeakRooFit(histoToFit_g1, false, j, Pi0EB);
		RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
		mean = mean_fitresult->getVal();
	      } else {
		TFitResultPtr fitresptr = FitEoverEtruePeak(histoToFit_g1, false, j, Pi0EB, false);
		mean = fitresptr->Parameter(1);
		if (mean >= 1.5) mean = 0.; 		
	      }

	    } else {

	      std::cout << "### g1 ### FIT_EPSILON: iR = " << j << ", integral() = " << integral << " , skipping the fit " << std::endl;
	      mean = 0.;
	      if (not fitEoverEtrueWithRooFit_) {
		EBmap_fitresptr_g1[j] = TFitResultPtr(-1); 
	      }

	    }

	    TH1F* histoToFit_g2 = (foldInSuperModule_ ? EoverEtrue_g2_EB_SM_hvec[crystalIndexInSM] : EoverEtrue_g2_EB_h[j]);
	    if (foldInSuperModule_ and (crystalIndexInSM == 1336 || crystalIndexInSM == 1611 || crystalIndexInSM == 1625)) {
	      histoToFit_g2 = EoverEtrue_g2_EB_SM_hvec[crystalIndexInSM-1];
	    }

	    // second photon 
	    // iMin = EoverEtrue_g2_EB_h[j]->GetXaxis()->FindFixBin(0.6); 
	    // iMax = EoverEtrue_g2_EB_h[j]->GetXaxis()->FindFixBin(1.1);
	    // integral = EoverEtrue_g2_EB_h[j]->Integral(iMin, iMax);  
	    integral = histoToFit_g2->Integral();  

	    if(integral > EoverEtrue_integralMin) {

	      if (fitEoverEtrueWithRooFit_) {
		Pi0FitResult fitres = FitEoverEtruePeakRooFit(histoToFit_g2, true, j, Pi0EB);
		RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
		mean_g2 = mean_fitresult->getVal();
	      } else {
		TFitResultPtr fitresptr = FitEoverEtruePeak(histoToFit_g2, true, j, Pi0EB, false);
		mean_g2 = fitresptr->Parameter(1);
		if (mean_g2 >= 1.5) mean_g2 = 0.; 
	      }
		    
	    } else {

	      std::cout << "### g2 ### FIT_EPSILON: iR = " << j << ", integral() = " << integral << " , skipping the fit " << std::endl;
	      mean_g2 = 0.;
	      if (not fitEoverEtrueWithRooFit_) {
		EBmap_fitresptr_g2[j] = TFitResultPtr(-1);
	      }

	    }
		  
	  } else {
		  
	    if(!useMassInsteadOfEpsilon_ && epsilon_EB_h[j]->Integral(epsilon_EB_h[j]->GetNbinsX()*(1./6.),epsilon_EB_h[j]->GetNbinsX()*0.5) > 20) 
	      {

		double Max = 0.;
		double Min = -0.5, bin = 0.0125;
		Max = Min+(bin*(double)epsilon_EB_h[j]->GetMaximumBin());
		double Bound1 = -0.15, Bound2 = 0.25;
		if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.1;}
		if ( Max+Bound2 > 0.34  ){ Bound2 = 0.15;}
		if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.075;}
		if ( Max+Bound2 > 0.34  ){ Bound2 = 0.1;}
		if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.03;}
		if ( Max+Bound2 > 0.34  ){ Bound2 = 0.05;}
		if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.009;}
		if ( Max+Bound2 > 0.34  ){ Bound2 = 0.01;}

		epsilon_EB_h[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
		if(ffit.GetNDF() != 0) {
		  double chi2 = ( ffit.GetChisquare()/ffit.GetNDF() );

		  if ( chi2  > 11 ){
		    ffit.SetParLimits(2,0.05,0.15);
		    ffit.SetParameters(100,0,0.1);
		    epsilon_EB_h[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
		    chi2 = (ffit.GetChisquare()/ffit.GetNDF());
		    if ( chi2  < 11 ){   cout<<"Saved 1 Level!!"<<endl;  }
		    else{
		      ffit.SetParameters(100,0,0.1);
		      ffit.SetParLimits(2,0.05,0.1);
		      epsilon_EB_h[j]->Fit(&ffit,"qB","",  Max+Bound1,Max+Bound2);
		      chi2 = (ffit.GetChisquare()/ffit.GetNDF());
		      if ( chi2  < 11 ){ cout<<"Saved 2 Level!!"<<endl; }
		      else{ cout<<"DAMN: High Chi square..."<<endl; }
		    }
		  }
		}
		else cout<<"DAMN: NDF == 0"<<endl;
		mean = ffit.GetParameter(1);
	      }
	    else if(useMassInsteadOfEpsilon_)
	      {

		int crystalIndexInSM = getArrayIndexOfFoldedSMfromDenseIndex(j);
		TH1F* histoToFit = (foldInSuperModule_ ? epsilon_EB_SM_hvec[crystalIndexInSM] : epsilon_EB_h[j]);

		int iMin = histoToFit->GetXaxis()->FindFixBin(Are_pi0_? 0.08:0.4 ); 
		int iMax = histoToFit->GetXaxis()->FindFixBin(Are_pi0_? 0.18:0.65 );
		double integral = histoToFit->Integral(iMin, iMax);  

		if(integral>60.) {

		  Pi0FitResult fitres = FitMassPeakRooFit( histoToFit, 
							   Are_pi0_? fitRange_low_pi0:fitRange_low_eta, 
							   Are_pi0_? fitRange_high_pi0:fitRange_high_eta, 
							   j, 1, Pi0EB, 0, isNot_2010_); //0.05-0.3
		  RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
		  mean = mean_fitresult->getVal();

		  float r2 = mean/(Are_pi0_? PI0MASS:ETAMASS);
		  r2 = r2*r2;
		  //cout<<"EBMEAN::"<<j<<":"<<mean<<" Saved if: "<<fitres.SoB<<">(isNot_2010_ ? 0.04:0.1) "<<(fitres.chi2/fitres.dof)<<" < 0.2 "<<fabs(mean-0.15)<<" >0.0000001) "<<endl;
		  //if( fitres.SoB>(isNot_2010_ ? 0.04:0.1) && (fitres.chi2/fitres.dof)< 0.5 && fabs(mean-0.15)>0.0000001) mean = 0.5 * ( r2 - 1. );
		  //if( fitres.chi2 < 5 && fabs(mean-(Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB))>0.0000001) mean = 0.5 * ( r2 - 1. );
		  if( fabs(mean-(Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB)) > 0.0000001 )
		    mean = 0.5 * ( r2 - 1. );
		  else 
		    mean = 0.;
		}
		else{
		  mean = 0.;
		}
	      }

	  }

	  if (isEoverEtrue_) {

	    if (foldInSuperModule_) {

	      // we assume j is a single xtal if we are here
	      // get SM 1 in EB+ to fill calibration map (it will be copied on all SM later).
	      // use EBDetId constructor that used SM number, EBDetId::ic() and EBDetId::SMCRYSTALMODE == 1
	      EBDetId ebid(1,j+1,1);
	      regionalCalibration_   ->getCalibMap()->coeff(ebid) *= (mean==0.)    ? 1. : 1./(mean);
	      regionalCalibration_g2_->getCalibMap()->coeff(ebid) *= (mean_g2==0.) ? 1. : 1./(mean_g2);
	      
	    } else {

	      std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
	      // actually it is just one crystal, unless we do a calibration based on trigger towers or etaring
	      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
		{
		    regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(mean);
		} // loop over DetId in regions

	      ids = regionalCalibration_g2_->allDetIdsInEBRegion(j);
	      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
		{
		  regionalCalibration_g2_->getCalibMap()->coeff(*iid) *= (mean_g2==0.) ? 1. : 1./(mean_g2);
		} // loop over DetId in regions	   

	    }	    

	  } else {

	    std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
	    // actually it is just one crystal, unless we do a calibration based on trigger towers or etaring
	    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	      {
		regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
	      } // loop over DetId in regions
	    
	  }
		  
	} // loop over regions

    }// if you have to fit barrel

    /// loop over EE crystals
    if( (EEoEB_ == "Endcap") && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){

      for(int jR = inRangeFit_; jR <=finRangeFit_ && jR < regionalCalibration_->getCalibMap()->getNRegionsEE(); jR++)
	{
	  cout << "FIT_EPSILON: Fitting EE Cristal--> " << jR << endl;
	  if(!(jR%1000))
	    cout << "FIT_EPSILON: fitting EE region " << jR << endl;

	  float mean = 0.;
	  float mean_g2 = 0.; // used only for E/Etrue with MC

	  if (isEoverEtrue_) {
		  
	    // int iMin = EoverEtrue_g1_EE_h[jR]->GetXaxis()->FindFixBin(0.6); 
	    // int iMax = EoverEtrue_g1_EE_h[jR]->GetXaxis()->FindFixBin(1.1);
	    // double integral = EoverEtrue_g1_EE_h[jR]->Integral(iMin, iMax);  
	    double integral = EoverEtrue_g1_EE_h[jR]->Integral();  

	    if(integral > EoverEtrue_integralMin) {

	      Pi0FitResult fitres = FitEoverEtruePeakRooFit(EoverEtrue_g1_EE_h[jR], false, jR, Pi0EE);
	      RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
	      mean = mean_fitresult->getVal();
		    
	    } else {

	      std::cout << "### g1 ### FIT_EPSILON: iR = " << jR << ", integral() = " << integral << " , skipping the fit " << std::endl;
	      mean = 0.;

	    }

	    // if(integral > EoverEtrue_integralMin) {

	    //   TFitResultPtr fitresptr = FitEoverEtruePeak( EoverEtrue_g1_EE_h[jR], false, jR, Pi0EE, false);
	    //   mean = fitresptr->Parameter(1);
	    //   if (mean >= 1.5) mean = 0.; 
		    
	    // } else {

	    //   std::cout << "### g1 ### FIT_EPSILON: iR = " << jR << ", integral() = " << integral << " , skipping the fit " << std::endl;
	    //   mean = 0.;
	    //   EEmap_fitresptr_g1[jR] = TFitResultPtr(-1);

	    // }

	    // iMin = EoverEtrue_g2_EE_h[jR]->GetXaxis()->FindFixBin(0.6); 
	    // iMax = EoverEtrue_g2_EE_h[jR]->GetXaxis()->FindFixBin(1.1);
	    // integral = EoverEtrue_g2_EE_h[jR]->Integral(iMin, iMax);  
	    integral = EoverEtrue_g2_EE_h[jR]->Integral();  

	    if(integral > EoverEtrue_integralMin) {

	      Pi0FitResult fitres = FitEoverEtruePeakRooFit(EoverEtrue_g1_EE_h[jR], true, jR, Pi0EE);
	      RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
	      mean_g2 = mean_fitresult->getVal();
		    
	    } else {

	      std::cout << "### g2 ### FIT_EPSILON: iR = " << jR << ", integral() = " << integral << " , skipping the fit " << std::endl;
	      mean_g2 = 0.;

	    }

	    // if(integral > EoverEtrue_integralMin) {

	    //   TFitResultPtr fitresptr = FitEoverEtruePeak( EoverEtrue_g2_EE_h[jR], true, jR, Pi0EE, false);
	    //   mean_g2 = fitresptr->Parameter(1);
	    //   if (mean_g2 >= 1.5) mean_g2 = 0.; 
		    
	    // } else {

	    //   std::cout << "### g2 ### FIT_EPSILON: iR = " << jR << ", integral() = " << integral << " , skipping the fit " << std::endl;
	    //   mean_g2 = 0.;
	    //   EEmap_fitresptr_g2[jR] = TFitResultPtr(-1);

	    // }
		  
	  } else {
		
	    if(!useMassInsteadOfEpsilon_ && epsilon_EE_h[jR]->Integral(epsilon_EE_h[jR]->GetNbinsX()*(1./6.),epsilon_EE_h[jR]->GetNbinsX()*0.5) > 20) 
	      {
		TF1 *ffit = new TF1("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
		ffit->SetParameters(100,0,0.1);
		ffit->SetParNames("Constant","Mean_value","Sigma","a","b");

		ffit->SetParLimits(0,0.,epsilon_EE_h[jR]->GetEntries()*1.1);
		ffit->SetParLimits(3,-500,500);
		ffit->SetParLimits(2,0.05,0.3);

		double Max = 0.;
		double Min = -0.5, bin = 0.0125;
		Max = Min+(bin*(double)epsilon_EE_h[jR]->GetMaximumBin());
		double Bound1 = -0.35, Bound2 = 0.35;
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.3;}
		if ( Max+Bound2 > 0.48  ){ Bound2 = 0.3;}
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.25;}
		if ( Max+Bound2 > 0.48  ){ Bound2 = 0.2;}
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.2;}
		if ( Max+Bound2 > 0.48  ){ Bound2 = 0.15;}
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.15;}
		if ( Max+Bound2 > 0.48  ){ Bound2 = 0.1;}
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.1;}
		if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.05;}
		//@@IterativeFit(epsilon_EE_h[jR], *ffit);
		//@@mean = ffit.GetParameter(1); 
		epsilon_EE_h[jR]->Fit(ffit,"qB","", Max+Bound1,Max+Bound2);

		if(ffit->GetNDF() != 0) {
		  double chi2 = ( ffit->GetChisquare()/ffit->GetNDF() );
		  if(chi2 > 11  ) { cout<<"DAMN:(EE) High Chi square..."<<endl; }
		}
		else cout<<"DAMN: NDF == 0"<<endl;
		mean = ffit->GetParameter(1);
	      }
	    else if(useMassInsteadOfEpsilon_)
	      {
		int iMin = epsilon_EE_h[jR]->GetXaxis()->FindFixBin(Are_pi0_? 0.08:0.4 ); 
		int iMax = epsilon_EE_h[jR]->GetXaxis()->FindFixBin(Are_pi0_? 0.18:0.65 );
		double integral = epsilon_EE_h[jR]->Integral(iMin, iMax);  

		if(integral>70.)
		  {
		    Pi0FitResult fitres = FitMassPeakRooFit( epsilon_EE_h[jR], 
							     Are_pi0_? fitRange_low_pi0:fitRange_low_etaEE, 
							     Are_pi0_? fitRange_high_pi0:fitRange_high_eta, 
							     jR, 1, Pi0EE, 0, isNot_2010_);//0.05-0.3
		    RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
		    mean = mean_fitresult->getVal();
		    float r2 = mean/(Are_pi0_? PI0MASS:ETAMASS);
		    r2 = r2*r2;
		    //cout<<"EEMEAN::"<<jR<<":"<<mean<<" Saved if: "<<fitres.SoB<<">0.3 "<<(fitres.chi2/fitres.dof)<<" < (isNot_2010_? 0.07:0.35) "<<fabs(mean-0.14)<<" >0.0000001) "<<endl;
		    //if( (fitres.chi2/fitres.dof)<0.3 && fitres.SoB>(isNot_2010_? 0.07:0.35) && fabs(mean-0.14)>0.0000001 ) mean = 0.5 * ( r2 - 1. );
		    //if( fitres.chi2 < 5 && fabs(mean-(Are_pi0_? upper_bound_pi0mass_EE:upper_bound_etamass_EE))>0.0000001 ) mean = 0.5 * ( r2 - 1. );
		    // do not use Chi2 for goodness of fit. If I have many events, then the chi2 will be huge because the model will not pass through all data points
		    // on the oter hand, if I have few events, the statistical uncertainty is large and the Chi2 tends to be little
		    // better not to use Chi2
		    if(fabs(mean-(Are_pi0_? upper_bound_pi0mass_EE:upper_bound_etamass_EE))>0.0000001 ) 
		      mean = 0.5 * ( r2 - 1. );
		    else
		      mean = 0.;
		  }
		else
		  {
		    mean = 0.; 
		  }
	      }

	  }

	  std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
	  for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	    {
	      if (isEoverEtrue_) regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(mean);
	      else               regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
	    }

	  // now loop on second photon if doing E/Etrue
	  if (isEoverEtrue_) {
	    ids = regionalCalibration_g2_->allDetIdsInEERegion(jR);
	    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	      {
		regionalCalibration_g2_->getCalibMap()->coeff(*iid) *= (mean_g2==0.) ? 1. : 1./(mean_g2);
	      } // loop over DetId in regions		  
	  }

	}//for EE

    }// if you have to fit Endcap

}



void FitEpsilonPlot::IterativeFit(TH1F* h, TF1 & ffit) 
{
    float chi2 = 1000.;
    const int iterMax = 10;

    h->Fit(&ffit,"q","",-0.4,0.4);

    float mean = (ffit.GetParameters())[1];
    float sigma = (ffit.GetParameters())[2];
    float xmin = mean-2.*sigma;
    float xmax = mean+2.*sigma;

    ffit.SetRange(xmin,xmax);

    double par[3] = { ffit.GetParameters()[0], mean, sigma };

    for(int iter=0; iter< iterMax && chi2>5.; ++iter) 
    {
	  ffit.SetParameters(par[0],par[1], par[2]);

	  h->Fit(&ffit,"q","",xmin,xmax);
	  par[0] = (ffit.GetParameters())[0];
	  par[1] = (ffit.GetParameters())[1];
	  par[2] = (ffit.GetParameters())[2];

	  if(ffit.GetNDF()!=0) {
		chi2 = ffit.GetChisquare()/ffit.GetNDF();
	  }

    }
    return;
}


//-----------------------------------------------------------------------------------

Pi0FitResult FitEpsilonPlot::FitMassPeakRooFit(TH1F* h, double xlo, double xhi,  uint32_t HistoIndex, int ngaus, FitMode mode, int niter, bool isNot_2010_) 
{
    //-----------------------------------------------------------------------------------

    std::stringstream ind;
    ind << (int) HistoIndex;
    TString nameHistofit = "Fit_n_" + ind.str() + Form("_attempt%d",niter);

    // add canvas to save rooplot on top (will save this in the file)
    TCanvas* canvas = new TCanvas((nameHistofit+Form("_c")).Data(),"",700,700);
    canvas->cd();
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->cd();
    canvas->SetRightMargin(0.06);
    canvas->SetLeftMargin(0.15);

    Double_t upMassBoundaryEB = Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB; 
    Double_t upMassBoundaryEE = Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EE; 
    Double_t upMassBoundary = (mode==Pi0EB) ? upMassBoundaryEB : upMassBoundaryEE;
    // need a patch for some crystals in EB that might have the peak around 160 MeV due to high laser corrections.
    // depending on the year, the containment corrections might also increase a bit the peak position
    // the problem is that the peak is expected to be below 150 MeV when defining the signal model, so we have to catch these exception
    Double_t xValMaxHisto = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()+1); // use value just above maximum, it will be used to set the mean of the gaussian
    // check the maximum is within xlo and xhi (the histogram range is larger than the fit range)
    Double_t maxMassForGaussianMean = 0.0; //upper_bound_pi0mass_EB;
    // first check if peak is in the fit range (for EB it will be, in EE the background rises up and the maximum might not coincide wth peak)
    if (xValMaxHisto < xhi) {

      if (xValMaxHisto > upMassBoundary) {
	maxMassForGaussianMean = xValMaxHisto;
	xhi = Are_pi0_? fitRange_high_pi0_ext : fitRange_high_eta_ext;      //xhi + 0.012; // increase a bit the fit range
      } else {
	maxMassForGaussianMean = upMassBoundary;
      }

    } else {

      // need to loop on bins in the fit window
      Double_t ymaxHisto = 0.0;
      Int_t binYmaxHisto = -1;
      for (Int_t ibin = h->GetXaxis()->FindFixBin(xlo); ibin <= h->GetXaxis()->FindFixBin(xhi); ibin++) {
	if (h->GetBinContent(ibin) > ymaxHisto) {
	  ymaxHisto = h->GetBinContent(ibin);
	  binYmaxHisto = ibin;
	}
      }
      // check if maximum was found and it was not the last-1 bin
      // in that case, use the next bin to get max value for mass, just to avoid biases (that's why we asked last-1)
      if (binYmaxHisto > 0 && binYmaxHisto < (h->GetXaxis()->FindFixBin(xhi)-1)) {
	maxMassForGaussianMean = h->GetXaxis()->GetBinCenter(binYmaxHisto+1);
	if (maxMassForGaussianMean > upMassBoundary) xhi = Are_pi0_? fitRange_high_pi0_ext : fitRange_high_eta_ext;  //xhi + 0.012; // increase a bit the fit range

      } else {
	maxMassForGaussianMean = upMassBoundary; // if all this mess didn't work, just use the value we would have used in the beginning
      }

    }
    
    RooRealVar x("x","#gamma#gamma invariant mass",xlo, xhi, "GeV/c^2");

    RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h);

    //RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.5, Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB,"GeV/c^{2}");
    RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.45, maxMassForGaussianMean,"GeV/c^{2}");
    RooRealVar sigma("sigma","#pi^{0} core #sigma",Are_pi0_ ? 0.011 : 0.02, Are_pi0_ ? 0.005 : 0.01,Are_pi0_ ? 0.015 : 0.035,"GeV/c^{2}");


    if(mode==Pi0EE)  {
      mean.setRange( Are_pi0_? 0.1:0.42, maxMassForGaussianMean);
      mean.setVal(Are_pi0_? 0.13:0.52);
      sigma.setRange(Are_pi0_ ? 0.005 : 0.01, Are_pi0_ ? 0.020 : 0.05);
    }
    if(mode==Pi0EB && niter==1){
      mean.setRange(Are_pi0_? 0.105:0.47, maxMassForGaussianMean);
      sigma.setRange(Are_pi0_ ? 0.003 : 0.016, Are_pi0_ ? 0.030 : 0.03);	  
    }

    //RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
    RooRealVar Nsig("Nsig","#pi^{0} yield",h->Integral()*0.15,0.,h->Integral()*10.0);
    //Nsig.setVal( h->GetSum()*0.1);

    RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

    RooRealVar sigmaTail("sigmaTail","#pi^{0} tail #sigma",0.040, 0.020,0.065,"GeV/c^{2}");

    RooGaussian gaus2("gaus2","Tail Gaussian",x, mean,sigmaTail);

    RooRealVar fcore("fcore","f_{core}",0.9,0.,1.);
    RooAddPdf  signal("signal","signal model",RooArgList(gaus,gaus2),fcore);

    RooRealVar p0("p0","p0", 1000.,-1.e5,1.e5);
    RooRealVar p1("p1","p1", -3000.,-1.e5,1.e5);
    RooRealVar p2("p2","p2", 10000.,-1.e5,1.e5);
    RooRealVar p3("p3","p3", -10000.,-1.e5,1.e5);
    RooRealVar p4("p4","p4",-4000.,-1.e5,1.e5);
    RooRealVar p5("p5","p5", 5.,-1.e5,1.e5);
    RooRealVar p6("p6","p6", 6.,-1.e5,1.e5);

    RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
    RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
    RooRealVar cb2("cb2","cb2", 0.1,  -1.,1.);
    RooRealVar cb3("cb3","cb3",-0.1, -0.5,0.5);
    RooRealVar cb4("cb4","cb4", 0.1, -1.,1.);
    RooRealVar cb5("cb5","cb5", 0.1, -1.,1.);
    RooRealVar cb6("cb6","cb6", 0.3, -1.,1.);


    //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2) );
    //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2,cb3) );

    RooArgList cbpars(cb0,cb1,cb2);
    //if(mode==Pi0EE) cbpars.add( cb4);
    //if(mode==Pi0EE) cbpars.add( cb5);

    // try to use a second order polynomial, if the fit is bad add other terms
    // if you start with many terms, the fit creates strange curvy shapes trying to fit the statistical fluctuations
    // 2nd order means a curve with no change of concavity
    
    if(niter==1){
      cbpars.add(cb3);
    }
    if(niter==2){
      cb3.setRange(-1,1.);
      cb4.setRange(-0.3,0.3);
      cbpars.add( cb3);
      cbpars.add( cb4 );     
    }
    if(niter==3){
      cb3.setRange(-1,1.);
      cb4.setRange(-1,1);
      cb5.setRange(-0.5, 0.5);
      cbpars.add( cb3);
      cbpars.add( cb4 );
      cbpars.add( cb5 );
    }

    RooChebychev bkg("bkg","bkg model", x, cbpars );

    //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3,p4,p5,p6) );
    //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3) );

    //RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,1.e8);
    RooRealVar Nbkg("Nbkg","background yield",h->Integral()*0.85,0.,h->Integral()*10.0);
    //Nbkg.setVal( h->GetSum()*0.8 );

    RooAbsPdf* model=0;

    RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
    RooAddPdf model2("model","sig+bkg",RooArgList(signal,bkg),RooArgList(Nsig,Nbkg));

    if(ngaus==1)      model = &model1;
    else if(ngaus==2) model = &model2;


    RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
    //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same

    RooFitResult* res = nullptr;
    RooMinuit m(nll);
    RooMinimizer mfit(nll);

    if (useFit_RooMinuit_) {

      // // original fit
      // // obsolete: see here --> https://root-forum.cern.ch/t/roominuit-and-roominimizer-difference/18230/8
      // // better to use RooMinimizer, but please read caveat below
      m.setVerbose(kFALSE);
      //m.setVerbose(kTRUE);
      m.migrad();
      m.hesse();  // sometimes it fails, caution
      res = m.save() ;

    } else {

      // alternative fit (results are pretty much the same)
      // IMPORTANT, READ CAREFULLY: sometimes this method fails.
      // This happens because at the boundaries of the fit range the pdf goea slightly below 0 (so it is negative). The fitter tries to cope wth it and should tipically
      // manage to converge. However, I noticed that after few attemps (even though the default number of attemps should be several hundreds or thousands of times) 
      // the job crashes, and this seems to be a feature of cmssw, not of RooFit
      // The reason why the pdf gets negative could be due to the fact that, regardless the chosen fit range given by xlo and xhi, the actual fit range goes from the 
      // lower edge of the leftmost bin containing xlo to the upper edge of the rightmost one containing xhi, but then the fit tries to "pass" across the bin centers
      // Therefore, for a sharply rising (or falling) distribution, the pdf can become negative
      // The consequence is that there are large areas in the calibration map of related 2D plots that are white (because the fit there was not done succesfully)
      // The previous method using RooMinuit seems to be more robust, so I suggest we should use that one even though it is said to be obsolete
      mfit.setVerbose(kFALSE);
      mfit.setPrintLevel(-1);
      mfit.setStrategy(2);  // 0,1,2:  MINUIT strategies for dealing most efficiently with fast FCNs (0), expensive FCNs (2) and 'intermediate' FCNs (1)
      //cout << "FIT_EPSILON: Minimize" << endl;
      mfit.minimize("Minuit2","minimize");
      //cout << "FIT_EPSILON: Minimize hesse " << endl;
      mfit.minimize("Minuit2","hesse");
      //cout<<"FIT_EPSILON: Estimate minos errors for all parameters"<<endl;
      mfit.minos(RooArgSet(Nsig,Nbkg,mean));
      res = mfit.save() ;

    }

    // cbasile [CMSSW_13_3_0_pre3] : with ROOT v6.24 and earlyer 
    //                                 RooChi2Var (const char *name, const char *title, RooAbsPdf &pdf, RooDataHist &data, Bool_t extended=kFALSE, ...)
    //                               with ROOT v6.26 
    //                                 RooChi2Var (const char *name, const char *title, RooAbsPdf &pdf, RooDataHist &data, RooAbsTestStatistic::Configuration const &cfg=RooAbsTestStatistic::Configuration{}, bool extended=false, RooDataHist::ErrorType=RooDataHist::SumW2) 
    //                                 > include RooAbsTestStatistic with deafult parameters (check https://root.cern/doc/v626/RooAbsTestStatistic_8h_source.html)
   
    RooChi2Var chi2("chi2","chi2 var",*model,dh, RooAbsTestStatistic::Configuration{}, true);
    // use only bins in fit range for ndof (dh is made with var x that already has the restricted range, but h is the full histogram)
    //int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();
    int ndof = h->FindFixBin(xhi) - h->FindFixBin(xlo) +1 - res->floatParsFinal().getSize(); 

    //compute S/B and chi2
    x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
    RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));

    RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

    float normSig = integralSig->getVal();
    float normBkg = integralBkg->getVal();

    Pi0FitResult pi0res; // this is the output value of this method
    pi0res.res = res;

    pi0res.S = normSig*Nsig.getVal();
    pi0res.Serr = normSig*Nsig.getError();

    pi0res.B = normBkg*Nbkg.getVal();
    pi0res.Berr = normBkg*Nbkg.getError();

    pi0res.SoB =  pi0res.S/pi0res.B;
    pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) + 
		pow(pi0res.Berr/pi0res.B,2) ) ;
    pi0res.dof = ndof;
    pi0res.nFitParam = res->floatParsFinal().getSize();


    RooPlot*  xframe = x.frame(h->GetNbinsX());
    //RooPlot*  xframe = x.frame(xlo, xhi);
    xframe->SetName((nameHistofit+Form("_rp")).Data());
    xframe->SetTitle(h->GetTitle());
    dh.plotOn(xframe, Name("data"));
    model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed), Name("bkgOnly"));
    model->plotOn(xframe,Components(gaus),LineStyle(kDashed), LineColor(kGreen+1), Name("sigOnly"));
    model->plotOn(xframe, Name("model"));

    // TMAth::Prob() uses Chi2, not reduced Chi2, while xframe->chiSquare() returns the reduced Chi2
    pi0res.chi2 = xframe->chiSquare("model","data",pi0res.nFitParam) * pi0res.dof;
    pi0res.probchi2 = TMath::Prob(pi0res.chi2, ndof);

    xframe->Draw();

    cout << "FIT_EPSILON: Nsig: " << Nsig.getVal() 
	 << " nsig 3sig: " << normSig*Nsig.getVal()
	 << " nbkg 3sig: " << normBkg*Nbkg.getVal()
	 << " S/B: " << pi0res.SoB << " +/- " << pi0res.SoBerr
	 << " chi2: " << pi0res.chi2
	 << " chi2 reduced: " << pi0res.chi2 / pi0res.dof
	 << " DOF: " << pi0res.dof
	 << " N(fit.param.): " << pi0res.nFitParam
	 << " prob(chi2): " << pi0res.probchi2
	 << endl;

    if(mode==Pi0EB){
	  EBmap_Signal[HistoIndex]=pi0res.S;
	  EBmap_Backgr[HistoIndex]=pi0res.B;
	  EBmap_Chisqu[HistoIndex]=xframe->chiSquare();
	  EBmap_ndof[HistoIndex]=ndof;
	  EBmap_mean[HistoIndex]=mean.getVal();
	  EBmap_mean_err[HistoIndex]=mean.getError();
	  EBmap_sigma[HistoIndex]=sigma.getVal();
	  EBmap_Snorm[HistoIndex]=normSig;
	  EBmap_b0[HistoIndex]=cb0.getVal();
	  EBmap_b1[HistoIndex]=cb1.getVal();
	  EBmap_b2[HistoIndex]=cb2.getVal();
	  EBmap_b3[HistoIndex]=cb3.getVal();
	  EBmap_Bnorm[HistoIndex]=normBkg;
    }
    if(mode==Pi0EE){
	  EEmap_Signal[HistoIndex]=pi0res.S;
	  EEmap_Backgr[HistoIndex]=pi0res.B;
	  EEmap_Chisqu[HistoIndex]=xframe->chiSquare();
	  EEmap_ndof[HistoIndex]=ndof;
	  EEmap_mean[HistoIndex]=mean.getVal();
	  EEmap_mean_err[HistoIndex]=mean.getError();
	  EEmap_sigma[HistoIndex]=sigma.getVal();
	  EEmap_Snorm[HistoIndex]=normSig;
	  EEmap_b0[HistoIndex]=cb0.getVal();
	  EEmap_b1[HistoIndex]=cb1.getVal();
	  EEmap_b2[HistoIndex]=cb2.getVal();
	  EEmap_b3[HistoIndex]=cb3.getVal();
	  EEmap_Bnorm[HistoIndex]=normBkg;
    }

    TLatex lat;
    std::string line = "";
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.SetTextColor(1);

    float xmin(0.2), yhi(0.80), ypass(0.05);
    if(mode==EtaEB) yhi=0.30;
    if(mode==Pi0EE) yhi=0.5;
    line = Form("Yield: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
    lat.DrawLatex(xmin,yhi, line.c_str());

    line = Form("m_{#gamma#gamma}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
    lat.DrawLatex(xmin,yhi-ypass, line.c_str());

    line = Form("#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
    lat.DrawLatex(xmin,yhi-2.*ypass, line.c_str());

    //sprintf(line,"S/B(3#sigma): %.2f #pm %.2f", pi0res.SoB, pi0res.SoBerr );
    line = Form("S/B(3#sigma): %.2f", pi0res.SoB );
    lat.DrawLatex(xmin,yhi-3.*ypass, line.c_str());

    line = Form("#Chi^{2}: %.2f (%d dof)", pi0res.chi2, pi0res.dof );
    lat.DrawLatex(xmin,yhi-4.*ypass, line.c_str());

    line = Form("B param. %d", cbpars.getSize() );
    lat.DrawLatex(xmin,yhi-5.*ypass, line.c_str());

    canvas->RedrawAxis("sameaxis");

    Pi0FitResult fitres = pi0res;
    //xframe->chiSquare() is the chi2 reduced, i.e., that whose expected value is 1
    // E[X^2]=v; Var[X^2]=2v --> fit is bad if |X^2-v|>5*sqrt(2v) 

    //if(mode==Pi0EB && ( xframe->chiSquare()/pi0res.dof>0.35 || pi0res.SoB<0.6 || fabs(mean.getVal()-(Are_pi0_? 0.150:0.62))<0.0000001 ) ){
    //bool badChi2 = fabs(xframe->chiSquare() - pi0res.dof) > 5.0 * sqrt(2. * pi0res.dof);

    if(mode==Pi0EB && ( fabs(mean.getVal()-maxMassForGaussianMean)<0.0000001 ) ){
	  if(niter==0) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 1, isNot_2010_);
	  if(niter==1) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 2, isNot_2010_);
	  if(niter==2) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 3, isNot_2010_);
    }
    if(mode==Pi0EE && ( fabs(mean.getVal()-maxMassForGaussianMean)<0.0000001 ) ){
	  if(niter==0) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 1, isNot_2010_);
	  if(niter==1) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 2, isNot_2010_);
	  if(niter==2) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 3, isNot_2010_);
    }

    // save last version of fit made
    // if(StoreForTest_ && niter==0){
    if(StoreForTest_){
      outfileTEST_->cd();
      xframe->Write();
      canvas->Write();
    }

    delete canvas;
    return fitres;
}


//------------------------------------------------
// method to fit E/Etrue

//-----------------------------------------------------------------------------------

//=====================================================================

Float_t my2sideCrystalBall(double* x, double* par) {

  // implementation of a 2-sided crystal ball
  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Float_t xcur = x[0];
  Float_t N = par[0];
  Float_t mu = par[1];
  Float_t sigma = par[2];
  Float_t alphaL = par[3];
  Float_t nL = par[4];
  Float_t alphaR = par[5];
  Float_t nR = par[6];
  Float_t t = (xcur-mu)/sigma;
  Float_t absAlphaL = fabs((Float_t)alphaL);
  Float_t invAbsAlphaL = 1./absAlphaL;
  Float_t absAlphaR = fabs((Float_t)alphaR);
  Float_t invAbsAlphaR = 1./absAlphaR;

  if ( t<-absAlphaL ) {
    Float_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Float_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    return N*exp(-0.5*t*t);
  } else {
    Float_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Float_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

//=====================================================================

Float_t myLeftTailCrystalBall(double* x, double* par) {

  // implementation of a left-tail crystal ball

  Float_t xcur = x[0];
  Float_t N = par[0];
  Float_t mu = par[1];
  Float_t sigma = par[2];
  Float_t alphaL = par[3];
  Float_t nL = par[4];
  Float_t t = (xcur-mu)/sigma;
  Float_t absAlphaL = fabs((Float_t)alphaL);
  Float_t invAbsAlphaL = 1./absAlphaL;

  if ( t<-absAlphaL ) {
    Float_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Float_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else {
    return N*exp(-0.5*t*t);
  }

}

//=====================================================================

Float_t myRightTailCrystalBall(double* x, double* par) {

  // implementation of a right-tail crystal ball

  Float_t xcur = x[0];
  Float_t N = par[0];
  Float_t mu = par[1];
  Float_t sigma = par[2];
  Float_t alphaR = par[3];
  Float_t nR = par[4];
  Float_t t = (xcur-mu)/sigma;
  Float_t absAlphaR = fabs((Float_t)alphaR);
  Float_t invAbsAlphaR = 1./absAlphaR;

  if ( t>absAlphaR ) {
    Float_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Float_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  } else {
    return N*exp(-0.5*t*t);
  }

}


//======================================================

TFitResultPtr FitEpsilonPlot::FitEoverEtruePeak(TH1F* h1, Bool_t isSecondGenPhoton, uint32_t HistoIndex, FitMode mode, Bool_t noDrawStatBox) 
{


  //////////////////////
  //
  // In release CMSSW_10_2_0 there seems to be an issue with the usage of the user define functions like my2sideCrystalBall
  // When compiling, I get an error like the following
  // Since this fucntion is not used anymore (we use FitEpsilonPlot::FitEoverEtruePeakRooFit), I just comment everything out and return 0
  return 0;

}

//=============================================

My_double_CB::My_double_CB(const char *name, const char *title, 
			   RooAbsReal& _x,
			   RooAbsReal& _mu,
			   RooAbsReal& _sig,
			   RooAbsReal& _a1,
			   RooAbsReal& _n1,
			   RooAbsReal& _a2,
			   RooAbsReal& _n2) :
RooAbsPdf(name,title), 
	    x("x","x",this,_x),
	    mu("mu","mu",this,_mu),
	    sig("sig","sig",this,_sig),
	    a1("a1","a1",this,_a1),
	    n1("n1","n1",this,_n1),
	    a2("a2","a2",this,_a2),
	    n2("n2","n2",this,_n2)
{ 
} 


My_double_CB::My_double_CB(const My_double_CB& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mu("mu",this,other.mu),
  sig("sig",this,other.sig),
  a1("a1",this,other.a1),
  n1("n1",this,other.n1),
  a2("a2",this,other.a2),
  n2("n2",this,other.n2)
{ 
} 



Double_t My_double_CB::evaluate() const 
{ 
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(n1/TMath::Abs(a1),n1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(n2/TMath::Abs(a2),n2)*TMath::Exp(-a2*a2/2);
  double B1  = n1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = n2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-n1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-n2);
  return result;
} 


//======================================================

Pi0FitResult FitEpsilonPlot::FitEoverEtruePeakRooFit(TH1F* h1, Bool_t isSecondGenPhoton, uint32_t HistoIndex, FitMode mode) 
{

  int ieta = -999;
  int iphi = -999;
    
  if (foldInSuperModule_) {
    //if (HistoIndex == 1336 || HistoIndex == 1611 || HistoIndex = 1625) HistoIndex = HistoIndex -1;
    EBDetId thisebid(1,HistoIndex+1,1);
    ieta = thisebid.ieta();
    iphi = thisebid.iphi();
  } else {
    EBDetId thisebid(EBDetId::detIdFromDenseIndex(HistoIndex) ); 
    ieta = thisebid.ieta();
    iphi = thisebid.iphi();
  }

  // some flags for the fit
  // using the double Crystal Ball seems the better choice, but might require some parameters tuning
  bool simpleFitTo = true; // will use neither RooMinuit nor RooMinimizer, but rather simple x.fitTo()
  bool noFitBkg = isSecondGenPhoton ? false : true;  // if true, use only signal model for the fit
  bool useRooCMSShapeAsBkg = isSecondGenPhoton ? false : false;
  bool useCBtoFit = isSecondGenPhoton ? true : false;  // use Crystal Ball (tail orientation depends on the parameter alpha given below)
  //  bool useCB2toFit = isSecondGenPhoton ? true : false;
  bool useCB2toFit = isSecondGenPhoton ? false : true;   // use double Crystal Ball (overrides useCBtoFit)
  //bool useCB2toFit = isSecondGenPhoton ? false : false;   // use double Crystal Ball (overrides useCBtoFit)
  bool usePol2 = false;  // when using background model as well, use pol2 (by default pol3 is used, if RooCMSShape is not specified)
  bool usePol1 = false;

  double hardCodedXmin = -1.0;
  double hardCodedXmax = -1.0;
  // for 2018, after bug fix
  if (isSecondGenPhoton) {  
    // ok except for some crystals, were it could actually work but it didn't
    hardCodedXmin = 0.25;
    hardCodedXmax = 1.45;
    useCB2toFit = false;
    noFitBkg = false;
    useRooCMSShapeAsBkg = false;
    usePol2 = true;
    // following is ok, but sometimes fit doesn't converge
    if (ieta >= 55) {
      // hardCodedXmin = 0.02;
      // hardCodedXmax = 1.45;
      // useCB2toFit = false;
      // noFitBkg = false;
      // useRooCMSShapeAsBkg = true;
      // usePol2 = false;
      hardCodedXmin = 0.02;
      hardCodedXmax = 1.45;
      useCB2toFit = false;
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = true;
    }
    if (ieta == 1 and (iphi == 1 || iphi == 19)) { 
      noFitBkg = true;      
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      useCB2toFit = true; 
      hardCodedXmin = 0.75;
      hardCodedXmax = 1.1;
    } else if (iphi == 1 and (ieta == 7 || ieta == 8 || ieta == 20 || ieta == 26)) {
      noFitBkg = true; 
      useRooCMSShapeAsBkg = false;   
      //useCB2toFit = true;
      useCB2toFit = true; 
      usePol2 = false;   
      hardCodedXmin = 0.75; 
      hardCodedXmax = 1.1;
    } else if ((iphi == 19 and (ieta == 19 || ieta == 60)) || (ieta == 25 and iphi == 3) || ((ieta == 32 || ieta == 54) and iphi == 12)) {
      noFitBkg = false; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      hardCodedXmin = 0.1;
      hardCodedXmax = 1.3; 
    } else if (ieta == 27 and iphi == 20) {
      noFitBkg = true; 
      useRooCMSShapeAsBkg = false;   
      usePol2 = false;   
      useCB2toFit = true;
      hardCodedXmin = 0.8;
      hardCodedXmax = 1.1; 
    } else if (iphi == 2 and ieta == 48) {
      noFitBkg = false; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      hardCodedXmin = 0.05;     
    } else if (iphi == 7 and (ieta == 45 || ieta == 53)) {
      noFitBkg = true; 
      useRooCMSShapeAsBkg = false;   
      usePol2 = false;   
      useCB2toFit = true;
      hardCodedXmin = 0.8;     
      hardCodedXmax = 1.1; 
    } else if ((iphi == 20 and ieta == 63) || (ieta == 64 and (iphi == 1 || iphi == 20))) {
      noFitBkg = true; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      useCB2toFit = true;
      hardCodedXmin = 0.8;
      hardCodedXmax = 1.1; 
    } else if (iphi == 4 and ieta == 65) {
      noFitBkg = false; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      hardCodedXmin = 0.1;
      hardCodedXmax = 1.3; 
    } else if (ieta == 66 and (iphi == 5 || iphi == 14)) {
      noFitBkg = true; 
      useRooCMSShapeAsBkg = false;   
      usePol2 = false;   
      hardCodedXmin = 0.8;
      hardCodedXmax = 1.1; 
    } else if ((iphi == 20 and ieta == 67) || (ieta == 69 and iphi == 10) || (ieta == 74 and iphi == 17) || (ieta == 76 and iphi == 9)) {
      noFitBkg = false; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      hardCodedXmin = 0.05; 
    } else if ((iphi == 20 and (ieta == 77 || ieta == 78)) || (ieta == 84 and iphi == 18)) {
      noFitBkg = false; 
      useRooCMSShapeAsBkg = true;   
      usePol2 = false;   
      hardCodedXmin = 0.05; 
    } else if (ieta ==78 and iphi == 1) { 
      noFitBkg = true;      
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      //useCB2toFit = true; 
      hardCodedXmin = 0.8;
      hardCodedXmax = 1.05;
    }

  } else {
    if (ieta == 83 and iphi == 18) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      hardCodedXmin = 0.8;
    } else if (ieta == 26 and iphi == 1) {
      noFitBkg = true;
      hardCodedXmin = 0.7;
    } else if (ieta == 82 and iphi == 18) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = true;
      usePol2 = false;
      hardCodedXmin = 0.2;
    } else if (ieta == 75 and iphi == 8) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      hardCodedXmin = 0.3;
    } else if (ieta == 78 and iphi == 20) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      hardCodedXmin = 0.3;
    // } else if (ieta == 64 and iphi == 5) {
    //   noFitBkg = false;
    //   useRooCMSShapeAsBkg = false;
    //   usePol2 = false;
    //   hardCodedXmin = 0.3;
    } else if (ieta == 46 and iphi == 5) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = false;
      hardCodedXmin = 0.8;
    } else if (ieta == 40 and iphi == 9) {
      noFitBkg = false;
      useRooCMSShapeAsBkg = false;
      usePol2 = true;
      hardCodedXmin = 0.4;
    } else if (ieta == 65 and iphi == 1) {
      // noFitBkg = false;
      // useRooCMSShapeAsBkg = false;
      // usePol2 = true;
      hardCodedXmin = 0.8;
    }
  }


  // hardcoded stuff for CC in 2018

  // // used with bugged CC for 2018 (bad second photon E/Etrue, where E was the one of photon 1)
  // if (isSecondGenPhoton) {

  //   if (ieta == 5 && iphi == 14) {
  //     usePol2 = true;
  //     hardCodedXmin = 0.4;
  //     hardCodedXmax = 1.45;      
  //   } else if ((ieta == 5 && iphi == 6) || (ieta == 8 && iphi == 14)) {
  //     usePol2 = true;
  //     hardCodedXmin = 0.3;
  //     hardCodedXmax = 1.45;      
  //   } else if ((ieta == 8 && iphi == 5) || (ieta >= 15 && ieta <= 16 && iphi == 6)) {
  //     usePol2 = true;
  //     hardCodedXmin = 0.4;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 36 && iphi == 18) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 40 && iphi == 7) {
  //     usePol2 = true;
  //     hardCodedXmin = 0.4;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 45 && iphi == 14) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 46 && iphi == 16) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 50 && iphi == 4) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta >= 50 && ieta <= 51 && iphi == 15) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 55 && iphi == 3) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 52 && iphi == 3) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } if (ieta == 57 && iphi == 1) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta > 55) {
  //     useRooCMSShapeAsBkg = true;
  //     noFitBkg = false;
  //     usePol2 = false;
  //     useCB2toFit = false;
  //     useCBtoFit = false;  // CB up to 1.4 was fine for many crystals     
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   }
  //   // override few xtals above ieta=55
  //   if (ieta == 74 && iphi >= 15 && iphi <= 16) {
  //     useRooCMSShapeAsBkg = true;
  //     useCBtoFit = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 72 && iphi == 2) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 73 && iphi == 5) {
  //     useRooCMSShapeAsBkg = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 76 && iphi == 8) {
  //     useRooCMSShapeAsBkg = true;
  //     useCBtoFit = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } else if (ieta == 83 && iphi == 7) {
  //     useRooCMSShapeAsBkg = true;
  //     //useCBtoFit = true;
  //     hardCodedXmin = 0.2;
  //     hardCodedXmax = 1.45;
  //   } 

  // } else {
  //   if (ieta == 83 and iphi == 18) {
  //     noFitBkg = false;
  //     useRooCMSShapeAsBkg = false;
  //     usePol2 = false;
  //     hardCodedXmin = 0.8;
  //   } else if (ieta == 82 and iphi == 18) {
  //     noFitBkg = false;
  //     useRooCMSShapeAsBkg = true;
  //     usePol2 = false;
  //     hardCodedXmin = 0.2;
  //   } else if (ieta == 46 and iphi == 5) {
  //     noFitBkg = false;
  //     useRooCMSShapeAsBkg = false;
  //     usePol2 = false;
  //     hardCodedXmin = 0.8;
  //   } else if (ieta == 40 and iphi == 9) {
  //     noFitBkg = false;
  //     useRooCMSShapeAsBkg = false;
  //     usePol2 = true;
  //     hardCodedXmin = 0.4;
  //   }
  // }

  // // hardcoded stuff for CC in 2017
  // if (isSecondGenPhoton) {
  //   noFitBkg = false;
  //   useRooCMSShapeAsBkg = true;
  //   useCB2toFit = true;
  //   hardCodedXmin = 0.07;
  // } else {
  //   if (iphi%20 == 1 || iphi%20 == 0 || ieta == 1 || ieta == 25 || ieta == 26 || ieta == 45 || ieta == 46 || ieta == 65 || ieta == 66 || ieta == 85) {
  //     noFitBkg = false;
  //     useRooCMSShapeAsBkg = true;    
  //   }
  // }
  

  // std::cout << "FitEpsilonPlot::FitEoverEtruePeak called " << std::endl;
  int nPhoton = isSecondGenPhoton ? 2 : 1;

  int niter = 0; // attempt of the fit, only 1 for the moment
  TString nameHistofit = Form("Fit_n_%u_attempt%d_g%d",HistoIndex,niter,nPhoton);

  // add canvas to save rooplot on top (will save this in the file)
  TCanvas* canvas = new TCanvas((nameHistofit+Form("_c")).Data(),"",700,600);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetRightMargin(0.06);

  // get RMS in narrow range around the peak
  TH1F* h1narrow = new TH1F("h1narrow","",
			    1 + h1->FindFixBin(1.1) - h1->FindFixBin(0.8), 
			    h1->GetBinLowEdge(h1->FindFixBin(0.8)), 
			    h1->GetBinLowEdge(1 + h1->FindFixBin(1.1)));
  for (int i = 0; i <= h1narrow->GetNbinsX(); i++) {
    h1narrow->SetBinContent(i,h1->GetBinContent(h1->FindFixBin(h1narrow->GetBinCenter(i))));
  }
  float rmsh1narrow = h1narrow->GetStdDev();
  float xmaxbin = h1narrow->GetBinCenter(h1narrow->GetMaximumBin());  
  delete h1narrow;  

  //cout << "FIT_EPSILON: photon " << (isSecondGenPhoton ?  2 : 1) << " --> xmaxbin :: rmsh1narrow = " << xmaxbin << " :: " << rmsh1narrow << endl;

  // fit range, allow for differences between two photons
  // float_t xlo = isSecondGenPhoton ? 0.82 : 0.82;
  // float_t xlo = std::max(0.82, isSecondGenPhoton ? (xmaxbin - 1.6 * rmsh1narrow) : (xmaxbin - 1.7 * rmsh1narrow));
  // //float_t xhi = isSecondGenPhoton ? 1.15 : 1.15;
  // float_t xhi = std::min(1.15, isSecondGenPhoton ? (xmaxbin + 2.0 * rmsh1narrow) : (xmaxbin + 2.2 * rmsh1narrow));

  float_t xlo = 0.5;
  float_t xhi = 1.4;

  if (noFitBkg) {
    if (useCB2toFit or useCBtoFit) {
      xlo = isSecondGenPhoton ? 0.87 : 0.87;
      xhi = 1.48;      
    } else {
      xlo = 0.8;
      xhi = 1.1;
    }
  } else {
    if (useRooCMSShapeAsBkg) {
      xlo = 0.15;
      xhi = 1.48;
    } else{
      xlo = 0.6;
      xhi = 1.48;
    }
  }

  if (hardCodedXmin > 0.0) xlo = hardCodedXmin;
  if (hardCodedXmax > 0.0) xhi = hardCodedXmax;

  RooRealVar x("x",Form("#gamma_{%d} E/E_{true}",nPhoton), 0.0, 1.5, "");
  RooDataHist dh("dh",Form("#gamma_{%d} E/E_{true}",nPhoton),RooArgList(x),h1);

  RooRealVar mean("mean","peak position", 
		  xmaxbin, 
		  std::min(xmaxbin - 0.01, std::max(xmaxbin-0.1,isSecondGenPhoton ? 0.88 : 0.89)), 
		  std::max(xmaxbin + 0.01,std::min(0.98,xmaxbin+0.1)), "");
  //Double_t g2sigma = noFitBkg ? 0.14 : 0.12;
  Double_t g2sigma = 0.14;
  RooRealVar sigma("sigma","core #sigma",rmsh1narrow, std::min(0.025,0.9*rmsh1narrow),isSecondGenPhoton ? g2sigma : 0.1,"");

  //  RooRealVar Nsig("Nsig","signal yield",h1->Integral()*0.7,0.,h1->Integral()*1.1); // signal represents the peak in E/Etrue (even though it is actually only signal)
  //Nsig.setVal( h->GetSum()*0.1);
  RooRealVar Nsig("Nsig","signal yield",
		  isSecondGenPhoton ? 0.7 : 0.8 ,
		  0.,
		  isSecondGenPhoton ? 0.95 : 1.0); // should use normalization
  if (noFitBkg) Nsig.setRange(0.0,1.0);  

  RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

  RooRealVar alphaCB  ("alphaCB","",
		       useCB2toFit ? 2.5  : -1, 
		       useCB2toFit ? 0.01 : -5.0,
		       useCB2toFit ? 5.0  : -0.01);  // for CB1, this alpha parameter is positive for left tail, negative for right one
  RooRealVar nCB ("nCB","",0.5,0.01,50);
  RooRealVar alphaCB2  ("alphaCB2","",1.07,0.05,5.0);  // the double CB I use requires positive parameter
  RooRealVar nCB2 ("nCB2","",2,0.01,50);
  RooCBShape cb_sig ("cb_sig","Crystal Ball",x, mean, sigma, alphaCB,nCB);
  My_double_CB cb2_sig = My_double_CB("cb2_sig", "cb2_sig", x, mean, sigma, alphaCB,nCB, alphaCB2,nCB2);
  
  RooRealVar alphaCMSshape("alphaCMSshape","alphaCMSshape", 1.0, -50.0,50.0);
  RooRealVar betaCMSshape("betaCMSshape","betaCMSshape", 3.0, 0.0, 50.0);
  RooRealVar gammaCMSshape("gammaCMSshape","gammaCMSshape", 3.0, -50.0,50.0);
  RooRealVar peakCMSshape("peakCMSshape","peakCMSshape", 0.2, 0.01,0.5);
  RooCMSShape cmsshape = RooCMSShape("rooCMSshape","cmsShape",x,alphaCMSshape,betaCMSshape,gammaCMSshape,peakCMSshape);  

  // cb0 is not used, because RooChebychev assumes we give parameters from the linear term
  RooRealVar cb0("cb0","cb0", 0.0, -10.0,200.0);
  RooRealVar cb1("cb1","cb1", 0.0, -50.,50);
  RooRealVar cb2("cb2","cb2",-1  ,  -5.,5.);
  RooRealVar cb3("cb3","cb3", 0.0,  -5.,5.);
  RooRealVar cb4("cb4","cb4", 0.0,  -5.,5.);
  RooRealVar cb5("cb5","cb5", 0.0,  -5.,5.);
  RooRealVar cb6("cb6","cb6", 0.0,  -5.,5.);
  //RooRealVar cb7("cb7","cb7", 0.0,  -5.,5.);

  //RooRealVar p0("p0","p0", 100.0, 0.0, 500.0);
  RooRealVar p1("p1","p1", -5, -500.0, 0.0);

  // define a background shape in addition for the bare gaussian or Crystal Ball for the peak
  //RooArgList cbpars(cb0,cb1,cb2);  
  RooArgList cbpars(cb1,cb2);
  //RooArgList cbparsMore(cb0,cb1,cb2,cb3,cb4,cb5,cb6);
  RooArgList cbparsMore(cb1,cb2,cb3);
  RooArgList pol1pars(p1);
  RooArgList *cbparsPtr = usePol2 ? &cbpars : &cbparsMore;
  if (usePol1) cbparsPtr = &pol1pars;
  RooChebychev bkg("bkg","bkg model", x, *cbparsPtr );
  //  RooRealVar Nbkg("Nbkg","background yield",h1->Integral()*0.3,0.,h1->Integral()*1.1);
  RooRealVar Nbkg("Nbkg","background yield",
		  isSecondGenPhoton ? 0.3 : 0.3 ,
		  0.,
		  isSecondGenPhoton ? 0.9 : 0.8);

  //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3,p4,p5,p6) );
  //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3) );

  RooAbsPdf* model=0;

  // pass only Nsig if it is a fraction of events (0 < x < 1). If both Nsig and Nbkg are passed, the extended likelihood case is assumed, and they would be
  // interpreted as actual number of events
  // RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
  // RooAddPdf model2("model","sig+bkg",RooArgList(cb_sig,bkg),RooArgList(Nsig,Nbkg));
  // RooAddPdf model3("model","sig+bkg",RooArgList(cb2_sig,bkg),RooArgList(Nsig,Nbkg));
  RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig));
  RooAddPdf model2("model","sig+bkg",RooArgList(cb_sig,bkg),RooArgList(Nsig));
  RooAddPdf model3("model","sig+bkg",RooArgList(cb2_sig,bkg),RooArgList(Nsig));

  RooAddPdf model4("model","sig",RooArgList(cb2_sig),RooArgList(Nsig));
  RooAddPdf model5("model","sig",RooArgList(cb_sig),RooArgList(Nsig));
  RooAddPdf model6("model","sig",RooArgList(gaus),RooArgList(Nsig));

  // RooAddPdf model7("model","sig+bkg",RooArgList(cb2_sig,cmsshape),RooArgList(Nsig,Nbkg));
  // RooAddPdf model8("model","sig+bkg",RooArgList(cb_sig,cmsshape),RooArgList(Nsig,Nbkg));
  // RooAddPdf model9("model","sig+bkg",RooArgList(gaus,cmsshape),RooArgList(Nsig,Nbkg));
  RooAddPdf model7("model","sig+bkg",RooArgList(cb2_sig,cmsshape),RooArgList(Nsig));
  RooAddPdf model8("model","sig+bkg",RooArgList(cb_sig,cmsshape),RooArgList(Nsig));
  RooAddPdf model9("model","sig+bkg",RooArgList(gaus,cmsshape),RooArgList(Nsig));

  if (noFitBkg) {
    if (useCB2toFit)     model = &model4;
    else if (useCBtoFit) model = &model5;
    else                 model = &model6;
  }
  else {
    if (useCB2toFit) {
      if (useRooCMSShapeAsBkg) model = &model7;
      else                     model = &model3;
    } else if (useCBtoFit) {
      if (useRooCMSShapeAsBkg) model = &model8;
      else                     model = &model2;
    } else {
      if (useRooCMSShapeAsBkg) model = &model9;
      else                     model = &model1;
    }
  }

  // cout << "===============================================" << endl;
  // cout << "===============================================" << endl;
  // cout << "===============================================" << endl;

  //RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::Range(xlo,xhi));
  //RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(0), RooFit::SumW2Error(kTRUE), RooFit::Range(xlo,xhi));
  //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same

  RooFitResult* res = nullptr;

  // if (useCB2toFit) res = cb2_sig.fitTo(dh,RooFit::SumW2Error(kTRUE),RooFit::Save(),RooFit::Range(xlo, xhi));
  // else if (useCBtoFit) res = cb_sig.fitTo(dh,RooFit::SumW2Error(kTRUE),RooFit::Save(),RooFit::Range(xlo, xhi));
  // else            res = gaus.fitTo(dh,RooFit::SumW2Error(kTRUE),RooFit::Save(),RooFit::Range(xlo, xhi));
  if (simpleFitTo) {
    //res = model->fitTo(dh,RooFit::Save(),RooFit::Range(xlo, xhi),RooFit::SumW2Error(kTRUE),RooFit::Strategy(2),RooFit::Minos(RooArgSet(mean)),RooFit::PrintLevel(-1));
    res = model->fitTo(dh,RooFit::Save(),RooFit::Range(xlo, xhi),RooFit::Strategy(2),RooFit::PrintLevel(-1),RooFit::SumW2Error(kFALSE));
    // if (useCB2toFit)     res = cb2_sig.fitTo(dh,RooFit::Save(),RooFit::Range(xlo, xhi),RooFit::SumW2Error(kTRUE));
    // else if (useCBtoFit) res = cb_sig.fitTo(dh,RooFit::Save(),RooFit::Range(xlo, xhi),RooFit::SumW2Error(kTRUE));
    // else                 res = gaus.fitTo(dh,RooFit::Save(),RooFit::Range(xlo, xhi),RooFit::SumW2Error(kTRUE));
  } else {

    // warning: I removed definition of nll and m and mfit from outside here, because I don't think I want to use them
    // in case I do, this might crash, because once res is returned, it might be destroyed outside this scope
    if (useFit_RooMinuit_) {

      RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(0), RooFit::SumW2Error(kTRUE), RooFit::Range(xlo,xhi));
      RooMinuit m(nll);

      // // original fit
      // // obsolete: see here --> https://root-forum.cern.ch/t/roominuit-and-roominimizer-difference/18230/8
      // // better to use RooMinimizer, but please read caveat below
      m.setVerbose(kFALSE);
      //m.setVerbose(kTRUE);
      m.migrad();
      m.hesse();  // sometimes it fails, caution
      res = m.save() ;

    } else {

      RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(0), RooFit::SumW2Error(kTRUE), RooFit::Range(xlo,xhi));
      RooMinimizer mfit(nll);
      // alternative fit (results are pretty much the same)
      // IMPORTANT, READ CAREFULLY: sometimes this method fails.
      // This happens because at the boundaries of the fit range the pdf goea slightly below 0 (so it is negative). The fitter tries to cope wth it and should tipically
      // manage to converge. However, I noticed that after few attemps (even though the default number of attemps should be several hundreds or thousands of times) 
      // the job crashes, and this seems to be a feature of cmssw, not of RooFit
      // The reason why the pdf gets negative could be due to the fact that, regardless the chosen fit range given by xlo and xhi, the actual fit range goes from the 
      // lower edge of the leftmost bin containing xlo to the upper edge of the rightmost one containing xhi, but then the fit tries to "pass" across the bin centers
      // Therefore, for a sharply rising (or falling) distribution, the pdf can become negative
      // The consequence is that there are large areas in the calibration map of related 2D plots that are white (because the fit there was not done succesfully)
      // The previous method using RooMinuit seems to be more robust, so I suggest we should use that one even though it is said to be obsolete
      mfit.setVerbose(kFALSE);
      mfit.setPrintLevel(-1);
      mfit.setStrategy(2);  // 0,1,2:  MINUIT strategies for dealing most efficiently with fast FCNs (0), expensive FCNs (2) and 'intermediate' FCNs (1)
      //cout << "FIT_EPSILON: Minimize" << endl;
      mfit.minimize("Minuit2","minimize");
      //cout << "FIT_EPSILON: Minimize hesse " << endl;
      mfit.minimize("Minuit2","hesse");
      //cout<<"FIT_EPSILON: Estimate minos errors for all parameters"<<endl;
      //mfit.minos(RooArgSet(Nsig,Nbkg,mean));
      //mfit.minos(RooArgSet(mean));
      res = mfit.save() ;

    }

  }

  // cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  // cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  // cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

  // use only bins in fit range for ndof (dh is made with var x that already has the restricted range, but h is the full histogram)
  //int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();
  int ndof = h1->FindFixBin(xhi) - h1->FindFixBin(xlo) + 1 - res->floatParsFinal().getSize(); 

  //compute S/B and chi2
  //x.setRange("sobRange",mean.getVal() - 2.0*sigma.getVal(), mean.getVal() + 2.*sigma.getVal());
  x.setRange("sobRange",xlo,xhi);
  //RooChi2Var chi2("chi2","chi2 var",*model,dh, true,"sobRange");
  RooChi2Var chi2("chi2","chi2 var",*model,dh, RooAbsTestStatistic::Configuration{"sobRange"}, false);
  // cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
  RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));
  // cout << "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY" << endl;
  RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

  // cout << "--------------------------------------------------------------" << endl;
  // cout << "--------------------------------------------------------------" << endl;
  // cout << "--------------------------------------------------------------" << endl;


  float normSig = integralSig->getVal();
  float normBkg = integralBkg->getVal();

  // here we do not have pi0, but let's be consistent with the original value (we do not really use these parameters)
  Pi0FitResult pi0res; // this is the output value of this method
  pi0res.res = res;

  pi0res.S = normSig*Nsig.getVal();
  pi0res.Serr = normSig*Nsig.getError();

  pi0res.B = normBkg*Nbkg.getVal();
  pi0res.Berr = normBkg*Nbkg.getError();

  pi0res.SoB =  pi0res.S/pi0res.B;
  pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) + 
				    pow(pi0res.Berr/pi0res.B,2) ) ;
  pi0res.dof = ndof;
  pi0res.nFitParam = res->floatParsFinal().getSize();


  RooPlot*  xframe = x.frame(h1->GetNbinsX());
  //RooPlot*  xframe = x.frame(xlo, xhi);
  xframe->SetName((nameHistofit+Form("_rp")).Data());
  xframe->SetTitle(h1->GetTitle());
  dh.plotOn(xframe, Name("data"));
  //model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed), RooFit::Range(xlo,xhi));
  //if (useCBtoFit and isSecondGenPhoton) model->plotOn(xframe,Components(cb_sig),LineStyle(kDashed), LineColor(kGreen+1), RooFit::Range(xlo,xhi));
  if (not noFitBkg) {
    if (useRooCMSShapeAsBkg) model->plotOn(xframe,Components(cmsshape),LineStyle(kDashed), LineColor(kRed), RooFit::Range(xlo,xhi));
    else                     model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed), RooFit::Range(xlo,xhi));
    if (useCB2toFit)     model->plotOn(xframe,Components(cb2_sig),LineStyle(kDashed), LineColor(kGreen+1), RooFit::Range(xlo,xhi));
    else if (useCBtoFit) model->plotOn(xframe,Components(cb_sig),LineStyle(kDashed), LineColor(kGreen+1), RooFit::Range(xlo,xhi));
    else                 model->plotOn(xframe,Components(gaus),LineStyle(kDashed), LineColor(kGreen+1), RooFit::Range(xlo,xhi));
  }
  model->plotOn(xframe, Name("model"),RooFit::Range(xlo,xhi));

  // TMAth::Prob() uses Chi2, not reduced Chi2, while xframe->chiSquare() returns the reduced Chi2
  pi0res.chi2 = xframe->chiSquare("model","data",pi0res.nFitParam) * pi0res.dof;
  pi0res.probchi2 = TMath::Prob(pi0res.chi2, ndof);

  // set better y axis range for photon 2 when plotting points, such that TLatex text do not largely overlap with peak at low E/Etrue
  if (isSecondGenPhoton) {
  // get RMS in narrow range around the peak
    TH1F* h1lowEoverEtrue = new TH1F("h1lowEoverEtrue","",
				     1 + h1->FindFixBin(0.701) - h1->FindFixBin(0.001), 
				     h1->GetBinLowEdge(h1->FindFixBin(0.001)), 
				     h1->GetBinLowEdge(1 + h1->FindFixBin(0.701))
				     );
    for (int i = 0; i <= h1lowEoverEtrue->GetNbinsX(); i++) {
      h1lowEoverEtrue->SetBinContent(i,h1->GetBinContent(h1->FindFixBin(h1lowEoverEtrue->GetBinCenter(i))));
    }
    double ymaxInLowerRange = h1lowEoverEtrue->GetBinContent(h1lowEoverEtrue->GetMaximumBin());
    cout << "ieta, iphi, ic() = " << ieta << ", " << iphi << ", " << HistoIndex+1 << endl;
    cout << ">>>> xmaxInLowerRange : xmaxbin = " << h1lowEoverEtrue->GetBinCenter(h1lowEoverEtrue->GetMaximumBin()) << " : " << xmaxbin << endl;
    cout << ">>>> ymaxInLowerRange : h1->GetBinContent(h1->FindFixBin(xmaxbin)) = " << ymaxInLowerRange << " : " << h1->GetBinContent(h1->FindFixBin(xmaxbin)) << endl;
    if (ymaxInLowerRange > 0.8 * h1->GetBinContent(h1->FindFixBin(xmaxbin))) {
      xframe->GetYaxis()->SetRangeUser(0, 1.5 * ymaxInLowerRange);
    }
    delete h1lowEoverEtrue;
  }

  xframe->Draw();

  cout << "FIT_EPSILON: "
       << "photon " << nPhoton << "  "
       << " mean " << mean.getVal() << " +/- " << mean.getError()
       << " sigma " << sigma.getVal() << " +/- " << sigma.getError()
    //<< " Nsig: " << Nsig.getVal() 
    //<< " nsig 2sig: " << normSig*Nsig.getVal()
    //<< " nbkg 2sig: " << normBkg*Nbkg.getVal()
    //<< " S/B: " << pi0res.SoB << " +/- " << pi0res.SoBerr
       << " chi2: " << pi0res.chi2 << "/" << pi0res.dof
    //<< " chi2 reduced: " << pi0res.chi2 / pi0res.dof
    //<< " DOF: " << pi0res.dof
    //<< " N(fit.param.): " << pi0res.nFitParam
    //<< " prob(chi2): " << pi0res.probchi2
       << endl;

  TLatex lat;
  std::string line = "";
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.SetTextColor(1);

  float xmin(0.15), yhi(0.82), ypass(0.05);
  if(mode==EtaEB) yhi=0.30;
  if(mode==Pi0EE) yhi=0.5;
  if(mode==Pi0EB) {
    if (foldInSuperModule_) {
      EBDetId thisebid(1,HistoIndex+1,1);
      line = Form("i#eta = %d, i#phi = %d, ic() = %d", ieta, iphi, thisebid.ic());
    } else {
      line = Form("i#eta = %d, i#phi = %d", ieta, iphi);
    }
  } else {
    line = Form("#gamma_{%d}", nPhoton);
  }
  lat.DrawLatex(xmin,yhi, line.c_str());

  line = Form("peak: %.3f #pm %.3f", mean.getVal(), mean.getError() );
  lat.DrawLatex(xmin,yhi-ypass, line.c_str());

  line = Form("#sigma: %.3f #pm %.3f", sigma.getVal(), sigma.getError());
  lat.DrawLatex(xmin,yhi-2.*ypass, line.c_str());

  line = Form("#Chi^{2}: %.1f / %d", pi0res.chi2, pi0res.dof );
  lat.DrawLatex(xmin,yhi-3.*ypass, line.c_str());

  line = Form("fit param. %d", pi0res.nFitParam );
  lat.DrawLatex(xmin,yhi-4.*ypass, line.c_str());

  canvas->RedrawAxis("sameaxis");


  //////////////////////////////////
  //////////////////////////////////  
  
  // some parameters do not make sense for the E/Etrue study, but for simplicity we keep the same structure as the mass fit
  // basically we just need the peak position
  if (isSecondGenPhoton) {

    if(mode==Pi0EB){
      EBmap_Signal_g2[HistoIndex]=pi0res.S;
      EBmap_Backgr_g2[HistoIndex]=pi0res.B;
      EBmap_Chisqu_g2[HistoIndex]=xframe->chiSquare();
      EBmap_ndof_g2[HistoIndex]=ndof;
      EBmap_mean_g2[HistoIndex]=mean.getVal();
      EBmap_mean_err_g2[HistoIndex]=mean.getError();
      EBmap_sigma_g2[HistoIndex]=sigma.getVal();
      EBmap_Snorm_g2[HistoIndex]=normSig;
      EBmap_b0_g2[HistoIndex]=cb0.getVal();
      EBmap_b1_g2[HistoIndex]=cb1.getVal();
      EBmap_b2_g2[HistoIndex]=cb2.getVal();
      EBmap_b3_g2[HistoIndex]=cb3.getVal();
      EBmap_Bnorm_g2[HistoIndex]=normBkg;
    }
    if(mode==Pi0EE){
      EEmap_Signal_g2[HistoIndex]=pi0res.S;
      EEmap_Backgr_g2[HistoIndex]=pi0res.B;
      EEmap_Chisqu_g2[HistoIndex]=xframe->chiSquare();
      EEmap_ndof_g2[HistoIndex]=ndof;
      EEmap_mean_g2[HistoIndex]=mean.getVal();
      EEmap_mean_err_g2[HistoIndex]=mean.getError();
      EEmap_sigma_g2[HistoIndex]=sigma.getVal();
      EEmap_Snorm_g2[HistoIndex]=normSig;
      EEmap_b0_g2[HistoIndex]=cb0.getVal();
      EEmap_b1_g2[HistoIndex]=cb1.getVal();
      EEmap_b2_g2[HistoIndex]=cb2.getVal();
      EEmap_b3_g2[HistoIndex]=cb3.getVal();
      EEmap_Bnorm_g2[HistoIndex]=normBkg;
    }


  } else {

    if(mode==Pi0EB){
      EBmap_Signal[HistoIndex]=pi0res.S;
      EBmap_Backgr[HistoIndex]=pi0res.B;
      EBmap_Chisqu[HistoIndex]=xframe->chiSquare();
      EBmap_ndof[HistoIndex]=ndof;
      EBmap_mean[HistoIndex]=mean.getVal();
      EBmap_mean_err[HistoIndex]=mean.getError();
      EBmap_sigma[HistoIndex]=sigma.getVal();
      EBmap_Snorm[HistoIndex]=normSig;
      EBmap_b0[HistoIndex]=cb0.getVal();
      EBmap_b1[HistoIndex]=cb1.getVal();
      EBmap_b2[HistoIndex]=cb2.getVal();
      EBmap_b3[HistoIndex]=cb3.getVal();
      EBmap_Bnorm[HistoIndex]=normBkg;
    }
    if(mode==Pi0EE){
      EEmap_Signal[HistoIndex]=pi0res.S;
      EEmap_Backgr[HistoIndex]=pi0res.B;
      EEmap_Chisqu[HistoIndex]=xframe->chiSquare();
      EEmap_ndof[HistoIndex]=ndof;
      EEmap_mean[HistoIndex]=mean.getVal();
      EEmap_mean_err[HistoIndex]=mean.getError();
      EEmap_sigma[HistoIndex]=sigma.getVal();
      EEmap_Snorm[HistoIndex]=normSig;
      EEmap_b0[HistoIndex]=cb0.getVal();
      EEmap_b1[HistoIndex]=cb1.getVal();
      EEmap_b2[HistoIndex]=cb2.getVal();
      EEmap_b3[HistoIndex]=cb3.getVal();
      EEmap_Bnorm[HistoIndex]=normBkg;
    }

  }

  Pi0FitResult fitres = pi0res;


  // if(StoreForTest_ && niter==0){
  if(StoreForTest_){
    outfileTEST_->cd();
    xframe->Write();
    canvas->Write();
  }

  delete canvas;
  return fitres;

}





// ------------ method called once each job just before starting event loop  ------------
    void 
FitEpsilonPlot::beginJob()
{
    if(StoreForTest_){
      outfileTEST_ = new TFile(fitFileName_.c_str(),"RECREATE");
      if(!outfileTEST_) cout << "WARNING: file " << fitFileName_ << " with fit not created." << endl;
      else outfileTEST_->cd();
    }
}

// ------------ method called once each job just after ending the event loop  ------------
    void 
FitEpsilonPlot::endJob() 
{

  if (foldInSuperModule_ and makeFoldedHistograms_) {
    cout << "FIT_EPSILON: not doing anything inside endJob(): we only wanted to fold histograms" << endl;
    return;
  }

  if (isEoverEtrue_) {
    // call it first with false to save first photon coefficients
    if (fitEoverEtrueWithRooFit_) {
      saveCoefficientsEoverEtrueRooFit(false);
      saveCoefficientsEoverEtrueRooFit(true);
    } else {
      saveCoefficientsEoverEtrue(false);
      saveCoefficientsEoverEtrue(true);
    }
  } else {
    saveCoefficients();
  }

  if(StoreForTest_){
    outfileTEST_->Write();
    outfileTEST_->Close();
    cout << "FIT_EPSILON: Fit stored in " << fitFileName_ << endl;
  }

}

// ------------ method called when starting to processes a run  ------------
    void 
FitEpsilonPlot::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
    void 
FitEpsilonPlot::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
    void 
FitEpsilonPlot::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
    void 
FitEpsilonPlot::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FitEpsilonPlot::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FitEpsilonPlot);
