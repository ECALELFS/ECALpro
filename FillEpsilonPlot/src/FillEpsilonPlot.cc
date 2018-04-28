// -*- C++ -*-
//
// Package:    FillEpsilonPlot
// Class:      FillEpsilonPlot
// 
/**\class FillEpsilonPlot FillEpsilonPlot.cc CalibCode/FillEpsilonPlot/src/FillEpsilonPlot.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Marco Grassi, CMS
//         Created:  Tue Sep 27 15:07:49 CEST 2011
// $Id: FillEpsilonPlot.cc,v 1.14 2013/06/19 13:42:19 lpernie Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <map>
#include <algorithm>

// user include files
#include "TFile.h"
#include "TRegexp.h"
//#include "TStopwatch.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapHardcodedTopology.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"

#include "CalibCode/FillEpsilonPlot/interface/FillEpsilonPlot.h"
#include "CalibCode/CalibTools/interface/GlobalFunctions.h"
#include "CalibCode/CalibTools/interface/EcalRecHitCompare.h"
#include "CalibCode/CalibTools/interface/PreshowerTools.h"
#include "CalibCode/CalibTools/interface/GeometryService.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
//Geom
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
//ES
#include "FWCore/Framework/interface/ESProducer.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
//HLT
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>
// for L1
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h" // included to get L1 info
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"

using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::max;

#include "CalibCode/FillEpsilonPlot/interface/JSON.h"
//MVA Stuff
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif
#include "CalibCode/GBRTrain/interface/GBRApply.h"
#include "CalibCode/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"

//#include "Cintex/Cintex.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#define DR_FOR_UNMERGED_GEN_PHOTONS 0.025 // if two gen photons are closer than this value, they will not be used for the gen-reco matching, because they are too close to be distinguished by the reco clustering algorithm (0.0175 in Dphi or Deta is ~1 ECAL cystal and the seeds must be farther than 1 crystal also on the diagonal)

using namespace TMVA;
using namespace edm;

//Function
double max_array(double *A, int n);
double max(double x, double y);
int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3);


FillEpsilonPlot::FillEpsilonPlot(const edm::ParameterSet& iConfig)
{

    /// parameters from python
    Are_pi0_                           = iConfig.getUntrackedParameter<bool>("Are_pi0",true);
    useContainmentCorrectionsFromEoverEtrue_ = iConfig.getUntrackedParameter<bool>("useContainmentCorrectionsFromEoverEtrue",true);
    useMVAContainmentCorrections_      = iConfig.getUntrackedParameter<bool>("useMVAContainmentCorrections",true);
    new_pi0ContainmentCorrections_     = iConfig.getUntrackedParameter<bool>("new_pi0ContainmentCorrections",false);

    EBRecHitCollectionToken_           = consumes<EBRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag"));
    EERecHitCollectionToken_           = consumes<EERecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag"));
    ESRecHitCollectionToken_           = consumes<ESRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag"));
    HLTResults_                        = iConfig.getUntrackedParameter<bool>("HLTResults",false);
    HLTResultsNameEB_                  = iConfig.getUntrackedParameter<std::string>("HLTResultsNameEB","AlCa_EcalPi0EB");
    HLTResultsNameEE_                  = iConfig.getUntrackedParameter<std::string>("HLTResultsNameEE","AlCa_EcalPi0EE");
    RemoveSeedsCloseToDeadXtal_        = iConfig.getUntrackedParameter<bool>("RemoveSeedsCloseToDeadXtal",false);
    RemoveDead_Flag_                   = iConfig.getUntrackedParameter<bool>("RemoveDead_Flag",false);
    RemoveDead_Map_                    = iConfig.getUntrackedParameter<std::string>("RemoveDead_Map");
    L1_Bit_Sele_                       = iConfig.getUntrackedParameter<std::string>("L1_Bit_Sele","");
    L1TriggerInfo_                     = iConfig.getUntrackedParameter<bool>("L1TriggerInfo",false);
    triggerResultsToken_               = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter("triggerTag",edm::InputTag("TriggerResults")));
    //L1GTobjmapToken_                   = consumes<GlobalAlgBlkBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("hltGtStage2Digis",edm::InputTag("hltGtStage2Digis"))); // for MC should use "gtStage2Digis",edm::InputTag("gtStage2Digis","","RECO")    
    L1GTobjmapToken_                   = consumes<GlobalAlgBlkBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("L1GTobjmapTag",edm::InputTag("hltGtStage2Digis"))); 
    GenPartCollectionToken_            = consumes<GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("GenPartCollectionTag",edm::InputTag("genParticles")));
    outfilename_                       = iConfig.getUntrackedParameter<std::string>("OutputFile");
    fileEoverEtrueContainmentCorrections_ = iConfig.getUntrackedParameter<std::string>("fileEoverEtrueContainmentCorrections");
    ebContainmentCorrections_          = iConfig.getUntrackedParameter<std::string>("EBContainmentCorrections");
    MVAEBContainmentCorrections_01_    = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_01");
    MVAEBContainmentCorrections_02_    = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_02");
    MVAEEContainmentCorrections_01_    = iConfig.getUntrackedParameter<std::string>("MVAEEContainmentCorrections_01");
    MVAEEContainmentCorrections_02_    = iConfig.getUntrackedParameter<std::string>("MVAEEContainmentCorrections_02");
    MVAEBContainmentCorrections_eta01_ = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_eta01");
    MVAEBContainmentCorrections_eta02_ = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_eta02");
    Endc_x_y_                          = iConfig.getUntrackedParameter<std::string>("Endc_x_y");
    EtaRingCalibEB_                    = iConfig.getUntrackedParameter<bool>("EtaRingCalibEB",false);
    SMCalibEB_                         = iConfig.getUntrackedParameter<bool>("SMCalibEB",false);
    EtaRingCalibEE_                    = iConfig.getUntrackedParameter<bool>("EtaRingCalibEE",false);
    SMCalibEE_                         = iConfig.getUntrackedParameter<bool>("SMCalibEE",false);
    CalibMapEtaRing_                   = iConfig.getUntrackedParameter<std::string>("CalibMapEtaRing","CalibCode/FillEpsilonPlot/data/calibMap.root");
    ebPHIContainmentCorrections_       = iConfig.getUntrackedParameter<std::string>("EBPHIContainmentCorrections");
    eeContainmentCorrections_          = iConfig.getUntrackedParameter<std::string>("EEContainmentCorrections");
    useEBContainmentCorrections_       = iConfig.getUntrackedParameter<bool>("useEBContainmentCorrections");
    useEEContainmentCorrections_       = iConfig.getUntrackedParameter<bool>("useEEContainmentCorrections");
    externalGeometry_                  = iConfig.getUntrackedParameter<std::string>("ExternalGeometry");
    currentIteration_                  = iConfig.getUntrackedParameter<int>("CurrentIteration");
    outputDir_                         = iConfig.getUntrackedParameter<std::string>("OutputDir");
    isCRAB_                            = iConfig.getUntrackedParameter<bool>("isCRAB",false);
    calibMapPath_                      = iConfig.getUntrackedParameter<std::string>("calibMapPath");
    Barrel_orEndcap_                   = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");
    EB_Seed_E_                         = iConfig.getUntrackedParameter<double>("EB_Seed_E",0.2);
    useEE_EtSeed_                      = iConfig.getUntrackedParameter<bool>("useEE_EtSeed",true);
    EE_Seed_E_                         = iConfig.getUntrackedParameter<double>("EE_Seed_E",0.5);
    EE_Seed_Et_                        = iConfig.getUntrackedParameter<double>("EE_Seed_Et",0.5);
    pi0PtCut_low_[EcalBarrel]          = iConfig.getUntrackedParameter<double>("Pi0PtCutEB_low");
    pi0PtCut_high_[EcalBarrel]         = iConfig.getUntrackedParameter<double>("Pi0PtCutEB_high");
    pi0PtCut_low_[EcalEndcap]          = iConfig.getUntrackedParameter<double>("Pi0PtCutEE_low");
    pi0PtCut_high_[EcalEndcap]         = iConfig.getUntrackedParameter<double>("Pi0PtCutEE_high");
    gPtCut_low_[EcalBarrel]            = iConfig.getUntrackedParameter<double>("gPtCutEB_low");
    gPtCut_high_[EcalBarrel]           = iConfig.getUntrackedParameter<double>("gPtCutEB_high");
    gPtCut_low_[EcalEndcap]            = iConfig.getUntrackedParameter<double>("gPtCutEE_low");
    gPtCut_high_[EcalEndcap]           = iConfig.getUntrackedParameter<double>("gPtCutEE_high");
    CutOnHLTIso_                       = iConfig.getUntrackedParameter<bool>("CutOnHLTIso",false);
    pi0HLTIsoCut_low_[EcalBarrel]      = iConfig.getUntrackedParameter<double>("Pi0HLTIsoCutEB_low");
    pi0HLTIsoCut_high_[EcalBarrel]     = iConfig.getUntrackedParameter<double>("Pi0HLTIsoCutEB_high");
    pi0HLTIsoCut_low_[EcalEndcap]      = iConfig.getUntrackedParameter<double>("Pi0HLTIsoCutEE_low");
    pi0HLTIsoCut_high_[EcalEndcap]     = iConfig.getUntrackedParameter<double>("Pi0HLTIsoCutEE_high");
    pi0IsoCut_low_[EcalBarrel]         = iConfig.getUntrackedParameter<double>("Pi0IsoCutEB_low");
    pi0IsoCut_high_[EcalBarrel]        = iConfig.getUntrackedParameter<double>("Pi0IsoCutEB_high");
    pi0IsoCut_low_[EcalEndcap]         = iConfig.getUntrackedParameter<double>("Pi0IsoCutEE_low");
    pi0IsoCut_high_[EcalEndcap]        = iConfig.getUntrackedParameter<double>("Pi0IsoCutEE_high");
    nXtal_1_cut_low_[EcalEndcap]       = iConfig.getUntrackedParameter<int>("nXtal_1_EE_low");
    nXtal_1_cut_high_[EcalEndcap]      = iConfig.getUntrackedParameter<int>("nXtal_1_EE_high");
    nXtal_2_cut_low_[EcalEndcap]       = iConfig.getUntrackedParameter<int>("nXtal_2_EE_low");
    nXtal_2_cut_high_[EcalEndcap]      = iConfig.getUntrackedParameter<int>("nXtal_2_EE_high");
    nXtal_1_cut_low_[EcalBarrel]       = iConfig.getUntrackedParameter<int>("nXtal_1_EB_low");
    nXtal_1_cut_high_[EcalBarrel]      = iConfig.getUntrackedParameter<int>("nXtal_1_EB_high");
    nXtal_2_cut_low_[EcalBarrel]       = iConfig.getUntrackedParameter<int>("nXtal_2_EB_low");
    nXtal_2_cut_high_[EcalBarrel]      = iConfig.getUntrackedParameter<int>("nXtal_2_EB_high");
    S4S9_cut_low_[EcalBarrel]          = iConfig.getUntrackedParameter<double>("S4S9_EB_low");
    S4S9_cut_high_[EcalBarrel]         = iConfig.getUntrackedParameter<double>("S4S9_EB_high");
    S4S9_cut_low_[EcalEndcap]          = iConfig.getUntrackedParameter<double>("S4S9_EE_low");
    S4S9_cut_high_[EcalEndcap]         = iConfig.getUntrackedParameter<double>("S4S9_EE_high");
    SystOrNot_                         = iConfig.getUntrackedParameter<int>("SystOrNot",0);
    useMassInsteadOfEpsilon_           = iConfig.getUntrackedParameter<bool>("useMassInsteadOfEpsilon",true);
    isMC_                              = iConfig.getUntrackedParameter<bool>("isMC",false);
    MC_Assoc_                          = iConfig.getUntrackedParameter<bool>("MC_Assoc",false);
    MC_Assoc_DeltaR                    = iConfig.getUntrackedParameter<double>("MC_Assoc_DeltaR",0.1);
    MakeNtuple4optimization_           = iConfig.getUntrackedParameter<bool>("MakeNtuple4optimization",false);
    GeometryFromFile_                  = iConfig.getUntrackedParameter<bool>("GeometryFromFile",false);
    JSONfile_                          = iConfig.getUntrackedParameter<std::string>("JSONfile","");

    L1SeedsPi0Stream_                  = iConfig.getUntrackedParameter<std::string>("L1SeedsPi0Stream");    
    nL1SeedsPi0Stream_                 = iConfig.getUntrackedParameter<int>("nL1SeedsPi0Stream");

    isDebug_                           = iConfig.getUntrackedParameter<bool>("isDebug",false);
    isEoverEtrue_                      = iConfig.getUntrackedParameter<bool>("isEoverEtrue",false);
    //pileupSummaryToken_                = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummaryTag",edm::InputTag("addPileupInfo")));
    if (isMC_)
      pileupSummaryToken_                = consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSummaryTag"));
    fillKinematicVariables_            = iConfig.getUntrackedParameter<bool>("fillKinematicVariables",false);

    // for MC-truth association
    // g4_simTk_Token_  = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
    // g4_simVtx_Token_ = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));

    if(useEE_EtSeed_) cout<<"SEEDS Used: EB "<<EB_Seed_E_<<" and EE "<<EE_Seed_Et_<<" (in Et) "<<endl;
    else              cout<<"SEEDS Used: EB "<<EB_Seed_E_<<" and EE "<<EE_Seed_E_<<" (in E) "<<endl;
    cout<<"Cut used: EB LOW)"<<endl;
    cout<<"Pt(pi0): "<<pi0PtCut_low_[EcalBarrel]
	<<", Pt(Clus): "<<gPtCut_low_[EcalBarrel]
	<<", Iso: "<<pi0IsoCut_low_[EcalBarrel]
	<<", Nxtal_1: "<<nXtal_1_cut_low_[EcalBarrel]
	<<", Nxtal_2: "<<nXtal_2_cut_low_[EcalBarrel]
	<<", S4S9: "<<S4S9_cut_low_[EcalBarrel]<<endl;
    cout<<"Cut used: EB HIGH)"<<endl;
    cout<<"Pt(pi0): "<<pi0PtCut_high_[EcalBarrel]
	<<", Pt(Clus): "<<gPtCut_high_[EcalBarrel]
	<<", Iso: "<<pi0IsoCut_high_[EcalBarrel]
	<<", Nxtal_1: "<<nXtal_1_cut_high_[EcalBarrel]
	<<", Nxtal_2: "<<nXtal_2_cut_high_[EcalBarrel]
	<<", S4S9: "<<S4S9_cut_high_[EcalBarrel]<<endl;
    cout<<"Cut used: EE LOW)"<<endl;
    cout<<"Pt(pi0): "<<pi0PtCut_low_[EcalEndcap]
	<<", Pt(Clus): "<<gPtCut_low_[EcalEndcap]
	<<", Iso: "<<pi0IsoCut_low_[EcalEndcap]
	<<", Nxtal_1: "<<nXtal_1_cut_low_[EcalEndcap]
	<<", Nxtal_2: "<<nXtal_2_cut_low_[EcalEndcap]
	<<", S4S9: "<<S4S9_cut_low_[EcalEndcap]<<endl;
    cout<<"Cut used: EE HIGH)"<<endl;
    cout<<"Pt(pi0): "<<pi0PtCut_high_[EcalEndcap]
	<<", Pt(Clus): "<<gPtCut_high_[EcalEndcap]
	<<", Iso: "<<pi0IsoCut_high_[EcalEndcap]
	<<", Nxtal_1: "<<nXtal_1_cut_high_[EcalEndcap]
	<<", Nxtal_2: "<<nXtal_2_cut_high_[EcalEndcap]
	<<", S4S9: "<<S4S9_cut_high_[EcalEndcap]<<endl;
    cout<<"The StatError option choose is: "<<SystOrNot_<<" [0= No error stat computation, 1 = yes only even events, 2 = yes only odd events]"<<endl;

    useOnlyEEClusterMatchedWithES_ = iConfig.getUntrackedParameter<bool>("useOnlyEEClusterMatchedWithES"); 
    //JSON
    if( JSONfile_!="" ) {
      if (JSONfile_.find("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/") != std::string::npos)
	// read from absolute /afs/ path
	myjson=new JSON( JSONfile_.c_str() );
      else 
	myjson=new JSON( edm::FileInPath( JSONfile_.c_str() ).fullPath().c_str() );      
    }
    // shower shape parameters
    PCparams_.param_LogWeighted_ = true;
    PCparams_.param_T0_barl_     = 7.4;
    PCparams_.param_T0_endc_     = 3.1;
    PCparams_.param_T0_endcES_   = 1.2;
    PCparams_.param_W0_          = 4.2;
    PCparams_.param_X0_          = 0.89;

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

    cout << "crosscheck: selected type: " << regionalCalibration_->printType() << endl;

    TH1::SetDefaultSumw2(); // all new histograms will automatically activate the storage of the sum of squares of errors (i.e, TH1::Sumw2 is automatically called).

    /// external hardcoded geometry

    externalGeometryFile_ = TFile::Open( edm::FileInPath( externalGeometry_.c_str() ).fullPath().c_str() );
    if(!externalGeometryFile_) cms::Exception("ExtGeom") << "External Geometry file (" << externalGeometry_ << ") not found" << endl;
    geom_ = ECALGeometry::getGeometry(externalGeometryFile_);
    GeometryService::setGeometryName(externalGeometry_);
    GeometryService::setGeometryPtr(geom_);

    // containment corrections
    if (useContainmentCorrectionsFromEoverEtrue_) loadEoverEtrueContainmentCorrections(fileEoverEtrueContainmentCorrections_);

#if (defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)) || defined(REGRESS_AND_PARAM_CONTCORR)
    if(useEEContainmentCorrections_)
	  containmentCorrections_.loadContainmentPointCorrectionsEE( edm::FileInPath( eeContainmentCorrections_.c_str() ).fullPath().c_str() );
    if(useEBContainmentCorrections_){
	  containmentCorrections_.loadContainmentCorrectionsEB( edm::FileInPath( ebContainmentCorrections_.c_str() ).fullPath().c_str());
	  EBPHI_Cont_Corr_load( edm::FileInPath( ebPHIContainmentCorrections_.c_str() ).fullPath() );
    }
#endif
    /// subdetector topology
    ebtopology_ = new CaloTopology();
    EcalBarrelHardcodedTopology* ebHTopology = new EcalBarrelHardcodedTopology();
    ebtopology_->setSubdetTopology(DetId::Ecal,EcalBarrel,ebHTopology);

    eetopology_ = new CaloTopology();  
    EcalEndcapHardcodedTopology* eeHTopology=new EcalEndcapHardcodedTopology();
    eetopology_->setSubdetTopology(DetId::Ecal,EcalEndcap,eeHTopology);

    /// retrieving calibration coefficients of the previous iteration
    // if currentIteration_ = 0, calibMapPath_ contains "iter_-1" unless the current set of ICs was started from another existing set (see parameters.py)
    // therefore, the case with extension is included below
    std::string stringToMatch = "iter_-1";  // used below: this string should not match to trigger true condition
    if (currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
    else if (currentIteration_ > 0 || (currentIteration_ == 0 && calibMapPath_.find(stringToMatch)==std::string::npos))
    {
      cout << "FillEpsilonPlot:: loading calibration map at " << calibMapPath_ << endl;
      regionalCalibration_->getCalibMap()->loadCalibMapFromFile(calibMapPath_.c_str(),false);
      if (isEoverEtrue_) regionalCalibration_g2_->getCalibMap()->loadCalibMapFromFile(calibMapPath_.c_str(),true);
    }

    /// epsilon histograms
    if(!MakeNtuple4optimization_){

      if (isEoverEtrue_) {

	if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  {
	  EoverEtrue_g1_EB_h = initializeEpsilonHistograms("EoverEtrue_g1_EB_iR_","reco/gen #gamma1 energy EB - iR", regionalCalibration_->getCalibMap()->getNRegionsEB() );
	  EoverEtrue_g2_EB_h = initializeEpsilonHistograms("EoverEtrue_g2_EB_iR_","reco/gen #gamma2 energy EB - iR", regionalCalibration_g2_->getCalibMap()->getNRegionsEB() );
	}
	if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  {
	  EoverEtrue_g1_EE_h = initializeEpsilonHistograms("EoverEtrue_g1_EE_iR_","reco/gen #gamma1 energy - iR", regionalCalibration_->getCalibMap()->getNRegionsEE() );
	  EoverEtrue_g2_EE_h = initializeEpsilonHistograms("EoverEtrue_g2_EE_iR_","reco/gen #gamma2 energy - iR", regionalCalibration_g2_->getCalibMap()->getNRegionsEE() );
	}

      } else {

	if(useMassInsteadOfEpsilon_ ) {

	  if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  
	    epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","#pi^{0} Mass distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
	  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  
	    epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","#pi^{0} Mass distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
	
	} else {
	
	  if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  
	    epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","Epsilon distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
	  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  
	    epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","Epsilon distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
	
	}

      }

    }

    EventFlow_EB  = new TH1F("EventFlow_EB", "EventFlow EB", 6, -0.5, 5.5 );
    EventFlow_EB->GetXaxis()->SetBinLabel(1,"All Events"); EventFlow_EB->GetXaxis()->SetBinLabel(2,"JSON"); EventFlow_EB->GetXaxis()->SetBinLabel(3,"Trigger Res");
    EventFlow_EB->GetXaxis()->SetBinLabel(4,"HLT"); EventFlow_EB->GetXaxis()->SetBinLabel(5,"Initial Comb."); EventFlow_EB->GetXaxis()->SetBinLabel(6,"Final Comb.");
    EventFlow_EE  = new TH1F("EventFlow_EE", "EventFlow EE", 6, -0.5, 5.5 );
    EventFlow_EE->GetXaxis()->SetBinLabel(1,"All Events"); EventFlow_EE->GetXaxis()->SetBinLabel(2,"JSON"); EventFlow_EE->GetXaxis()->SetBinLabel(3,"Trigger Res");
    EventFlow_EE->GetXaxis()->SetBinLabel(4,"HLT"); EventFlow_EE->GetXaxis()->SetBinLabel(5,"Initial Comb."); EventFlow_EE->GetXaxis()->SetBinLabel(6,"Final Comb.");
    if (isDebug_) {
      EventFlow_EB_debug  = new TH1F("EventFlow_EB_debug", "EventFlow EB", 5, -0.5, 4.5 );
      EventFlow_EB_debug->GetXaxis()->SetBinLabel(1,"Initial Comb.");
      EventFlow_EB_debug->GetXaxis()->SetBinLabel(2,"pi0pt"); 
      EventFlow_EB_debug->GetXaxis()->SetBinLabel(3,"isocut"); 
      EventFlow_EB_debug->GetXaxis()->SetBinLabel(4,"hltiso");
      EventFlow_EB_debug->GetXaxis()->SetBinLabel(5,"nxtal");
      EventFlow_EE_debug  = new TH1F("EventFlow_EE_debug", "EventFlow EE", 5, -0.5, 4.5 );
      EventFlow_EE_debug->GetXaxis()->SetBinLabel(1,"Initial Comb.");
      EventFlow_EE_debug->GetXaxis()->SetBinLabel(2,"pi0pt"); 
      EventFlow_EE_debug->GetXaxis()->SetBinLabel(3,"isocut"); 
      EventFlow_EE_debug->GetXaxis()->SetBinLabel(4,"hltiso");
      EventFlow_EE_debug->GetXaxis()->SetBinLabel(5,"nxtal");
    }

    if (isEoverEtrue_) {
      allEoverEtrue_g1_EB = new TH1F("allEoverEtrue_g1_EB", "allEoverEtrue_g1_EB", 150, 0.0, 1.5);
      allEoverEtrue_g1_EBnw = new TH1F("allEoverEtrue_g1_EBnw", "allEoverEtrue_g1_EBnw", 150, 0.0, 1.5);
      allEoverEtrue_g1_EE = new TH1F("allEoverEtrue_g1_EE", "allEoverEtrue_g1_EE", 150, 0.0, 1.5);
      allEoverEtrue_g1_EEnw = new TH1F("allEoverEtrue_g1_EEnw", "allEoverEtrue_g1_EEnw", 150, 0.0, 1.5);
      allEoverEtrue_g2_EB = new TH1F("allEoverEtrue_g2_EB", "allEoverEtrue_g2_EB", 150, 0.0, 1.5);
      allEoverEtrue_g2_EBnw = new TH1F("allEoverEtrue_g2_EBnw", "allEoverEtrue_g2_EBnw", 150, 0.0, 1.5);
      allEoverEtrue_g2_EE = new TH1F("allEoverEtrue_g2_EE", "allEoverEtrue_g2_EE", 150, 0.0, 1.5);
      allEoverEtrue_g2_EEnw = new TH1F("allEoverEtrue_g2_EEnw", "allEoverEtrue_g2_EEnw", 150, 0.0, 1.5);
    } else { 
      allEpsilon_EB = new TH1F("allEpsilon_EB", "allEpsilon_EB",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
      allEpsilon_EBnw = new TH1F("allEpsilon_EBnw", "allEpsilon_EBnw",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
      allEpsilon_EE = new TH1F("allEpsilon_EE", "allEpsilon_EE",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
      allEpsilon_EEnw = new TH1F("allEpsilon_EEnw", "allEpsilon_EEnw",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
    }

    if (fillKinematicVariables_) {

      // 4 regions: 2 in EB, 2 in EE
      std::vector< std::string> regionStreamPi0; // will contain a name for the region (e.g.: region1EB, region2EE, ecc...)
      regionStreamPi0.clear();
      regionStreamPi0.push_back("region1EB");
      regionStreamPi0.push_back("region2EB");
      regionStreamPi0.push_back("region1EE");
      regionStreamPi0.push_back("region2EE");

      for (int iregEcal = 0; iregEcal < 4; ++iregEcal) {

	pi0pt_afterCuts.push_back( new TH1F(Form("pi0pt_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"#pi^{0} p_{T} after cuts",60,0.0,15.0) );
	g1pt_afterCuts.push_back( new TH1F(Form("g1pt_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"leading (seed) #gamma p_{T} after cuts",60,0.0,10.0) );
	g2pt_afterCuts.push_back( new TH1F(Form("g2pt_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"trailing (seed) #gamma p_{T} after cuts",60,0.0,10.0) );
	g1Nxtal_afterCuts.push_back( new TH1F(Form("g1Nxtal_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"leading (seed) #gamma number of crystals after cuts",9,0.5,9.5) );
	g2Nxtal_afterCuts.push_back( new TH1F(Form("g2Nxtal_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"trailing (seed) #gamma number of crystals after cuts",9,0.5,9.5) );
	pi0PhotonsNoverlappingXtals_afterCuts.push_back( new TH1F(Form("pi0PhotonsNoverlappingXtals_afterCuts_%s",regionStreamPi0[iregEcal].c_str()),"number of overlapping crystals in #pi^{0}->#gamma#gamma after cuts",10,-0.5,9.5) );
	if (isMC_) {
	  pi0MassVsPU.push_back( new TH2F(Form("pi0MassVsPU_%s",regionStreamPi0[iregEcal].c_str()),"#pi^{0} mass vs PU",100,0.05,0.25,50,0.5,50.5) );
	}

      }

      regionStreamPi0.clear();

    }

    entries_EEp   = new TH2F("entries_EEp","entries_EEp",101,-0.5,100.5,101,-0.5,100.5);
    entries_EEm   = new TH2F("entries_EEm","entries_EEm",101,-0.5,100.5,101,-0.5,100.5);
    entries_EB    = new TH2F("entries_EB","entries_EB",2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
    Occupancy_EEp = new TH2F("Occupancy_EEp","Occupancy_EEp",101,-0.5,100.5,101,-0.5,100.5);
    Occupancy_EEm = new TH2F("Occupancy_EEm","Occupancy_EEm",101,-0.5,100.5,101,-0.5,100.5);
    Occupancy_EB  = new TH2F("Occupancy_EB","Occupancy_EB",2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );


    pi0MassVsIetaEB = new TH2F("pi0MassVsIetaEB","#pi^{0} mass vs i#eta",85,0.5,85.5,120,Are_pi0_? 0.:0.3, Are_pi0_? 0.3:0.8);
    pi0MassVsIetaEB->GetXaxis()->SetTitle("i#eta");
    pi0MassVsIetaEB->GetYaxis()->SetTitle("#pi^{0} mass");
    pi0MassVsETEB = new TH2F("pi0MassVsETEB", "#pi^{0} mass vs E_{T}(pi^{0})",120,0.,20.,120,Are_pi0_? 0.:0.3, Are_pi0_? 0.3:0.8);
    pi0MassVsETEB->GetXaxis()->SetTitle("E_{T}(pi^{0})");
    pi0MassVsETEB->GetYaxis()->SetTitle("#pi^{0} mass");

    //DeadXtal from Map
    if( RemoveDead_Map_!="" ){
      DeadMap         = TFile::Open( RemoveDead_Map_.Data() );
      EBMap_DeadXtal  = (TH2F*) DeadMap->Get("rms_EB");
      EEmMap_DeadXtal = (TH2F*) DeadMap->Get("rms_EEm");
      EEpMap_DeadXtal = (TH2F*) DeadMap->Get("rms_EEp");
    }
    // output file
    string fileName = "";
    fileName = outputDir_ + outfilename_;
    outfile_ = new TFile(fileName.c_str(),"RECREATE");
    if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << fileName << "\n";
#ifdef SELECTION_TREE
    CutVariables_EB = new TTree("CutVariables_EB","(EB) Variables used at first cuts");
    CutVariables_EB->Branch("NSeeds_EB", &NSeeds_EB, "NSeeds_EB/F");
    CutVariables_EB->Branch("Xclus_EB", &Xclus_EB, "Xclus_EB/F");
    CutVariables_EB->Branch("Yclus_EB", &Yclus_EB, "Yclus_EB/F");
    CutVariables_EB->Branch("Zclus_EB", &Zclus_EB, "Zclus_EB/F");
    CutVariables_EB->Branch("e3x3_EB", &e3x3_EB, "e3x3_EB/F");
    CutVariables_EB->Branch("S4S9_EB", &S4S9_EB, "S4S9_EB/F");
    CutVariables_EB->Branch("PTClus_EB", &PTClus_EB, "PTClus_EB/F");
    CutVariables_EE = new TTree("CutVariables_EE","(EE) Variables used at first cuts");
    CutVariables_EE->Branch("NSeeds_EE", &NSeeds_EB, "NSeeds_EB/F");
    CutVariables_EE->Branch("Xclus_EE", &Xclus_EE, "Xclus_EE/F");
    CutVariables_EE->Branch("Yclus_EE", &Yclus_EE, "Yclus_EE/F");
    CutVariables_EE->Branch("Zclus_EE", &Zclus_EE, "Zclus_EE/F");
    CutVariables_EE->Branch("e3x3_EE", &e3x3_EE, "e3x3_EE/F");
    CutVariables_EE->Branch("S4S9_EE", &S4S9_EE, "S4S9_EE/F");
    CutVariables_EE->Branch("PTClus_EE", &PTClus_EE, "PTClus_EE/F");

    Pi0Info_EB= new TTree("Pi0Info_EB","(EB) Pi0 informations");
    Pi0Info_EB->Branch("PtPi0_EB", &PtPi0_EB, "PtPi0_EB/F");
    Pi0Info_EB->Branch("mpi0_EB", &mpi0_EB, "mpi0_EB/F");
    Pi0Info_EB->Branch("Etapi0_EB", &Etapi0_EB, "Etapi0_EB/F");
    Pi0Info_EB->Branch("Phipi0_EB", &Phipi0_EB, "Phipi0_EB/F");
    Pi0Info_EB->Branch("Epsilon_EB", &Epsilon_EB, "Epsilon_EB/F");
    Pi0Info_EE= new TTree("Pi0Info_EE","(EE) Pi0 informations");
    Pi0Info_EE->Branch("PtPi0_EE", &PtPi0_EE, "PtPi0_EE/F");
    Pi0Info_EE->Branch("mpi0_EE", &mpi0_EE, "mpi0_EE/F");
    Pi0Info_EE->Branch("Etapi0_EE", &Etapi0_EE, "Etapi0_EE/F");
    Pi0Info_EE->Branch("Phipi0_EE", &Phipi0_EE, "Phipi0_EE/F");
    Pi0Info_EE->Branch("Epsilon_EE", &Epsilon_EE, "Epsilon_EE/F");
#endif
    if(MakeNtuple4optimization_){
	Tree_Optim = new TTree("Tree_Optim","Output TTree");
	// event info for data
	Tree_Optim->Branch( "Event",     &myEvent,     "Event/l"); // l is for ULong64_t
	Tree_Optim->Branch( "LumiBlock", &myLumiBlock, "LumiBlock/I");
	Tree_Optim->Branch( "Run",       &myRun,       "Run/I");
	Tree_Optim->Branch( "BunchCrossing",       &myBunchCrossing,       "BunchCrossing/I");
	if (HLTResults_) {
	  if (Are_pi0_) {
	    Tree_Optim->Branch( "AlCa_EcalPi0EBonly", &EB_HLT, "AlCa_EcalPi0EBonly/O"); // O (capital letter o, not zero) is for a Bool_t
	    Tree_Optim->Branch( "AlCa_EcalPi0EEonly", &EE_HLT, "AlCa_EcalPi0EEonly/O");	  
	  } else {
	    Tree_Optim->Branch( "AlCa_EcalEtaEBonly", &EB_HLT, "AlCa_EcalEtaEBonly/O"); 
	    Tree_Optim->Branch( "AlCa_EcalEtaEEonly", &EE_HLT, "AlCa_EcalEtaEEonly/O");	  
	  }
	}
	//Tree_Optim->Branch( "STr2_L1Seed",        &Op_L1Seed,           Form("STr2_L1Seed[%d]/I",NL1SEED));
	Tree_Optim->Branch( "nPi0",          &Op_NPi0,             "nPi0/I");
	Tree_Optim->Branch( "Pi0recIsEB",    &Op_Pi0recIsEB,       "Pi0recIsEB/I");
	Tree_Optim->Branch( "ClusIsoPi0",    &Op_ClusIsoPi0,       "ClusIsoPi0/F");
	Tree_Optim->Branch( "HLTIsoPi0",     &Op_HLTIsoPi0,    "HLTIsoPi0/F");
	Tree_Optim->Branch( "n1CrisPi0",     &Op_nCrisG1,    "nCrisG1Pi0/I");
	Tree_Optim->Branch( "n2CrisPi0",     &Op_nCrisG2,    "nCrisG2Pi0/I");
	Tree_Optim->Branch( "mPi0_cor",      &Op_mPi0_cor,         "mPi0_cor/F");
	Tree_Optim->Branch( "enG1_cor",      &Op_enG1_cor,         "enG1_cor/F");
	Tree_Optim->Branch( "enG2_cor",      &Op_enG2_cor,         "enG2_cor/F");
	Tree_Optim->Branch( "etaG1_cor",     &Op_etaG1_cor,         "etaG1_cor/F");
	Tree_Optim->Branch( "etaG2_cor",     &Op_etaG2_cor,         "etaG2_cor/F");
	Tree_Optim->Branch( "phiG1_cor",     &Op_phiG1_cor,         "phiG1_cor/F");
	Tree_Optim->Branch( "phiG2_cor",     &Op_phiG2_cor,         "phiG2_cor/F");
	Tree_Optim->Branch( "etaPi0_cor",    &Op_etaPi0_cor,       "etaPi0_cor/F");
	Tree_Optim->Branch( "ptPi0_cor",     &Op_ptPi0_cor,        "ptPi0_cor/F");
	Tree_Optim->Branch( "ptPi0_nocor",   &Op_ptPi0_nocor,      "ptPi0_nocor/F");
	Tree_Optim->Branch( "enG1_nocor",    &Op_enG1_nocor,       "enG1_nocor/F");
	Tree_Optim->Branch( "enG2_nocor",    &Op_enG2_nocor,       "enG2_nocor/F");
	Tree_Optim->Branch( "etaG1_nocor",   &Op_etaG1_nocor,       "etaG1_nocor/F");
	Tree_Optim->Branch( "etaG2_nocor",   &Op_etaG2_nocor,       "etaG2_nocor/F");
	Tree_Optim->Branch( "phiG1_nocor",   &Op_phiG1_nocor,       "phiG1_nocor/F");
	Tree_Optim->Branch( "phiG2_nocor",   &Op_phiG2_nocor,       "phiG2_nocor/F");
	Tree_Optim->Branch( "mPi0_nocor",    &Op_mPi0_nocor,       "mPi0_nocor/F");
	Tree_Optim->Branch( "DeltaRG1G2",    &Op_DeltaRG1G2,       "DeltaRG1G2/F");
	Tree_Optim->Branch( "Es_e1_1",       &Op_Es_e1_1,          "Es_e1_1/F");
	Tree_Optim->Branch( "Es_e1_2",       &Op_Es_e1_2,          "Es_e1_2/F");
	Tree_Optim->Branch( "Es_e2_1",       &Op_Es_e2_1,          "Es_e2_1/F");
	Tree_Optim->Branch( "Es_e2_2",       &Op_Es_e2_2,          "Es_e2_2/F");
	Tree_Optim->Branch( "S4S9_1",        &Op_S4S9_1,           "S4S9_1/F");
	Tree_Optim->Branch( "S4S9_2",        &Op_S4S9_2,           "S4S9_2/F");
	Tree_Optim->Branch( "S2S9_1",        &Op_S2S9_1,           "S2S9_1/F");
	Tree_Optim->Branch( "S2S9_2",        &Op_S2S9_2,           "S2S9_2/F");
	Tree_Optim->Branch( "S1S9_1",        &Op_S1S9_1,           "S1S9_1/F");
	Tree_Optim->Branch( "S1S9_2",        &Op_S1S9_2,           "S1S9_2/F");
	Tree_Optim->Branch( "Time_1",        &Op_Time_1,           "Time_1/F");
	Tree_Optim->Branch( "Time_2",        &Op_Time_2,           "Time_2/F");
	Tree_Optim->Branch( "iEtaiX_1",      &Op_iEtaiX_1,         "iEtaiX_1/I");
	Tree_Optim->Branch( "iEtaiX_2",      &Op_iEtaiX_2,         "iEtaiX_2/I");
	Tree_Optim->Branch( "iPhiiY_1",      &Op_iPhiiY_1,         "iPhiiY_1/I");
	Tree_Optim->Branch( "iPhiiY_2",      &Op_iPhiiY_2,         "iPhiiY_2/I");
	Tree_Optim->Branch( "iEta_1on5",     &Op_iEta_1on5,        "iEta_1on5/I");
	Tree_Optim->Branch( "iEta_2on5",     &Op_iEta_2on5,        "iEta_2on5/I");
	Tree_Optim->Branch( "iPhi_1on2",     &Op_iPhi_1on2,        "iPhi_1on2/I");
	Tree_Optim->Branch( "iPhi_2on2",     &Op_iPhi_2on2,        "iPhi_2on2/I");
	Tree_Optim->Branch( "iEta_1on2520",  &Op_iEta_1on2520,     "iEta_1on2520/I");
	Tree_Optim->Branch( "iEta_2on2520",  &Op_iEta_2on2520,     "iEta_2on2520/I");
	Tree_Optim->Branch( "iPhi_1on20",    &Op_iPhi_1on20,       "iPhi_1on20/I");
	Tree_Optim->Branch( "iPhi_2on20",    &Op_iPhi_2on20,       "iPhi_2on20/I");
	if( isMC_ && MC_Assoc_ ) {
	  Tree_Optim->Branch( "enG1_true",     &Op_enG1_true,        "enG1_true/F");
	  Tree_Optim->Branch( "enG2_true",     &Op_enG2_true,        "enG2_true/F");
	  Tree_Optim->Branch( "DeltaR_1",      &Op_DeltaR_1,         "DeltaR_1/F");
	  Tree_Optim->Branch( "DeltaR_2",      &Op_DeltaR_2,         "DeltaR_2/F");
	}
    }

    /// trigger histo
    if (L1TriggerInfo_) {
      triggerComposition = new TH1F("triggerComposition", "Trigger Composition", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);    
      triggerComposition_EB = new TH1F("triggerComposition_EB", "Trigger Composition in EB", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);
      triggerComposition_EE = new TH1F("triggerComposition_EE", "Trigger Composition in EE", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);
      // for L1
      seedIsInStream = new int[GlobalAlgBlk::maxPhysicsTriggers];
      algoBitToName = new TString[GlobalAlgBlk::maxPhysicsTriggers];
      l1flag = new short[GlobalAlgBlk::maxPhysicsTriggers];
    }
    areLabelsSet_ = false;
    //L1_nameAndNumb.clear();
    //for(unsigned int i=0; i<NL1SEED; i++) L1BitCollection_[i]=-1;

    if (isMC_ and MC_Assoc_) {
      // since we have 20 gen pi0, use 21,-0.5,20.5 as range if counting integer number
      h_numberUnmergedGenPhotonPairs_EB = new TH1F("h_numberUnmergedGenPhotonPairs_EB",Form("fraction of gen photon pairs in EB with #DeltaR > %.3f",DR_FOR_UNMERGED_GEN_PHOTONS),20,0,1.01);
      h_numberMatchedGenPhotonPairs_EB = new TH1F("h_numberMatchedGenPhotonPairs_EB","fraction of gen photon pairs in EB matched to reco clusters",20,0,1.01);
      h_numberUnmergedGenPhotonPairs_EE = new TH1F("h_numberUnmergedGenPhotonPairs_EE",Form("fraction of gen photon pairs in EE with #DeltaR > %.3f",DR_FOR_UNMERGED_GEN_PHOTONS),20,0,1.01);
      h_numberMatchedGenPhotonPairs_EE = new TH1F("h_numberMatchedGenPhotonPairs_EE","fraction of gen photon pairs in EE matched to reco clusters",20,0,1.01);
      h_numberUnmergedGenPhotonPairs = new TH1F("h_numberUnmergedGenPhotonPairs",Form("gen photon pairs with #DeltaR > %.3f",DR_FOR_UNMERGED_GEN_PHOTONS),21,-0.5,20.5);
      h_numberMatchedGenPhotonPairs = new TH1F("h_numberMatchedGenPhotonPairs","gen photon pairs succesfully matched to reco clusters",21,-0.5,20.5);
    }

#ifdef MVA_REGRESSIO
    EBweight_file_1 = TFile::Open( Are_pi0_? edm::FileInPath( MVAEBContainmentCorrections_01_.c_str() ).fullPath().c_str() : edm::FileInPath( MVAEBContainmentCorrections_eta01_.c_str() ).fullPath().c_str() );
    EBweight_file_2 = TFile::Open( Are_pi0_? edm::FileInPath( MVAEBContainmentCorrections_02_.c_str() ).fullPath().c_str() : edm::FileInPath( MVAEBContainmentCorrections_eta02_.c_str() ).fullPath().c_str() );
    if(new_pi0ContainmentCorrections_)
        {
        forestD_EB_1 = (GBRForestD *)EBweight_file_1->Get("Correction");
        forestD_EB_2 = (GBRForestD *)EBweight_file_2->Get("Correction");
        }
    else
        {
        forest_EB_1 = (GBRForest *)EBweight_file_1->Get("Correction");
        forest_EB_2 = (GBRForest *)EBweight_file_2->Get("Correction");
        }
#endif
#ifdef MVA_REGRESSIO_EE
    EEweight_file_pi01 = TFile::Open( edm::FileInPath( MVAEEContainmentCorrections_01_.c_str() ).fullPath().c_str() );
    EEweight_file_pi02 = TFile::Open( edm::FileInPath( MVAEEContainmentCorrections_02_.c_str() ).fullPath().c_str() );
    if(new_pi0ContainmentCorrections_)
        {
        forestD_EE_pi01 = (GBRForestD *)EEweight_file_pi01->Get("Correction");
        forestD_EE_pi02 = (GBRForestD *)EEweight_file_pi02->Get("Correction");
        }
    else
        {
        forest_EE_pi01 = (GBRForest *)EEweight_file_pi01->Get("Correction");
        forest_EE_pi02 = (GBRForest *)EEweight_file_pi02->Get("Correction");
        }
#endif

}

FillEpsilonPlot::~FillEpsilonPlot()
{
  externalGeometryFile_->Close();
  outfile_->Write();
  outfile_->Close();

  if (useContainmentCorrectionsFromEoverEtrue_) {
    delete hCC_EoverEtrue_g1;
    delete hCC_EoverEtrue_g2;
  }

  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    if (isEoverEtrue_) {
      deleteEpsilonPlot(EoverEtrue_g1_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
      deleteEpsilonPlot(EoverEtrue_g2_EB_h, regionalCalibration_g2_->getCalibMap()->getNRegionsEB() );
    } else {
      deleteEpsilonPlot(epsilon_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
    }
  }

  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    if (isEoverEtrue_) {
      deleteEpsilonPlot(EoverEtrue_g1_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
      deleteEpsilonPlot(EoverEtrue_g2_EE_h, regionalCalibration_g2_->getCalibMap()->getNRegionsEE() );
    } else {
      deleteEpsilonPlot(epsilon_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
    }
  }

  if (isEoverEtrue_) {
    delete allEoverEtrue_g1_EB;
    delete allEoverEtrue_g1_EBnw;
    delete allEoverEtrue_g1_EE;
    delete allEoverEtrue_g1_EEnw;
    delete allEoverEtrue_g2_EB;
    delete allEoverEtrue_g2_EBnw;
    delete allEoverEtrue_g2_EE;
    delete allEoverEtrue_g2_EEnw;
  } else {
    delete allEpsilon_EB;
    delete allEpsilon_EBnw;
    delete allEpsilon_EE;
    delete allEpsilon_EEnw;
  }

  if (fillKinematicVariables_) {
    for (uint32_t i = 0; i < pi0pt_afterCuts.size(); ++i) {
      delete pi0pt_afterCuts[i];
      delete g1pt_afterCuts[i];
      delete g2pt_afterCuts[i];
      delete g1Nxtal_afterCuts[i];
      delete g2Nxtal_afterCuts[i];
      delete pi0PhotonsNoverlappingXtals_afterCuts[i];
      if (isMC_) {
	delete pi0MassVsPU[i];
      }
    }
    pi0pt_afterCuts.clear();
    g1pt_afterCuts.clear();
    g2pt_afterCuts.clear();
    g1Nxtal_afterCuts.clear();
    g2Nxtal_afterCuts.clear();
    pi0PhotonsNoverlappingXtals_afterCuts.clear();
  }

  delete entries_EEp;
  delete entries_EEm;
  delete entries_EB;
  delete Occupancy_EEp;
  delete Occupancy_EEm;
  delete Occupancy_EB;
  delete pi0MassVsIetaEB;
  delete pi0MassVsETEB;
  if (L1TriggerInfo_) {
    delete triggerComposition;
    delete triggerComposition_EB;
    delete triggerComposition_EE;
  }

  if (isMC_ and MC_Assoc_) {
    delete h_numberUnmergedGenPhotonPairs_EB;
    delete h_numberMatchedGenPhotonPairs_EB;
    delete h_numberUnmergedGenPhotonPairs_EE;
    delete h_numberMatchedGenPhotonPairs_EE;
    delete h_numberUnmergedGenPhotonPairs;
    delete h_numberMatchedGenPhotonPairs;
  }


#ifdef SELECTION_TREE
  delete CutVariables_EB;
  delete CutVariables_EE;
  delete Pi0Info_EB;
  delete Pi0Info_EE;
#endif

  EndcapTools::freeMemory();
  delete geom_;
  delete ebtopology_;
  delete eetopology_;

#if (defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)) || defined(REGRESS_AND_PARAM_CONTCORR)
  delete EBPHI_ConCorr_p;
  delete EBPHI_ConCorr_m;
#endif
  //JSON
  if(JSONfile_ != "") delete myjson;
  //#ifdef MVA_REGRESSIO
  //  // if the analyzer did not run it crash because you do not create it. Better never delete it
  //  if(!isMC_){
  //    delete forest_EB_1;
  //    delete forest_EB_2;
  //  }
  //#endif
  //#ifdef MVA_REGRESSIO_EE
  //  delete forest_EE_pi01;
  //  delete forest_EE_pi02;
  //#endif
  //if( calibMapPath_.find("iter_-1")!=std::string::npos ){
  //Write the PassPreselection Map
  //cout<<"Preselection:: Siamo al primo iter: Scrivo le correzioni"<<endl;
  //PassPreselection
  //}

  // for L1
  if (L1TriggerInfo_) {
    delete[] seedIsInStream;
    delete[] algoBitToName;
    delete[] l1flag;
  }

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
FillEpsilonPlot::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // std::cout<<"Event: "<<iEvent.id().event()<<" Run "<<iEvent.id().run()<<" LS "<<iEvent.id().luminosityBlock()<<endl;
  //JSON
  EventFlow_EB->Fill(0.); EventFlow_EE->Fill(0.);
  if ( JSONfile_!="" && !myjson->isGoodLS(iEvent.id().run(),iEvent.id().luminosityBlock()) ) return;
  EventFlow_EB->Fill(1.); EventFlow_EE->Fill(1.);

  myEvent = iEvent.id().event();
  myLumiBlock = iEvent.id().luminosityBlock();
  myRun = iEvent.id().run();
  myBunchCrossing = iEvent.bunchCrossing();

  // std::cout << "iEvent.bunchCrossing() = " << iEvent.bunchCrossing() << std::endl;

  if (MakeNtuple4optimization_) {

    Op_Pi0recIsEB.clear();
    Op_ClusIsoPi0.clear();
    Op_HLTIsoPi0.clear();
    Op_nCrisG1.clear();
    Op_nCrisG2.clear();
    Op_enG1_cor.clear();
    Op_enG2_cor.clear();
    Op_etaG1_cor.clear();
    Op_etaG2_cor.clear();
    Op_phiG1_cor.clear();
    Op_phiG2_cor.clear();
    Op_mPi0_cor.clear();
    Op_etaPi0_cor.clear();
    Op_ptPi0_cor.clear();
    Op_DeltaRG1G2.clear();
    Op_Es_e1_1.clear();
    Op_Es_e1_2.clear();
    Op_Es_e2_1.clear();
    Op_Es_e2_2.clear();
    Op_S4S9_1.clear();
    Op_S4S9_2.clear();
    Op_S1S9_1.clear();
    Op_S1S9_2.clear();
    Op_S2S9_1.clear();
    Op_S2S9_2.clear();
    Op_Time_1.clear();
    Op_Time_2.clear();
    Op_DeltaR_1.clear();
    Op_DeltaR_2.clear();
    Op_enG1_nocor.clear();
    Op_enG2_nocor.clear();
    Op_etaG1_nocor.clear();
    Op_etaG2_nocor.clear();
    Op_phiG1_nocor.clear();
    Op_phiG2_nocor.clear();
    Op_ptPi0_nocor.clear();
    Op_mPi0_nocor.clear();
    Op_enG1_true.clear();
    Op_enG2_true.clear();
    Op_iEtaiX_1.clear();
    Op_iEtaiX_2.clear();
    Op_iPhiiY_1.clear();
    Op_iPhiiY_2.clear();
    Op_iEta_1on5.clear();
    Op_iEta_2on5.clear();
    Op_iPhi_1on2.clear();
    Op_iPhi_2on2.clear();
    Op_iEta_1on2520.clear();
    Op_iEta_2on2520.clear();
    Op_iPhi_1on20.clear();
    Op_iPhi_2on20.clear();

  }


  if( !areLabelsSet_ && L1TriggerInfo_ ){
    
    edm::Handle< GlobalAlgBlkBxCollection > gtReadoutRecord;
    iEvent.getByToken( L1GTobjmapToken_, gtReadoutRecord);
 
    if (gtReadoutRecord.isValid()) { 

      const GlobalAlgBlkBxCollection *l1results = gtReadoutRecord.product(); 
      if (l1results->size() == 0) std::cout << "%L1Results -- No trigger name given in TriggerResults of the input " << std::endl;

 	
      edm::ESHandle<L1TUtmTriggerMenu> menu;
      iSetup.get<L1TUtmTriggerMenuRcd>().get(menu);

      // get the bit/name association         
      for (auto const & keyval: menu->getAlgorithmMap()) { 
	std::string const & trigName  = keyval.second.getName(); 
	unsigned int iTrigIndex = keyval.second.getIndex(); 
	std::cerr << "bit: " << iTrigIndex << "\tname: " << trigName << std::endl;                                                         
	algoBitToName[iTrigIndex] = TString( trigName );

      } // end algo Map
 
      int trigCompBin = 1;
      GlobalAlgBlk const &result = l1results->at(0, 0);

      for (unsigned int itrig = 0; itrig < result.maxPhysicsTriggers; ++itrig) {
	//      std::cerr << "bit: " << itrig << "\tresult: " << results.getAlgoDecisionFinal(itrig) << std::endl;

	// some indices are empty: name them appropriately
	if (std::string(algoBitToName[itrig]) == "") algoBitToName[itrig] = Form("EMPTY_%d",itrig);

	//std::string l1triggername = std::string(algoBitToName[itrig]);
	//L1_nameAndNumb[std::string(algoBitToName[itrig])] = itrig;

	// check if index is valid
	if ( std::string(algoBitToName[itrig]).find("EMPTY") != std::string::npos ) {

	  // -1 for non valid index
	  seedIsInStream[itrig] = -1;
	  l1flag[itrig] = -2; 

	} else {

	  // check if seed is used by the stream: seed expression is "seed1 OR seed2 OR seed3 ... "
	  // the space at the end is important: see parameters.py
	  if ( L1SeedsPi0Stream_.find((algoBitToName[itrig]+" ")) != std::string::npos ) { 

	    seedIsInStream[itrig] = 1;
	    bool myflag = result.getAlgoDecisionFinal(itrig) ; 
	    if (myflag ) { l1flag[itrig] = 1; }
	    else {l1flag[itrig] = 0 ; }
	    // save bit for any seed n the stream in the tree
	    std::string trigName = std::string(algoBitToName[itrig]);
	    cout << "trigCompBin = " << trigCompBin << "    Seed name = " << algoBitToName[itrig] << endl;
	    if (trigCompBin <= triggerComposition->GetNbinsX()) triggerComposition->GetXaxis()->SetBinLabel(trigCompBin,trigName.c_str());
	    else cout << "Warning: trigCompBin is exceeding the allowed number of bins. Check! " << endl;
	    trigCompBin++;
	    if(MakeNtuple4optimization_) Tree_Optim->Branch(trigName.c_str(),l1flag+itrig,(trigName+"/S").c_str());   // l1flag+(int)itrig is the pointer to the itrig-th object of l1flag
	    //Tree_Optim->Branch((trigName+"_Prescl").c_str(),l1Prescl+(int)itrig,(trigName+"_Prescl/I").c_str());  // not implemented yet
	  } else {

	    seedIsInStream[itrig] = 0;
	    l1flag[itrig] = -1;

	  }

	}
	
	// std::cout << "L1 TD: "<<itrig<<" "<<algoBitToName[itrig]<<" " 
	// 	  << l1flag[itrig] <<" " 
	// 	  << std::endl;           
	

      }

      if(!areLabelsSet_){
	areLabelsSet_ = true;
	cout << "setting labels of triggerComposition histogram" << endl;
      }

    }

  }
  // end of --> if (!areLabelsSet_ && L1TriggerInfo_)

  //MC Photons (they will be associated to the clusters later)

  if( isMC_ ) {

    // get the PU collection
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    nPUtrue_ = PupInfo->begin()->getTrueNumInteractions(); // it is the same for each PVI and it is a float

    nPUobs_BX0_ = -1;
    Int_t nBX = 0;
    nPUobs_.clear();
    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

      nPUobs_[PVI->getBunchCrossing()] = PVI->getPU_NumInteractions();
      if (PVI->getBunchCrossing() == 0) nPUobs_BX0_ = PVI->getPU_NumInteractions();     
      // std::cout << "   BX[" << nBX << "]: " << PVI->getBunchCrossing()
      // 		<< "   nPUtrue: " << PVI->getTrueNumInteractions()
      // 		<< "   nPUobs[" << PVI->getBunchCrossing() << "]: " << nPUobs_[PVI->getBunchCrossing()]
      // 		<< std::endl;
      nBX++;

    }

    // vecGamma1MC.clear();
    // vecGamma2MC.clear();

    vecGamma1MC_EB.clear();
    vecGamma2MC_EB.clear();

    vecGamma1MC_EE.clear();
    vecGamma2MC_EE.clear();

    // from 2017 we have pi0 and photons in MC
    iEvent.getByToken( GenPartCollectionToken_,genParticles);

    // taken from Zhicai's private ntuplizer: https://github.com/RazorCMS/Pi0Tuplizer/blob/master/plugins/Pi0Tuplizer.cc#L1240-L1298
    for (size_t iG = 0; iG < genParticles->size(); ++iG ) {

      if((*genParticles)[iG].status()!=2) continue;

      unsigned int ndau = (*genParticles)[iG].numberOfDaughters();
      if( (*genParticles)[iG].pdgId() != (Are_pi0_ ? 111 : 221)  ) continue;

      if((*genParticles)[iG].pdgId() == 111) {
	// pi0
	  // ptPi0_genall[N_Pi0_genall] = (*genParticles)[iG].pt(); 
	  // etaPi0_genall[N_Pi0_genall] = (*genParticles)[iG].p4().Eta(); 
	  // phiPi0_genall[N_Pi0_genall] = (*genParticles)[iG].p4().Phi(); 
	  // N_Pi0_genall ++;
      }
      
      // if((*genParticles)[iG].pdgId() == 221) {
      // 	//eta
      // 	  // ptEta_genall[N_Eta_genall] = (*genParticles)[iG].pt(); 
      // 	  // etaEta_genall[N_Eta_genall] = (*genParticles)[iG].p4().Eta(); 
      // 	  // phiEta_genall[N_Eta_genall] = (*genParticles)[iG].p4().Phi(); 
      // 	  // N_Eta_genall ++;
      // }

	
      if(ndau != 2 ) continue;
      bool isDiphoton = true;

      for (unsigned int jD=0; jD < ndau; ++jD) {
	const reco::Candidate *dau = (*genParticles)[iG].daughter(jD);
	if(dau->pdgId() != 22) isDiphoton=false;
      }

      if(!isDiphoton) continue;

      //fill GEN pi0
      if ((*genParticles)[iG].pdgId() == 111) {

	TLorentzVector gamma1_temp, gamma2_temp;
	//  g1_tmp.SetXYZ(*genParticles)[iG].daughter(0)->momentum().x(),(*genParticles)[iG].daughter(0)->momentum().y(),(*genParticles)[iG].daughter(0)->momentum().z());
	//  g2_tmp.SetXYZ(*genParticles)[iG].daughter(1)->momentum().x(),(*genParticles)[iG].daughter(1)->momentum().y(),(*genParticles)[iG].daughter(1)->momentum().z());

	const reco::Candidate *dau1 = (*genParticles)[iG].daughter(0);
	const reco::Candidate *dau2 = (*genParticles)[iG].daughter(1);

	if(dau1->pt() > dau2->pt()) 
	  {
	    gamma1_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
	    gamma2_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
	  }
	else
	  {
	    gamma2_temp.SetPtEtaPhiE(dau1->pt(), dau1->p4().Eta(), dau1->p4().Phi(), dau1->p4().E());
	    gamma1_temp.SetPtEtaPhiE(dau2->pt(), dau2->p4().Eta(), dau2->p4().Phi(), dau2->p4().E());
	  }
	
	// keep only gen photons that do not merge too much. Threshold is chosen in order to have the seed crystals not in the other photon's 3x3 matrix
	// if commented, keep them until the gen-reco matching to count how many are merged
	//if (gamma1_temp.DeltaR(gamma2_temp) > DR_FOR_UNMERGED_GEN_PHOTONS) {
	if (fabs( (*genParticles)[iG].p4().Eta()) <= 1.479 ) {
	  vecGamma1MC_EB.push_back(gamma1_temp);
	  vecGamma2MC_EB.push_back(gamma2_temp);
	} else {
	  vecGamma1MC_EE.push_back(gamma1_temp);
	  vecGamma2MC_EE.push_back(gamma2_temp);
	}
	// vecGamma1MC.push_back(gamma1_temp);
	// vecGamma2MC.push_back(gamma2_temp);
	//}

      }
      
    }  // end of loop on genParticles

    //std::cout << "There are " << vecGamma1MC.size() << " diphoton pairs" << std::endl;

    // OLD PART FOR GEN LEVEL, WILL BE REMOVED AT SOME POINT

    // // const reco::GenParticleCollection *GenPars = 0;
    // // std::cout << "MC truth taken" << std::endl;
    // // //if ( ! GenParProd.isValid() )  edm::LogWarning("GenParSummary") << "GenPars not found";
    // // if ( ! GenParProd.isValid() )  std::cout << "GenPars not found" << std::endl;
    // // GenPars = GenParProd.product();


    // // GUN sample made with PYTHIA6 doesn't decay the pi0, need to look at simtracks by GEANT
    // // get GEANT sim tracks and vertices (includes conversions)
    // Handle<SimTrackContainer> simTracks_h;
    // const SimTrackContainer* simTracks;
    // iEvent.getByToken(g4_simTk_Token_, simTracks_h);
    // simTracks = (simTracks_h.isValid()) ? simTracks_h.product() : 0;

    // Handle<SimVertexContainer> simVert_h;
    // const SimVertexContainer* simVertices;
    // iEvent.getByToken(g4_simVtx_Token_, simVert_h);
    // simVertices = (simVert_h.isValid()) ? simVert_h.product() : 0;


    // // Vertices only return trackID of their parent SimTrack
    // // Figure out the mapping from trackID to SimTrack
    // map<unsigned int, const SimTrack*> trackMap;
    // for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {
    //   if (!iSim->noVertex()) {
    //     assert(trackMap.find(iSim->trackId())==trackMap.end());
    //     trackMap[iSim->trackId()] = &(*iSim);
    //   }
    // }

    // // Find all SimTracks that come from decays before the ECAL
    // // and find their parent SimTracks
    // map<const SimTrack*, const SimTrack*> promptParent; // daughter->mother
    // map<const SimTrack*, set<const SimTrack*> > promptDecays; // m->ds
    // map<const SimTrack*, const SimVertex*> promptVertex; // daughter->vertex
    // map<const SimTrack*, const SimVertex*> promptALLVertex; // daughter->vertex in Any Occasion
    // map<const SimTrack*, const SimTrack*> promptALLParent; // daughter->mother in Any Occasion

    // int num=0;
    // for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim, num++) {
    //   if (!iSim->noVertex()) {

    // 	// Find the parent vertex and see if it classifies as an early decay
    // 	// Exclude the primary vertex (noParent)
    // 	SimVertex const& vtx = (*simVertices)[iSim->vertIndex()];
    // 	if ( !vtx.noParent() ) {
	  
    // 	  assert(trackMap.find(vtx.parentIndex())!=trackMap.end());
    // 	  const SimTrack* p = trackMap[vtx.parentIndex()];
	  
    // 	  if ( vtx.position().Rho() < 129 && fabs(vtx.position().z()) < 304) {
    // 	    // Find parent SimParticle that produced this vertex
    // 	    // vtx->parentIndex is NOT a vector index :( so use trackMap
    // 	    promptParent[&(*iSim)] = p; // nel Pi0Gun: ->genpartIndex() e' -1
    // 	    promptDecays[p].insert(&(*iSim));
    // 	    promptVertex[&(*iSim)] = &vtx;
    // 	  } // early decay
	  
    // 	  promptALLVertex[&(*iSim)] = &vtx;
    // 	  promptALLParent[&(*iSim)] = p; 
    // 	  // cout<<num<<" PDG: "<<iSim->type()<<" id: "<<iSim->trackId()<<" Son of: "<<p->type()<<" id: "<<p->trackId()
    // 	  // <<" x vtx:" <<vtx.position().x()<<" Z Vtx: "<<vtx.position().z()<<"Px "<<iSim->momentum().x()<<" py "<<iSim->momentum().y()<<" pz "<<iSim->momentum().z()<<
    // 	  // "pt "<<sqrt(pow(iSim->momentum().x(),2)+pow(iSim->momentum().y(),2)+pow(iSim->momentum().z(),2))<<endl;
    // 	}
	
    //   } // has vertex
      
    // } // for simTracks
    
    // //cout<<"Event:"<<endl;
    // //Store Pi0 & gamma
    // //    unsigned int IdGamma1=0, IdGamma2=0;
    // TVector3 pi0_pos, ga1, ga2;
    // num=0;
    // for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim, num++) {
    //   if (!iSim->noVertex() ){
    //     SimVertex const& vtx = (*simVertices)[iSim->vertIndex()];
    //     if( !vtx.noParent() ) {
          
    //       if( iSim->type()==22 && promptALLParent[&(*iSim)]->type()==(Are_pi0_ ? 111:221) && num==1){
    //         pi0_pos.SetXYZ(promptALLVertex[&(*iSim)]->position().x(),promptALLVertex[&(*iSim)]->position().y(),promptALLVertex[&(*iSim)]->position().z());
    //         Gamma1MC.SetXYZ(iSim->momentum().x(),iSim->momentum().y(),iSim->momentum().z());
    //         //            IdGamma1 = iSim->trackId();
    //         //            std::cout << "Photon1 eta,phi = " << Gamma1MC.eta() << "  " << Gamma1MC.phi() << std::endl;
    // 	  }
    //       if( iSim->type()==22 && promptALLParent[&(*iSim)]->type()==(Are_pi0_ ? 111:221) && num==2){
    //         Gamma2MC.SetXYZ(iSim->momentum().x(),iSim->momentum().y(),iSim->momentum().z());
    //         //            IdGamma2 = iSim->trackId();
    //         //            std::cout << "Photon2 eta,phi = " << Gamma2MC.eta() << "  " << Gamma2MC.phi() << std::endl;
    //       }
    //     }
    //   }
    // }

    // //Find MC photons
    // /*
    // std::cout << "coll size = " << GenPars->size() << std::endl;
    // bool firstnotfound = true;
    // //    for (auto& GenPar : *GenPars){
    // for (reco::GenParticleCollection::const_iterator GenPar = GenPars->begin(); GenPar != GenPars->end(); ++GenPar) {
    //   std::cout << "id = " << GenPar->pdgId() << std::endl;
    //   if(GenPar->mother()!=0) std::cout << " mothId = " << GenPar->mother()->pdgId() << std::endl;
      
    //   int motherID = Are_pi0_ ? 111:221;
    //   if( GenPar->pdgId()==22 && GenPar->mother()->pdgId()==motherID && firstnotfound ){
    //     std::cout << "Found 1st photon, pt = " << GenPar->pt() << "  " << GenPar->p4().Eta() << "  " << GenPar->p4().Phi() << std::endl;
    //     std::cout << "dentro id = " << GenPar->pdgId() << std::endl;
    //     Gamma1MC.SetPtEtaPhiE( GenPar->pt(), GenPar->p4().Eta(), GenPar->p4().Phi(), GenPar->p4().E() );
    //     firstnotfound = false;
    //   }
    //   if( GenPar->pdgId()==22 && GenPar->mother()->pdgId()==motherID && GenPar->p4().Eta() != Gamma1MC.Eta() ){
    //     std::cout << "Found 2nd photon, pt = " << GenPar->pt() << "  " << GenPar->p4().Eta() << "  " << GenPar->p4().Phi() << std::endl;
    //     Gamma2MC.SetPtEtaPhiE( GenPar->pt(), GenPar->p4().Eta(), GenPar->p4().Phi(), GenPar->p4().E() );
    //   }
    //   std::cout << "running PT1,PT2 = " << Gamma1MC.Pt() << " , " << Gamma2MC.Pt() << std::endl;
    // }
    // std::cout << "==> final PT1,PT2 = " << Gamma1MC.Pt() << " , " << Gamma2MC.Pt() << std::endl;
    // */

  }


  if (isDebug_) cout << "\n --------------- [DEBUG] Beginning New Event ------------------"<< endl;

  using namespace edm;
  nPi0=0;
  //For Syst error SystOrNot_=1 or 2, for normal calib is 0
  if(SystOrNot_==1. && int(iEvent.id().event())%2!=0 ) return;
  else if(SystOrNot_==2. && int(iEvent.id().event())%2==0 ) return;

  iEvent.getByToken ( EBRecHitCollectionToken_, ebHandle);
  iEvent.getByToken ( EERecHitCollectionToken_, eeHandle);
  iEvent.getByToken ( ESRecHitCollectionToken_, esHandle);


  //Internal Geometry
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  geometry = geoHandle.product();
  estopology_ = new EcalPreshowerTopology(geoHandle);
  esGeometry_ = (dynamic_cast<const EcalPreshowerGeometry*>( (CaloSubdetectorGeometry*) geometry->getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ));

  ///////////////////////
  // I moved the evaluation of HLT before that of the L1 seeds because the triggerComposition histogram is filled inside getTriggerResult() method
  // Since I have also the triggerComposition for EB or EE only, I need to know in advance if EB_HLT or EE_HLT fired (one of them should have fired)
  ///////////////////////////////
  // I put definition of these variables in FillEpsilonPlot.h, so they are accessible from any method of fillEpsilonPlot
  EB_HLT=true, EE_HLT=true;

  // Warning: when you are filling ntuples for data, GetHLTResults() should be used, otherwise when entering fillEEClusters()
  // the code crushes saying
 
 // ----- Begin Fatal Exception 12-Apr-2017 09:58:12 CEST-----------------------
  //   An exception of category 'ProductNotFound' occurred while
  //   [0] Processing run: 282814 lumi: 797 event: 1519368689
  //   [1] Running path 'p'
  //   [2] Calling event method for module FillEpsilonPlot/'analyzerFillEpsilon'
  // 	Exception Message:
  //    Principal::getByToken: Found zero products matching all criteria
  //    Looking for type: edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >
  //    Looking for module label: hltAlCaPi0RecHitsFilterEEonlyRegional
  //    Looking for productInstanceName: pi0EcalRecHitsES

  // For MC this is not necessary probably (I didn't check)
  // if( HLTResults_ && (!MakeNtuple4optimization_) ){
  if( HLTResults_){
    EB_HLT = GetHLTResults(iEvent, HLTResultsNameEB_); //Adding * at the end of the sentence make always true the "->Contains" method. So do not use it.
    EE_HLT = GetHLTResults(iEvent, HLTResultsNameEE_);
  }
  //std::cout << "EB_HLT,EE_HLT = " << EB_HLT << "," << EE_HLT << endl;

  //L1 Trigget bit list (and cut if L1_Bit_Sele_ is not empty)
  // this function is not meant to apply the global L1 seed expression to accept or not the event
  // it can reject events if you decide to select one or more specific seeds (this must still be implemented)
  if( L1TriggerInfo_ ){ if( !getTriggerResult(iEvent, iSetup) ) return; }


  //Vectors
  std::vector< CaloCluster > ebclusters;
  ebclusters.clear();
  vs4s9.clear(); vs2s9.clear(); vs2s9.clear(); vSeedTime.clear();
  vs4s9EE.clear(); Es_1.clear(); Es_2.clear(); vSeedTimeEE.clear();
  vs2s9EE.clear(); vs2s9EE.clear(); ESratio.clear();
  //cout << "I'm before std::vector< CaloCluster > eseeclusters; eseeclusters.clear(); " << endl;
  std::vector< CaloCluster > eseeclusters; eseeclusters.clear();
  std::vector< CaloCluster > eseeclusters_tot; eseeclusters_tot.clear();
  Ncristal_EB.clear(); Ncristal_EE.clear();
  //cout << "I'm after Ncristal_EB.clear(); Ncristal_EE.clear(); " << endl;

  //get status from DB
  edm::ESHandle<EcalChannelStatus> csHandle; 
  iSetup.get<EcalChannelStatusRcd>().get(csHandle);
  const EcalChannelStatus &channelStatus = *csHandle; 
  ////cout << "I'm after const EcalChannelStatus &channelStatus = *csHandle; " << endl;

  EventFlow_EB->Fill(2.); EventFlow_EE->Fill(2.);
  if ( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) && EB_HLT ) { 
    EventFlow_EB->Fill(3.); 
    fillEBClusters(ebclusters, iEvent, channelStatus);
  }
  ////cout << "I'm after fillEBClusters(ebclusters, iEvent, channelStatus) " << endl;
  if ( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) && EE_HLT ) { 
    EventFlow_EE->Fill(3.); 
    fillEEClusters(eseeclusters, eseeclusters_tot, iEvent, channelStatus);
  }
  // std::cout << "ebclusters.size() = " << ebclusters.size() << std::endl;
  // std::cout << "eseeclusters.size() = " << eseeclusters.size() << std::endl;
  ////cout << "I'm after fillEEClusters(eseeclusters, eseeclusters_tot, iEvent, ...) " << endl;

  std::vector< CaloCluster > ebclusters_used, eeclusters_used;
  ebclusters_used.clear(); eeclusters_used.clear();
  // not used anymore, keep track of the whole gen photons through ebclusters_matchedGenPhoton and eeclusters_matchedGenPhoton 
  // std::vector< float > ebclusters_matchedGenPhotonEnergy;  // will store the energy of the gen photon corresponding to a given reco cluster
  // std::vector< float > eeclusters_matchedGenPhotonEnergy;  // will store the energy of the gen photon corresponding to a given reco cluster
  // ebclusters_matchedGenPhotonEnergy.clear();
  // eeclusters_matchedGenPhotonEnergy.clear();
  // moved in FillEpsilonPlot.h
  // std::vector< TLorentzVector* > ebclusters_matchedGenPhoton;  // will store the gen photon corresponding to a given reco cluster
  // std::vector< TLorentzVector* > eeclusters_matchedGenPhoton;  // will store the gen photon corresponding to a given reco cluster
  ebclusters_matchedGenPhoton.clear();
  eeclusters_matchedGenPhoton.clear();

  if(isMC_ && MC_Assoc_) {

    // ebclusters_used = MCTruthAssociate(ebclusters,MC_Assoc_DeltaR,true);
    // eeclusters_used = MCTruthAssociate(eseeclusters_tot,MC_Assoc_DeltaR,false);
    int unmergedGenPairs_EB = 0;
    int unmergedGenPairs_EE = 0;
    int matchedGenPairs_EB = 0;
    int matchedGenPairs_EE = 0;
    // based on the function implementation, the vector of clusters is such that clusters originating from the same gen pi0 are consecutive
    // therefore, we can consider each pair of clusters (i-th and (i+1)-th object, i=0,2,4,...) to be from a pi0  
    ebclusters_used = MCTruthAssociateMultiPi0(ebclusters,
					       unmergedGenPairs_EB,matchedGenPairs_EB,
					       ebclusters_matchedGenPhoton,
					       MC_Assoc_DeltaR,true);
    eeclusters_used = MCTruthAssociateMultiPi0(eseeclusters_tot,
					       unmergedGenPairs_EE,matchedGenPairs_EE,
					       eeclusters_matchedGenPhoton,
					       MC_Assoc_DeltaR,false);
    h_numberUnmergedGenPhotonPairs->Fill(unmergedGenPairs_EB+unmergedGenPairs_EE);
    h_numberMatchedGenPhotonPairs->Fill(matchedGenPairs_EB+matchedGenPairs_EE);

  } else {

    // Things named *_used are meaningful for MC because they are obtained filtering on MC truth to match the reco photons to the gen ones.
    // For data or in case we decide to disable MC_Assoc_ when using MC (but currently it is basically the same flag as isMC_) we just assign 
    // the collection to the "*_used" so that we can always use the latter without loss of generality

    Ncristal_EB_used = Ncristal_EB;
    Ncristal_EE_used = Ncristal_EE;
    ebclusters_used = ebclusters;
    eeclusters_used = eseeclusters_tot;

  }

  if(isMC_ && MC_Assoc_ && isEoverEtrue_) {  // asking just for isEoverEtrue_ would be enough because it can be true only for MC (see parameters.py)

    if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEoverEtrue(ebclusters_used, ebclusters_matchedGenPhoton, EcalBarrel);
    if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEoverEtrue(eeclusters_used, eeclusters_matchedGenPhoton, EcalEndcap);

  } else {

    if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(ebclusters_used, ebclusters_matchedGenPhoton, EcalBarrel);
    if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(eeclusters_used, ebclusters_matchedGenPhoton, EcalEndcap);
  
  }

  delete estopology_;

}


/*===============================================================*/
void FillEpsilonPlot::fillEBClusters(std::vector< CaloCluster > & ebclusters, const edm::Event& iEvent, const EcalChannelStatus &channelStatus)
  /*===============================================================*/
{

  std::vector<EcalRecHit> ebseeds;

  typedef std::set<EBDetId> XtalInUse;
  XtalInUse isUsed; // map of which xtals have been used

  int dc = 0;

  // sort by energy and find the seeds
  //bool founded=false;
  for(EBRecHitCollection::const_iterator itb= ebHandle->begin(); itb != ebHandle->end(); ++itb, ++dc) 
  {
    EBDetId tmp_id(itb->id());
    Occupancy_EB->Fill(tmp_id.ieta(), tmp_id.iphi());
    if(itb->energy() > EB_Seed_E_)  ebseeds.push_back( *itb );
    ////Preselection
    //if(itb->energy() > 0.200-0.200*(28.3/100)) founded=true;
  }
  //if(founded) FailPreselEB=false;
  //else        FailPreselEB=true;

#ifdef SELECTION_TREE
  Fill_NSeeds_EB(ebseeds.size());
#endif

  sort(ebseeds.begin(), ebseeds.end(), ecalRecHitLess());
  int seed_c = 0;
  // loop over seeds and make clusters
  for (std::vector<EcalRecHit>::iterator itseed=ebseeds.begin(); itseed!=ebseeds.end(); itseed++, seed_c++) 
  {
    EBDetId seed_id( itseed->id() );
    float SeedTime = itseed->time();

    // check if seed already in use. If so go to next seed
    if(isUsed.count(seed_id)!=0) continue;

    // find 3x3 matrix of xtals
    std::vector<DetId> clus_v = ebtopology_->getWindow(seed_id,3,3);       
    // needed for position calculator
    std::vector<std::pair<DetId,float> > clus_used;

    // xtals actually used after removing those already used
    vector<const EcalRecHit*> RecHitsInWindow;
    vector<const EcalRecHit*> RecHitsInWindow5x5;

    float simple_energy = 0; 
    float posTotalEnergy(0.); // need for position calculation

    // make 3x3  cluster - reject overlaps
    int i_clus=0;
    for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++, i_clus++) 
    {
	EBDetId thisId( *det );
	// skip this xtal if already used
	if(isUsed.count(thisId)!=0) continue; //already used

	// find the rec hit
	EBRecHitCollection::const_iterator ixtal = ebHandle->find( thisId );
	if( ixtal == ebHandle->end() ) continue; // xtal not found

	RecHitsInWindow.push_back( &(*ixtal) );
	clus_used.push_back(std::make_pair(*det,1.));  // it seems it is not used anywhere

	simple_energy +=  ixtal->energy();
	if(ixtal->energy()>0.) posTotalEnergy += ixtal->energy(); // use only pos energy for position
    } // loop over xtals in the region
    ///debug
    //cout << "seed #" << seed_c << "RecHitsInWindow.size() = " << RecHitsInWindow.size() << endl;

    if(simple_energy <= 0) { 
	//cout << "skipping cluster with negative energy " << simple_energy << endl; 
	continue;
    }

    float s4s9_tmp[4]={0.,0.,0.,0.};

    int seed_ieta = seed_id.ieta();
    int seed_iphi = seed_id.iphi();

    // following method defined in http://cmslxr.fnal.gov/source/HLTrigger/special/src/HLTRegionalEcalResonanceFilter.cc?v=CMSSW_9_4_1#0903
    // it is used to use diff_neta_s and diff_nphi_s to get distances between crystals in number of crystals
    // iphi is ported to 0,..,359, and ieta from -85,...,-1,0,1,2,...,84
    // note that ieta and iphi can be 0 and iphi= 360 and ieta = 85 are "lost", so do no tuse them to access histogram content   
    convxtalid( seed_iphi,seed_ieta);

    // energy of 3x3 cluster
    float e3x3(0.);
    std::vector<std::pair<DetId,float> > enFracs;

    // variables for position caculation
    float xclu(0.), yclu(0.), zclu(0.); // temp var to compute weighted average
    float total_weight(0.);// to compute position

    // Calculate shower depth
    float T0 = PCparams_.param_T0_barl_;
    float maxDepth = PCparams_.param_X0_ * ( T0 + log( posTotalEnergy ) );
    float maxToFront;
    if( GeometryFromFile_ ) maxToFront = geom_->getPosition(seed_id).mag(); // to front face
    else                  {
      const CaloCellGeometry* cell = geometry->getGeometry( seed_id ).get();
      GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
      maxToFront = posit.mag();
    }

    double EnergyCristals[9] = {0.};

    bool All_rechit_good=true;

    // loop over xtals and compute energy and position
    for(unsigned int j=0; j<RecHitsInWindow.size();j++)
    {

	EBDetId det(RecHitsInWindow[j]->id());

	if( RemoveDead_Flag_){
	  if(!checkStatusOfEcalRecHit(channelStatus, *RecHitsInWindow[j] )  )  {
	    All_rechit_good = false; 
	    break;  // exit this loop as soon as one crystal is dead, because at the end of this loop there is a continue if All_rechit_good = false
	    // at the moment I see the rechits are added in the list of used ones, which is probably wrong if I reject this cluster
	  }
	}
	if( RemoveDead_Map_!="" ){
	  if( isInDeadMap( true, *RecHitsInWindow[j] ) ) {
	    All_rechit_good = false;
	    break;  // exit this loop as soon as one crystal is dead
	  }
	}

	int ieta = det.ieta();
	int iphi = det.iphi();
	// if (useContainmentCorrectionsFromEoverEtrue_) {
	//   std::cout << "ieta,iphi,IC,CC1,CC2 " << ieta << "  " << iphi << "  " 
	// 	    << regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id()) << "  "
	// 	    << hCC_EoverEtrue_g1->GetBinContent(hCC_EoverEtrue_g1->FindFixBin(ieta,iphi)) << "  "
	// 	    << hCC_EoverEtrue_g2->GetBinContent(hCC_EoverEtrue_g2->FindFixBin(ieta,iphi)) << std::endl;
	  
	// }

	// following method defined in http://cmslxr.fnal.gov/source/HLTrigger/special/src/HLTRegionalEcalResonanceFilter.cc?v=CMSSW_9_4_1#0903
	// it is used to use diff_neta_s and diff_nphi_s to get distances between crystals in number of crystals
	// iphi is ported to 0,..,359, and ieta from -85,...,-1,0,1,2,...,84
	// note that ieta and iphi can be 0 and iphi= 360 and ieta = 85 are "lost", so do no tuse them to access histogram content   
	convxtalid(iphi,ieta);

	// use calibration coeff for energy and position
	// FIXME: if isEoverEtrue_ is true, then we are not using the pi0 IC, but the photon dependent correction based on E/Etrue
	// this means we should know which photon we are looking at
	// at this moment, we are not able to assess which photon will be the first one (we would know it if we could rely on the fact that the clusters 
	// will always be the same so that the first and second clusters are always the same (therefore we could rely on their DetId).
	float en = RecHitsInWindow[j]->energy();
	if (not isEoverEtrue_) {
	  en *= regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id());
	}

	int dx = diff_neta_s(seed_ieta,ieta);
	int dy = diff_nphi_s(seed_iphi,iphi);
	EnergyCristals[j] = en;

	if(abs(dx)<=1 && abs(dy)<=1) 
	{
	  e3x3 += en;
	  if(dx <= 0 && dy <=0){ s4s9_tmp[0] += en; }
	  if(dx >= 0 && dy <=0){ s4s9_tmp[1] += en; }
	  if(dx <= 0 && dy >=0){ s4s9_tmp[2] += en; }
	  if(dx >= 0 && dy >=0){ s4s9_tmp[3] += en; }
	  enFracs.push_back( std::make_pair( RecHitsInWindow[j]->id(), en ) );
	  // Note: I'm using fractions to save rechit energy
	  // isUsed.insert(RecHitsInWindow[j]->id());  // these crystals should be readded if the matrix is discarded, I add this line later

	}

	// compute position
	if(en>0.) 
	{
	  float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );
	  float pos_geo;
	  if( GeometryFromFile_ ) pos_geo = geom_->getPosition(det).mag(); // to front face
	  else                  {
	    const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	    GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    pos_geo = posit.mag();
	  }
	  float depth = maxDepth + maxToFront - pos_geo;
	  GlobalPoint posThis;
	  if( GeometryFromFile_ ) posThis = geom_->getPosition(det,depth);
	  else{
	    const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	    posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( depth );
	  }

	  xclu += weight*posThis.x(); 
	  yclu += weight*posThis.y(); 
	  zclu += weight*posThis.z(); 
	  total_weight += weight;
	}

    } // loop over 3x3 rechits

    if(!All_rechit_good) continue;

    float e2x2 = *max_element( s4s9_tmp,s4s9_tmp+4);
    float s4s9 = e2x2/e3x3;
    math::XYZPoint clusPos( xclu/total_weight, 
	  yclu/total_weight,
	  zclu/total_weight ); 

    //cout << "seed #" << seed_c << " ptClus(before): " << e3x3/cosh(clusPos.eta()) << endl;

#ifdef SELECTION_TREE
    Fill_xClus_EB(xclu/total_weight);
    Fill_yClus_EB(yclu/total_weight);
    Fill_zClus_EB(zclu/total_weight);
    Fill_S4S9_EB(s4s9);
    Fill_e3x3_EB(e3x3);
    Fill_PtClus_EB(e3x3/cosh(clusPos.eta()));
    CutVariables_EB->Fill();
#endif
    ////Preselection
    //if(i_clus==0){
    //  if(s4s9<S4S9_cut_[EcalBarrel]-S4S9_cut_[EcalBarrel]*(28.3/100) || e3x3/cosh(clusPos.eta())<gPtCut_[EcalBarrel]-gPtCut_[EcalBarrel]*(28.3/100) ) FailPreselEB=true;
    //  else FailPreselEB=false;
    //}
    //else{
    //  if(!FailPreselEB){
    //    if(s4s9<S4S9_cut_[EcalBarrel]-S4S9_cut_[EcalBarrel]*(28.3/100) || e3x3/cosh(clusPos.eta())<gPtCut_[EcalBarrel]-gPtCut_[EcalBarrel]*(28.3/100)) FailPreselEB=true;
    //  }
    //}

    // adding cut on number of crystals
    if( fabs( clusPos.eta() )<1. ) {
      if (s4s9<S4S9_cut_low_[EcalBarrel]) continue;
      if ( RecHitsInWindow.size() < min(nXtal_1_cut_low_[EcalBarrel],nXtal_2_cut_low_[EcalBarrel]) ) continue;
    } else { 
      if (s4s9<S4S9_cut_high_[EcalBarrel]) continue;
      if ( RecHitsInWindow.size() <  min(nXtal_1_cut_high_[EcalBarrel],nXtal_2_cut_high_[EcalBarrel]) ) continue;
    }


#if (defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)) || defined(REGRESS_AND_PARAM_CONTCORR)
    if(useEBContainmentCorrections_) 
    {
	e3x3 *=  containmentCorrections_.getContainmentCorrectionsEB(e3x3, seed_id.ieta() );
	e3x3 *=  EBPHI_Cont_Corr(e3x3/cosh(clusPos.eta()), seed_id.iphi()%20, seed_id.ieta() );
    }
#endif

    // compute pt of gamma and cut
    float ptClus = e3x3/cosh(clusPos.eta());

    if( fabs( clusPos.eta() )<1. ){ if( ptClus<gPtCut_low_[EcalBarrel]) continue; }
    else                          { if( ptClus<gPtCut_high_[EcalBarrel]) continue; }


    // reloop on recHits to add them in the set of used recHits
    for(unsigned int j=0; j<RecHitsInWindow.size();j++)
    {
      isUsed.insert(RecHitsInWindow[j]->id());
    }

    // make calo clusters
    vs4s9.push_back( s4s9 ); 
    vs1s9.push_back( itseed->energy()/e3x3 );
    double maxEne = max_array( EnergyCristals, 9 );
    for(int i=0; i<9; i++){ if( EnergyCristals[i]==maxEne ) EnergyCristals[i]=0.; }
    double maxEne2 = max_array( EnergyCristals, 9);
    vs2s9.push_back( (maxEne+maxEne2)/e3x3 );
    Ncristal_EB.push_back(RecHitsInWindow.size() );
    ebclusters.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_BARREL), enFracs, CaloCluster::undefined, seed_id ) );
    vSeedTime.push_back( SeedTime ); 
  } //loop over seeds to make EB clusters

  //std::cout << "### FillEpsilonPlot::fillEBClusters():   dc = " << dc << std::endl;


}

/*===============================================================*/
void FillEpsilonPlot::fillEEClusters(std::vector< CaloCluster > & eseeclusters, std::vector< CaloCluster > & eseeclusters_tot, const edm::Event& iEvent, const EcalChannelStatus &channelStatus)
  /*===============================================================*/
{

  //cout << "I'm before PreshowerTools esClusteringAlgo(geometry, estopology_, esHandle); " << endl;
  PreshowerTools esClusteringAlgo(geometry, estopology_, esHandle);
  ////cout << "I'm after PreshowerTools esClusteringAlgo(geometry, estopology_, esHandle); " << endl;

  std::vector<EcalRecHit> eeseeds;

  vector <double> eeclusterS4S9; eeclusterS4S9.clear();
  vector <double> SeedTime_v;    SeedTime_v.clear();
  vector <double> eeclusterS1S9; eeclusterS1S9.clear();
  vector <double> eeclusterS2S9; eeclusterS2S9.clear();

  std::vector< CaloCluster > eeclusters; // contains the output eeclusters
  eeclusters.clear();

  int dc = 0;

  // sort by energy and find the eeseeds
  //bool found=false;
  for(EERecHitCollection::const_iterator ite= eeHandle->begin(); ite != eeHandle->end(); ++ite, ++dc) {
    EEDetId idXtal( ite->id() );
    if(idXtal.zside()<0) Occupancy_EEm->Fill(idXtal.ix(),idXtal.iy()); 
    if(idXtal.zside()>0) Occupancy_EEp->Fill(idXtal.ix(),idXtal.iy()); 
    GlobalPoint posThis;
    if( GeometryFromFile_ ) posThis = geom_->getPosition(idXtal,0.);
    else{
      const CaloCellGeometry* cell = geometry->getGeometry(idXtal).get();
      posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
    }
    if( useEE_EtSeed_ ){ if(ite->energy()/cosh(posThis.eta()) > EE_Seed_Et_ )              eeseeds.push_back( *ite ); }
    else               { if(ite->energy()                     > EE_Seed_E_  )              eeseeds.push_back( *ite ); }
  } // loop over xtals

#ifdef SELECTION_TREE
  Fill_NSeeds_EE(eeseeds.size());
#endif

  sort(eeseeds.begin(), eeseeds.end(), ecalRecHitLess());

  typedef std::map< EEDetId, bool > EEXtalInUse;
  EEXtalInUse EEXisUsed;  //map of which eextals have been used

  //loop over seeds to make eeclusters
  for (std::vector<EcalRecHit>::iterator eeitseed=eeseeds.begin(); eeitseed!=eeseeds.end(); eeitseed++) 
  {
    EEDetId eeseed_id( eeitseed->id() );
    float SeedTimeEE = eeitseed->time();
    // check if seed already in use. If so go to next seed
    EEXtalInUse::const_iterator mapit = EEXisUsed.find( eeseed_id );
    if( mapit != EEXisUsed.end() ) continue; // seed already in use

    // find 3x3 matrix of xtals
    int clusEtaSize_(3), clusPhiSize_(3);
    std::vector<DetId> clus_v = eetopology_->getWindow(eeseed_id,clusEtaSize_,clusPhiSize_); 

    // needed for position calculator
    std::vector<std::pair<DetId,float> > clus_used;

    // xtals actually used after removing those already used
    vector<const EcalRecHit*> RecHitsInWindow;
    vector<const EcalRecHit*> RecHitsInWindow5x5;

    float simple_energy = 0.; 
    float posTotalEnergy(0.); // need for position calculation

    // make 3x3  cluster - reject overlaps
    int i_clus=0;
    for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++,i_clus++) 
    {
	EEDetId thisId( *det );
	// skip this xtal if already used
	EEXtalInUse::const_iterator mapit = EEXisUsed.find( thisId );
	if( mapit != EEXisUsed.end() ) continue; // xtal already used

	// find the rec hit
	EERecHitCollection::const_iterator ixtal = eeHandle->find( thisId );

	//cout<<"ixtal output: "<< ixtal->energy() <<endl;

	if( ixtal == eeHandle->end() ) continue; // xtal not found

	RecHitsInWindow.push_back( &(*ixtal) );
	clus_used.push_back(std::make_pair(*det,1.)); // it seems it is not used anywhereisUsed
	simple_energy +=  ixtal->energy();
	if(ixtal->energy()>0.) posTotalEnergy += ixtal->energy(); // use only pos energy for position
    }  // loop over xtals in the region
    if(simple_energy <= 0) { 
	//cout << "skipping cluster with negative energy " << simple_energy << endl; 
	continue;
    }

    float s4s9_tmp[4];
    for(int i=0;i<4;i++){ 
	s4s9_tmp[i]= 0;
    }

    int seed_ix = eeseed_id.ix(); 
    int seed_iy = eeseed_id.iy();   

    // energy of 3x3 cluster
    float e3x3(0.);

    std::vector<std::pair<DetId,float> > enFracs;

    // variables for position caculation
    float xclu(0.), yclu(0.), zclu(0.); // temp var to compute weighted average
    float total_weight(0.);// to compute position

    // Calculate shower depth
    float T0 = PCparams_.param_T0_endc_;
    float maxDepth = PCparams_.param_X0_ * ( T0 + log( posTotalEnergy ) );
    float maxToFront;
    if( GeometryFromFile_ ) maxToFront = geom_->getPosition(eeseed_id).mag(); // to front face
    else                   {
      const CaloCellGeometry* cell = geometry->getGeometry( eeseed_id ).get();
      GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
      maxToFront = posit.mag();
    }
    double EnergyCristals[9] = {0.};
    bool All_rechit_good=true;
    // loop over xtals and compute energy and position
    for(unsigned int j=0; j<RecHitsInWindow.size();j++)
    { 
	EEDetId det(RecHitsInWindow[j]->id());

	if( RemoveDead_Flag_ ){
	  if( !checkStatusOfEcalRecHit(channelStatus, *RecHitsInWindow[j] )  ) {
	    All_rechit_good = false;
	    break;  // exit this loop as soon as one crystal is dead, because at the end of this loop there is a continue if All_rechit_good = false
	    // at the moment I see the rechits are added in the list of used ones, which is probably wrong if I reject this cluster
	  }
	}
	if( RemoveDead_Map_!="" ){
	  if( isInDeadMap( false, *RecHitsInWindow[j] ) ) {
	    All_rechit_good = false;
	    break;  // exit this loop as soon as one crystal is dead
	  }
	}

	int ix = det.ix();
	int iy = det.iy();

	// use calibration coeff for energy and position
	// FIXME: if isEoverEtrue_ is true, then we are not using the pi0 IC, but the photon dependent correction based on E/Etrue
	// this means we should know which photon we are looking at
	float en = RecHitsInWindow[j]->energy();
	if (not isEoverEtrue_) {
	  en *= regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id());
	} 
	int dx = seed_ix-ix;
	int dy = seed_iy-iy;
	EnergyCristals[j] = en;
	if(abs(dx)<=1 && abs(dy)<=1) 
	{
	  e3x3 += en;
	  if(dx <= 0 && dy <=0){ s4s9_tmp[0] += en; }
	  if(dx >= 0 && dy <=0){ s4s9_tmp[1] += en; }
	  if(dx <= 0 && dy >=0){ s4s9_tmp[2] += en; }
	  if(dx >= 0 && dy >=0){ s4s9_tmp[3] += en; }
	  enFracs.push_back( std::make_pair( RecHitsInWindow[j]->id(), en ) );
	}

	// compute position
	if(en>0.) 
	{
	  float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );
	  float pos_geo;
	  if( GeometryFromFile_ ) pos_geo = geom_->getPosition(det).mag();
	  else                   {
	    const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	    GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    pos_geo = posit.mag();
	  }
	  float depth = maxDepth + maxToFront - pos_geo;
	  GlobalPoint posThis;
	  if( GeometryFromFile_ ) posThis = geom_->getPosition(det,depth);
	  else{
	    const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	    posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( depth );
	  }
	  xclu += weight*posThis.x(); 
	  yclu += weight*posThis.y(); 
	  zclu += weight*posThis.z(); 
	  total_weight += weight;
	}
    } // loop over 3x3 eerechits

    if(!All_rechit_good) continue;
    float e2x2 = *max_element( s4s9_tmp,s4s9_tmp+4);
    float s4s9 = e2x2/e3x3;
    math::XYZPoint clusPos( xclu/total_weight, yclu/total_weight, zclu/total_weight ); 

#ifdef SELECTION_TREE
    Fill_xClus_EE(xclu/total_weight);
    Fill_yClus_EE(yclu/total_weight);
    Fill_zClus_EE(zclu/total_weight);
    Fill_S4S9_EE(s4s9);
    Fill_e3x3_EE(e3x3);
    Fill_PtClus_EE(e3x3/cosh(clusPos.eta()));
    CutVariables_EE->Fill();
#endif

    //Preselection
    //if(i_clus==0){
    //  if(s4s9<S4S9_cut_[EcalEndcap]-S4S9_cut_[EcalEndcap]*(42.5/100) || e3x3/cosh(clusPos.eta())<gPtCut_[EcalEndcap]-gPtCut_[EcalEndcap]*(42.5/100) ) FailPreselEE=true;
    //  else FailPreselEE=false;
    //}
    //else{
    //  if(!FailPreselEE){
    //    if(s4s9<S4S9_cut_[EcalEndcap]-S4S9_cut_[EcalEndcap]*(42.5/100) || e3x3/cosh(clusPos.eta())<gPtCut_[EcalEndcap]-gPtCut_[EcalEndcap]*(42.5/100)) FailPreselEE=true;
    //  }
    //}
    if ( fabs( clusPos.eta() )<1.8 ) { 
      if (s4s9<S4S9_cut_low_[EcalEndcap]) continue; 
      if ( RecHitsInWindow.size() < min(nXtal_1_cut_low_[EcalEndcap],nXtal_2_cut_low_[EcalEndcap])) continue;
    } else { 
      if (s4s9<S4S9_cut_high_[EcalEndcap]) continue;
      if ( RecHitsInWindow.size() < min(nXtal_1_cut_high_[EcalEndcap],nXtal_2_cut_high_[EcalEndcap])) continue;
    }

    float ptClus = e3x3/cosh(clusPos.eta());

    if( fabs( clusPos.eta() )<1.8 ){ if(ptClus<gPtCut_low_[EcalEndcap]) continue; }
    else                           { if(ptClus<gPtCut_high_[EcalEndcap]) continue; }
    // make calo clusters
    for(unsigned int j=0; j<RecHitsInWindow.size();j++){
	EEXisUsed [RecHitsInWindow[j]->id()] = true;
    }
    Ncristal_EE.push_back( RecHitsInWindow.size() );
    eeclusters.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP),
	    enFracs, CaloCluster::undefined, eeseed_id ) );

    eeclusterS4S9.push_back(s4s9);
    SeedTime_v.push_back(SeedTimeEE);
    eeclusterS1S9.push_back(eeitseed->energy()/e3x3);
    double maxEne = max_array( EnergyCristals, 9 );
    for(int i=0; i<9; i++){ if( EnergyCristals[i]==maxEne ) EnergyCristals[i]=0.; }
    double maxEne2 = max_array( EnergyCristals, 9);
    eeclusterS2S9.push_back( (maxEne+maxEne2)/e3x3 );

  } //loop over seeds to make eeclusters


  /************************** ENDCAP-PRESHOWER MATCHING ************************/

  //loop over eecluster to find matches with preshower
  int ind=0;
  std::vector<int> Nxtal; Nxtal.clear();
  std::vector<int> Nxtal_tot; Nxtal_tot.clear();

  for( std::vector<CaloCluster>::const_iterator eeclus_iter  = eeclusters.begin(); eeclus_iter != eeclusters.end(); ++eeclus_iter, ++ind)
  {

    if(fabs(eeclus_iter->position().Eta())>1.7 && fabs(eeclus_iter->position().Eta())<2.55){
	double X = eeclus_iter->x();
	double Y = eeclus_iter->y(); 
	double Z = eeclus_iter->z();
	const GlobalPoint point(X,Y,Z);

	DetId tmp1 = esGeometry_->getClosestCellInPlane(point,1);
	DetId tmp2 = esGeometry_->getClosestCellInPlane(point,2);

	if ((tmp1.rawId()!=0) && (tmp2.rawId()!=0)) 
	{

	  ESDetId tmp1_conversion (tmp1);
	  ESDetId tmp2_conversion (tmp2);

          // replace the std PreshowerTools::clusterwindowsize_ = 15 with 5, smaller for 3x3 clusters
          float es_clusterwindowsize = 5;
	  //cout << "I'm before 	  PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp1_conversion); " << endl;
	  PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp1_conversion);
	  ////cout << "I'm after 	  PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp1_conversion); " << endl;
	  PreshowerCluster preshowerclusterp2 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp2_conversion);
	  ////cout << "I'm after 	  PreshowerCluster preshowerclusterp2 = esClusteringAlgo.makeOnePreshowerCluster( es_clusterwindowsize, &tmp2_conversion); " << endl;


	  double e1 = preshowerclusterp1.energy();
	  double e2 = preshowerclusterp2.energy();
	  // GeV to #MIPs
	  e1 = e1 / PreshowerTools::mip_;
	  e2 = e2 / PreshowerTools::mip_;
	  double tempenergy = eeclus_iter->energy();

	  //if(e1+e2 > 1.0e-10) 
	  if(e1 > 2.0 && e2 > 2.0) /// cut @ 2 MIPs as suggested by Ming @ DPG/EGamma Joint Meeting 19.03.2012 
	  {
	    double deltaE = PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1 + PreshowerTools::calib_planeY_*e2);

	    tempenergy = deltaE + eeclus_iter->energy();
#if (defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)) || defined(REGRESS_AND_PARAM_CONTCORR) 
	    if(useEEContainmentCorrections_) tempenergy *= containmentCorrections_.getContainmentPointCorrectionsEE( tempenergy , (eeclus_iter->position()).eta() );
#endif

	    eseeclusters.push_back( CaloCluster( tempenergy, eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
	    Nxtal.push_back(Ncristal_EE[ind]);

	    double DZ2 = (preshowerclusterp2.z()-preshowerclusterp1.z() )/2.;
	    GlobalPoint posClu(preshowerclusterp1.x()*(1.+DZ2/preshowerclusterp1.z() ),preshowerclusterp2.y()*(1.-DZ2/preshowerclusterp2.z()),(preshowerclusterp1.z()+preshowerclusterp2.z() )/2. );

	    if( fabs(preshowerclusterp1.z())>30  && fabs(preshowerclusterp2.z())>30){

		math::XYZPoint posit(posClu.x(),posClu.y(),posClu.z());
		eseeclusters_tot.push_back( CaloCluster( tempenergy, posit, CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
		Nxtal_tot.push_back(Ncristal_EE[ind]);
		vs4s9EE.push_back( eeclusterS4S9[ind] );
		vSeedTimeEE.push_back( SeedTime_v[ind] );
		Es_1.push_back( e1 ); Es_2.push_back( e2 );
#ifdef MVA_REGRESSIO_EE
		vs1s9EE.push_back( eeclusterS1S9[ind] );
		vs2s9EE.push_back( eeclusterS2S9[ind] );
		ESratio.push_back( deltaE/eeclus_iter->energy() );
#endif
	    }
	  }
	}
    }
    else{
	eseeclusters_tot.push_back( CaloCluster( eeclus_iter->energy(), eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
	Nxtal_tot.push_back(Ncristal_EE[ind]);
	vs4s9EE.push_back( eeclusterS4S9[ind] );
	vSeedTimeEE.push_back( SeedTime_v[ind] );
	Es_1.push_back( -999. ); Es_2.push_back( -999. );
	vs1s9EE.push_back( eeclusterS1S9[ind] );
	vs2s9EE.push_back( eeclusterS2S9[ind] );
	ESratio.push_back( (-1998.)/eeclus_iter->energy() );
    }
  }//end of the matching loop

  Ncristal_EE.clear();
  Ncristal_EE = Nxtal_tot; 
  Nxtal.clear();
  //cout << "Exiting fillEEClusters(...) " << endl;


}


TH1F** FillEpsilonPlot::initializeEpsilonHistograms(const char *name, const char *title, int size )
{

  TH1::SetDefaultSumw2(); // all new histograms will automatically activate the storage of the sum of squares of errors (i.e, TH1::Sumw2 is automatically called).

  TH1F **h = new TH1F*[size];
  std::string name_c = "";
  std::string title_c = "";
  if (isEoverEtrue_) std::cout << "FillEpsilonPlot::initializeEpsilonHistograms::isEoverEtrue_ = "            << isEoverEtrue_            << std::endl;
  else               std::cout << "FillEpsilonPlot::initializeEpsilonHistograms::useMassInsteadOfEpsilon_ = " << useMassInsteadOfEpsilon_ << std::endl;

  for(int jR=0; jR<size; jR++)
  {
    name_c = Form("%s%d", name, jR);
    title_c = Form("%s%d", title, jR);

    if (isEoverEtrue_) {

      h[jR] = new TH1F(name_c.c_str(), title_c.c_str(), 75, 0.0, 1.5); 
      h[jR]->GetXaxis()->SetTitle("photon E/E_{true}");

    } else {

      if(useMassInsteadOfEpsilon_)
	{
	  h[jR] = new TH1F(name_c.c_str(), title_c.c_str(), Are_pi0_? 70 : 80, Are_pi0_? 0.02:0.38, Are_pi0_? 0.30:0.70); // let's keep 0.004 GeV/bin
	  h[jR]->GetXaxis()->SetTitle("Mass(#gamma#gamma)");
	}
      else
	{
	  h[jR] = new TH1F(name_c.c_str(), title_c.c_str(), 120,-0.5,1);
	  h[jR]->GetXaxis()->SetTitle("Epsilon");
	}

    }

  }

  return h;

}


void  FillEpsilonPlot::deleteEpsilonPlot(TH1F **h, int size)
{
  for(int jR=0; jR<size; jR++)
    delete h[jR];

  delete h;
}


void  FillEpsilonPlot::writeEpsilonPlot(TH1F **h, const char *folder, int size)
{
  if (not outfile_->GetKey(folder)) outfile_->mkdir(folder);
  outfile_->cd(folder);
  for(int jR=0; jR<size; jR++)
    h[jR]->Write();
}

// std::vector< CaloCluster > FillEpsilonPlot::MCTruthAssociate(std::vector< CaloCluster > & clusters, double deltaR, bool isEB) {

//   // obsolete function used with old MC
//   // now using FillEpsilonPlot::MCTruthAssociateMultiPi0

//   std::vector< CaloCluster > ret;
//   ret.clear();
//   int n_tmp1 = -1;    int n_tmp2 = -1;
//   double deltaR1_tmp = 999;    double deltaR2_tmp = 999;
//   if(isMC_ && MC_Assoc_) {
//     //    std::cout << "Association with MC: initial collection size = " << clusters.size() << std::endl;
//     for(unsigned int i=0; i<clusters.size(); i++){
//       const CaloCluster g = clusters[i];
//       double deltaR1 = reco::deltaR(g.eta(),g.phi(), Gamma1MC.Eta(),Gamma1MC.Phi());
//       double deltaR2 = reco::deltaR(g.eta(),g.phi(), Gamma2MC.Eta(),Gamma2MC.Phi());
//       //      std::cout << "DR1,2 = " << deltaR1 << "  " << deltaR2 << std::endl;
//       if(deltaR1<deltaR1_tmp) {
//         deltaR1_tmp = deltaR1; n_tmp1 = i;
//       }
//       if(deltaR2<deltaR2_tmp) {
//         deltaR2_tmp = deltaR2; n_tmp2 = i;
//       }
//     }
//   } else return ret;

//   // std::cout << "deltaR1,2_tmp = " << deltaR1_tmp << "  " << deltaR2_tmp << std::endl;
//   // std::cout << "ntmp = " << n_tmp1 << "  " << n_tmp2 << std::endl;

//   // one could look to another close cluster, but if the two MC photons are closer than 0.05 to the same cluster, 
//   // that is merged in one cluster most probably, so another one would be fake
//   if(n_tmp1==n_tmp2) return ret;
  
//   if(deltaR1_tmp < deltaR) {
//     if(isEB) Ncristal_EB_used.push_back(Ncristal_EB[n_tmp1]);
//     else Ncristal_EE_used.push_back(Ncristal_EE[n_tmp1]);
//     ret.push_back(clusters[n_tmp1]);
//   }
//   if(deltaR2_tmp < deltaR) {
//     if(isEB) Ncristal_EB_used.push_back(Ncristal_EB[n_tmp2]);
//     else Ncristal_EE_used.push_back(Ncristal_EE[n_tmp2]);
//     ret.push_back(clusters[n_tmp2]);
//   }
//   return ret;
// }

int FillEpsilonPlot::getNumberOverlappingCrystals(std::vector<CaloCluster>::const_iterator g1, std::vector<CaloCluster>::const_iterator g2, const bool isEB = true) {

  // To count how many crystals overlap between the two 3x3 photon clusters we simply
  // open a geometric 3x3 matrix around the 2 seeds and see how many DetId are in common.
  // we also check that there is a RecHit in those overlapping crystals, otherwise there is no real overlap

  int nOverlapXtals = 0;
  std::vector<DetId> clus3x3_g1;       
  std::vector<DetId> clus3x3_g2;       
  std::set<DetId> DetIdUsed;

  if (isEB) {
  
    EBDetId  id_g1(g1->seed());
    EBDetId  id_g2(g2->seed());
    clus3x3_g1 = ebtopology_->getWindow(id_g1,3,3);       
    clus3x3_g2 = ebtopology_->getWindow(id_g2,3,3);       

  } else {

    EEDetId  id_g1(g1->seed());
    EEDetId  id_g2(g2->seed());
    clus3x3_g1 = eetopology_->getWindow(id_g1,3,3);       
    clus3x3_g2 = eetopology_->getWindow(id_g2,3,3);       

  }

  for (std::vector<DetId>::const_iterator det = clus3x3_g1.begin(); det != clus3x3_g1.end(); ++det) {
    DetIdUsed.insert(*det);
  }

  for (std::vector<DetId>::const_iterator det = clus3x3_g2.begin(); det != clus3x3_g2.end(); ++det) {
    EcalRecHitCollection::const_iterator rechit = isEB ? ebHandle->find( *det ) : eeHandle->find( *det );
    if ( (rechit != ebHandle->end()) || (rechit != eeHandle->end()) ) { 
      if (DetIdUsed.count(*det) != 0) nOverlapXtals++;
    }
  }

  return nOverlapXtals;

}


std::vector< CaloCluster > FillEpsilonPlot::MCTruthAssociateMultiPi0(std::vector< CaloCluster > & clusters, 
								     int& retNumberUnmergedGen,
								     int& retNumberMatchedGen,
								     std::vector<TLorentzVector*>& retClusters_matchedGenPhoton,
								     const double deltaR = 0.1, 
								     const bool isEB = true
								 ) 
{

  // std::cout << "### Entering FillEpsilonPlot::MCTruthAssociateMultiPi0" << std::endl;

  vector<TLorentzVector>* vecGamma1MC_ptr = (isEB ? &vecGamma1MC_EB : &vecGamma1MC_EE);
  vector<TLorentzVector>* vecGamma2MC_ptr = (isEB ? &vecGamma2MC_EB : &vecGamma2MC_EE);

  int numberUnmatchedGenPhotonPairs = 0;
  int numberGenPhotonsNotMerged = 0;

  // current strategy for the matching:
  // 1) loop on gen photons (use only the first one, the second has the same index)
  // 2) for each, find the closest reco CaloCluster in DR (loop twice to match both gen photons)
  // 3) the 2 clusters found in this way will be associated to the gen photons when we have to compute E/Etrue

  // for the moment, if we find the same reco cluster to be closer to both gen photons wrt other clusters, we just skip this event

  std::set<unsigned int> clusterIndexAlreadyUsed;
  std::vector< CaloCluster > ret;
  ret.clear();

  //    std::cout << "Association with MC: initial collection size = " << clusters.size() << std::endl;

  // gen quantities for DR computation: avoid recomputing them for each reco cluster
  double etaGen1 = 999;
  double etaGen2 = 999;
  double phiGen1 = 999;
  double phiGen2 = 999;
  double deltaR_gen = 999;

  //std::cout << "clusters.size() = " << clusters.size() << endl;

  // we loop on first gen photon because the index on the vector is the same for both photons from the same pi0
  // therefore, no need to loop on both gen photons
  //std::cout << "vecGamma1MC_ptr->size() = " << vecGamma1MC_ptr->size() << std::endl;
  for (unsigned int ig = 0; ig < vecGamma1MC_ptr->size(); ++ig) {    // size should be 20 when using MC made by Zhicai, because we have 20 pi0/event

    unsigned int n_tmp1 = -1;    unsigned int n_tmp2 = -1;
    double deltaR1_tmp = 999;    double deltaR2_tmp = 999;

    // avoid recomputing them for each reco cluster
    etaGen1 = vecGamma1MC_ptr->at(ig).Eta();
    etaGen2 = vecGamma2MC_ptr->at(ig).Eta();
    phiGen1 = vecGamma1MC_ptr->at(ig).Phi();
    phiGen2 = vecGamma2MC_ptr->at(ig).Phi();
    deltaR_gen = GetDeltaR(etaGen1,etaGen2,phiGen1,phiGen2);
    if (deltaR_gen < DR_FOR_UNMERGED_GEN_PHOTONS) continue;
    else numberGenPhotonsNotMerged++;

    // std::cout << std::endl;
    // std::cout << "ig = " << ig << "   gen eta1,2,phi1,2 " << etaGen1 << "," << etaGen2 << "," << phiGen1 << "," << phiGen2;
    // std::cout << "    DR = " << deltaR_gen << std::endl;

    // loop on reco cluster the first time to match the first gen photon
    for (unsigned int iclus = 0; iclus < clusters.size(); ++iclus) {

      if (clusterIndexAlreadyUsed.count(iclus) != 0 ) continue;
      const CaloCluster g = clusters[iclus];
      double deltaR1 = GetDeltaR(g.eta(), etaGen1, g.phi(), phiGen1);
      // std::cout << "iclus1 = " << iclus << "     deltaR1 = " << deltaR1 << std::endl;
      if (deltaR1 < deltaR1_tmp ) {
	deltaR1_tmp = deltaR1; 
	n_tmp1 = iclus;
      }
  
    } // loop on reco cluster the first time

    // now, if we find a reco cluster with DR lower than the chosen threshold, we remove the cluster index from the list of available cluster and
    // go on to match the second one
    // then, if the second gen photon cannot be matched to any reco photon, we will add the cluster index back in the list and remove the cluster from ret
    // in this way we will build only the list of cluster matched unambiguously to the gen level photon 
    // we also have to order the reco pair based on seed energy to make the flow consistent with pi0 calibration
    
    if(deltaR1_tmp < deltaR) {
      clusterIndexAlreadyUsed.insert(n_tmp1);
    } else {
      numberUnmatchedGenPhotonPairs++;
      continue;
    }

    // loop again on reco clusters to match the second gen photon
    for (unsigned int iclus = 0; iclus < clusters.size(); ++iclus) {

      if (clusterIndexAlreadyUsed.count(iclus) != 0 ) continue;
      const CaloCluster g = clusters[iclus];
      double deltaR2 = GetDeltaR(g.eta(), etaGen2, g.phi(), phiGen2);
      // std::cout << "iclus2 = " << iclus << "     deltaR2 = " << deltaR2 << std::endl;
      if( deltaR2 < deltaR2_tmp) {
	deltaR2_tmp = deltaR2; 
	n_tmp2 = iclus;
      }

    } // loop on reco cluster the second time 

    // if we arrived here, it means the first gen photon was matched to a reco cluster within a given DR
    // now, if we cannot matched the second photon, we consider the photon pair to be unmatched and we add back the n_tmp1 cluster index
    // otherwise we evaluate the ordering based on seed energy and fill the containers 

    if(deltaR2_tmp < deltaR) {

      bool g1_seed_energy_is_bigger = true;
      clusterIndexAlreadyUsed.insert(n_tmp2);

      if(isEB) {

	// find the rec hit of the seed and then get the energy
	EBRecHitCollection::const_iterator rechit_seed_g1 = ebHandle->find( clusters[n_tmp1].seed() );
	EBRecHitCollection::const_iterator rechit_seed_g2 = ebHandle->find( clusters[n_tmp2].seed() );
	if (rechit_seed_g1->energy() > rechit_seed_g2->energy()) {
	  g1_seed_energy_is_bigger = true;
	  Ncristal_EB_used.push_back(Ncristal_EB[n_tmp1]);
	  Ncristal_EB_used.push_back(Ncristal_EB[n_tmp2]);
	  ret.push_back(clusters[n_tmp1]);
	  ret.push_back(clusters[n_tmp2]);
	} else {
	  g1_seed_energy_is_bigger = false;
	  Ncristal_EB_used.push_back(Ncristal_EB[n_tmp2]);
	  Ncristal_EB_used.push_back(Ncristal_EB[n_tmp1]);
	  ret.push_back(clusters[n_tmp2]);
	  ret.push_back(clusters[n_tmp1]);
	}

      } else {

	// find the rec hit of the seed and then get the energy
	EERecHitCollection::const_iterator rechit_seed_g1 = eeHandle->find( clusters[n_tmp1].seed() );
	EERecHitCollection::const_iterator rechit_seed_g2 = eeHandle->find( clusters[n_tmp2].seed() );
	if (rechit_seed_g1->energy() > rechit_seed_g2->energy()) {
	  g1_seed_energy_is_bigger = true;
	  Ncristal_EE_used.push_back(Ncristal_EE[n_tmp1]);
	  Ncristal_EE_used.push_back(Ncristal_EE[n_tmp2]);
	  ret.push_back(clusters[n_tmp1]);
	  ret.push_back(clusters[n_tmp2]);
	} else {
	  g1_seed_energy_is_bigger = false;
	  Ncristal_EE_used.push_back(Ncristal_EE[n_tmp2]);
	  Ncristal_EE_used.push_back(Ncristal_EE[n_tmp1]);
	  ret.push_back(clusters[n_tmp2]);
	  ret.push_back(clusters[n_tmp1]);
	}

      }

      if (g1_seed_energy_is_bigger) {
	retClusters_matchedGenPhoton.push_back(&(vecGamma1MC_ptr->at(ig)));
	retClusters_matchedGenPhoton.push_back(&(vecGamma2MC_ptr->at(ig)));
	// retClusters_matchedGenPhotonEnergy.push_back(vecGamma1MC_ptr->at(ig).Energy());
	// retClusters_matchedGenPhotonEnergy.push_back(vecGamma2MC_ptr->at(ig).Energy());
      } else {
	retClusters_matchedGenPhoton.push_back(&(vecGamma2MC_ptr->at(ig)));
	retClusters_matchedGenPhoton.push_back(&(vecGamma1MC_ptr->at(ig)));
	// retClusters_matchedGenPhotonEnergy.push_back(vecGamma2MC_ptr->at(ig).Energy());
	// retClusters_matchedGenPhotonEnergy.push_back(vecGamma1MC_ptr->at(ig).Energy());
      }

    } else {
      clusterIndexAlreadyUsed.erase(n_tmp1); 
      numberUnmatchedGenPhotonPairs++;
      continue;
    }

    // now, if we arrived here it means we matched both photons, and according to our algorithm the two cluster indices should be different
    // anyway, we do the check and skip the gen pair if the same reco cluster is matched to both gen photons (it can happen if the two are really close to each other)
    // we must also remove the last two elements from the vectors and add the indices back to the set
    if(n_tmp1 == n_tmp2) {
      clusterIndexAlreadyUsed.erase(n_tmp1); // can erase the element just once because it was the same and a std::set only stores one element with the same value 
      if(isEB) {
	Ncristal_EB_used.pop_back();
	Ncristal_EB_used.pop_back();
      } else {
	Ncristal_EE_used.pop_back();
	Ncristal_EE_used.pop_back();
      }
      ret.pop_back();
      ret.pop_back();
      numberUnmatchedGenPhotonPairs++;
      continue;
    }

    // std::cout << "deltaR1,2_tmp = " << deltaR1_tmp << "  " << deltaR2_tmp << std::endl;
    // std::cout << "ntmp = " << n_tmp1 << "  " << n_tmp2 << std::endl;

    //cout << "Warning in FillEpsilonPlot::MCTruthAssociateMultiPi0: found the same CaloCluster closer to both gen photons. Will skip this gen photon pair." << endl;

  } // loop on gen photons

  // std::cout << "### Exiting FillEpsilonPlot::MCTruthAssociateMultiPi0: number of unmatched gen photon pairs is " << numberUnmatchedGenPhotonPairs << std::endl;
  // std::cout << "### Exiting FillEpsilonPlot::MCTruthAssociateMultiPi0: number of unmerged  gen photon pairs is " << numberGenPhotonsNotMerged << std::endl;
  retNumberUnmergedGen = numberGenPhotonsNotMerged;
  retNumberMatchedGen = numberGenPhotonsNotMerged-numberUnmatchedGenPhotonPairs;

  if (isEB) {
    h_numberUnmergedGenPhotonPairs_EB->Fill(((double)retNumberUnmergedGen)/vecGamma1MC_ptr->size());
    h_numberMatchedGenPhotonPairs_EB->Fill(((double)retNumberMatchedGen)/vecGamma1MC_ptr->size());
  } else {
    h_numberUnmergedGenPhotonPairs_EE->Fill(((double)retNumberUnmergedGen)/vecGamma1MC_ptr->size());
    h_numberMatchedGenPhotonPairs_EE->Fill(((double)retNumberMatchedGen)/vecGamma1MC_ptr->size());
  }

  return ret;

}


CaloCluster FillEpsilonPlot::getClusterAfterContainmentCorrections(std::vector<CaloCluster>::const_iterator gam, const bool isSecondPhoton = false, const bool isEB = true) {

  // this method is used to correct photon recHits energy based on containment corrections derived with E/Etrue in MC
  // we need to correct recHits and recompute energy and position

  // only EB for the moment!
  if (not isEB) return *gam;

  // corrected energy is obtained by correcting energy in each RecHit of photon
  float totalCorrectedClusterEnergy = 0.0;

  if (!hCC_EoverEtrue_g1 || !hCC_EoverEtrue_g2)
    throw cms::Exception("FillEpsilonPlot::getClusterAfterContainmentCorrections") << "Pointers to histograms with CC are null\n";
  TH2F* hContainmentCorrection = (isSecondPhoton ? hCC_EoverEtrue_g2 : hCC_EoverEtrue_g1);

  std::vector< std::pair<DetId, float> > hitsAndFrac = gam->hitsAndFractions();
  std::vector< std::pair<DetId, float> > correctedHitsAndFrac; // as hitsAndFrac but with corrections
  float correctedEnergy_it = 0;

  //  int ind = 1; // needed for cout below, only for debugging

  // FIXME: hardcoded correction: CC was obtained from V1 MC. We had a V2 MC but statistics was not enough
  // since the ratio is nearly flat, we divided CC from V1 MC by the mean of the ratio
  // it would be better to use a histogram for CC which is already scaled, but this is more straightforward if we want to add other corrections
  float CC_meanRatioV1overV2MC = isSecondPhoton ? 1.006 : 1.01;
  //float CC_meanRatioV1overV2MC = 1.0;


  for (std::vector< std::pair<DetId, float> >::const_iterator it  = hitsAndFrac.begin(); it != hitsAndFrac.end(); ++it) {	  

    EBDetId ebId(it->first);
    correctedEnergy_it = it->second * hContainmentCorrection->GetBinContent(hContainmentCorrection->FindFixBin(ebId.ieta(),ebId.iphi())) / CC_meanRatioV1overV2MC;

    // DEBUG
    // only for debugging
    // EBRecHitCollection::const_iterator ixtal = ebHandle->find( ebId );
    // if( ixtal == ebHandle->end() ) continue; // xtal not found
    // if (ixtal->energy() < 0) continue; // should not happen
    // std::cout << ind 
    // 	      << ":  ieta,iphi = " << ebId.ieta() << "," << ebId.iphi() 
    // 	      << "   RecHit_energy = " << ixtal->energy()
    // 	      << "   ChannelStatus = " <<
    // 	      << "   IC = " << regionalCalibration_->getCalibMap()->coeff(it->first)
    // 	      << "   RecHit_energy_IC = " << it->second
    // 	      << "   CC = " << hContainmentCorrection->GetBinContent(hContainmentCorrection->FindFixBin(ebId.ieta(),ebId.iphi()))
    // 	      << "   RecHit_energy_IC_CC = " << correctedEnergy_it << std::endl;

    totalCorrectedClusterEnergy += correctedEnergy_it;
    // in the rest of the code the fraction was defined as the energy of the RecHit, not the ratio with the total one
    correctedHitsAndFrac.push_back( std::make_pair(it->first, correctedEnergy_it));  
    //ind++;

  }

  // variables for position caculation
  float xclu(0.), yclu(0.), zclu(0.); // temp var to compute weighted average
  float total_weight(0.);// to compute position
  EBDetId seed_id(gam->seed());

  // Calculate shower depth
  float T0 = PCparams_.param_T0_barl_;
  float maxDepth = PCparams_.param_X0_ * ( T0 + log( totalCorrectedClusterEnergy ) ); 
  float maxToFront;
  if( GeometryFromFile_ ) maxToFront = geom_->getPosition(seed_id).mag(); // to front face
  else {
    const CaloCellGeometry* cell = geometry->getGeometry( seed_id ).get();
    GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
    maxToFront = posit.mag();
  }

  // loop over xtals (only those with positive energy) and compute energy and position
  for(unsigned int j = 0; j < correctedHitsAndFrac.size(); ++j) {

      EBDetId det(correctedHitsAndFrac[j].first);

      // compute position
      float weight = std::max( float(0.), PCparams_.param_W0_ + log(correctedHitsAndFrac[j].second/totalCorrectedClusterEnergy) );  // here it requires the fraction Ei/Etot
      float pos_geo;
      if( GeometryFromFile_ ) pos_geo = geom_->getPosition(det).mag(); // to front face
      else                  {
	const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	pos_geo = posit.mag();
      }
      float depth = maxDepth + maxToFront - pos_geo;
      GlobalPoint posThis;
      if( GeometryFromFile_ ) posThis = geom_->getPosition(det,depth);
      else{
	const CaloCellGeometry* cell = geometry->getGeometry(det).get();
	posThis = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( depth );
      }

      xclu += weight*posThis.x(); 
      yclu += weight*posThis.y(); 
      zclu += weight*posThis.z(); 
      total_weight += weight;

  } // loop over 3x3 rechits

  math::XYZPoint clusPos( xclu/total_weight, 
			  yclu/total_weight,
			  zclu/total_weight ); 

  CaloCluster correctedCaloCluster( totalCorrectedClusterEnergy, clusPos, CaloID(CaloID::DET_ECAL_BARREL), correctedHitsAndFrac, CaloCluster::undefined, seed_id );
  // cout << "before correction: " << *gam << endl;
  // cout << "after correction: " << correctedCaloCluster << endl;

  return correctedCaloCluster;

}

//========================================================


void FillEpsilonPlot::computeEpsilon(std::vector< CaloCluster > & clusters, std::vector<TLorentzVector*>& clusters_matchedGenPhoton, int subDetId ) 
{
  if(subDetId!=EcalBarrel && subDetId != EcalEndcap) 
    throw cms::Exception("FillEpsilonPlot::computeEpsilon") << "Subdetector Id not recognized\n";

  if (isDebug_) cout << "[DEBUG] Beginning cluster loop.."<< endl;

  // loop over clusters to make Pi0
  size_t i=0;
  for(std::vector<CaloCluster>::const_iterator g1  = clusters.begin(); g1 != clusters.end(); ++g1, ++i) 
  {
    size_t j=i+1;
    for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != clusters.end(); ++g2, ++j ) {

        if (isDebug_) cout << "\n[DEBUG] New Pair of Clusters"<< endl;

	if( subDetId==EcalBarrel ) {EventFlow_EB->Fill(4.); if (isDebug_) EventFlow_EB_debug->Fill(0.);}
	else                       {EventFlow_EE->Fill(4.); if (isDebug_) EventFlow_EE_debug->Fill(0.);}
	float Corr1 = 1., Corr2 = 1.;

	// g1 and g2 are ordered with the energy of the seed, but their respective clusters don't necessarily follow the same order
	// also, their pTs are not necessarily ordered 
	// Defining few variables to save photon quantities that are used more than once, to avoid recomputing them every time
	Double_t g1eta = g1->eta();
	Double_t g2eta = g2->eta();
	Double_t g1phi = g1->phi();
	Double_t g2phi = g2->phi();
	Double_t g1pt = g1->energy()/cosh(g1eta);
	Double_t g2pt = g2->energy()/cosh(g2eta);
	// following two object store the two photons ordered by pt
	TLorentzVector G_Sort_1, G_Sort_2, GSort1plus2;

	if( g1pt > g2pt ){
	  G_Sort_1.SetPtEtaPhiE( g1pt, g1eta, g1phi, g1->energy() );
	  G_Sort_2.SetPtEtaPhiE( g2pt, g2eta, g2phi, g2->energy() );
	}
	else{
	  G_Sort_1.SetPtEtaPhiE( g2pt, g2eta, g2phi, g2->energy() );
	  G_Sort_2.SetPtEtaPhiE( g1pt, g1eta, g1phi, g1->energy() );
	}

	GSort1plus2 = G_Sort_1 + G_Sort_2;
	  
#if !defined(NEW_CONTCORR) && defined(MVA_REGRESSIO) || defined(REGRESS_AND_PARAM_CONTCORR)
	if( subDetId==EcalBarrel && (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ){

	  // cout << "################################" << endl;
	  // cout << "### We are in the barrel! ###" << endl;
	  // cout << "################################" << endl;
	  
	  // following variable should be equivalent to transverse energy of the photon pair (for massless object it is equal to Pt() )
	  // this will store G.Energy()/cosh(G.Eta()), in order to compute it only once
	  //	  double GSort1plus2_EoverCoshEta = GSort1plus2.Energy()/cosh(GSort1plus2.Eta());  // currently not used here for EB

	  int ind1 = i, ind2 = j;
	  EBDetId  id_1(g1->seed()); int iEta1 = id_1.ieta(); int iPhi1 = id_1.iphi();
	  EBDetId  id_2(g2->seed()); int iEta2 = id_2.ieta(); int iPhi2 = id_2.iphi();
#ifdef MVA_REGRESSIO_Tree
	  int iSMod_1 = id_1.ism(); int iSMod_2 = id_2.ism();
#endif

	  bool Inverted=false;

          if( g1pt < g2pt ){
            iEta1=id_2.ieta(); iEta2 = id_1.ieta();
            iPhi1=id_2.iphi(); iPhi2 = id_1.iphi();
#ifdef MVA_REGRESSIO_Tree
            iSMod_1=id_2.ism(); iSMod_2=id_1.ism();
#endif
            ind1=j; ind2=i;
            Inverted=true;
          }

	  float Correct1(1.), Correct2(1.);
	  if(Are_pi0_){
	    float value_pi01[14];
	    //input list for regression in 2017 EB:
	    //enG_rec
	    //Nxtal
	    //S4S9
	    //S2S9
	    //iEta
	    //iPhi
	    //SM_dist: 
	    //M_dist:
	    float new_value_pi01[8]; 

	    value_pi01[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	    value_pi01[1] = ( G_Sort_1.Pt() );
	    value_pi01[2] = ( Ncristal_EB_used[ind1] );
	    value_pi01[3] = ( Ncristal_EB_used[ind2] );
	    value_pi01[4] = ( vs4s9[ind1] );
	    value_pi01[5] = ( vs1s9[ind1] );
	    value_pi01[6] = ( vs2s9[ind1] );
	    value_pi01[7] = ( iEta1 );
	    value_pi01[8] = ( iPhi1 );
	    value_pi01[9] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi01[10] = ( iEta1%5 );
	    value_pi01[11] = ( iPhi1%2 );
	    value_pi01[12] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
	    value_pi01[13] = ( iPhi1%20 );
	    
	    new_value_pi01[0] = ( g1->energy() );
            new_value_pi01[1] = ( Ncristal_EB_used[i] );
            new_value_pi01[2] = ( vs4s9[i] );
            new_value_pi01[3] = ( vs2s9[i] );
            new_value_pi01[4] = ( id_1.ieta() );
            new_value_pi01[5] = ( id_1.iphi() );
	    float temp_SM_dist_1 = ((id_1.iphi()-1)%20<10)*((id_1.iphi()-1)%20) + (((id_1.iphi()-1)%20)>=10)*(19-(id_1.iphi()-1)%20);
	    float temp_M_dist_1  = (abs(id_1.ieta())<=25)*(((abs(id_1.ieta())-1)%25<12)*((abs(id_1.ieta())-1)%25) + (((abs(id_1.ieta())-1)%25)>=12)*(24-(abs(id_1.ieta())-1)%25))
                               +(abs(id_1.ieta())>25) * (((abs(id_1.ieta())-26)%20<10)*((abs(id_1.ieta())-26)%20) + (((abs(id_1.ieta())-26)%20)>=10)*(19-(abs(id_1.ieta())-26)%20));
            new_value_pi01[6] = ( temp_SM_dist_1 );
            new_value_pi01[7] = ( temp_M_dist_1 );

	    //if( fabs((G_Sort_1+G_Sort_2).Eta())>1 ) value_pi01[14] = true ;
	    //else                                    value_pi01[14] = false ;
	    if(useMVAContainmentCorrections_)
            {
	      if(new_pi0ContainmentCorrections_)
                {
		  float Correct1_tmp = forestD_EB_1->GetResponse(new_value_pi01);
		  Correct1 = meanoffset + meanscale*TMath::Sin(Correct1_tmp);
		  // cout<<"DEBUG in FillEpsilonPlot.cc... computeEpsilon... new regression Correct1 = "<<Correct1<<endl;
                }
	      else
                {
		  Correct1 = forest_EB_1->GetResponse(value_pi01);
		  // cout<<"DEBUG in FillEpsilonPlot.cc... computeEpsilon... old regression Correct1 = "<<Correct1<<endl;
                }
	    }

	    float value_pi02[14];//#
	    float new_value_pi02[8];

	    value_pi02[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	    value_pi02[1] = ( G_Sort_2.Pt() );
	    value_pi02[2] = ( Ncristal_EB_used[ind1] );
	    value_pi02[3] = ( Ncristal_EB_used[ind2] );
	    value_pi02[4] = ( vs4s9[ind2] );
	    value_pi02[5] = ( vs1s9[ind2] );
	    value_pi02[6] = ( vs2s9[ind2] );
	    value_pi02[7] = ( iEta2 );
	    value_pi02[8] = ( iPhi2 );
	    value_pi02[9] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi02[10] = ( iEta2%5 );
	    value_pi02[11] = ( iPhi2%2 );
	    value_pi02[12] = ( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
	    value_pi02[13] = ( iPhi2%20 );
	    //if( fabs((G_Sort_1+G_Sort_2).Eta())>1 ) value_pi02[14] = true ;
	    //else                                    value_pi02[14] = false ;
	   
	    new_value_pi02[0] = ( g2->energy() );
            new_value_pi02[1] = ( Ncristal_EB_used[j] );
            new_value_pi02[2] = ( vs4s9[j] );
            new_value_pi02[3] = ( vs2s9[j] );
            new_value_pi02[4] = ( id_2.ieta() );
            new_value_pi02[5] = ( id_2.iphi() );
	    float temp_SM_dist_2 = ((id_2.iphi()-1)%20<10)*((id_2.iphi()-1)%20) + (((id_2.iphi()-1)%20)>=10)*(19-(id_2.iphi()-1)%20);
	    float temp_M_dist_2  = (abs(id_2.ieta())<=25)*(((abs(id_2.ieta())-1)%25<12)*((abs(id_2.ieta())-1)%25) + (((abs(id_2.ieta())-1)%25)>=12)*(24-(abs(id_2.ieta())-1)%25))
                               +(abs(id_2.ieta())>25) * (((abs(id_2.ieta())-26)%20<10)*((abs(id_2.ieta())-26)%20) + (((abs(id_2.ieta())-26)%20)>=10)*(19-(abs(id_2.ieta())-26)%20));
            new_value_pi02[6] = ( temp_SM_dist_2 );
            new_value_pi02[7] = ( temp_M_dist_2 );


	    if(useMVAContainmentCorrections_)
            {
            if(new_pi0ContainmentCorrections_)
                {
                float Correct2_tmp = forestD_EB_2->GetResponse(new_value_pi02);
                Correct2 = meanoffset + meanscale*TMath::Sin(Correct2_tmp);
                }
            else
                {
                Correct2 = forest_EB_2->GetResponse(value_pi02);
                }
             }
	  }
	  else{
	    float value_pi01[10];
 	    value_pi01[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	    value_pi01[1] = ( G_Sort_1.Pt() );
	    value_pi01[2] = ( Ncristal_EB_used[ind1] );
	    value_pi01[3] = ( iEta1 );
	    value_pi01[4] = ( iPhi1 );
	    value_pi01[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi01[6] = ( iEta1%5 );
	    value_pi01[7] = ( iPhi1%2 );
	    value_pi01[8] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
	    value_pi01[9] = ( iPhi1%20 );
	    Correct1 = forest_EB_1->GetResponse(value_pi01);
	    float value_pi02[10];
	    value_pi02[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	    value_pi02[1] = ( G_Sort_2.Pt() );
	    value_pi02[2] = ( Ncristal_EB_used[ind2] );
	    value_pi02[3] = ( iEta2 );
	    value_pi02[4] = ( iPhi2 );
	    value_pi02[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi02[6] = ( iEta2%5 );
	    value_pi02[7] = ( iPhi2%2 );
	    value_pi02[8] = ( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
	    value_pi02[9] = ( iPhi2%20 );
	    Correct2 = forest_EB_2->GetResponse(value_pi02);
	  }
	
	  // new regression trained in 2017: photon order is based on seed energy (energy in overlapped crystal taken by the first photon)
	  // old regression trained in 2012: photon order is based on photon Pt
	  if( (!Inverted) || new_pi0ContainmentCorrections_) { Corr1 = Correct1; Corr2 = Correct2; }
	  else           { Corr1 = Correct2; Corr2 = Correct1; }

	  //WARNIGN no CC for now! Put back in CKM 20/10/2015
	  //	  Corr1 = 1.; Corr2 = 1.; 
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
	  //In case ES give same posizion for different clusters
	  Correction1_mva = Correct1; Correction2_mva = Correct2;
	  iEta1_mva = iEta1; iEta2_mva = iEta2; iPhi1_mva = iPhi1; iPhi2_mva = iPhi2; Pt1_mva = G_Sort_1.Pt(); Pt2_mva = G_Sort_2.Pt();
	  iSM1_mva = iSMod_1; iSM2_mva = iSMod_2;

	  // Define just the original vector, the corrected vector has same direction and magnitude scaled by correction
	  // actually we just need G_sort_1 and G_sort_2, no need to define another object

	  // TLorentzVector mvag1P4; mvag1P4.SetPtEtaPhiE( Correct1*G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), Correct1*G_Sort_1.E() );
	  // TLorentzVector mvag2P4; mvag2P4.SetPtEtaPhiE( Correct2*G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), Correct2*G_Sort_2.E() );
	  // TLorentzVector mvaOrg1P4; mvaOrg1P4.SetPtEtaPhiE( G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), G_Sort_1.E() );
	  // TLorentzVector mvaOrg2P4; mvaOrg2P4.SetPtEtaPhiE( G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), G_Sort_2.E() );

	  MassOr_mva = GSort1plus2.M();
	  pi0Eta     = GSort1plus2.Eta();
	  // get corrected mass from sum of corrected photons
	  Mass_mva = (Correct1 * G_Sort_1 + Correct2 * G_Sort_2).M();
#endif
	}
#endif

#ifdef MVA_REGRESSIO_EE

	if( subDetId==EcalEndcap && (g1->seed().subdetId()==2) && (g2->seed().subdetId()==2) ){

	  // cout << "################################" << endl;
	  // cout << "### We are in the endcap! ###" << endl;
	  // cout << "################################" << endl;

	  // following variable should be equivalent to transverse energy of the photon pair (for massless object it is equal to Pt() )
	  // this will store G.E()/cosh(G.Eta()), in order to compute it only once
	  double GSort1plus2_EoverCoshEta = GSort1plus2.Energy()/cosh(GSort1plus2.Eta());

	  int ind1 = i, ind2 = j;
	  EEDetId  id_1(g1->seed()); int iX1 = id_1.ix(); int iY1 = id_1.iy();
	  EEDetId  id_2(g2->seed()); int iX2 = id_2.ix(); int iY2 = id_2.iy();

          if( g1pt < g2pt ){
            iX1=id_2.ix(); iX2 = id_1.ix();
            iY1=id_2.iy(); iY2 = id_1.iy();
            ind1=j; ind2=i;
          }

	  int EtaRing_1=GetRing( iX1, iY1, VectRing, false), EtaRing_2=GetRing( iX2, iY2, VectRing, false);
	  float value_pi01[10];
	    //input list for regression in 2017 EE:
	    //enG_rec
	    //Nxtal
	    //S4S9
	    //S2S9
	    //iX
	    //iY
	  float new_value_pi01[6];
	  value_pi01[0] = ( GSort1plus2_EoverCoshEta );
	  value_pi01[1] = ( G_Sort_1.Energy()/ GSort1plus2_EoverCoshEta );
	  value_pi01[2] = ( G_Sort_1.Pt() );
	  value_pi01[3] = ( Ncristal_EE_used[ind1] );
	  value_pi01[4] = ( Ncristal_EE_used[ind2] );
	  value_pi01[5] = ( vs4s9EE[ind1] );
	  value_pi01[6] = ( vs1s9EE[ind1] );
	  value_pi01[7] = ( vs2s9EE[ind1] );
	  value_pi01[8] = ( ESratio[ind1] );
	  value_pi01[9] = ( EtaRing_1 );
	  
	  new_value_pi01[0] = ( g1->energy() );
          new_value_pi01[1] = ( Ncristal_EE_used[i] );
          new_value_pi01[2] = ( vs4s9[i] );
          new_value_pi01[3] = ( vs2s9[i] );
          new_value_pi01[4] = ( id_1.ix() );
          new_value_pi01[5] = ( id_1.iy() );
	   
           float  Correct1 = 1.0;
          if(Are_pi0_ && useMVAContainmentCorrections_)
                {
                if(new_pi0ContainmentCorrections_)
                {
                float Correct1_tmp = forestD_EE_pi01->GetResponse(new_value_pi01);
                Correct1 = meanoffset + meanscale*TMath::Sin(Correct1_tmp);
                }
                else
                {
                Correct1 = forest_EE_pi01->GetResponse(value_pi01);
                }

                }

	  cout<<"Correction1: "<<Correct1<<" iX: "<<iX1<<" iY "<<iY1<<" Epi0 "<<GSort1plus2_EoverCoshEta
	    <<" ratio E "<< G_Sort_1.Energy()/GSort1plus2_EoverCoshEta<<" Pt "<<G_Sort_1.Pt()
	    <<" xtal "<<Ncristal_EE_used[ind1]<<" vs4s9EE "<<vs4s9EE[ind1]<<" vs1s9EE "<<vs1s9EE[ind1]<<" vs2s9EE "<<vs2s9EE[ind1]
	    <<" ESratio "<<ESratio[ind1]<<" EtaRing_1 "<<EtaRing_1<<endl;

	  float value_pi02[10];
	  float new_value_pi02[6];
	  value_pi02[0] = ( GSort1plus2_EoverCoshEta );
	  value_pi02[1] = ( G_Sort_2.Energy()/GSort1plus2_EoverCoshEta );
	  value_pi02[2] = ( G_Sort_2.Pt() );
	  value_pi02[3] = ( Ncristal_EE_used[ind1] );
	  value_pi02[4] = ( Ncristal_EE_used[ind2] );
	  value_pi02[5] = ( vs4s9EE[ind2] );
	  value_pi02[6] = ( vs1s9EE[ind2] );
	  value_pi02[7] = ( vs2s9EE[ind2] );
	  value_pi02[8] = ( ESratio[ind2] );
	  value_pi02[9] = ( EtaRing_2 );

	  new_value_pi02[0] = ( g2->energy() );
          new_value_pi02[1] = ( Ncristal_EE_used[j] );
          new_value_pi02[2] = ( vs4s9[j] );
          new_value_pi02[3] = ( vs2s9[j] );
          new_value_pi02[4] = ( id_2.ix() );
          new_value_pi02[5] = ( id_2.iy() );

	 float  Correct2 = 1.0;
          if(Are_pi0_ && useMVAContainmentCorrections_)
          {
                if(new_pi0ContainmentCorrections_)
                {
                float Correct2_tmp = forestD_EE_pi02->GetResponse(new_value_pi02);
                Correct2 = meanoffset + meanscale*TMath::Sin(Correct2_tmp);
                }
                else
                {
                Correct2 = forest_EE_pi02->GetResponse(value_pi02);
                }

          }

	  cout<<"Correction2: "<<Correct2<<" iX: "<<iX2<<" iY "<<iY2<<" Epi0 "<<GSort1plus2_EoverCoshEta
	    <<" ratio E "<< G_Sort_2.Energy()/GSort1plus2_EoverCoshEta<<" Pt "<<G_Sort_2.Pt()
	    <<" xtal "<<Ncristal_EE_used[ind2]<<" vs4s9EE "<<vs4s9EE[ind2]<<" vs1s9EE "<<vs1s9EE[ind2]<<" vs2s9EE "<<vs2s9EE[ind2]
	    <<" ESratio "<<ESratio[ind2]<<" EtaRing_1 "<<EtaRing_2<<endl;

	  // FIXME: should we uncomment these lines as for the barrel?
	  //  if( !Inverted ){ Corr1 = Correct1; Corr2 = Correct2; }
	  //  else           { Corr1 = Correct2; Corr2 = Correct1; }

	  //at least for 2017 regression:
	  if(new_pi0ContainmentCorrections_) { Corr1 = Correct1; Corr2 = Correct2; }

	  Correction1EE_mva = Correct1; Correction2EE_mva = Correct2;
	  iX1_mva = iX1; iX2_mva = iX2; iY1_mva = iY1; iY1_mva = iY2; Pt1EE_mva = G_Sort_1.Pt(); Pt2EE_mva = G_Sort_2.Pt();
	  EtaRing1_mva = EtaRing_1; EtaRing2_mva = EtaRing_2;

	  // Define just the original vector, the corrected vector has same direction and magnitude scaled by correction
	  // actually we just need G_sort_1 and G_sort_2, no need to define another object

	  // TLorentzVector mvag1P4; mvag1P4.SetPtEtaPhiE( Correct1*G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), Correct1*G_Sort_1.E() );
	  // TLorentzVector mvag2P4; mvag2P4.SetPtEtaPhiE( Correct2*G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), Correct2*G_Sort_2.E() );
	  // TLorentzVector mvaOrg1P4; mvaOrg1P4.SetPtEtaPhiE( G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), G_Sort_1.E() );
	  // TLorentzVector mvaOrg2P4; mvaOrg2P4.SetPtEtaPhiE( G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), G_Sort_2.E() );

	  // get corrected mass from sum of corrected photons
	  MassEEOr_mva = (GSort1plus2).M();
	  MassEE_mva = (Correct1 * G_Sort_1 + Correct2 * G_Sort_2).M();
	  TTree_JoshMva_EE->Fill();   
	}
#endif
	
	// uncorrected versions of photons
	// math::PtEtaPhiMLorentzVector g1P4_nocor( g1pt, g1eta, g1phi, 0. );
	// math::PtEtaPhiMLorentzVector g2P4_nocor( g2pt, g2eta, g2phi, 0. );
	//math::PtEtaPhiMLorentzVector pi0P4_nocor = g1P4_nocor + g2P4_nocor;
	// here the order is not important (while for the correction we apply corr1 to first cluster in list (not necessarily the leading in pt))
	double pi0P4_nocor_pt = GSort1plus2.Pt();
	double pi0P4_nocor_mass = GSort1plus2.M();
	//corrected version; note that Corr1 and Corr2 refers to first and second photon as selected looping on CaloCluster 
	// this means g1 is not necessarily the leading photon
	TLorentzVector pi0P4;
	TLorentzVector g1_contCorr_tlv;
	TLorentzVector g2_contCorr_tlv;

	if (useContainmentCorrectionsFromEoverEtrue_) {

	  // if (subDetId==EcalBarrel)                                                                                           
	  //   std::cout << "Barrel: npi0 = " << nPi0 << "   pi0,g1,g2 eta = " << GSort1plus2.Eta(),g1.eta(),g2.eta() << std::endl;                
	  // else                                                                                                                          
	  //   std::cout << "   Endcap: npi0 = " << nPi0 << "   pi0,g1,g2 eta = " << GSort1plus2.Eta(),g1.eta(),g2.eta() << std::endl;       
	  CaloCluster g1_contCorr(getClusterAfterContainmentCorrections(g1,false,subDetId==EcalBarrel));
	  CaloCluster g2_contCorr(getClusterAfterContainmentCorrections(g2,true,subDetId==EcalBarrel));

	  g1_contCorr_tlv.SetPtEtaPhiE(g1_contCorr.energy()/cosh(g1_contCorr.eta()),g1_contCorr.eta(),g1_contCorr.phi(),g1_contCorr.energy());
	  g2_contCorr_tlv.SetPtEtaPhiE(g2_contCorr.energy()/cosh(g2_contCorr.eta()),g2_contCorr.eta(),g2_contCorr.phi(),g2_contCorr.energy());                            

	} else {

	  if (g1pt > g2pt) {
	    g1_contCorr_tlv = Corr1 * G_Sort_1;
	    g2_contCorr_tlv = Corr2 * G_Sort_2;
	  } else {
	    // when g1pt < g2pt, G_Sort_1 is made with g2, and Corr2 must be applied to it                      
	    g1_contCorr_tlv = Corr1 * G_Sort_2;
	    g2_contCorr_tlv = Corr2 * G_Sort_1; 
	  }
	}

	pi0P4 = g1_contCorr_tlv + g2_contCorr_tlv;                                                                                                                      

	// eta, pt, phi of corrected photons are used many times. Since their computation is tipically time consuming, store them in doubles for later usage
	double pi0P4_pt = pi0P4.Pt();
	double pi0P4_eta = pi0P4.Eta();
	double pi0P4_phi = pi0P4.Phi();
	double pi0P4_mass = pi0P4.M();
	//std::cout << "pio mass: w/o corr, w/ corr --> " << pi0P4_nocor_mass << "," << pi0P4_mass << endl;  

	// note that photon eta and phi are not modified by correction (only pT) since Corr * vector modifies the cartesian coordinates of the vector (pT and pZ)
	// double g1P4_eta = g1eta; 
	// double g1P4_phi = g1phi;
	// double g2P4_eta = g2eta;
	// double g2P4_phi = g2phi;

	//cout << "pi0P4_nocor_mass, pi0P4_mass " << pi0P4_nocor_mass << "  " << pi0P4_mass << endl;

	//In case ES give same posizion for different clusters
	if( pi0P4_nocor_mass<0.03 && pi0P4_mass < 0.03 ) continue;

#ifdef SELECTION_TREE
	if( subDetId == EcalBarrel ){ 
	  Fill_PtPi0_EB( pi0P4_pt );
	  Fill_mpi0_EB( pi0P4_mass );
	  Fill_etapi0_EB( pi0P4_eta );
	  Fill_phipi0_EB( pi0P4_phi );
	  //adding other variables  WARNING: MUST STILL ADD TO TTREE DEFINITION 
	  Fill_PtGamma_EB( g1pt * Corr1, g2pt * Corr2 );
	  Fill_EtaGamma_EB( g1eta, g2eta );
	  // to be implemented
	  if(isMC_ && MC_Assoc_) Fill_NcrystalUsedGamma_EB(Ncristal_EB_used[0], Ncristal_EB_used[1]);
	  Fill_S4S9Gamma_EB(vs4s9[0], vs4s9[1]);
	  //Fill_NxtalEnergGamma_EB(Nxtal_EnergGamma);  
	  //Fill_NxtalEnergGamma2_EB(Nxtal_EnergGamma2);
	  //
	  //the difference of Nxtal_EnergGamma wrt Ncristal_EB_used is that Nxtal_EnergGamma is set equal to Ncristal_EB[i]. For real data it is the same because Ncristal_EB_used is set equal to Ncristal_EB, but for MC it can be different due to the MC truth matching, which selects a subset of Ncristal_EB
	  //
	  //
	  Fill_Epsilon_EB( 0.5 * ( pow(pi0P4_mass/PI0MASS,2)  - 1. ) );
	  Pi0Info_EB->Fill();
	}
	if( subDetId == EcalEndcap ){
	  Fill_PtPi0_EE( pi0P4_pt );
	  Fill_mpi0_EE( pi0P4_mass );
	  Fill_etapi0_EE( pi0P4_eta );
	  Fill_phipi0_EE( pi0P4_phi );
	  //adding other variables  WARNING: MUST STILL ADD TO TTREE DEFINITION 
	  Fill_PtGamma_EE( g1pt * Corr1, g2pt * Corr2 );
	  Fill_EtaGamma_EE( g1eta, g2eta );
	  // to be implemented
	  if(isMC_ && MC_Assoc_) Fill_NcrystalUsedGamma_EE(Ncristal_EE_used[0], Ncristal_EE_used[1]);
	  Fill_S4S9Gamma_EE(vs4s9[0], vs4s9[1]);
	  //Fill_NxtalEnergGamma_EE(Nxtal_EnergGamma);  //which is the difference wrt Ncristal_EE_used ?!? 
	  //Fill_NxtalEnergGamma2_EE(Nxtal_EnergGamma2);
	  //
	  Fill_Epsilon_EE( 0.5 * ( pow(pi0P4_mass/PI0MASS,2)  - 1. ) );
	  Pi0Info_EE->Fill();
	}
#endif

	

	if (isDebug_) cout << "[DEBUG] Apply kinematic selection cuts" << endl;

	if( g1eta == g2eta && g1phi == g2phi ) continue;

	// pi0/eta pT cut
	if (subDetId == EcalBarrel) {

	  if (fabs(pi0P4_eta)<1.0)       { if( pi0P4_nocor_pt < pi0PtCut_low_[subDetId]) continue; }
	  else if (fabs(pi0P4_eta)<1.479) { if( pi0P4_nocor_pt < pi0PtCut_high_[subDetId]) continue; }
	  if (isDebug_) EventFlow_EB_debug->Fill(1.);

	} else {
	  
	  if (fabs(pi0P4_eta)<1.8 )     { if( pi0P4_nocor_pt < pi0PtCut_low_[subDetId]) continue; }	  
	  else                          { if( pi0P4_nocor_pt < pi0PtCut_high_[subDetId]) continue; }
	  if (isDebug_) EventFlow_EE_debug->Fill(1.);

	}

	float nextClu = 999., Drtmp = 999.;
	for(size_t ind=0; ind<clusters.size(); ++ind){
	  const CaloCluster* Gtmp = &(clusters[ind]);
	  double deltaR1 = GetDeltaR(Gtmp->eta(),g1eta,Gtmp->phi(),g1phi);
	  double deltaR2 = GetDeltaR(Gtmp->eta(),g2eta,Gtmp->phi(),g2phi);
	  if( ind!=i && ind!=j && (deltaR1<Drtmp || deltaR2<Drtmp ) ){
	    nextClu = min(deltaR1,deltaR2);
	    Drtmp = nextClu;
	  }
	}

	// pi0/eta isolation cut (distance to other clusters)
	if (subDetId == EcalBarrel) {

	  if (fabs(pi0P4_eta)<1.0)       { if( nextClu<pi0IsoCut_low_[subDetId] ) continue; }
	  else if (fabs(pi0P4_eta)<1.479) { if( nextClu<pi0IsoCut_high_[subDetId] ) continue; }
	  if (isDebug_) EventFlow_EB_debug->Fill(2.);

	} else {
	  
	  if (fabs(pi0P4_eta)<1.8 )     { if( nextClu<pi0IsoCut_low_[subDetId] ) continue; }	  
	  else                          { if( nextClu<pi0IsoCut_high_[subDetId] ) continue; }
	  if (isDebug_) EventFlow_EE_debug->Fill(2.);

	}

	// Implementation of HLT Filter Isolation - Eta Band Isolation 
	// implemented in HLT: CMSSW_7_1_0/src/HLTrigger/special/src/HLTEcalResonanceFilter.cc
	// see Yong Yang's  Thesis: http://thesis.library.caltech.edu/7345/

	if (isDebug_) cout << "[DEBUG] Running HLT Isolation" << endl;

	float hlt_iso = 0;
	for(size_t ind=0; ind < clusters.size(); ++ind){
	  if( clusters[ind].seed() == clusters[i].seed() || clusters[ind].seed() == clusters[j].seed()) continue;
	  const CaloCluster* Gtmp = &(clusters[ind]);
	  TLorentzVector GtmpP4;  
	  GtmpP4.SetPtEtaPhiE(Gtmp->energy()/cosh(Gtmp->eta()), Gtmp->eta(), Gtmp->phi(), Gtmp->energy());
	  if (GtmpP4.Pt() < 0.5) continue;  // FIXME: based on the stream, it should represent "ptMinForIsolation*"
	  // delta R from the pi0 candidates
	  double deltaR0 = GetDeltaR(Gtmp->eta(), pi0P4_eta, Gtmp->phi(), pi0P4_phi);
	  if (deltaR0  > ((Are_pi0_) ? 0.2:0.3)) continue;
	  // cluster must be inside of an eta strip 
	  double deta = fabs(Gtmp->eta() - pi0P4_eta); 
	  if (deta > ((Are_pi0_) ? 0.05:0.1)) continue;
	  hlt_iso += GtmpP4.Pt();
	}
	// the cut is taken relative to the pi0 pt
	hlt_iso /= pi0P4_nocor_pt;
	//category break down of cuts
	// pi0/eta isolation cut
	if (subDetId == EcalBarrel) {

	  if (fabs(pi0P4_eta)<1.0)       { if( hlt_iso > pi0HLTIsoCut_low_[subDetId]  && CutOnHLTIso_ ) continue; }
	  else if (fabs(pi0P4_eta)<1.479) { if( hlt_iso > pi0HLTIsoCut_high_[subDetId] && CutOnHLTIso_ ) continue; }
	  if (isDebug_) EventFlow_EB_debug->Fill(3.);

	} else {
	  
	  if (fabs(pi0P4_eta)<1.8 )     { if( hlt_iso > pi0HLTIsoCut_low_[subDetId]  && CutOnHLTIso_ ) continue; }	  
	  else                          { if( hlt_iso > pi0HLTIsoCut_high_[subDetId] && CutOnHLTIso_ ) continue; }
	  if (isDebug_) EventFlow_EE_debug->Fill(3.);

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////

	if (isDebug_) cout << "[DEBUG] N Cristal Cuts" << endl;


	int Nxtal_EnergGamma = 0;
	int Nxtal_EnergGamma2 = 0;
	int Nxtal_g1 = 0;
	int Nxtal_g2 = 0;
	// for the selection we use the more energetic photon, but in general we want to know how many crystals thw two clusters have.
	// Now, the clusters are ordered based on the energy of their seed, so that the second photon (with less energetic seed) should tipically have less crystals
	// index i refers to leading seed photon, j to the oher one
	if(subDetId==EcalEndcap){

	  Nxtal_g1 = Ncristal_EE_used[i];
	  Nxtal_g2 = Ncristal_EE_used[j];
	  if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EE_used[i]; Nxtal_EnergGamma2 = Ncristal_EE_used[j]; }
	  else                           {  Nxtal_EnergGamma = Ncristal_EE_used[j]; Nxtal_EnergGamma2 = Ncristal_EE_used[i]; }

	  if( fabs(pi0P4_eta)<1.8 ) { 
	    if( Nxtal_EnergGamma < nXtal_1_cut_low_[subDetId] ) continue; 
	    if( Nxtal_EnergGamma2 < nXtal_2_cut_low_[subDetId] ) continue;
	  } else {
	    if( Nxtal_EnergGamma < nXtal_1_cut_high_[subDetId] ) continue; 
	    if( Nxtal_EnergGamma2 < nXtal_2_cut_high_[subDetId] ) continue;
	  }
	  if (isDebug_) EventFlow_EE_debug->Fill(4.);
	  EventFlow_EE->Fill(5.);

	}
	else{

	  Nxtal_g1 = Ncristal_EB_used[i];
	  Nxtal_g2 = Ncristal_EB_used[j];
	  if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EB_used[i]; Nxtal_EnergGamma2 = Ncristal_EB_used[j]; }
	  else                           {  Nxtal_EnergGamma = Ncristal_EB_used[j]; Nxtal_EnergGamma2 = Ncristal_EB_used[i]; }

	  if( fabs(pi0P4_eta)<1.0 ) { 
	    if( Nxtal_EnergGamma < nXtal_1_cut_low_[subDetId] ) continue; 
	    if( Nxtal_EnergGamma2 < nXtal_2_cut_low_[subDetId] ) continue;
	  } else if( fabs(pi0P4_eta)<1.479 )  { 
	    if( Nxtal_EnergGamma < nXtal_1_cut_high_[subDetId] ) continue; 
	    if( Nxtal_EnergGamma2 < nXtal_2_cut_high_[subDetId] ) continue;
	  } 

	  pi0MassVsIetaEB->Fill( fabs(pi0P4_eta)/0.0174, pi0P4_mass);
	  pi0MassVsETEB->Fill(pi0P4_pt, pi0P4_mass);
	  if (isDebug_) EventFlow_EB_debug->Fill(4.);
	  EventFlow_EB->Fill(5.);

	}

	if (isDebug_) cout << "[DEBUG] Fill Optimization Variables..." << endl;

	if (fillKinematicVariables_) {

	  if (fabs(pi0P4_eta)<1.0) whichRegionEcalStreamPi0 = 0;
	  else if (fabs(pi0P4_eta)<1.479) whichRegionEcalStreamPi0 = 1;
	  else if (fabs(pi0P4_eta)<1.8) whichRegionEcalStreamPi0 = 2;
	  else whichRegionEcalStreamPi0 = 3;

	  pi0pt_afterCuts[whichRegionEcalStreamPi0]->Fill(pi0P4_nocor_pt);
	  g1pt_afterCuts[whichRegionEcalStreamPi0]->Fill(g1pt);
	  g2pt_afterCuts[whichRegionEcalStreamPi0]->Fill(g2pt);
	  g1Nxtal_afterCuts[whichRegionEcalStreamPi0]->Fill(Nxtal_g1);
	  g2Nxtal_afterCuts[whichRegionEcalStreamPi0]->Fill(Nxtal_g2);
	  pi0PhotonsNoverlappingXtals_afterCuts[whichRegionEcalStreamPi0]->Fill(getNumberOverlappingCrystals(g1,g2,subDetId==EcalBarrel));
	  if (isMC_) {
	    pi0MassVsPU[whichRegionEcalStreamPi0]->Fill(pi0P4_nocor_mass,nPUobs_BX0_);
	  }

	}


	//Fill Optimization
	if( MakeNtuple4optimization_ && pi0P4_mass > ((Are_pi0_)?0.03:0.2) && pi0P4_mass < ((Are_pi0_)?0.25:1.) ) {

	  //FIXME: check how the tree is filled when using E/Etrue corrections (evaluate to save both corrected and uncorrected photons, because also position is changed)
	  // add in case a flag saying which corrections are used
	  if( nPi0>NPI0MAX-2 ) { 
	    cout<<"nPi0::TOO MANY PI0: ("<<nPi0<<")!!!"<<endl; 
	  } else{
	    Op_Pi0recIsEB.push_back(  (subDetId==EcalBarrel)? 1:0);
	    Op_ClusIsoPi0.push_back(  nextClu);  
	    Op_HLTIsoPi0.push_back(   hlt_iso);
	    Op_nCrisG1.push_back(     Nxtal_g1); 
	    Op_nCrisG2.push_back(     Nxtal_g2);
	    Op_enG1_cor.push_back(    g1_contCorr_tlv.Energy());
	    Op_enG2_cor.push_back(    g2_contCorr_tlv.Energy());
	    Op_etaG1_cor.push_back(   g1_contCorr_tlv.Eta());
	    Op_etaG2_cor.push_back(   g2_contCorr_tlv.Eta());
	    Op_phiG1_cor.push_back(   g1_contCorr_tlv.Phi());
	    Op_phiG2_cor.push_back(   g2_contCorr_tlv.Phi());
	    Op_mPi0_cor.push_back(    pi0P4_mass);
	    Op_etaPi0_cor.push_back(  pi0P4_eta);
	    Op_ptPi0_cor.push_back(   pi0P4_pt);
	    Op_DeltaRG1G2.push_back(  GetDeltaR( g1eta, g2eta, g1phi, g2phi ));
	    Op_ptPi0_nocor.push_back( pi0P4_nocor_pt);
	    Op_mPi0_nocor.push_back(  pi0P4_nocor_mass);
	    Op_Es_e1_1.push_back(     (subDetId==EcalBarrel) ? 0. : Es_1[i]);
	    Op_Es_e1_2.push_back(     (subDetId==EcalBarrel) ? 0. : Es_1[j]);
	    Op_Es_e2_1.push_back(     (subDetId==EcalBarrel) ? 0. : Es_2[i]);
	    Op_Es_e2_2.push_back(     (subDetId==EcalBarrel) ? 0. : Es_2[j]);
	    Op_S4S9_1.push_back(      (subDetId==EcalBarrel) ? vs4s9[i] : vs4s9EE[i]);
	    Op_S4S9_2.push_back(      (subDetId==EcalBarrel) ? vs4s9[j] : vs4s9EE[j]);
	    Op_S2S9_1.push_back(      (subDetId==EcalBarrel) ? vs2s9[i] : vs2s9EE[i]);
	    Op_S2S9_2.push_back(      (subDetId==EcalBarrel) ? vs2s9[j] : vs2s9EE[j]);
	    Op_S1S9_1.push_back(      (subDetId==EcalBarrel) ? vs1s9[i] : vs1s9EE[i]);
	    Op_S1S9_2.push_back(      (subDetId==EcalBarrel) ? vs1s9[j] : vs1s9EE[j]);
	    Op_enG1_nocor.push_back(  g1->energy()); // g1P4_nocor.E();
	    Op_enG2_nocor.push_back(  g2->energy()); // g2P4_nocor.E();
	    Op_etaG1_nocor.push_back( g1eta);
            Op_etaG2_nocor.push_back( g2eta);
	    Op_phiG1_nocor.push_back( g1phi);
            Op_phiG2_nocor.push_back( g2phi);
	    Op_Time_1.push_back(      (subDetId==EcalBarrel) ? vSeedTime[i] : vSeedTimeEE[i]);
	    Op_Time_2.push_back(      (subDetId==EcalBarrel) ? vSeedTime[j] : vSeedTimeEE[j]);
	    if( isMC_ && MC_Assoc_ ) {
	      Op_enG1_true.push_back(   clusters_matchedGenPhoton[i]->Energy());
	      Op_enG2_true.push_back(   clusters_matchedGenPhoton[j]->Energy());
	      Op_DeltaR_1.push_back(    GetDeltaR(g1eta, clusters_matchedGenPhoton[i]->Eta(), g1phi, clusters_matchedGenPhoton[i]->Phi()));
	      Op_DeltaR_2.push_back(    GetDeltaR(g2eta, clusters_matchedGenPhoton[j]->Eta(), g2phi, clusters_matchedGenPhoton[j]->Phi()));
	    }

            if( (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ) {

              EBDetId  id_1(g1->seed()); int iEta1 = id_1.ieta(); int iPhi1 = id_1.iphi();
              EBDetId  id_2(g2->seed()); int iEta2 = id_2.ieta(); int iPhi2 = id_2.iphi();

              Op_iEtaiX_1.push_back(   iEta1);
              Op_iEtaiX_2.push_back(   iEta2);
              Op_iPhiiY_1.push_back(   iPhi1);
              Op_iPhiiY_2.push_back(   iPhi2);
              Op_iEta_1on5.push_back(  iEta1%5);
              Op_iEta_2on5.push_back(  iEta2%5);
              Op_iPhi_1on2.push_back(  iPhi1%2);
              Op_iPhi_2on2.push_back(  iPhi2%2);
              Op_iEta_1on2520.push_back( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20)); //Distance in xtal from module boundaries
              Op_iEta_2on2520.push_back((TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20));
              Op_iPhi_1on20.push_back(  iPhi1%20);
              Op_iPhi_2on20.push_back(  iPhi2%20);

            } else if( (g1->seed().subdetId()==2) && (g2->seed().subdetId()==2) ) {
              
              EEDetId  id_1(g1->seed()); int iX1 = id_1.ix(); int iY1 = id_1.iy();
              EEDetId  id_2(g2->seed()); int iX2 = id_2.ix(); int iY2 = id_2.iy();

              Op_iEtaiX_1.push_back(      (iX1 < 50) ? iX1 : 100-iX1);
              Op_iEtaiX_2.push_back(      (iX2 < 50) ? iX2 : 100-iX2);
              Op_iPhiiY_1.push_back(      (iY1 < 50) ? iY1 : 100-iY1);
              Op_iPhiiY_2.push_back(      (iY2 < 50) ? iY2 : 100-iY2);
              Op_iEta_1on5.push_back(     999);
              Op_iEta_2on5.push_back(     999);
              Op_iPhi_1on2.push_back(     999);
              Op_iPhi_2on2.push_back(     999);
              Op_iEta_1on2520.push_back(  999);
              Op_iEta_2on2520.push_back(  999);
              Op_iPhi_1on20.push_back(    999);
              Op_iPhi_2on20.push_back(    999);
              
            } else {
              Op_iEtaiX_1.push_back(      -999);
              Op_iEtaiX_2.push_back(      -999);
              Op_iPhiiY_1.push_back(      -999);
              Op_iPhiiY_2.push_back(      -999);
              Op_iEta_1on5.push_back(     -999);
              Op_iEta_2on5.push_back(     -999);
              Op_iPhi_1on2.push_back(     -999);
              Op_iPhi_2on2.push_back(     -999);
              Op_iEta_1on2520.push_back(  -999);
              Op_iEta_2on2520.push_back(  -999);
              Op_iPhi_1on20.push_back(    -999);
              Op_iPhi_2on20.push_back(    -999);
            }
	    nPi0++;
	  }

	}
	
	if (isDebug_) cout << "[DEBUG] End Accessing Optmization Variables..." << endl;

	//Check the Conteinment correction for Barrel
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
	if( pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.28:0.75) ){
	  if( subDetId==EcalBarrel && (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ) TTree_JoshMva->Fill();
	}
#endif

	if (!MakeNtuple4optimization_) {

	  if (isDebug_) cout << "[DEBUG] computing region weights" << endl; 

	  // compute region weights
	  RegionWeightVector w1 = regionalCalibration_->getWeights( &(*g1), subDetId ); // region weights W_j^k for clu1
	  RegionWeightVector w2 = regionalCalibration_->getWeights( &(*g2), subDetId ); // region weights W_j^k for clu2

	  // append w2 to w1
	  w1.insert( w1.end(), w2.begin(), w2.end() );

	  float r2 = pi0P4_mass/PI0MASS;
	  r2 = r2*r2;
	  //average <eps> for cand k
	  float eps_k = 0.5 * ( r2 - 1. );
	  // compute quantities needed for <eps>_j in each region j
	  if(subDetId==EcalBarrel) allEpsilon_EBnw->Fill( pi0P4_mass );
	  else                     allEpsilon_EEnw->Fill( pi0P4_mass );
	  for(RegionWeightVector::const_iterator it = w1.begin(); it != w1.end(); ++it) {
	    const uint32_t& iR = (*it).iRegion;
	    const float& w = (*it).value;

	    if(subDetId==EcalBarrel){
		if( pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.23:0.7) ){
		  if( !EtaRingCalibEB_ && !SMCalibEB_ ) epsilon_EB_h[iR]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w );
		  allEpsilon_EB->Fill( pi0P4_mass, w );
		  std::vector<DetId> mioId(regionalCalibration_->allDetIdsInEERegion(iR));
		  //allDetIdsInEERegion is not reliable for EB and probably wrong. Getting iEta and iPhi elsewhere
		  std::map<int,vector<int>>::iterator it; 
		  int iEta = List_IR_EtaPhi.find(iR)->second[0]; int iPhi = List_IR_EtaPhi.find(iR)->second[1]; int iSM = List_IR_EtaPhi.find(iR)->second[2];
		  entries_EB->Fill( iEta, iPhi, w );
		  //If Low Statistic fill all the Eta Ring
		  if( EtaRingCalibEB_ ){
		    for(auto const &iterator : ListEtaFix_xtalEB){
			if( iterator.first == iEta ){ 
			  for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EB_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			}
		    }
		  }
		  if( SMCalibEB_ ){
		    for(auto const &iterator : ListSMFix_xtalEB){
			if( iterator.first == iSM ){ 
			  for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EB_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			}
		    }
		  }
		  //		  for(unsigned int i=0; i<mioId.size(); ++i){//Actually size is 1 for this loop, it is just to access the recHit
		  //
		  //		    EBDetId tmp_id(mioId.at(i));
		  //		    entries_EB->Fill( tmp_id.ieta(), tmp_id.iphi(), w );
		  //cout<<"His iEta and iPhi is : "<<tmp_id.ieta()<<" "<<tmp_id.iphi()<<" iR "<<iR<<" "<<mioId.at(i).rawId()<<endl;
		  //		    //If Low Statistic fill all the Eta Ring
		  //		    if( EtaRingCalib_ ){
		  //		      for(auto const &iterator : ListEtaFix_xtalEB){
		  //			  if( iterator.first == tmp_id.ieta() ){ 
		  //			    for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EB_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
		  //			  }
		  //			}
		  //		    }
		  //		  }
		}
	    }
	    else {
		if( pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.28:0.75) ){
		  if( !EtaRingCalibEE_ && !SMCalibEE_ ) epsilon_EE_h[iR]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w );
		  allEpsilon_EE->Fill( pi0P4_mass, w );
		  std::vector<DetId> mioId(regionalCalibration_->allDetIdsInEERegion(iR));
		  //allDetIdsInEERegion is not reliable for EE. Getting ix and iy elsewhere
		  std::map<int,vector<int>>::iterator it; 
		  int iX = List_IR_XYZ.find(iR)->second[0]; int iY = List_IR_XYZ.find(iR)->second[1]; int iZ = List_IR_XYZ.find(iR)->second[2]; int Quad = List_IR_XYZ.find(iR)->second[3];
		  if( iZ==-1 ){
		    entries_EEm->Fill( iX, iY, w );
		    //If Low Statistic fill all the Eta Ring
		    if( EtaRingCalibEE_ ){
			for(auto const &iterator : ListEtaFix_xtalEEm){
			  if( iterator.first == GetRing( iX, iY, VectRing,false) ){ 
			    for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			  }
			}
		    }
		    if( SMCalibEE_ ){
			for(auto const &iterator : ListQuadFix_xtalEEm){
			  if( iterator.first == Quad ){ 
			    for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			  }
			}
		    }
		  }
		  else{
		    entries_EEp->Fill( iX, iY, w );
		    //If Low Statistic fill all the Eta Ring
		    if( EtaRingCalibEE_ ){
			for(auto const &iterator : ListEtaFix_xtalEEp){
			  if( iterator.first == GetRing( iX, iY, VectRing,false) ){
			    for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			  }
			}
		    }
		    if( SMCalibEE_ ){
			for(auto const &iterator : ListQuadFix_xtalEEp){
			  if( iterator.first == Quad ){ 
			    for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
			  }
			}
		    }
		  }
		  //		  for(unsigned int i=0; i<mioId.size(); ++i){//Actually size is 1 for this loop, it is just to access the recHit
		  //		    EEDetId tmp_id(mioId.at(i));
		  //		    if( tmp_id.zside()==-1 ){
		  //			entries_EEm->Fill( tmp_id.ix(), tmp_id.iy(), w );
		  //			//If Low Statistic fill all the Eta Ring
		  //			if( EtaRingCalib_ ){
		  //			  for(auto const &iterator : ListEtaFix_xtalEEm){
		  //			    if( iterator.first == GetRing( tmp_id.ix(),tmp_id.iy(),VectRing,false) ){ 
		  //				for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
		  //			    }
		  //			  }
		  //			}
		  //		    }
		  //		    else{
		  //			entries_EEp->Fill( tmp_id.ix(), tmp_id.iy(), w );
		  //			//If Low Statistic fill all the Eta Ring
		  //			if( EtaRingCalib_ ){
		  //			  for(auto const &iterator : ListEtaFix_xtalEEp){
		  //			    if( iterator.first == GetRing( tmp_id.ix(),tmp_id.iy(),VectRing,false) ){ 
		  //				for(unsigned int iRtmp=0; iRtmp<iterator.second.size(); iRtmp++){ epsilon_EE_h[ iterator.second[iRtmp] ]->Fill( useMassInsteadOfEpsilon_? pi0P4_mass : eps_k, w ); }
		  //			    }
		  //			  }
		  //			}
		  //		    }
		  //		  }
		}
	    }
	  }  
	} // end filling histograms with mass

	if (isDebug_) cout << "[DEBUG] End of Cluster Loop" << endl;


    } // loop over clusters (g2)
  } // loop over clusters to make pi0 

  if(MakeNtuple4optimization_){
    if (isDebug_) cout << "[DEBUG] Filling Tree" << endl; 
    //for(unsigned int i=0; i<NL1SEED; i++) Op_L1Seed[i] = L1BitCollection_[i];
    //for(unsigned int i=0; i<NL1SEED; i++) Op_L1Seed[i] = l1flag[i];
    Op_NPi0 = nPi0; 
    if(nPi0>0) Tree_Optim->Fill();
  }
 
}


///////======================================


void FillEpsilonPlot::computeEoverEtrue(std::vector< CaloCluster > & clusters, std::vector<TLorentzVector*>& clusters_matchedGenPhoton, int subDetId ) 
{

  // this method is meant to be used with MC using MC truth
  // the cluster vector is such that two consecutive clusters belong to the same gen pi0
  // therefore, we loop on the vector using a step of 2, and assign g2 = g1+1

  if(subDetId!=EcalBarrel && subDetId != EcalEndcap) 
    throw cms::Exception("FillEpsilonPlot::computeEoverEtrue") << "Subdetector Id not recognized\n";

  if (isDebug_) cout << "[DEBUG] Beginning cluster loop for E/Etrue ..."<< endl;

  // loop over clusters to make Pi0
  size_t i = 0;

  for(std::vector<CaloCluster>::const_iterator g1  = clusters.begin(); g1 != clusters.end(); g1+=2, i+=2) 
  {

    std::vector<CaloCluster>::const_iterator g2 = g1 + 1;
    size_t j= i + 1;

    if (isDebug_) cout << "\n[DEBUG] New Pair of Clusters"<< endl;

    if( subDetId==EcalBarrel ) {EventFlow_EB->Fill(4.); if (isDebug_) EventFlow_EB_debug->Fill(0.);}
    else                       {EventFlow_EE->Fill(4.); if (isDebug_) EventFlow_EE_debug->Fill(0.);}

    // the following correction parameters are actually useless for MC. They would be a correction coming from a regression when using data
    // but this method is meant to be used with MC exactly to compute the containment corrections
    // for the moment I keep them because they might be useful for something later
    float Corr1 = 1., Corr2 = 1.;

    // g1 and g2 are ordered with the energy of the seed, but their respective clusters don't necessarily follow the same order
    // also, their pTs are not necessarily ordered 
    // Defining few variables to save photon quantities that are used more than once, to avoid recomputing them every time
    Double_t g1eta = g1->eta();
    Double_t g2eta = g2->eta();
    Double_t g1phi = g1->phi();
    Double_t g2phi = g2->phi();
    Double_t g1pt = g1->energy()/cosh(g1eta);
    Double_t g2pt = g2->energy()/cosh(g2eta);
    // following two object store the two photons ordered by pt
    TLorentzVector G_Sort_1, G_Sort_2, GSort1plus2;

    if( g1pt > g2pt ){
      G_Sort_1.SetPtEtaPhiE( g1pt, g1eta, g1phi, g1->energy() );
      G_Sort_2.SetPtEtaPhiE( g2pt, g2eta, g2phi, g2->energy() );
    }
    else{
      G_Sort_1.SetPtEtaPhiE( g2pt, g2eta, g2phi, g2->energy() );
      G_Sort_2.SetPtEtaPhiE( g1pt, g1eta, g1phi, g1->energy() );
    }

    GSort1plus2 = G_Sort_1 + G_Sort_2;


#if !defined(NEW_CONTCORR) && defined(MVA_REGRESSIO) || defined(REGRESS_AND_PARAM_CONTCORR)
    if( subDetId==EcalBarrel && (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ){

      // cout << "################################" << endl;
      // cout << "### We are in the barrel! ###" << endl;
      // cout << "################################" << endl;
	  
      // following variable should be equivalent to transverse energy of the photon pair (for massless object it is equal to Pt() )
      // this will store G.E()/cosh(G.Eta()), in order to compute it only once
      //	  double GSort1plus2_EoverCoshEta = GSort1plus2.E()/cosh(GSort1plus2.Eta());  // currently not used here for EB

      int ind1 = i, ind2 = j;
      EBDetId  id_1(g1->seed()); int iEta1 = id_1.ieta(); int iPhi1 = id_1.iphi();
      EBDetId  id_2(g2->seed()); int iEta2 = id_2.ieta(); int iPhi2 = id_2.iphi();
#ifdef MVA_REGRESSIO_Tree
      int iSMod_1 = id_1.ism(); int iSMod_2 = id_2.ism();
#endif

      bool Inverted=false;

      if( g1pt < g2pt ){
	iEta1=id_2.ieta(); iEta2 = id_1.ieta();
	iPhi1=id_2.iphi(); iPhi2 = id_1.iphi();
#ifdef MVA_REGRESSIO_Tree
	iSMod_1=id_2.ism(); iSMod_2=id_1.ism();
#endif
	ind1=j; ind2=i;
	Inverted=true;
      }

      float Correct1(1.), Correct2(1.);
      if(Are_pi0_){
	float value_pi01[14];
	//input list for regression in 2017 EB:
	//enG_rec
	//Nxtal
	//S4S9
	//S2S9
	//iEta
	//iPhi
	//SM_dist: 
	//M_dist:
	float new_value_pi01[8]; 

	value_pi01[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	value_pi01[1] = ( G_Sort_1.Pt() );
	value_pi01[2] = ( Ncristal_EB_used[ind1] );
	value_pi01[3] = ( Ncristal_EB_used[ind2] );
	value_pi01[4] = ( vs4s9[ind1] );
	value_pi01[5] = ( vs1s9[ind1] );
	value_pi01[6] = ( vs2s9[ind1] );
	value_pi01[7] = ( iEta1 );
	value_pi01[8] = ( iPhi1 );
	value_pi01[9] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	value_pi01[10] = ( iEta1%5 );
	value_pi01[11] = ( iPhi1%2 );
	value_pi01[12] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
	value_pi01[13] = ( iPhi1%20 );
	    
	new_value_pi01[0] = ( g1->energy() );
	new_value_pi01[1] = ( Ncristal_EB_used[i] );
	new_value_pi01[2] = ( vs4s9[i] );
	new_value_pi01[3] = ( vs2s9[i] );
	new_value_pi01[4] = ( id_1.ieta() );
	new_value_pi01[5] = ( id_1.iphi() );
	float temp_SM_dist_1 = ((id_1.iphi()-1)%20<10)*((id_1.iphi()-1)%20) + (((id_1.iphi()-1)%20)>=10)*(19-(id_1.iphi()-1)%20);
	float temp_M_dist_1  = (abs(id_1.ieta())<=25)*(((abs(id_1.ieta())-1)%25<12)*((abs(id_1.ieta())-1)%25) + (((abs(id_1.ieta())-1)%25)>=12)*(24-(abs(id_1.ieta())-1)%25))
	  +(abs(id_1.ieta())>25) * (((abs(id_1.ieta())-26)%20<10)*((abs(id_1.ieta())-26)%20) + (((abs(id_1.ieta())-26)%20)>=10)*(19-(abs(id_1.ieta())-26)%20));
	new_value_pi01[6] = ( temp_SM_dist_1 );
	new_value_pi01[7] = ( temp_M_dist_1 );

	//if( fabs((G_Sort_1+G_Sort_2).Eta())>1 ) value_pi01[14] = true ;
	//else                                    value_pi01[14] = false ;
	if(useMVAContainmentCorrections_)
	  {
	    if(new_pi0ContainmentCorrections_)
	      {
		float Correct1_tmp = forestD_EB_1->GetResponse(new_value_pi01);
		Correct1 = meanoffset + meanscale*TMath::Sin(Correct1_tmp);
		// cout<<"DEBUG in FillEpsilonPlot.cc... computeEpsilon... new regression Correct1 = "<<Correct1<<endl;
	      }
	    else
	      {
		Correct1 = forest_EB_1->GetResponse(value_pi01);
		// cout<<"DEBUG in FillEpsilonPlot.cc... computeEpsilon... old regression Correct1 = "<<Correct1<<endl;
	      }
	  }

	float value_pi02[14];//#
	float new_value_pi02[8];

	value_pi02[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	value_pi02[1] = ( G_Sort_2.Pt() );
	value_pi02[2] = ( Ncristal_EB_used[ind1] );
	value_pi02[3] = ( Ncristal_EB_used[ind2] );
	value_pi02[4] = ( vs4s9[ind2] );
	value_pi02[5] = ( vs1s9[ind2] );
	value_pi02[6] = ( vs2s9[ind2] );
	value_pi02[7] = ( iEta2 );
	value_pi02[8] = ( iPhi2 );
	value_pi02[9] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	value_pi02[10] = ( iEta2%5 );
	value_pi02[11] = ( iPhi2%2 );
	value_pi02[12] = ( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
	value_pi02[13] = ( iPhi2%20 );
	//if( fabs((G_Sort_1+G_Sort_2).Eta())>1 ) value_pi02[14] = true ;
	//else                                    value_pi02[14] = false ;
	   
	new_value_pi02[0] = ( g2->energy() );
	new_value_pi02[1] = ( Ncristal_EB_used[j] );
	new_value_pi02[2] = ( vs4s9[j] );
	new_value_pi02[3] = ( vs2s9[j] );
	new_value_pi02[4] = ( id_2.ieta() );
	new_value_pi02[5] = ( id_2.iphi() );
	float temp_SM_dist_2 = ((id_2.iphi()-1)%20<10)*((id_2.iphi()-1)%20) + (((id_2.iphi()-1)%20)>=10)*(19-(id_2.iphi()-1)%20);
	float temp_M_dist_2  = (abs(id_2.ieta())<=25)*(((abs(id_2.ieta())-1)%25<12)*((abs(id_2.ieta())-1)%25) + (((abs(id_2.ieta())-1)%25)>=12)*(24-(abs(id_2.ieta())-1)%25))
	  +(abs(id_2.ieta())>25) * (((abs(id_2.ieta())-26)%20<10)*((abs(id_2.ieta())-26)%20) + (((abs(id_2.ieta())-26)%20)>=10)*(19-(abs(id_2.ieta())-26)%20));
	new_value_pi02[6] = ( temp_SM_dist_2 );
	new_value_pi02[7] = ( temp_M_dist_2 );


	if(useMVAContainmentCorrections_)
	  {
            if(new_pi0ContainmentCorrections_)
	      {
                float Correct2_tmp = forestD_EB_2->GetResponse(new_value_pi02);
                Correct2 = meanoffset + meanscale*TMath::Sin(Correct2_tmp);
	      }
            else
	      {
                Correct2 = forest_EB_2->GetResponse(value_pi02);
	      }
	  }
      }
      else{
	float value_pi01[10];
	value_pi01[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	value_pi01[1] = ( G_Sort_1.Pt() );
	value_pi01[2] = ( Ncristal_EB_used[ind1] );
	value_pi01[3] = ( iEta1 );
	value_pi01[4] = ( iPhi1 );
	value_pi01[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	value_pi01[6] = ( iEta1%5 );
	value_pi01[7] = ( iPhi1%2 );
	value_pi01[8] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
	value_pi01[9] = ( iPhi1%20 );
	Correct1 = forest_EB_1->GetResponse(value_pi01);
	float value_pi02[10];
	value_pi02[0] = ( G_Sort_1.Energy()/G_Sort_2.Energy() );
	value_pi02[1] = ( G_Sort_2.Pt() );
	value_pi02[2] = ( Ncristal_EB_used[ind2] );
	value_pi02[3] = ( iEta2 );
	value_pi02[4] = ( iPhi2 );
	value_pi02[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	value_pi02[6] = ( iEta2%5 );
	value_pi02[7] = ( iPhi2%2 );
	value_pi02[8] = ( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
	value_pi02[9] = ( iPhi2%20 );
	Correct2 = forest_EB_2->GetResponse(value_pi02);
      }
	
      // new regression trained in 2017: photon order is based on seed energy (energy in overlapped crystal taken by the first photon)
      // old regression trained in 2012: photon order is based on photon Pt
      if( (!Inverted) || new_pi0ContainmentCorrections_) { Corr1 = Correct1; Corr2 = Correct2; }
      else           { Corr1 = Correct2; Corr2 = Correct1; }

      //WARNIGN no CC for now! Put back in CKM 20/10/2015
      //	  Corr1 = 1.; Corr2 = 1.; 
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
      //In case ES give same posizion for different clusters
      Correction1_mva = Correct1; Correction2_mva = Correct2;
      iEta1_mva = iEta1; iEta2_mva = iEta2; iPhi1_mva = iPhi1; iPhi2_mva = iPhi2; Pt1_mva = G_Sort_1.Pt(); Pt2_mva = G_Sort_2.Pt();
      iSM1_mva = iSMod_1; iSM2_mva = iSMod_2;

      // Define just the original vector, the corrected vector has same direction and magnitude scaled by correction
      // actually we just need G_sort_1 and G_sort_2, no need to define another object

      // TLorentzVector mvag1P4; mvag1P4.SetPtEtaPhiE( Correct1*G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), Correct1*G_Sort_1.E() );
      // TLorentzVector mvag2P4; mvag2P4.SetPtEtaPhiE( Correct2*G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), Correct2*G_Sort_2.E() );
      // TLorentzVector mvaOrg1P4; mvaOrg1P4.SetPtEtaPhiE( G_Sort_1.Pt(), G_Sort_1.Eta(), G_Sort_1.Phi(), G_Sort_1.E() );
      // TLorentzVector mvaOrg2P4; mvaOrg2P4.SetPtEtaPhiE( G_Sort_2.Pt(), G_Sort_2.Eta(), G_Sort_2.Phi(), G_Sort_2.E() );

      MassOr_mva = GSort1plus2.M();
      pi0Eta     = GSort1plus2.Eta();
      // get corrected mass from sum of corrected photons
      Mass_mva = (Correct1 * G_Sort_1 + Correct2 * G_Sort_2).M();
#endif
    }
#endif

    // end of part using regression for EB (forget EE, we don't use any correction there at the moment
    ////////////////////////////////////
    /////////////////////////////////////
	  	
    // uncorrected versions of photons
    // math::PtEtaPhiMLorentzVector g1P4_nocor( g1pt, g1eta, g1phi, 0. );
    // math::PtEtaPhiMLorentzVector g2P4_nocor( g2pt, g2eta, g2phi, 0. );
    //math::PtEtaPhiMLorentzVector pi0P4_nocor = g1P4_nocor + g2P4_nocor;
    // here the order is not important (while for the correction we apply corr1 to first cluster in list (not necessarily the leading in pt))
    double pi0P4_nocor_pt = GSort1plus2.Pt();
    double pi0P4_nocor_mass = GSort1plus2.M();

    //corrected version; note that Corr1 and Corr2 refers to first and second photon as selected looping on CaloCluster 
    // this means g1 is not necessarily the leading photon
    TLorentzVector pi0P4;
    TLorentzVector g1_contCorr_tlv;
    TLorentzVector g2_contCorr_tlv;

    if (useContainmentCorrectionsFromEoverEtrue_) {

      // if (subDetId==EcalBarrel)                                                                                                                                       
      //   std::cout << "Barrel: npi0 = " << nPi0 << "   pi0,g1,g2 eta = " << GSort1plus2.Eta(),g1.eta(),g2.eta() << std::endl;                                          
      // else                                                                                                                                                            
      //   std::cout << "   Endcap: npi0 = " << nPi0 << "   pi0,g1,g2 eta = " << GSort1plus2.Eta(),g1.eta(),g2.eta() << std::endl;                                       

      CaloCluster g1_contCorr(getClusterAfterContainmentCorrections(g1,false,subDetId==EcalBarrel));
      CaloCluster g2_contCorr(getClusterAfterContainmentCorrections(g2,true,subDetId==EcalBarrel));

      g1_contCorr_tlv.SetPtEtaPhiE(g1_contCorr.energy()/cosh(g1_contCorr.eta()),g1_contCorr.eta(),g1_contCorr.phi(),g1_contCorr.energy());
      g2_contCorr_tlv.SetPtEtaPhiE(g2_contCorr.energy()/cosh(g2_contCorr.eta()),g2_contCorr.eta(),g2_contCorr.phi(),g2_contCorr.energy());                            

    } else {

      if (g1pt > g2pt) {
	g1_contCorr_tlv = Corr1 * G_Sort_1;
	g2_contCorr_tlv = Corr2 * G_Sort_2;
      } else {
	// when g1pt < g2pt, G_Sort_1 is made with g2, and Corr2 must be applied to it                      
	g1_contCorr_tlv = Corr1 * G_Sort_2;
	g2_contCorr_tlv = Corr2 * G_Sort_1; 
      }
    }

    pi0P4 = g1_contCorr_tlv + g2_contCorr_tlv;                                                                                                                      

    // eta, pt, phi of corrected photons are used many times. Since their computation is tipically time consuming, store them in doubles for later usage
    double pi0P4_pt = pi0P4.Pt();
    double pi0P4_eta = pi0P4.Eta();
    double pi0P4_phi = pi0P4.Phi();
    double pi0P4_mass = pi0P4.M();
    // note that photon eta and phi are not modified by correction (only pT) since Corr * vector modifies the cartesian coordinates of the vector (pT and pZ)
    // double g1P4_eta = g1eta; 
    // double g1P4_phi = g1phi;
    // double g2P4_eta = g2eta;
    // double g2P4_phi = g2phi;

    //cout << "pi0P4_nocor_mass, pi0P4_mass " << pi0P4_nocor_mass << "  " << pi0P4_mass << endl;


    // g1 and g2 are already ordered, the same for the vectors Ncristal_EE_used ad Ncristal_EB_used (index i corresponds to g1, j to g2)
    int Nxtal_EnergGamma = 0;
    int Nxtal_EnergGamma2 = 0;
    if(subDetId==EcalEndcap){
      Nxtal_EnergGamma = Ncristal_EE_used[i]; 
      Nxtal_EnergGamma2 = Ncristal_EE_used[j];
    }
    else{
      Nxtal_EnergGamma = Ncristal_EB_used[i]; 
      Nxtal_EnergGamma2 = Ncristal_EB_used[j]; 
    }


    // I'll have to decide if I want to apply the reco pi0  selection to fill the E/Etrue
    // However, in case we are matching the clusters to gen photons, it is implicit that the clusters belong to a pi0 (even though they might fail the reco selection)
    // The reco selection is already applied at the singlecluster level in fillEBcluster()
    //
    ///////////////////
    // BEGIN SELECTION
    ///////////////////

    bool applySelectionForEoverEtrue = true;

    if (applySelectionForEoverEtrue) {

      //In case ES give same posizion for different clusters
      if( pi0P4_nocor_mass<0.03 && pi0P4_mass < 0.03 ) continue;
	

      if (isDebug_) cout << "[DEBUG] Apply kinematic selection cuts" << endl;

      if( g1eta == g2eta && g1phi == g2phi ) continue; // this should already be impossible ...

      // pi0/eta pT cut
      if (subDetId == EcalBarrel) {

	if (fabs(pi0P4_eta)<1.0)       { if( pi0P4_nocor_pt < pi0PtCut_low_[subDetId]) continue; }
	else if (fabs(pi0P4_eta)<1.479) { if( pi0P4_nocor_pt < pi0PtCut_high_[subDetId]) continue; }
	if (isDebug_) EventFlow_EB_debug->Fill(1.);

      } else {
	  
	if (fabs(pi0P4_eta)<1.8 )     { if( pi0P4_nocor_pt < pi0PtCut_low_[subDetId]) continue; }	  
	else                          { if( pi0P4_nocor_pt < pi0PtCut_high_[subDetId]) continue; }
	if (isDebug_) EventFlow_EE_debug->Fill(1.);

      }

      float nextClu = 999., Drtmp = 999.;
      for(size_t ind=0; ind<clusters.size(); ++ind){
	const CaloCluster* Gtmp = &(clusters[ind]);
	double deltaR1 = GetDeltaR(Gtmp->eta(),g1eta,Gtmp->phi(),g1phi);
	double deltaR2 = GetDeltaR(Gtmp->eta(),g2eta,Gtmp->phi(),g2phi);
	if( ind!=i && ind!=j && (deltaR1<Drtmp || deltaR2<Drtmp ) ){
	  nextClu = min(deltaR1,deltaR2);
	  Drtmp = nextClu;
	}
      }

      // pi0/eta isolation cut (distance to other clusters)
      if (subDetId == EcalBarrel) {

	if (fabs(pi0P4_eta)<1.0)       { if( nextClu<pi0IsoCut_low_[subDetId] ) continue; }
	else if (fabs(pi0P4_eta)<1.479) { if( nextClu<pi0IsoCut_high_[subDetId] ) continue; }
	if (isDebug_) EventFlow_EB_debug->Fill(2.);

      } else {
	  
	if (fabs(pi0P4_eta)<1.8 )     { if( nextClu<pi0IsoCut_low_[subDetId] ) continue; }	  
	else                          { if( nextClu<pi0IsoCut_high_[subDetId] ) continue; }
	if (isDebug_) EventFlow_EE_debug->Fill(2.);

      }

      // Implementation of HLT Filter Isolation - Eta Band Isolation 
      // implemented in HLT: CMSSW_7_1_0/src/HLTrigger/special/src/HLTEcalResonanceFilter.cc
      // see Yong Yang's  Thesis: http://thesis.library.caltech.edu/7345/

      if (isDebug_) cout << "[DEBUG] Running HLT Isolation" << endl;

      float hlt_iso = 0;
      for(size_t ind=0; ind < clusters.size(); ++ind){
	if( clusters[ind].seed() == clusters[i].seed() || clusters[ind].seed() == clusters[j].seed()) continue;
	const CaloCluster* Gtmp = &(clusters[ind]);
	TLorentzVector GtmpP4;  
	GtmpP4.SetPtEtaPhiE(Gtmp->energy()/cosh(Gtmp->eta()), Gtmp->eta(), Gtmp->phi(), Gtmp->energy());
	if (GtmpP4.Pt() < 0.5) continue;  // FIXME: based on the stream, it should represent "ptMinForIsolation*"
	// delta R from the pi0 candidates
	double deltaR0 = GetDeltaR(Gtmp->eta(), pi0P4_eta, Gtmp->phi(), pi0P4_phi);
	if (deltaR0  > ((Are_pi0_) ? 0.2:0.3)) continue;
	// cluster must be inside of an eta strip 
	double deta = fabs(Gtmp->eta() - pi0P4_eta); 
	if (deta > ((Are_pi0_) ? 0.05:0.1)) continue;
	hlt_iso += GtmpP4.Pt();
      }
      // the cut is taken relative to the pi0 pt
      hlt_iso /= pi0P4_nocor_pt;
      //category break down of cuts
      // pi0/eta isolation cut
      if (subDetId == EcalBarrel) {

	if (fabs(pi0P4_eta)<1.0)       { if( hlt_iso > pi0HLTIsoCut_low_[subDetId]  && CutOnHLTIso_ ) continue; }
	else if (fabs(pi0P4_eta)<1.479) { if( hlt_iso > pi0HLTIsoCut_high_[subDetId] && CutOnHLTIso_ ) continue; }
	if (isDebug_) EventFlow_EB_debug->Fill(3.);

      } else {
	  
	if (fabs(pi0P4_eta)<1.8 )     { if( hlt_iso > pi0HLTIsoCut_low_[subDetId]  && CutOnHLTIso_ ) continue; }	  
	else                          { if( hlt_iso > pi0HLTIsoCut_high_[subDetId] && CutOnHLTIso_ ) continue; }
	if (isDebug_) EventFlow_EE_debug->Fill(3.);

      }
      //////////////////////////////////////////////////////////////////////////////////////////////////

      if (isDebug_) cout << "[DEBUG] N Cristal Cuts" << endl;

      if (subDetId == EcalBarrel) {

	if( fabs(pi0P4_eta)<1.0 ) { 
	  if( Nxtal_EnergGamma < nXtal_1_cut_low_[subDetId] ) continue; 
	  if( Nxtal_EnergGamma2 < nXtal_2_cut_low_[subDetId] ) continue;
	} else if( fabs(pi0P4_eta)<1.479 )  { 
	  if( Nxtal_EnergGamma < nXtal_1_cut_high_[subDetId] ) continue; 
	  if( Nxtal_EnergGamma2 < nXtal_2_cut_high_[subDetId] ) continue;
	} 

	pi0MassVsIetaEB->Fill( fabs(pi0P4_eta)/0.0174, pi0P4_mass);
	pi0MassVsETEB->Fill(pi0P4_pt, pi0P4_mass);
	if (isDebug_) EventFlow_EB_debug->Fill(4.);
	EventFlow_EB->Fill(5.);

      } else {

	if( fabs(pi0P4_eta)<1.8 ) { 
	  if( Nxtal_EnergGamma < nXtal_1_cut_low_[subDetId] ) continue; 
	  if( Nxtal_EnergGamma2 < nXtal_2_cut_low_[subDetId] ) continue;
	} else {
	  if( Nxtal_EnergGamma < nXtal_1_cut_high_[subDetId] ) continue; 
	  if( Nxtal_EnergGamma2 < nXtal_2_cut_high_[subDetId] ) continue;
	}
	if (isDebug_) EventFlow_EE_debug->Fill(4.);
	EventFlow_EE->Fill(5.);

      }

    } // end of if (applySelectionForEoverEtrue)

    ///////////////////
    // END SELECTION
    ///////////////////
    ////////////////////////////////////////////

    if (fillKinematicVariables_) {

      if (fabs(pi0P4_eta)<1.0) whichRegionEcalStreamPi0 = 0;
      else if (fabs(pi0P4_eta)<1.479) whichRegionEcalStreamPi0 = 1;
      else if (fabs(pi0P4_eta)<1.8) whichRegionEcalStreamPi0 = 2;
      else whichRegionEcalStreamPi0 = 3;

      pi0pt_afterCuts[whichRegionEcalStreamPi0]->Fill(pi0P4_nocor_pt);
      g1pt_afterCuts[whichRegionEcalStreamPi0]->Fill(g1pt);
      g2pt_afterCuts[whichRegionEcalStreamPi0]->Fill(g2pt);
      g1Nxtal_afterCuts[whichRegionEcalStreamPi0]->Fill(Nxtal_EnergGamma);
      g2Nxtal_afterCuts[whichRegionEcalStreamPi0]->Fill(Nxtal_EnergGamma2);
      pi0PhotonsNoverlappingXtals_afterCuts[whichRegionEcalStreamPi0]->Fill(getNumberOverlappingCrystals(g1,g2,subDetId==EcalBarrel));
      if (isMC_) {
	pi0MassVsPU[whichRegionEcalStreamPi0]->Fill(pi0P4_nocor_mass,nPUobs_BX0_);
      }
    }

    if (!MakeNtuple4optimization_) {

      if (isDebug_) cout << "[DEBUG] computing region weights" << endl; 

      // compute region weights
      RegionWeightVector w1 = regionalCalibration_->getWeights( &(*g1), subDetId ); // region weights W_j^k for clu1
      RegionWeightVector w2 = regionalCalibration_g2_->getWeights( &(*g2), subDetId ); // region weights W_j^k for clu2

      // do not append w2 to w1, keep them separate to disentangle the 2 photons
      // w1.insert( w1.end(), w2.begin(), w2.end() );

      double EoverEtrue_g1 = 0.0;
      double EoverEtrue_g2 = 0.0;

      // if (useContainmentCorrectionsFromEoverEtrue_) {
      // 	EoverEtrue_g1 = g1_contCorr.energy()/clusters_matchedGenPhoton[i]->Energy();
      // 	EoverEtrue_g2 = g2_contCorr.energy()/clusters_matchedGenPhoton[j]->Energy();
      // } else {
      // 	// if no corrections are used, Corr1 and Corr2 are 1.0
      // 	EoverEtrue_g1 = Corr1 * g1->energy()/clusters_matchedGenPhoton[i]->Energy();
      // 	EoverEtrue_g2 = Corr2 * g2->energy()/clusters_matchedGenPhoton[j]->Energy();
      // }

      EoverEtrue_g1 = g1_contCorr_tlv.Energy()/clusters_matchedGenPhoton[i]->Energy();
      EoverEtrue_g2 = g1_contCorr_tlv.Energy()/clusters_matchedGenPhoton[j]->Energy();


      // compute quantities needed for <eps>_j in each region j
      if (subDetId==EcalBarrel) {
	allEoverEtrue_g1_EBnw->Fill( EoverEtrue_g1 );
	allEoverEtrue_g2_EBnw->Fill( EoverEtrue_g2 );
      } else {
	allEoverEtrue_g1_EEnw->Fill( EoverEtrue_g1 );
	allEoverEtrue_g2_EEnw->Fill( EoverEtrue_g2 );
      }

      for (RegionWeightVector::const_iterator it = w1.begin(); it != w1.end(); ++it) {

	const uint32_t& iR = (*it).iRegion;
	const float& w = (*it).value;

	if (subDetId==EcalBarrel) {

	  // if using the reco selection, select only a window in pi0 mass, otherwise just skip the selection
	  if ( not applySelectionForEoverEtrue || (pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.23:0.7)) ) {

	    EoverEtrue_g1_EB_h[iR]->Fill( EoverEtrue_g1, w );
	    allEoverEtrue_g1_EB->Fill( EoverEtrue_g1, w );
	    int iEta = List_IR_EtaPhi.find(iR)->second[0]; 
	    int iPhi = List_IR_EtaPhi.find(iR)->second[1]; 
	    //int iSM = List_IR_EtaPhi.find(iR)->second[2];
	    entries_EB->Fill( iEta, iPhi, w );
	  
	  }

	} else {

	  if (  not applySelectionForEoverEtrue || (pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.28:0.75)) ) {

	    EoverEtrue_g1_EE_h[iR]->Fill( EoverEtrue_g1, w );
	    allEoverEtrue_g1_EE->Fill( EoverEtrue_g1, w );
	    int iX = List_IR_XYZ.find(iR)->second[0]; 
	    int iY = List_IR_XYZ.find(iR)->second[1]; 
	    int iZ = List_IR_XYZ.find(iR)->second[2]; 
	    //int Quad = List_IR_XYZ.find(iR)->second[3];
	    if ( iZ==-1 ) entries_EEm->Fill( iX, iY, w );
	    else entries_EEp->Fill( iX, iY, w );
	 
	  }  // closes condition on mass boundary

	}   // if subDetId == Endcap (closes else)

      }   // for (RegionWeightVector::const_iterator it = w1.begin(); it != w1.end(); ++it) {  


      for (RegionWeightVector::const_iterator it = w2.begin(); it != w2.end(); ++it) {

	const uint32_t& iR = (*it).iRegion;
	const float& w = (*it).value;

	if (subDetId==EcalBarrel) {

	  // if using the reco selection, select only a window in pi0 mass, otherwise just skip the selection
	  if ( not applySelectionForEoverEtrue || (pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.23:0.7)) ) {

	    EoverEtrue_g2_EB_h[iR]->Fill( EoverEtrue_g2, w );
	    allEoverEtrue_g2_EB->Fill( EoverEtrue_g2, w );
	    int iEta = List_IR_EtaPhi.find(iR)->second[0]; 
	    int iPhi = List_IR_EtaPhi.find(iR)->second[1]; 
	    //int iSM = List_IR_EtaPhi.find(iR)->second[2];
	    entries_EB->Fill( iEta, iPhi, w );
	  
	  }

	} else {

	  if (  not applySelectionForEoverEtrue || (pi0P4_mass>((Are_pi0_)?0.03:0.35) && pi0P4_mass<((Are_pi0_)?0.28:0.75)) ) {

	    EoverEtrue_g2_EE_h[iR]->Fill( EoverEtrue_g2, w );
	    allEoverEtrue_g2_EE->Fill( EoverEtrue_g2, w );
	    int iX = List_IR_XYZ.find(iR)->second[0]; 
	    int iY = List_IR_XYZ.find(iR)->second[1]; 
	    int iZ = List_IR_XYZ.find(iR)->second[2]; 
	    //int Quad = List_IR_XYZ.find(iR)->second[3];
	    if ( iZ==-1 ) entries_EEm->Fill( iX, iY, w );
	    else entries_EEp->Fill( iX, iY, w );
	 
	  }  // closes condition on mass boundary

	}   // if subDetId == Endcap (closes else)

      }   // for (RegionWeightVector::const_iterator it = w2.begin(); it != w2.end(); ++it) {  

    } // end filling histograms with mass

    if (isDebug_) cout << "[DEBUG] End of Cluster Loop for E/Etrue" << endl;


  } // loop over clusters to make pi0 
 
}


///////======================================


// ------------ method called once each job just before starting event loop  ------------
  void 
FillEpsilonPlot::beginJob()
{

  if (isDebug_) cout << "[DEBUG] beginJob" << endl;

  /// testing the EE eta ring
  TH2F eep("eep","EE+",102,0.5,101.5,102,-0.5,101.5);
  TH2F eem("eem","EE-",102,0.5,101.5,102,-0.5,101.5);

  for(int etaring=0; etaring < EndcapTools::N_RING_ENDCAP; ++etaring)
  {

    float fillValue = (etaring%2)==0 ? 1. : 2.;
    std::vector<DetId> allDetIds = EcalCalibType::EtaRing::allDetIdsInEERegion(etaring);
    for(int ixtal=0; ixtal<int(allDetIds.size()); ixtal++)
    {
	EEDetId eeid(allDetIds.at(ixtal));
	if(eeid.zside()==-1)
	  eem.SetBinContent(eeid.ix(),eeid.iy(),fillValue);
	else
	  eep.SetBinContent(eeid.ix(),eeid.iy(),fillValue);

    }
  }
  outfile_->cd();
  eep.Write();
  eem.Write();

  ifstream file;
  file.open( edm::FileInPath ( Endc_x_y_.c_str() ).fullPath().c_str(), ifstream::in);
  VectRing.clear();
  while ( !file.eof() ) {
    string Line;
    getline( file, Line);
    string value;
    std::stringstream MyLine(Line);

    char * cstr, *p;
    cstr = new char [Line.size()+1];
    strcpy (cstr, Line.c_str());
    p=strtok (cstr," ");
    int i(0);
    while (p!=NULL){
	if(i==0)  GiveRing.iX = atoi(p);
	if(i==1)  GiveRing.iY = atoi(p);
	if(i==2)  GiveRing.sign = atoi(p);
	if(i==3){
	  GiveRing.Ring = atoi(p);
	  VectRing.push_back(GiveRing);
	}
	p=strtok(NULL," ");
	i++;
    }
    delete[] cstr;
  }
  //Initialize Map iR vs Eta
  if( (SMCalibEB_ && EtaRingCalibEB_) || (SMCalibEE_ && EtaRingCalibEE_) ) cout<<"WARNING: Intercalibrating with EtaRing and SM!!!"<<endl; 
  std::vector<int> InitV; InitV.clear();
  for(Long64_t i=-85; i<86; i++) ListEtaFix_xtalEB[i]   = InitV;
  for(Long64_t i=0; i<37; i++)   ListSMFix_xtalEB[i]    = InitV;
  for(Long64_t i=0; i<40; i++)   ListEtaFix_xtalEEm[i]  = InitV;
  for(Long64_t i=0; i<40; i++)   ListEtaFix_xtalEEp[i]  = InitV;
  for(Long64_t i=0; i<10; i++)   ListQuadFix_xtalEEm[i] = InitV;
  for(Long64_t i=0; i<10; i++)   ListQuadFix_xtalEEp[i] = InitV;
  //Open File where to take iR vs Eta
  TFile *CalibMapEtaRingF = TFile::Open( edm::FileInPath( CalibMapEtaRing_.c_str() ).fullPath().c_str() );
  TTree *calibMap_EB      = (TTree*) CalibMapEtaRingF->Get("calibEB");
  TTree *calibMap_EE      = (TTree*) CalibMapEtaRingF->Get("calibEE");
  Int_t hashedIndexEB_, hashedIndexEE_, ieta_, iphi_, iSM_, ix_, iy_, zside_, iquadrant_;
  calibMap_EB->SetBranchAddress( "hashedIndex_", &hashedIndexEB_);
  calibMap_EB->SetBranchAddress( "ieta_", &ieta_);
  calibMap_EB->SetBranchAddress( "iphi_", &iphi_);
  calibMap_EB->SetBranchAddress( "iSM_", &iSM_);
  calibMap_EE->SetBranchAddress( "hashedIndex_", &hashedIndexEE_);
  calibMap_EE->SetBranchAddress( "ix_", &ix_);
  calibMap_EE->SetBranchAddress( "iy_", &iy_);
  calibMap_EE->SetBranchAddress( "zside_", &zside_);
  calibMap_EE->SetBranchAddress( "iquadrant_", &iquadrant_);
  //Loop to Fill the map in EB
  Long64_t nentries = calibMap_EB->GetEntriesFast();
  for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
    calibMap_EB->GetEntry(iEntry);
    ListEtaFix_xtalEB[ieta_].push_back( hashedIndexEB_ );
    ListSMFix_xtalEB[iSM_].push_back( hashedIndexEB_ );
    std::vector<int> EtaPhi; EtaPhi.clear(); EtaPhi.push_back( ieta_ ); EtaPhi.push_back( iphi_ ); EtaPhi.push_back( iSM_ );
    List_IR_EtaPhi[hashedIndexEB_] = EtaPhi;
  }
  //Loop to Fill the map in EE
  nentries = calibMap_EE->GetEntriesFast();
  for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
    calibMap_EE->GetEntry(iEntry);
    if(zside_<0.) ListEtaFix_xtalEEm[ GetRing( ix_,iy_, VectRing, false) ].push_back( hashedIndexEE_ );
    if(zside_>0.) ListEtaFix_xtalEEp[ GetRing( ix_,iy_, VectRing, false) ].push_back( hashedIndexEE_ );
    if(zside_<0.) ListQuadFix_xtalEEm[ iquadrant_ ].push_back( hashedIndexEE_ );
    if(zside_>0.) ListQuadFix_xtalEEp[ iquadrant_ ].push_back( hashedIndexEE_ );
    std::vector<int> iXYZ; iXYZ.clear(); iXYZ.push_back( ix_ ); iXYZ.push_back( iy_ ); iXYZ.push_back( zside_ ); iXYZ.push_back( iquadrant_ );
    List_IR_XYZ[hashedIndexEE_] = iXYZ;
  }
  //  //###########
  //  fstream  file_Ix;
  //  file_Ix.open( "/afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_5_3_6/src/CalibCode/submit/common/ix_iy_iz_EtaRing_Eta.txt", ios::out);
  //  for(int x=0; x<100;x++){
  //    for(int y=0; y<100;y++){
  //	int ring = GetRing( x, y, VectRing, false);
  //	if(ring!=-1){
  //	  EEDetId EE_id(x, y, 1, 0);
  //	  file_Ix << x << " "<< y << " " << ring << " " <<endl;
  //	}
  //    }
  //  }
  //  file_Ix.close();

#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
  TTree_JoshMva = new TTree("TTree_JoshMva","MVA corrections");
  TTree_JoshMva->Branch("Correction1_mva", &Correction1_mva, "Correction1_mva/F");
  TTree_JoshMva->Branch("Correction2_mva", &Correction2_mva, "Correction2_mva/F");
  TTree_JoshMva->Branch("iEta1_mva", &iEta1_mva, "iEta1_mva/I");
  TTree_JoshMva->Branch("iEta2_mva", &iEta2_mva, "iEta2_mva/I");
  TTree_JoshMva->Branch("iPhi1_mva", &iPhi1_mva, "iPhi1_mva/I");
  TTree_JoshMva->Branch("iPhi2_mva", &iPhi2_mva, "iPhi2_mva/I");
  TTree_JoshMva->Branch("iSM1_mva", &iSM1_mva, "iSM1_mva/I");
  TTree_JoshMva->Branch("iSM2_mva", &iSM2_mva, "iSM2_mva/I");
  TTree_JoshMva->Branch("Pt1_mva", &Pt1_mva, "Pt1_mva/F");
  TTree_JoshMva->Branch("Pt2_mva", &Pt2_mva, "Pt2_mva/F");
  TTree_JoshMva->Branch("Mass_mva", &Mass_mva, "Mass_mva/F");
  TTree_JoshMva->Branch("MassOr_mva", &MassOr_mva, "MassOr_mva/F");
  TTree_JoshMva->Branch("pi0Eta", &pi0Eta, "pi0Eta/F");
#endif
#ifdef MVA_REGRESSIO_EE
  TTree_JoshMva_EE = new TTree("TTree_JoshMva_EE","EE MVA corrections");
  TTree_JoshMva_EE->Branch("Correction1EE_mva", &Correction1EE_mva, "Correction1EE_mva/F");
  TTree_JoshMva_EE->Branch("Correction2EE_mva", &Correction2EE_mva, "Correction2EE_mva/F");
  TTree_JoshMva_EE->Branch("iX1_mva", &iX1_mva, "iX1_mva/I");
  TTree_JoshMva_EE->Branch("iX2_mva", &iX2_mva, "iX2_mva/I");
  TTree_JoshMva_EE->Branch("iY1_mva", &iY1_mva, "iY1_mva/I");
  TTree_JoshMva_EE->Branch("iY2_mva", &iY2_mva, "iY2_mva/I");
  TTree_JoshMva_EE->Branch("EtaRing1_mva", &EtaRing1_mva, "EtaRing1_mva/I");
  TTree_JoshMva_EE->Branch("EtaRing2_mva", &EtaRing2_mva, "EtaRing2_mva/I");
  TTree_JoshMva_EE->Branch("Pt1EE_mva", &Pt1EE_mva, "Pt1EE_mva/F");
  TTree_JoshMva_EE->Branch("Pt2EE_mva", &Pt2EE_mva, "Pt2EE_mva/F");
  TTree_JoshMva_EE->Branch("MassEE_mva", &MassEE_mva, "MassEE_mva/F");
  TTree_JoshMva_EE->Branch("MassEEOr_mva", &MassEEOr_mva, "MassEEOr_mva/F");
#endif
}



float 
FillEpsilonPlot::GetDeltaR(float eta1, float eta2, float phi1, float phi2){

  return sqrt( (eta1-eta2)*(eta1-eta2) 
	+ DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );

}


float 
FillEpsilonPlot::DeltaPhi(float phi1, float phi2){

  //float diff = fabs(phi2 - phi1);
  // while (diff >acos(-1)) diff -= 2*acos(-1);
  // while (diff <= -acos(-1)) diff += 2*acos(-1);

  float diff = phi2 - phi1;
  while (diff > TMath::Pi()) diff -= 2.*TMath::Pi();
  while (diff <= -TMath::Pi()) diff += 2.*TMath::Pi();

  return diff;

}

bool FillEpsilonPlot::GetHLTResults(const edm::Event& iEvent, std::string s){

  edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
  iEvent.getByToken(triggerResultsToken_, hltTriggerResultHandle);

  edm::TriggerNames HLTNames;

  HLTNames = iEvent.triggerNames(*hltTriggerResultHandle);
  std::string tempnames;
  int hltCount = hltTriggerResultHandle->size();
  TRegexp reg(TString( s.c_str()) );
  //std::cout << "hlCount = " << hltCount << std::endl;
  for (int i = 0 ; i != hltCount; ++i) {
    TString hltName_tstr(HLTNames.triggerName(i));
    //std::cout << "hltName_tstr is: " << hltName_tstr << " and reg is: " << s << std::endl;
    if ( hltName_tstr.Contains(reg) ){          // If reg contains * ir will say always True. So you ask for ->accept(i) to the first HLTName always.
      //std::cout << "hltName_tstr.Contains(reg) give: " << hltName_tstr << " --> " << hltTriggerResultHandle->accept(i) << std::endl;
      return hltTriggerResultHandle->accept(i); // False or True depending if it fired.
    }
  }
  return false;
} // HLT isValid

// not used anymore
// bool FillEpsilonPlot::getTriggerByName( std::string s ) {
//   std::map< std::string, int >::iterator currentTrigger;
//   currentTrigger = l1TrigNames_.find(s);
//   if(currentTrigger != l1TrigNames_.end())
//     return l1TrigBit_[currentTrigger->second];
//   else 
//     std::cout << "Trigger Name not found" << std::endl;
//   return false;
// }

bool FillEpsilonPlot::getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle< GlobalAlgBlkBxCollection > gtReadoutRecord;
  iEvent.getByToken( L1GTobjmapToken_, gtReadoutRecord);

  if (gtReadoutRecord.isValid()) {

    const GlobalAlgBlkBxCollection *l1results = gtReadoutRecord.product(); 

    ////////////////////
    ////////////////////
    // here we redo the association bit number <--> bit name
    // the reason is that this is not a constant 
    //e.g. during data taking in 2017 I noticed the number associated to a name changed, for instance SingleJet16 was 130 and then it became 131)
    edm::ESHandle<L1TUtmTriggerMenu> menu;
    iSetup.get<L1TUtmTriggerMenuRcd>().get(menu);
    // get the bit/name association         
    for (auto const & keyval: menu->getAlgorithmMap()) { 
      std::string const & trigName  = keyval.second.getName(); 
      unsigned int iTrigIndex = keyval.second.getIndex(); 
      algoBitToName[iTrigIndex] = TString( trigName );  
    } // end algo Map
    ////////////////////
    ////////////////////

    GlobalAlgBlk const &result = l1results->at(0, 0);

    for (unsigned int itrig = 0; itrig < result.maxPhysicsTriggers; ++itrig) {
      //      std::cerr << "bit: " << itrig << "\tresult: " << results.getAlgoDecisionFinal(itrig) << std::endl;

      // L1 decision below is 1 if seed fired, 0 if it didn't. 
      // seedIsInStream[itrig] = -1 if index is not valid, = 0 if index is valid but the seed is not used by the stream
      if (seedIsInStream[itrig] > 0) { 
    
	bool myflag = result.getAlgoDecisionFinal(itrig) ; 
	if (myflag ) { 
	  l1flag[itrig] = 1; 
	  triggerComposition->Fill(algoBitToName[itrig], l1flag[itrig]); 
	  if (EB_HLT) triggerComposition_EB->Fill(algoBitToName[itrig], l1flag[itrig]); 
	  if (EE_HLT) triggerComposition_EE->Fill(algoBitToName[itrig], l1flag[itrig]); 
	  // cout << " itrig = " << itrig << "    ";
	  // cout << " l1flag[itrig] = " << l1flag[itrig] << "    ";
	  // cout << " algoBitToName[itrig] = " << algoBitToName[itrig] << endl;	  
	} else {
	  l1flag[itrig] = 0 ; 
	}

      } 
 
    }

  }

  // edm::Handle< L1GlobalTriggerObjectMapRecord > gtReadoutRecord;
  // iEvent.getByToken( L1GTobjmapToken_, gtReadoutRecord);
  // const L1GlobalTriggerObjectMapRecord *l1trig = gtReadoutRecord.product();
  // for( int i=0; i<NL1SEED; i++ ){
  //   const L1GlobalTriggerObjectMap* trg = l1trig->getObjectMap(i);
  //   if(trg){
  // 	L1BitCollection_[trg->algoBitNumber()] = trg->algoGtlResult();
  // 	if( trg->algoGtlResult() ){
  // 	  triggerComposition->Fill( trg->algoBitNumber() );
  // 	}
  //   }
  // }
  // if( L1_Bit_Sele_!="" ){
  //   if ( L1_nameAndNumb.find(L1_Bit_Sele_.Data()) != L1_nameAndNumb.end() ){
  // 	const L1GlobalTriggerObjectMap* trg = l1trig->getObjectMap( L1_nameAndNumb[L1_Bit_Sele_.Data()] );
  // 	return trg->algoGtlResult();
  //   }
  //   else{
  // 	cout<<"WARNING!! L1_Bit_Sele_ is not in the list, I will return true!"<<endl;
  // 	return true;
  //   }
  // }
  // else{ return true;}  
  //cout << "Going out of FillEpsilonPlot::getTriggerResult()" << endl;

  return true;

  //  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  //  iEvent.getByLabel( l1TriggerTag_, gtReadoutRecord);
  //  const DecisionWord& gtDecisionWord = gtReadoutRecord->decisionWord();
  //  int thisBit =0;
  //  for (std::vector<bool>::const_iterator itBit = gtDecisionWord.begin(); itBit != gtDecisionWord.end(); ++itBit, ++thisBit) {
  //    L1BitCollection_[thisBit] = gtDecisionWord.at(thisBit);
  //    if( gtDecisionWord.at(thisBit) ) triggerComposition->Fill(thisBit);
  //  }
  //  if( !L1_Bit_Sele_.Contains("") ){
  //    edm::ESHandle<L1GtTriggerMenu> menuRcd;
  //    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  //    return gtDecisionWord.at(l1TrigNames_[L1_Bit_Sele_.Data()]);
  //  }
  //  else{ return true;}
}

void FillEpsilonPlot::endJob(){

  outfile_->cd();
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
  TTree_JoshMva->Write();
#endif
#ifdef MVA_REGRESSIO_EE
  TTree_JoshMva_EE->Write();
#endif
#ifdef SELECTION_TREE
  CutVariables_EB->Write();
  CutVariables_EE->Write();
  Pi0Info_EB->Write();
  Pi0Info_EE->Write();
#endif
  if(MakeNtuple4optimization_){
    Tree_Optim->Write();
  }
  EventFlow_EB->Write();  if (isDebug_) EventFlow_EB_debug->Write();
  EventFlow_EE->Write();  if (isDebug_) EventFlow_EE_debug->Write();

  if (isEoverEtrue_) {
    allEoverEtrue_g1_EB->Write();
    allEoverEtrue_g1_EBnw->Write();
    allEoverEtrue_g1_EE->Write();
    allEoverEtrue_g1_EEnw->Write();
    allEoverEtrue_g2_EB->Write();
    allEoverEtrue_g2_EBnw->Write();
    allEoverEtrue_g2_EE->Write();
    allEoverEtrue_g2_EEnw->Write();
  } else {
    allEpsilon_EB->Write();
    allEpsilon_EBnw->Write();
    allEpsilon_EE->Write();
    allEpsilon_EEnw->Write();
  }
  entries_EEp->Write();
  entries_EEm->Write();
  entries_EB->Write();
  Occupancy_EEp->Write();
  Occupancy_EEm->Write();
  Occupancy_EB->Write();
  pi0MassVsIetaEB->Write();
  pi0MassVsETEB->Write();
  if (L1TriggerInfo_) {
    triggerComposition->Write();
    triggerComposition_EB->Write();
    triggerComposition_EE->Write();
  }
  if (isMC_ and MC_Assoc_) {
    h_numberUnmergedGenPhotonPairs_EB->Write();
    h_numberMatchedGenPhotonPairs_EB->Write();
    h_numberUnmergedGenPhotonPairs_EE->Write();
    h_numberMatchedGenPhotonPairs_EE->Write();
    h_numberUnmergedGenPhotonPairs->Write();
    h_numberMatchedGenPhotonPairs->Write();
  }

  if (fillKinematicVariables_) {

    for (uint32_t i = 0; i < pi0pt_afterCuts.size(); ++i) {
      pi0pt_afterCuts[i]->Write();
      g1pt_afterCuts[i]->Write();
      g2pt_afterCuts[i]->Write();
      g1Nxtal_afterCuts[i]->Write();
      g2Nxtal_afterCuts[i]->Write();
      pi0PhotonsNoverlappingXtals_afterCuts[i]->Write();
      if (isMC_) {
	pi0MassVsPU[i]->Write();
      }

    }

  }


  if( !MakeNtuple4optimization_ &&(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    if (isEoverEtrue_) {
      writeEpsilonPlot(EoverEtrue_g1_EB_h, "Barrel" ,  regionalCalibration_->getCalibMap()->getNRegionsEB() );
      writeEpsilonPlot(EoverEtrue_g2_EB_h, "Barrel" ,  regionalCalibration_g2_->getCalibMap()->getNRegionsEB() );
    } else {
      writeEpsilonPlot(epsilon_EB_h, "Barrel" ,  regionalCalibration_->getCalibMap()->getNRegionsEB() );
    }
  }

  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {
    if (isEoverEtrue_) {
      writeEpsilonPlot(EoverEtrue_g1_EE_h, "Endcap" ,  regionalCalibration_->getCalibMap()->getNRegionsEE() );
      writeEpsilonPlot(EoverEtrue_g2_EE_h, "Endcap" ,  regionalCalibration_g2_->getCalibMap()->getNRegionsEE() );
    } else {
      writeEpsilonPlot(epsilon_EE_h, "Endcap" ,  regionalCalibration_->getCalibMap()->getNRegionsEE() );
    }
  }

#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
  delete TTree_JoshMva;
#endif
#ifdef MVA_REGRESSIO_EE
  delete TTree_JoshMva_EE;
#endif

  std::cout << "### FillEpsilonPlot::endJob()" << std::endl;

}

// ------------ EBPHI LOAD Containment correction  ------------
void FillEpsilonPlot::EBPHI_Cont_Corr_load(std::string FileName)
{
  cout << "FillEpsilonPlot:: loading phi contaiment corrections from " << FileName << endl;

  TFile* f = TFile::Open(FileName.c_str());

  if(!f)     cout << "Invalid file " << FileName << " .. try again" << endl;
  else{
    EBPHI_ConCorr_p = (TH1F*) f->Get("EBp_PHIFitContCorr");
    EBPHI_ConCorr_m = (TH1F*) f->Get("EBm_PHIFitContCorr");
  }
  f->Close();
}
 
// ------------ EBPHI Containment correction  ------------
float FillEpsilonPlot::EBPHI_Cont_Corr(float PT, int giPhi, int ieta)
{

  // Choos PT bin
  int ien=0;
  double PtBinBoundEB[7];
  PtBinBoundEB[0]=0.; PtBinBoundEB[1]=0.9; PtBinBoundEB[2]=1.5; PtBinBoundEB[3]=2.1; PtBinBoundEB[4]=3.; PtBinBoundEB[5]=5.; PtBinBoundEB[6]=8.;

  for(ien=0; ien < 7; ++ien) {
    if(PT <= PtBinBoundEB[ien+1]) break;
  }
  if(ien==7) ien=6;
  if(giPhi==0) giPhi=20;
  int nBin = 20*ien+giPhi;

  float Correction = 1.;
  if(ieta>0) Correction = EBPHI_ConCorr_p->GetBinContent(nBin);    
  else       Correction = EBPHI_ConCorr_m->GetBinContent(nBin);    

  if(Correction > 0.85){ return 1./Correction;}
  else{               
    //cout<<"Cont. Correction too low... I'm using 1. Check if all is right please. (nBin = "<<nBin<<" )"<<endl;  
    return 1.;
  }
}

// ------------ EB E/Etrue Containment correction  ------------
void FillEpsilonPlot::loadEoverEtrueContainmentCorrections(const std::string& fileName = "")
{

  cout << "FillEpsilonPlot:: loading E/Etrue containment corrections from " << fileName << endl;

  hCC_EoverEtrue_g1 = new TH2F("hCC_EoverEtrue_g1","",171,-85.5,85.5,360,0.5,360.5);
  hCC_EoverEtrue_g2 = new TH2F("hCC_EoverEtrue_g2","",171,-85.5,85.5,360,0.5,360.5);

  TH2F* hCC_tmp1 = nullptr;
  TH2F* hCC_tmp2 = nullptr;

  TFile* f = TFile::Open(fileName.c_str());

  if (!f) throw cms::Exception("loadEoverEtrueCC") << "Could not open file with containment corrections\n";
  else {

    // hCC_tmp1 = (TH2F*) ((TH2F*) f->Get("calibMap_EB"))->Clone();
    // hCC_tmp2 = (TH2F*) ((TH2F*) f->Get("calibMap_EB_g2"))->Clone();
    // make sure to change this name to avoid possible conflicts with the IC histograms, which have the same name
    hCC_tmp1 = (TH2F*) ((TH2F*) f->Get("calibMap_EB"))->Clone("hCC_tmp1_CC_EoverEtrue");
    hCC_tmp2 = (TH2F*) ((TH2F*) f->Get("calibMap_EB_g2"))->Clone("hCC_tmp2_CC_EoverEtrue");
    if (!hCC_tmp1) throw cms::Exception("loadEoverEtrueCC") << "Could not get histograms with containment corrections for photon 1\n";
    if (!hCC_tmp2) throw cms::Exception("loadEoverEtrueCC") << "Could not get histograms with containment corrections for photon 2\n";
  }
  // detach histogram from file so that we can safely close the file
  // hCC_tmp1->SetDirectory(0);
  // hCC_tmp2->SetDirectory(0);
  // f->Close();
  // better to close the file directly at the end

  if (hCC_tmp1->GetNbinsX() != hCC_tmp2->GetNbinsX()) 
    throw cms::Exception("loadEoverEtrueCC") << "Histograms with containment corrections for photon 1 and 2 have different number of X bins\n";
  if (hCC_tmp1->GetNbinsY() != hCC_tmp2->GetNbinsY()) 
    throw cms::Exception("loadEoverEtrueCC") << "Histograms with containment corrections for photon 1 and 2 have different number of Y bins\n";

  // map with CC from E/Etrue has ieta on X axis, but do check if you are not sure (someone might have changed it)
  for (Int_t ix = 1; ix <= hCC_tmp1->GetNbinsX(); ++ix) {
    for (Int_t iy = 1; iy <= hCC_tmp1->GetNbinsY(); ++iy) {
      hCC_EoverEtrue_g1->SetBinContent(ix,iy, hCC_tmp1->GetBinContent(ix,iy));
      hCC_EoverEtrue_g2->SetBinContent(ix,iy, hCC_tmp2->GetBinContent(ix,iy));
    }
  }

  f->Close();

}



// ------------ method called when starting to processes a run  ------------
void FillEpsilonPlot::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {

  //    edm::ESHandle<L1GtTriggerMenu> menuRcd;
  //    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  //    const L1GtTriggerMenu* menu = menuRcd.product();
  //    std::map< std::string, int >::iterator currentTrigger;
  //
  //    if(l1TrigNames_.size()>0) {
  //	bool triggerChanged = false;
  //	for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
  //	  currentTrigger = l1TrigNames_.find((algo->second).algoName());
  //	  if (currentTrigger == l1TrigNames_.end() || currentTrigger->second != (algo->second).algoBitNumber()) {
  //	    triggerChanged = true;
  //	    break;
  //	  }
  //	}
  //	if(!triggerChanged) return;
  //	cout << "beginRun:: Trigger names / ordering changed" << endl;
  //    }
  //    cout << "beginRun:: Filling trigger names" << endl;
  //    // filling trigger map
  //    l1TrigNames_.clear();
  //    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
  //	l1TrigNames_[(algo->second).algoName()] = (algo->second).algoBitNumber();
  //	// using same loop to set trigger histogram labels
  //	if(!areLabelsSet_)
  ////cout<<"NAME "<<(algo->second).algoBitNumber()+1<<" "<<(algo->second).algoName().c_str()<<endl;
  //	  triggerComposition->GetXaxis()->SetBinLabel((algo->second).algoBitNumber()+1,(algo->second).algoName().c_str());
  //    }
  //    if(!areLabelsSet_){
  //	areLabelsSet_ = true;
  //	cout << "beginRun:: setting labels of triggerComposition histogram" << endl;
  //    }
}

bool FillEpsilonPlot::checkStatusOfEcalRecHit(const EcalChannelStatus &channelStatus,const EcalRecHit &rh){
  int status =  int(channelStatus[rh.id().rawId()].getStatusCode()); 
  if ( status > 0/*statusLevelRecHitsToUsea_*/ ) return false; 
  return true; 
}

bool FillEpsilonPlot::isInDeadMap( bool isEB, const EcalRecHit &rh ){
  bool isBad=false;
  if(isEB){
    EBDetId det(rh.id());
    int ieta = det.ieta();
    int iphi = det.iphi();
    if( EBMap_DeadXtal->GetBinContent( iphi+1, ieta+86 ) == 1 ) isBad=true;
  }
  else{
    EEDetId det(rh.id());
    int ix = det.ix();
    int iy = det.iy();
    int iz = det.zside();
    if( iz==-1 ){
	if( EEmMap_DeadXtal->GetBinContent( ix+1, iy+1 ) == 1 ) isBad=true;
    }
    else{
	if( EEpMap_DeadXtal->GetBinContent( ix+1, iy+1 ) == 1 ) isBad=true;
    }
  }
  return isBad;
}

// ------------ method called when ending the processing of a run  ------------
  void 
FillEpsilonPlot::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
  void 
FillEpsilonPlot::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
  void 
FillEpsilonPlot::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FillEpsilonPlot::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double FillEpsilonPlot::min( double a, double b)
{
  if( a<=b ) return a;
  else       return b;
}

double max_array(double *A, int n){
  if(n==1)  return A[0];
  else      return max(A[0],max_array(A+1, n-1));
}

double max (double x, double y){
  if(x>=y)  return x;
  else      return y;
}

int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3){
  int index(0);
  bool found = false;
  if(debug3) cout<<"--> Looking for "<<x<<" "<<y<<endl;
  for( size_t i=0; i<VectRing.size(); i++){
    if(  VectRing[i].iX != x || VectRing[i].iY != y ){index++;    if(debug3){cout<<"Is not "<<VectRing[i].iX<<" and "<<VectRing[i].iY<<" index is: "<<index<<endl;}}
    if(  VectRing[i].iX == x && VectRing[i].iY == y ){found=true; if(debug3){cout<<"FOUND! "<<VectRing[i].iX<<" "<<VectRing[i].iY<<endl;}                   break;}
  }
  if(found){ if(debug3){cout<<"Returning: "<<VectRing[index].iX<<" "<<VectRing[index].iY<<" "<<VectRing[index].Ring<<endl;} return VectRing[index].Ring;}
  else{      if(debug3){cout<<"NOT found: "<<x<<" "<<y<<endl;}             return -1;}
}


//define this as a plug-in
DEFINE_FWK_MODULE(FillEpsilonPlot);
