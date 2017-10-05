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

using namespace TMVA;
using namespace edm;

//Function
double max_array(double *A, int n);
double max(double x, double y);
int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3);


FillEpsilonPlot::FillEpsilonPlot(const edm::ParameterSet& iConfig)
{
    /// to be moved in parameters.py
    useMassInsteadOfEpsilon_ = 1;

    /// parameters from python
    Are_pi0_                           = iConfig.getUntrackedParameter<bool>("Are_pi0",true);
    useMVAContainmentCorrections_      = iConfig.getUntrackedParameter<bool>("useMVAContainmentCorrections",true);
    new_pi0ContainmentCorrections_     = iConfig.getUntrackedParameter<bool>("new_pi0ContainmentCorrections",false);

    EBRecHitCollectionToken_           = consumes<EBRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag"));
    EERecHitCollectionToken_           = consumes<EERecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag"));
    ESRecHitCollectionToken_           = consumes<ESRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag"));
    //ESRecHitCollectionToken_           = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag"));
    HLTResults_                        = iConfig.getUntrackedParameter<bool>("HLTResults",false);
    HLTResultsNameEB_                  = iConfig.getUntrackedParameter<std::string>("HLTResultsNameEB","AlCa_EcalPi0EB");
    HLTResultsNameEE_                  = iConfig.getUntrackedParameter<std::string>("HLTResultsNameEE","AlCa_EcalPi0EE");
    RemoveDead_Flag_                   = iConfig.getUntrackedParameter<bool>("RemoveDead_Flag",false);
    RemoveDead_Map_                    = iConfig.getUntrackedParameter<std::string>("RemoveDead_Map");
    L1_Bit_Sele_                       = iConfig.getUntrackedParameter<std::string>("L1_Bit_Sele","");
    L1TriggerInfo_                     = iConfig.getUntrackedParameter<bool>("L1TriggerInfo",false);
    l1TriggerTag_                      = iConfig.getUntrackedParameter<edm::InputTag>("L1TriggerTag");
    triggerResultsToken_               = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter("triggerTag",edm::InputTag("TriggerResults")));
    //L1GTobjmapToken_                   = consumes<L1GlobalTriggerObjectMapRecord>(iConfig.getUntrackedParameter<edm::InputTag>("hltL1GtObjectMap",edm::InputTag("hltL1GtObjectMap")));
    // edm::Handle<GlobalAlgBlkBxCollection> & l1results
    L1GTobjmapToken_                   = consumes<GlobalAlgBlkBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("hltGtStage2Digis",edm::InputTag("hltGtStage2Digis")));
    GenPartCollectionToken_            = consumes<GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("GenPartCollectionTag",edm::InputTag("genParticles")));
    outfilename_                       = iConfig.getUntrackedParameter<std::string>("OutputFile");
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
    isMC_                              = iConfig.getUntrackedParameter<bool>("isMC",false);
    MC_Asssoc_                         = iConfig.getUntrackedParameter<bool>("MC_Asssoc",false);
    MC_Asssoc_DeltaR                   = iConfig.getUntrackedParameter<double>("MC_Asssoc_DeltaR",0.1);
    MakeNtuple4optimization_           = iConfig.getUntrackedParameter<bool>("MakeNtuple4optimization",false);
    GeometryFromFile_                  = iConfig.getUntrackedParameter<bool>("GeometryFromFile",false);
    JSONfile_                          = iConfig.getUntrackedParameter<std::string>("JSONfile","");

    L1SeedsPi0Stream_                  = iConfig.getUntrackedParameter<std::string>("L1SeedsPi0Stream");    
    nL1SeedsPi0Stream_                 = iConfig.getUntrackedParameter<int>("nL1SeedsPi0Stream");

    isDebug_                           = iConfig.getUntrackedParameter<bool>("isDebug",false);

    // for MC-truth association
    g4_simTk_Token_  = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
    g4_simVtx_Token_ = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));

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
    if( JSONfile_!="" ) myjson=new JSON( edm::FileInPath( JSONfile_.c_str() ).fullPath().c_str() );
    // shower shape parameters
    PCparams_.param_LogWeighted_ = true;
    PCparams_.param_T0_barl_     = 7.4;
    PCparams_.param_T0_endc_     = 3.1;
    PCparams_.param_T0_endcES_   = 1.2;
    PCparams_.param_W0_          = 4.2;
    PCparams_.param_X0_          = 0.89;

    /// setting calibration type
    calibTypeString_ = iConfig.getUntrackedParameter<std::string>("CalibType");
    if(     calibTypeString_.compare("xtal")    == 0 ) { calibTypeNumber_ = xtal;    regionalCalibration_ = &xtalCalib; } 
    else if(calibTypeString_.compare("tt")      == 0 ) { calibTypeNumber_ = tt;      regionalCalibration_ = &TTCalib;   }
    else if(calibTypeString_.compare("etaring") == 0 ) { calibTypeNumber_ = etaring; regionalCalibration_ = &etaCalib;  }
    else throw cms::Exception("CalibType") << "Calib type not recognized\n";
    cout << "crosscheck: selected type: " << regionalCalibration_->printType() << endl;

    /// external hardcoded geometry

    externalGeometryFile_ = TFile::Open( edm::FileInPath( externalGeometry_.c_str() ).fullPath().c_str() );
    if(!externalGeometryFile_) cms::Exception("ExtGeom") << "External Geometry file (" << externalGeometry_ << ") not found" << endl;
    geom_ = ECALGeometry::getGeometry(externalGeometryFile_);
    GeometryService::setGeometryName(externalGeometry_);
    GeometryService::setGeometryPtr(geom_);
    // containment corrections
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
	  char fileName[200];
	  cout << "FillEpsilonPlot:: loading calibraion map at " << calibMapPath_ << endl;
	  if( isCRAB_ ) sprintf(fileName,"%s",  edm::FileInPath( calibMapPath_.c_str() ).fullPath().c_str() );
	  else          sprintf(fileName,"%s", calibMapPath_.c_str());
	  regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
    }

    /// epsilon histograms
    if(!MakeNtuple4optimization_){
      if(useMassInsteadOfEpsilon_ ){
	  if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","#pi^{0} Mass distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
	  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","#pi^{0} Mass distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
	}
      else{
	  if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","Epsilon distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
	  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","Epsilon distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
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

    allEpsilon_EB = new TH1F("allEpsilon_EB", "allEpsilon_EB",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
    allEpsilon_EBnw = new TH1F("allEpsilon_EBnw", "allEpsilon_EBnw",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
    allEpsilon_EE = new TH1F("allEpsilon_EE", "allEpsilon_EE",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
    allEpsilon_EEnw = new TH1F("allEpsilon_EEnw", "allEpsilon_EEnw",240, Are_pi0_? 0.:0.3 , Are_pi0_? 0.5:0.8 );
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
    char fileName[200];
    sprintf(fileName,"%s%s", outputDir_.c_str(), outfilename_.c_str());
    outfile_ = new TFile(fileName,"RECREATE");
    if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << string(fileName) << "\n";
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
	Tree_Optim->Branch( "STr2_NPi0_rec",      &Op_NPi0_rec,         "STr2_NPi0_rec/I");
	Tree_Optim->Branch( "STr2_Pi0recIsEB",    &Op_Pi0recIsEB,       "STr2_Pi0recIsEB[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_IsoPi0_rec",    &Op_IsoPi0_rec,       "STr2_IsoPi0_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_HLTIsoPi0_rec", &Op_HLTIsoPi0_rec,    "STr2_HLTIsoPi0_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_n1CrisPi0_rec", &Op_n1CrisPi0_rec,    "STr2_n1CrisPi0_rec[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_n2CrisPi0_rec", &Op_n2CrisPi0_rec,    "STr2_n2CrisPi0_rec[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_mPi0_rec",      &Op_mPi0_rec,         "STr2_mPi0_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_enG1_rec",      &Op_enG1_rec,         "STr2_enG1_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_enG2_rec",      &Op_enG2_rec,         "STr2_enG2_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_etaPi0_rec",    &Op_etaPi0_rec,       "STr2_etaPi0_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_ptPi0_rec",     &Op_ptPi0_rec,        "STr2_ptPi0_rec[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_ptPi0_nocor",   &Op_ptPi0_nocor,      "STr2_ptPi0_nocor[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_enG1_nocor",    &Op_enG1_nocor,       "STr2_enG1_nocor[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_enG2_nocor",    &Op_enG2_nocor,       "STr2_enG2_nocor[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_mPi0_nocor",    &Op_mPi0_nocor,       "STr2_mPi0_nocor[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_DeltaRG1G2",    &Op_DeltaRG1G2,       "STr2_DeltaRG1G2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Es_e1_1",       &Op_Es_e1_1,          "STr2_Es_e1_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Es_e1_2",       &Op_Es_e1_2,          "STr2_Es_e1_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Es_e2_1",       &Op_Es_e2_1,          "STr2_Es_e2_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Es_e2_2",       &Op_Es_e2_2,          "STr2_Es_e2_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S4S9_1",        &Op_S4S9_1,           "STr2_S4S9_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S4S9_2",        &Op_S4S9_2,           "STr2_S4S9_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S2S9_1",        &Op_S2S9_1,           "STr2_S2S9_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S2S9_2",        &Op_S2S9_2,           "STr2_S2S9_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S1S9_1",        &Op_S1S9_1,           "STr2_S1S9_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_S1S9_2",        &Op_S1S9_2,           "STr2_S1S9_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Eta_1",         &Op_Eta_1,            "STr2_Eta_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Eta_2",         &Op_Eta_2,            "STr2_Eta_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Phi_1",         &Op_Phi_1,            "STr2_Phi_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Phi_2",         &Op_Phi_2,            "STr2_Phi_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Time_1",        &Op_Time_1,           "STr2_Time_1[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_Time_2",        &Op_Time_2,           "STr2_Time_2[STr2_NPi0_rec]/F");
	Tree_Optim->Branch( "STr2_iEtaiX_1",      &Op_iEtaiX_1,         "STr2_iEtaiX_1[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iEtaiX_2",      &Op_iEtaiX_2,         "STr2_iEtaiX_2[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhiiY_1",      &Op_iPhiiY_1,         "STr2_iPhiiY_1[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhiiY_2",      &Op_iPhiiY_2,         "STr2_iPhiiY_2[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iEta_1on5",     &Op_iEta_1on5,        "STr2_iEta_1on5[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iEta_2on5",     &Op_iEta_2on5,        "STr2_iEta_2on5[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhi_1on2",     &Op_iPhi_1on2,        "STr2_iPhi_1on2[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhi_2on2",     &Op_iPhi_2on2,        "STr2_iPhi_2on2[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iEta_1on2520",  &Op_iEta_1on2520,     "STr2_iEta_1on2520[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iEta_2on2520",  &Op_iEta_2on2520,     "STr2_iEta_2on2520[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhi_1on20",    &Op_iPhi_1on20,       "STr2_iPhi_1on20[STr2_NPi0_rec]/I");
	Tree_Optim->Branch( "STr2_iPhi_2on20",    &Op_iPhi_2on20,       "STr2_iPhi_2on20[STr2_NPi0_rec]/I");
	if( isMC_ && MC_Asssoc_ ) {
	  Tree_Optim->Branch( "STr2_enG1_true",     &Op_enG1_true,        "STr2_enG1_true[STr2_NPi0_rec]/F");
	  Tree_Optim->Branch( "STr2_enG2_true",     &Op_enG2_true,        "STr2_enG2_true[STr2_NPi0_rec]/F");
	  Tree_Optim->Branch( "STr2_DeltaR_1",      &Op_DeltaR_1,         "STr2_DeltaR_1[STr2_NPi0_rec]/F");
	  Tree_Optim->Branch( "STr2_DeltaR_2",      &Op_DeltaR_2,         "STr2_DeltaR_2[STr2_NPi0_rec]/F");
	  Tree_Optim->Branch( "STr2_Nxtal_1",       &Op_Nxtal_1,          "STr2_Nxtal_1[STr2_NPi0_rec]/I");
	  Tree_Optim->Branch( "STr2_Nxtal_2",       &Op_Nxtal_2,          "STr2_Nxtal_2[STr2_NPi0_rec]/I");
	}
    }

    /// trigger histo
    if (L1TriggerInfo_) {
      triggerComposition = new TH1F("triggerComposition", "Trigger Composition", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);    
      triggerComposition_EB = new TH1F("triggerComposition_EB", "Trigger Composition in EB", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);
      triggerComposition_EE = new TH1F("triggerComposition_EE", "Trigger Composition in EE", nL1SeedsPi0Stream_, -0.5, (double)nL1SeedsPi0Stream_ -0.5);
    }
    areLabelsSet_ = false;
    //L1_nameAndNumb.clear();
    //for(unsigned int i=0; i<NL1SEED; i++) L1BitCollection_[i]=-1;

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
    // for L1
    seedIsInStream = new int[GlobalAlgBlk::maxPhysicsTriggers];
    algoBitToName = new TString[GlobalAlgBlk::maxPhysicsTriggers];
    l1flag = new short[GlobalAlgBlk::maxPhysicsTriggers];

}

FillEpsilonPlot::~FillEpsilonPlot()
{
  externalGeometryFile_->Close();
  outfile_->Write();
  outfile_->Close();

  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) deleteEpsilonPlot(epsilon_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) deleteEpsilonPlot(epsilon_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
  delete allEpsilon_EB;
  delete allEpsilon_EBnw;
  delete allEpsilon_EE;
  delete allEpsilon_EEnw;
  delete entries_EEp;
  delete entries_EEm;
  delete entries_EB;
  delete Occupancy_EEp;
  delete Occupancy_EEm;
  delete Occupancy_EB;
  delete pi0MassVsIetaEB;
  delete pi0MassVsETEB;
  delete triggerComposition;
  delete triggerComposition_EB;
  delete triggerComposition_EE;

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
  delete myjson;
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
  delete[] seedIsInStream;
  delete[] algoBitToName;
  delete[] l1flag;


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
  //Trigger Histo

  myEvent = iEvent.id().event();
  myLumiBlock = iEvent.id().luminosityBlock();
  myRun = iEvent.id().run();
  myBunchCrossing = iEvent.bunchCrossing();

  if( !areLabelsSet_ && L1TriggerInfo_ ){
    // edm::Handle< L1GlobalTriggerObjectMapRecord > gtReadoutRecord;
    // iEvent.getByToken( L1GTobjmapToken_, gtReadoutRecord);
    // const L1GlobalTriggerObjectMapRecord *l1trig = gtReadoutRecord.product();
    // for( int i=0; i<NL1SEED; i++ ){
    // 	const L1GlobalTriggerObjectMap* trg = l1trig->getObjectMap(i);
    // 	if(trg){
    // 	  L1_nameAndNumb[trg->algoName()] = trg->algoBitNumber();
    // 	  triggerComposition->GetXaxis()->SetBinLabel(trg->algoBitNumber()+1,trg->algoName().c_str());
    // 	}
    //  if(!areLabelsSet_){
    //    areLabelsSet_ = true;
    // 	  cout << "setting labels of triggerComposition histogram" << endl;
    //  }
    // }
    
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
  if( isMC_ && MC_Asssoc_ ){
    edm::Handle<std::vector<reco::GenParticle>> GenParProd;
    iEvent.getByToken( GenPartCollectionToken_, GenParProd);//Fatal Root Error: @SUB=TBufferFile::CheckByteCount object of class edm::RefCore read too many bytes: 10 instead of 8
    // const reco::GenParticleCollection *GenPars = 0;
    // std::cout << "MC truth taken" << std::endl;
    // //if ( ! GenParProd.isValid() )  edm::LogWarning("GenParSummary") << "GenPars not found";
    // if ( ! GenParProd.isValid() )  std::cout << "GenPars not found" << std::endl;
    // GenPars = GenParProd.product();


    // GUN sample made with PYTHIA6 doesn't decay the pi0, need to look at simtracks by GEANT
    // get GEANT sim tracks and vertices (includes conversions)
    Handle<SimTrackContainer> simTracks_h;
    const SimTrackContainer* simTracks;
    iEvent.getByToken(g4_simTk_Token_, simTracks_h);
    simTracks = (simTracks_h.isValid()) ? simTracks_h.product() : 0;

    Handle<SimVertexContainer> simVert_h;
    const SimVertexContainer* simVertices;
    iEvent.getByToken(g4_simVtx_Token_, simVert_h);
    simVertices = (simVert_h.isValid()) ? simVert_h.product() : 0;


    // Vertices only return trackID of their parent SimTrack
    // Figure out the mapping from trackID to SimTrack
    map<unsigned int, const SimTrack*> trackMap;
    for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {
      if (!iSim->noVertex()) {
        assert(trackMap.find(iSim->trackId())==trackMap.end());
        trackMap[iSim->trackId()] = &(*iSim);
      }
    }

    // Find all SimTracks that come from decays before the ECAL
    // and find their parent SimTracks
    map<const SimTrack*, const SimTrack*> promptParent; // daughter->mother
    map<const SimTrack*, set<const SimTrack*> > promptDecays; // m->ds
    map<const SimTrack*, const SimVertex*> promptVertex; // daughter->vertex
    map<const SimTrack*, const SimVertex*> promptALLVertex; // daughter->vertex in Any Occasion
    map<const SimTrack*, const SimTrack*> promptALLParent; // daughter->mother in Any Occasion

    int num=0;
    for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim, num++) 
      {
        if (!iSim->noVertex()) 
          {
            // Find the parent vertex and see if it classifies as an early decay
            // Exclude the primary vertex (noParent)
            SimVertex const& vtx = (*simVertices)[iSim->vertIndex()];
            if (!vtx.noParent() && vtx.position().Rho() < 129 && fabs(vtx.position().z()) < 304) 
              {
                // Find parent SimParticle that produced this vertex
                // vtx->parentIndex is NOT a vector index :( so use trackMap
                assert(trackMap.find(vtx.parentIndex())!=trackMap.end());
                const SimTrack* p = trackMap[vtx.parentIndex()];
                promptParent[&(*iSim)] = p; // nel Pi0Gun: ->genpartIndex() e' -1
                promptDecays[p].insert(&(*iSim));
                promptVertex[&(*iSim)] = &vtx;
              } // early decay
            if (!vtx.noParent() ){
              promptALLVertex[&(*iSim)] = &vtx;
              const SimTrack* p = trackMap[vtx.parentIndex()];
              promptALLParent[&(*iSim)] = p; 
              // cout<<num<<" PDG: "<<iSim->type()<<" id: "<<iSim->trackId()<<" Son of: "<<p->type()<<" id: "<<p->trackId()
              // <<" x vtx:" <<vtx.position().x()<<" Z Vtx: "<<vtx.position().z()<<"Px "<<iSim->momentum().x()<<" py "<<iSim->momentum().y()<<" pz "<<iSim->momentum().z()<<
              // "pt "<<sqrt(pow(iSim->momentum().x(),2)+pow(iSim->momentum().y(),2)+pow(iSim->momentum().z(),2))<<endl;
            }
          } // has vertex
      } // for simTracks

    //cout<<"Event:"<<endl;
    //Store Pi0 & gamma
    //    unsigned int IdGamma1=0, IdGamma2=0;
    TVector3 pi0_pos, ga1, ga2;
    num=0;
    for (SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim, num++){
      if (!iSim->noVertex() ){
        SimVertex const& vtx = (*simVertices)[iSim->vertIndex()];
        if( !vtx.noParent() ) {
          
          if( iSim->type()==22 && promptALLParent[&(*iSim)]->type()==(Are_pi0_ ? 111:221) && num==1){
            pi0_pos.SetXYZ(promptALLVertex[&(*iSim)]->position().x(),promptALLVertex[&(*iSim)]->position().y(),promptALLVertex[&(*iSim)]->position().z());
            Gamma1MC.SetXYZ(iSim->momentum().x(),iSim->momentum().y(),iSim->momentum().z());
            //            IdGamma1 = iSim->trackId();
            //            std::cout << "Photon1 eta,phi = " << Gamma1MC.eta() << "  " << Gamma1MC.phi() << std::endl;
         }
          if( iSim->type()==22 && promptALLParent[&(*iSim)]->type()==(Are_pi0_ ? 111:221) && num==2){
            Gamma2MC.SetXYZ(iSim->momentum().x(),iSim->momentum().y(),iSim->momentum().z());
            //            IdGamma2 = iSim->trackId();
            //            std::cout << "Photon2 eta,phi = " << Gamma2MC.eta() << "  " << Gamma2MC.phi() << std::endl;
          }
        }
      }
    }

    //Find MC photons
    /*
    std::cout << "coll size = " << GenPars->size() << std::endl;
    bool firstnotfound = true;
    //    for (auto& GenPar : *GenPars){
    for (reco::GenParticleCollection::const_iterator GenPar = GenPars->begin(); GenPar != GenPars->end(); ++GenPar) {
      std::cout << "id = " << GenPar->pdgId() << std::endl;
      if(GenPar->mother()!=0) std::cout << " mothId = " << GenPar->mother()->pdgId() << std::endl;
      
      int motherID = Are_pi0_ ? 111:221;
      if( GenPar->pdgId()==22 && GenPar->mother()->pdgId()==motherID && firstnotfound ){
        std::cout << "Found 1st photon, pt = " << GenPar->pt() << "  " << GenPar->p4().Eta() << "  " << GenPar->p4().Phi() << std::endl;
        std::cout << "dentro id = " << GenPar->pdgId() << std::endl;
        Gamma1MC.SetPtEtaPhiE( GenPar->pt(), GenPar->p4().Eta(), GenPar->p4().Phi(), GenPar->p4().E() );
        firstnotfound = false;
      }
      if( GenPar->pdgId()==22 && GenPar->mother()->pdgId()==motherID && GenPar->p4().Eta() != Gamma1MC.Eta() ){
        std::cout << "Found 2nd photon, pt = " << GenPar->pt() << "  " << GenPar->p4().Eta() << "  " << GenPar->p4().Phi() << std::endl;
        Gamma2MC.SetPtEtaPhiE( GenPar->pt(), GenPar->p4().Eta(), GenPar->p4().Phi(), GenPar->p4().E() );
      }
      std::cout << "running PT1,PT2 = " << Gamma1MC.Pt() << " , " << Gamma2MC.Pt() << std::endl;
    }
    std::cout << "==> final PT1,PT2 = " << Gamma1MC.Pt() << " , " << Gamma2MC.Pt() << std::endl;
    */
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
  //if( HLTResults_ && (!MakeNtuple4optimization_) ){
  if( HLTResults_ ){
    EB_HLT = GetHLTResults(iEvent, HLTResultsNameEB_); //Adding * at the end of the sentence make always true the "->Contains" method. So do not use it.
    EE_HLT = GetHLTResults(iEvent, HLTResultsNameEE_);
  }


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
  if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) && EB_HLT ){ EventFlow_EB->Fill(3.); fillEBClusters(ebclusters, iEvent, channelStatus);}
  ////cout << "I'm after fillEBClusters(ebclusters, iEvent, channelStatus) " << endl;
  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) && EE_HLT ){ EventFlow_EE->Fill(3.); fillEEClusters(eseeclusters, eseeclusters_tot, iEvent, channelStatus);}
  ////cout << "I'm after fillEEClusters(eseeclusters, eseeclusters_tot, iEvent, ...) " << endl;

  std::vector< CaloCluster > ebclusters_used, eeclusters_used;
  if(isMC_ && MC_Asssoc_) {
    ebclusters_used = MCTruthAssociate(ebclusters,MC_Asssoc_DeltaR,true);
    eeclusters_used = MCTruthAssociate(eseeclusters_tot,MC_Asssoc_DeltaR,false);
  } else {
    ebclusters_used = ebclusters;
    eeclusters_used = eseeclusters_tot;
  }

  if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(ebclusters_used, EcalBarrel);
  if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(eeclusters_used, EcalEndcap);
  
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
	clus_used.push_back(std::make_pair(*det,1.));

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
	const CaloCellGeometry* cell=geometry->getGeometry( seed_id );
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
	  if(!checkStatusOfEcalRecHit(channelStatus, *RecHitsInWindow[j] )  ) All_rechit_good = false; 
	}
	if( RemoveDead_Map_!="" ){
	  if( isInDeadMap( true, *RecHitsInWindow[j] ) ) All_rechit_good = false;
	}

	int ieta = det.ieta();
	int iphi = det.iphi();
	convxtalid(iphi,ieta);

	// use calibration coeff for energy and position
	float en = RecHitsInWindow[j]->energy() * regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id());
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
	  isUsed.insert(RecHitsInWindow[j]->id());

	}

	// compute position
	if(en>0.) 
	{
	  float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );
	  float pos_geo;
	  if( GeometryFromFile_ ) pos_geo = geom_->getPosition(det).mag(); // to front face
	  else                  {
	    const CaloCellGeometry* cell=geometry->getGeometry(det);
	    GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    pos_geo = posit.mag();
	  }
	  float depth = maxDepth + maxToFront - pos_geo;
	  GlobalPoint posThis;
	  if( GeometryFromFile_ ) posThis = geom_->getPosition(det,depth);
	  else{
	    const CaloCellGeometry* cell=geometry->getGeometry(det);
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
    if( fabs( clusPos.eta() )<1. ){ if(s4s9<S4S9_cut_low_[EcalBarrel]) continue;}
    else                          { if(s4s9<S4S9_cut_high_[EcalBarrel]) continue;}


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
	const CaloCellGeometry* cell=geometry->getGeometry(idXtal);
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
	clus_used.push_back(std::make_pair(*det,1.));
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
	const CaloCellGeometry* cell=geometry->getGeometry( eeseed_id );
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
	  if( !checkStatusOfEcalRecHit(channelStatus, *RecHitsInWindow[j] )  ) All_rechit_good = false;
	}
	if( RemoveDead_Map_!="" ){
	  if( isInDeadMap( false, *RecHitsInWindow[j] ) ) All_rechit_good = false;
	}

	int ix = det.ix();
	int iy = det.iy();

	// use calibration coeff for energy and position
	float en = RecHitsInWindow[j]->energy() * regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id());
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
	    const CaloCellGeometry* cell=geometry->getGeometry(det);
	    GlobalPoint posit = ( dynamic_cast<const TruncatedPyramid*>(cell) )->getPosition( 0. );
	    pos_geo = posit.mag();
	  }
	  float depth = maxDepth + maxToFront - pos_geo;
	  GlobalPoint posThis;
	  if( GeometryFromFile_ ) posThis = geom_->getPosition(det,depth);
	  else{
	    const CaloCellGeometry* cell=geometry->getGeometry(det);
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
    if( fabs( clusPos.eta() )<1.8 ){ if(s4s9<S4S9_cut_low_[EcalEndcap]) continue; }
    else                           { if(s4s9<S4S9_cut_high_[EcalEndcap]) continue; }

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
  TH1F **h = new TH1F*[size];
  char name_c[100];
  char title_c[200];

  cout << "FillEpsilonPlot::initializeEpsilonHistograms::useMassInsteadOfEpsilon_ = " << useMassInsteadOfEpsilon_ << endl;

  for(int jR=0; jR<size; jR++)
  {
    sprintf(name_c, "%s%d", name, jR);
    sprintf(title_c, "%s%d", title, jR);

    if(useMassInsteadOfEpsilon_)
    {
	h[jR] = new TH1F(name_c, title_c, 120, Are_pi0_? 0.:0.3, Are_pi0_? 0.5:0.8);
	h[jR]->GetXaxis()->SetTitle("Mass(#gamma#gamma)");
    }
    else
    {
	h[jR] = new TH1F(name_c, title_c, 120,-0.5,1);
	h[jR]->GetXaxis()->SetTitle("Epsilon");
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
  outfile_->mkdir(folder);
  outfile_->cd(folder);
  for(int jR=0; jR<size; jR++)
    h[jR]->Write();
}

std::vector< CaloCluster > FillEpsilonPlot::MCTruthAssociate(std::vector< CaloCluster > & clusters, double deltaR, bool isEB) {

  std::vector< CaloCluster > ret;
  ret.clear();
  int n_tmp1 = -1;    int n_tmp2 = -1;
  double deltaR1_tmp = 999;    double deltaR2_tmp = 999;
  if(isMC_ && MC_Asssoc_) {
    //    std::cout << "Association with MC: initial collection size = " << clusters.size() << std::endl;
    for(unsigned int i=0; i<clusters.size(); i++){
      const CaloCluster g = clusters[i];
      double deltaR1 = reco::deltaR(g.eta(),g.phi(), Gamma1MC.Eta(),Gamma1MC.Phi());
      double deltaR2 = reco::deltaR(g.eta(),g.phi(), Gamma2MC.Eta(),Gamma2MC.Phi());
      //      std::cout << "DR1,2 = " << deltaR1 << "  " << deltaR2 << std::endl;
      if(deltaR1<deltaR1_tmp) {
        deltaR1_tmp = deltaR1; n_tmp1 = i;
      }
      if(deltaR2<deltaR2_tmp) {
        deltaR2_tmp = deltaR2; n_tmp2 = i;
      }
    }
  } else return ret;

  // std::cout << "deltaR1,2_tmp = " << deltaR1_tmp << "  " << deltaR2_tmp << std::endl;
  // std::cout << "ntmp = " << n_tmp1 << "  " << n_tmp2 << std::endl;

  // one could look to another close cluster, but if the two MC photons are closer than 0.05 to the same cluster, 
  // that is merged in one cluster most probably, so another one would be fake
  if(n_tmp1==n_tmp2) return ret;
  
  if(deltaR1_tmp < deltaR) {
    if(isEB) Ncristal_EB_used.push_back(Ncristal_EB[n_tmp1]);
    else Ncristal_EE_used.push_back(Ncristal_EE[n_tmp1]);
    ret.push_back(clusters[n_tmp1]);
  }
  if(deltaR2_tmp < deltaR) {
    if(isEB) Ncristal_EB_used.push_back(Ncristal_EB[n_tmp2]);
    else Ncristal_EE_used.push_back(Ncristal_EE[n_tmp2]);
    ret.push_back(clusters[n_tmp2]);
  }
  return ret;
}

void FillEpsilonPlot::computeEpsilon(std::vector< CaloCluster > & clusters, int subDetId ) 
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
	    float new_value_pi01[12];

	    value_pi01[0] = ( G_Sort_1.E()/G_Sort_2.E() );
	    value_pi01[1] = ( G_Sort_1.Pt() );
	    value_pi01[2] = ( Ncristal_EB[ind1] );
	    value_pi01[3] = ( Ncristal_EB[ind2] );
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
	    
	    new_value_pi01[0] = ( G_Sort_1.E() );
            new_value_pi01[1] = ( G_Sort_1.Eta() );
            new_value_pi01[2] = ( G_Sort_1.Phi() );
            new_value_pi01[3] = ( Ncristal_EB[ind1] );
            new_value_pi01[4] = ( vs4s9[ind1] );
            new_value_pi01[5] = ( vs2s9[ind1] );
            new_value_pi01[6] = ( iEta1 );
            new_value_pi01[7] = ( iPhi1 );
            new_value_pi01[8] = ( iEta1%5 );
            new_value_pi01[9] = ( iEta1%2 );
            new_value_pi01[10] = ( iEta1%20 );
            new_value_pi01[11] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );

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
	    float new_value_pi02[12];

	    value_pi02[0] = ( G_Sort_1.E()/G_Sort_2.E() );
	    value_pi02[1] = ( G_Sort_2.Pt() );
	    value_pi02[2] = ( Ncristal_EB[ind1] );
	    value_pi02[3] = ( Ncristal_EB[ind2] );
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
	   
	    new_value_pi02[0] = ( G_Sort_2.E() );
            new_value_pi02[1] = ( G_Sort_2.Eta() );
            new_value_pi02[2] = ( G_Sort_2.Phi() );
            new_value_pi02[3] = ( Ncristal_EB[ind2] );
            new_value_pi02[4] = ( vs4s9[ind2] );
            new_value_pi02[5] = ( vs2s9[ind2] );
            new_value_pi02[6] = ( iEta2 );
            new_value_pi02[7] = ( iPhi2 );
            new_value_pi02[8] = ( iEta2%5 );
            new_value_pi02[9] = ( iEta2%2 );
            new_value_pi02[10] = ( iEta2%20 );
            new_value_pi02[11] = ((TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::     Abs(iEta2)/iEta2)%20) );

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
 	    value_pi01[0] = ( G_Sort_1.E()/G_Sort_2.E() );
	    value_pi01[1] = ( G_Sort_1.Pt() );
	    value_pi01[2] = ( Ncristal_EB[ind1] );
	    value_pi01[3] = ( iEta1 );
	    value_pi01[4] = ( iPhi1 );
	    value_pi01[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi01[6] = ( iEta1%5 );
	    value_pi01[7] = ( iPhi1%2 );
	    value_pi01[8] = ( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
	    value_pi01[9] = ( iPhi1%20 );
	    Correct1 = forest_EB_1->GetResponse(value_pi01);
	    float value_pi02[10];
	    value_pi02[0] = ( G_Sort_1.E()/G_Sort_2.E() );
	    value_pi02[1] = ( G_Sort_2.Pt() );
	    value_pi02[2] = ( Ncristal_EB[ind2] );
	    value_pi02[3] = ( iEta2 );
	    value_pi02[4] = ( iPhi2 );
	    value_pi02[5] = ( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
	    value_pi02[6] = ( iEta2%5 );
	    value_pi02[7] = ( iPhi2%2 );
	    value_pi02[8] = ( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
	    value_pi02[9] = ( iPhi2%20 );
	    Correct2 = forest_EB_2->GetResponse(value_pi02);
	  }

	  if( !Inverted ){ Corr1 = Correct1; Corr2 = Correct2; }
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
	  double GSort1plus2_EoverCoshEta = GSort1plus2.E()/cosh(GSort1plus2.Eta());

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
	  float new_value_pi01[8];
	  value_pi01[0] = ( GSort1plus2_EoverCoshEta );
	  value_pi01[1] = ( G_Sort_1.E()/ GSort1plus2_EoverCoshEta );
	  value_pi01[2] = ( G_Sort_1.Pt() );
	  value_pi01[3] = ( Ncristal_EE[ind1] );
	  value_pi01[4] = ( Ncristal_EE[ind2] );
	  value_pi01[5] = ( vs4s9EE[ind1] );
	  value_pi01[6] = ( vs1s9EE[ind1] );
	  value_pi01[7] = ( vs2s9EE[ind1] );
	  value_pi01[8] = ( ESratio[ind1] );
	  value_pi01[9] = ( EtaRing_1 );
	    new_value_pi01[0] = ( G_Sort_1.E() );
            new_value_pi01[1] = ( G_Sort_1.Eta() );
            new_value_pi01[2] = ( G_Sort_1.Phi() );
            new_value_pi01[3] = ( Ncristal_EE[ind1] );
            new_value_pi01[4] = ( vs4s9[ind1] );
            new_value_pi01[5] = ( vs2s9[ind1] );
            new_value_pi01[6] = ( iX1 );
            new_value_pi01[7] = ( iY1 );
	   
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
	    <<" ratio E "<< G_Sort_1.E()/GSort1plus2_EoverCoshEta<<" Pt "<<G_Sort_1.Pt()
	    <<" xtal "<<Ncristal_EE[ind1]<<" vs4s9EE "<<vs4s9EE[ind1]<<" vs1s9EE "<<vs1s9EE[ind1]<<" vs2s9EE "<<vs2s9EE[ind1]
	    <<" ESratio "<<ESratio[ind1]<<" EtaRing_1 "<<EtaRing_1<<endl;

	  float value_pi02[10];
	  float new_value_pi02[8];
	  value_pi02[0] = ( GSort1plus2_EoverCoshEta );
	  value_pi02[1] = ( G_Sort_2.E()/GSort1plus2_EoverCoshEta );
	  value_pi02[2] = ( G_Sort_2.Pt() );
	  value_pi02[3] = ( Ncristal_EE[ind1] );
	  value_pi02[4] = ( Ncristal_EE[ind2] );
	  value_pi02[5] = ( vs4s9EE[ind2] );
	  value_pi02[6] = ( vs1s9EE[ind2] );
	  value_pi02[7] = ( vs2s9EE[ind2] );
	  value_pi02[8] = ( ESratio[ind2] );
	  value_pi02[9] = ( EtaRing_2 );
	    new_value_pi02[0] = ( G_Sort_2.E() );
            new_value_pi02[1] = ( G_Sort_2.Eta() );
            new_value_pi02[2] = ( G_Sort_2.Phi() );
            new_value_pi02[3] = ( Ncristal_EE[ind2] );
            new_value_pi02[4] = ( vs4s9[ind2] );
            new_value_pi02[5] = ( vs2s9[ind2] );
            new_value_pi02[6] = ( iX2 );
            new_value_pi02[7] = ( iY2 );

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
	    <<" ratio E "<< G_Sort_2.E()/GSort1plus2_EoverCoshEta<<" Pt "<<G_Sort_2.Pt()
	    <<" xtal "<<Ncristal_EE[ind2]<<" vs4s9EE "<<vs4s9EE[ind2]<<" vs1s9EE "<<vs1s9EE[ind2]<<" vs2s9EE "<<vs2s9EE[ind2]
	    <<" ESratio "<<ESratio[ind2]<<" EtaRing_1 "<<EtaRing_2<<endl;

	  // FIXME: should we uncomment these lines as for the barrel?
	  //  if( !Inverted ){ Corr1 = Correct1; Corr2 = Correct2; }
	  //  else           { Corr1 = Correct2; Corr2 = Correct1; }

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
	if (g1pt > g2pt) pi0P4 = Corr1 * G_Sort_1 + Corr2 * G_Sort_2; 
	else             pi0P4 = Corr1 * G_Sort_2 + Corr2 * G_Sort_1;  // when g1pt < g2pt, G_Sort_1 is made with g2, and Corr2 must be applied to it
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
	  if(isMC_ && MC_Asssoc_) Fill_NcrystalUsedGamma_EB(Ncristal_EB_used[0], Ncristal_EB_used[1]);
	  Fill_S4S9Gamma_EB(vs4s9[0], vs4s9[1]);
	  //Fill_NxtalEnergGamma_EB(Nxtal_EnergGamma);  //which is the difference wrt Ncristal_EB_used ?!? 
	  //Fill_NxtalEnergGamma2_EB(Nxtal_EnergGamma2);
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
	  if(isMC_ && MC_Asssoc_) Fill_NcrystalUsedGamma_EE(Ncristal_EE_used[0], Ncristal_EE_used[1]);
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

	  if (fabs(pi0P4_eta)<.1)       { if( pi0P4_nocor_pt < pi0PtCut_low_[subDetId]) continue; }
	  else if (fabs(pi0P4_eta)<1.5) { if( pi0P4_nocor_pt < pi0PtCut_high_[subDetId]) continue; }
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

	  if (fabs(pi0P4_eta)<1.)       { if( nextClu<pi0IsoCut_low_[subDetId] ) continue; }
	  else if (fabs(pi0P4_eta)<1.5) { if( nextClu<pi0IsoCut_high_[subDetId] ) continue; }
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

	  if (fabs(pi0P4_eta)<1.)       { if( hlt_iso > pi0HLTIsoCut_low_[subDetId]  && CutOnHLTIso_ ) continue; }
	  else if (fabs(pi0P4_eta)<1.5) { if( hlt_iso > pi0HLTIsoCut_high_[subDetId] && CutOnHLTIso_ ) continue; }
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
	if(subDetId==EcalEndcap){
	  if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EE[i]; Nxtal_EnergGamma2 = Ncristal_EE[j]; }
	  else                           {  Nxtal_EnergGamma = Ncristal_EE[j]; Nxtal_EnergGamma2 = Ncristal_EE[i]; }
	}
	else{
	  if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EB[i]; Nxtal_EnergGamma2 = Ncristal_EB[j]; }
	  else                           {  Nxtal_EnergGamma = Ncristal_EB[j]; Nxtal_EnergGamma2 = Ncristal_EB[i]; }
	}

	if (subDetId == EcalBarrel) {

	  if( fabs(pi0P4_eta)<1. ) { 
	    if( Nxtal_EnergGamma < nXtal_1_cut_low_[subDetId] ) continue; 
	    if( Nxtal_EnergGamma2 < nXtal_2_cut_low_[subDetId] ) continue;
	  } else if( fabs(pi0P4_eta)<1.5 )  { 
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

	if (isDebug_) cout << "[DEBUG] Fill Optimization Variables..." << endl;

	//Fill Optimization
	if( MakeNtuple4optimization_ && pi0P4_mass > ((Are_pi0_)?0.03:0.2) && pi0P4_mass < ((Are_pi0_)?0.25:1.) ){
	  if( nPi0>NPI0MAX-2 ){ cout<<"nPi0::TOO MANY PI0: ("<<nPi0<<")!!!"<<endl; }
	  else{
	    Op_Pi0recIsEB[nPi0]    = (subDetId==EcalBarrel)? 1:0;
	    Op_IsoPi0_rec[nPi0]    = nextClu;  
	    Op_HLTIsoPi0_rec[nPi0] = hlt_iso;
	    Op_n1CrisPi0_rec[nPi0] = Nxtal_EnergGamma; 
	    Op_n2CrisPi0_rec[nPi0] = Nxtal_EnergGamma2;
	    Op_mPi0_rec[nPi0]      = pi0P4_mass;
	    Op_enG1_rec[nPi0]      = Corr1 * g1->energy();
	    Op_enG2_rec[nPi0]      = Corr2 * g2->energy();
	    Op_etaPi0_rec[nPi0]    = pi0P4_eta;
	    Op_ptPi0_rec[nPi0]     = pi0P4_pt;
	    Op_DeltaRG1G2[nPi0]    = GetDeltaR( g1eta, g2eta, g1phi, g2phi );
	    Op_enG1_nocor[nPi0]    = g1->energy(); // g1P4_nocor.E();
	    Op_enG2_nocor[nPi0]    = g2->energy(); // g2P4_nocor.E();
	    Op_ptPi0_nocor[nPi0]   = pi0P4_nocor_pt;
	    Op_mPi0_nocor[nPi0]    = pi0P4_nocor_mass;
	    Op_Es_e1_1[nPi0]       = (subDetId==EcalBarrel) ? 0. : Es_1[0];
	    Op_Es_e1_2[nPi0]       = (subDetId==EcalBarrel) ? 0. : Es_1[1];
	    Op_Es_e2_1[nPi0]       = (subDetId==EcalBarrel) ? 0. : Es_2[0];
	    Op_Es_e2_2[nPi0]       = (subDetId==EcalBarrel) ? 0. : Es_2[1];
	    Op_S4S9_1[nPi0]        = (subDetId==EcalBarrel) ? vs4s9[0] : vs4s9EE[0];
	    Op_S4S9_2[nPi0]        = (subDetId==EcalBarrel) ? vs4s9[1] : vs4s9EE[1];
	    Op_S2S9_1[nPi0]        = (subDetId==EcalBarrel) ? vs2s9[0] : vs2s9EE[0];
	    Op_S2S9_2[nPi0]        = (subDetId==EcalBarrel) ? vs2s9[1] : vs2s9EE[1];
	    Op_S1S9_1[nPi0]        = (subDetId==EcalBarrel) ? vs1s9[0] : vs1s9EE[0];
	    Op_S1S9_2[nPi0]        = (subDetId==EcalBarrel) ? vs1s9[1] : vs1s9EE[1];
	    Op_Eta_1[nPi0]         = g1eta;
            Op_Eta_2[nPi0]         = g2eta;
	    Op_Phi_1[nPi0]         = g1phi;
	    Op_Phi_2[nPi0]         = g2phi;
	    Op_Time_1[nPi0]        = (subDetId==EcalBarrel) ? vSeedTime[0] : vSeedTimeEE[0];
	    Op_Time_2[nPi0]        = (subDetId==EcalBarrel) ? vSeedTime[1] : vSeedTimeEE[1];
	    if( isMC_ && MC_Asssoc_ ) {
	      Op_Nxtal_1[nPi0]       = (subDetId==EcalBarrel) ? Ncristal_EB_used[0] : Ncristal_EE_used[0];
	      Op_Nxtal_2[nPi0]       = (subDetId==EcalBarrel) ? Ncristal_EB_used[1] : Ncristal_EE_used[1];
	      Op_enG1_true[nPi0]     = Gamma1MC.R();
	      Op_enG2_true[nPi0]     = Gamma2MC.R();
	      Op_DeltaR_1[nPi0]      = reco::deltaR(g1eta,g1phi, Gamma1MC.Eta(),Gamma1MC.Phi());
	      Op_DeltaR_2[nPi0]      = reco::deltaR(g2eta,g2phi, Gamma2MC.Eta(),Gamma2MC.Phi());
	    }

            if( (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ) {

              EBDetId  id_1(g1->seed()); int iEta1 = id_1.ieta(); int iPhi1 = id_1.iphi();
              EBDetId  id_2(g2->seed()); int iEta2 = id_2.ieta(); int iPhi2 = id_2.iphi();

              Op_iEtaiX_1[nPi0]     = iEta1;
              Op_iEtaiX_2[nPi0]     = iEta2;
              Op_iPhiiY_1[nPi0]     = iPhi1;
              Op_iPhiiY_2[nPi0]     = iPhi2;
              Op_iEta_1on5[nPi0]    = iEta1%5;
              Op_iEta_2on5[nPi0]    = iEta2%5;
              Op_iPhi_1on2[nPi0]    = iPhi1%2;
              Op_iPhi_2on2[nPi0]    = iPhi2%2;
              Op_iEta_1on2520[nPi0] =  (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20); //Distance in xtal from module boundaries
              Op_iEta_2on2520[nPi0] =  (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20);
              Op_iPhi_1on20[nPi0]   =  iPhi1%20;
              Op_iPhi_2on20[nPi0]   =  iPhi2%20;
            } else if( (g1->seed().subdetId()==2) && (g2->seed().subdetId()==2) ) {
              
              EEDetId  id_1(g1->seed()); int iX1 = id_1.ix(); int iY1 = id_1.iy();
              EEDetId  id_2(g2->seed()); int iX2 = id_2.ix(); int iY2 = id_2.iy();

              Op_iEtaiX_1[nPi0]        = (iX1 < 50) ? iX1 : 100-iX1;
              Op_iEtaiX_2[nPi0]        = (iX2 < 50) ? iX2 : 100-iX2;
              Op_iPhiiY_1[nPi0]        = (iY1 < 50) ? iY1 : 100-iY1;
              Op_iPhiiY_2[nPi0]        = (iY2 < 50) ? iY2 : 100-iY2;
              Op_iEta_1on5[nPi0]       = 999;
              Op_iEta_2on5[nPi0]       = 999;
              Op_iPhi_1on2[nPi0]       = 999;
              Op_iPhi_2on2[nPi0]       = 999;
              Op_iEta_1on2520[nPi0]    = 999;
              Op_iEta_2on2520[nPi0]    = 999;
              Op_iPhi_1on20[nPi0]      = 999;
              Op_iPhi_2on20[nPi0]      = 999;
              
            } else {
              Op_iEtaiX_1[nPi0]        = -999;
              Op_iEtaiX_2[nPi0]        = -999;
              Op_iPhiiY_1[nPi0]        = -999;
              Op_iPhiiY_2[nPi0]        = -999;
              Op_iEta_1on5[nPi0]       = -999;
              Op_iEta_2on5[nPi0]       = -999;
              Op_iPhi_1on2[nPi0]       = -999;
              Op_iPhi_2on2[nPi0]       = -999;
              Op_iEta_1on2520[nPi0]    = -999;
              Op_iEta_2on2520[nPi0]    = -999;
              Op_iPhi_1on20[nPi0]      = -999;
              Op_iPhi_2on20[nPi0]      = -999;
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
	  if(subDetId!=EcalBarrel) allEpsilon_EEnw->Fill( pi0P4_mass );
	  if(subDetId==EcalBarrel) allEpsilon_EBnw->Fill( pi0P4_mass );
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

  if (isDebug_) cout << "[DEBUG] Filling Tree" << endl; 

  if(MakeNtuple4optimization_){
    //for(unsigned int i=0; i<NL1SEED; i++) Op_L1Seed[i] = L1BitCollection_[i];
    //for(unsigned int i=0; i<NL1SEED; i++) Op_L1Seed[i] = l1flag[i];
    Op_NPi0_rec = nPi0; 
    if(nPi0>0) Tree_Optim->Fill();
  }
 
}


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
  for (int i = 0 ; i != hltCount; ++i) {
    TString hltName_tstr(HLTNames.triggerName(i));
    //std::string hltName_str(HLTNames.triggerName(i));
    //cout<<"hltName_tstr is: "<<hltName_tstr<<" and reg is: "<<s<<endl;
    if ( hltName_tstr.Contains(reg) ){          // If reg contains * ir will say always True. So you ask for ->accept(i) to the first HLTName always.
	//cout<<"hltName_tstr.Contains(reg) give: "<<hltTriggerResultHandle->accept(i)<<endl;
	return hltTriggerResultHandle->accept(i); // False or True depending if it fired.
    }
  }
  return false;
} // HLT isValid

bool FillEpsilonPlot::getTriggerByName( std::string s ) {
  std::map< std::string, int >::iterator currentTrigger;
  currentTrigger = l1TrigNames_.find(s);
  if(currentTrigger != l1TrigNames_.end())
    return l1TrigBit_[currentTrigger->second];
  else 
    std::cout << "Trigger Name not found" << std::endl;
  return false;
}

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
  allEpsilon_EB->Write();
  allEpsilon_EBnw->Write();
  allEpsilon_EE->Write();
  allEpsilon_EEnw->Write();
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
  if( !MakeNtuple4optimization_ &&(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EB_h, "Barrel" ,  regionalCalibration_->getCalibMap()->getNRegionsEB() );
  if( !MakeNtuple4optimization_ && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EE_h, "Endcap" ,  regionalCalibration_->getCalibMap()->getNRegionsEE() );
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
  delete TTree_JoshMva;
#endif
#ifdef MVA_REGRESSIO_EE
  delete TTree_JoshMva_EE;
#endif
}

// ------------ EBPHI LOAD Containment correction  ------------
void FillEpsilonPlot::EBPHI_Cont_Corr_load(std::string FileName )
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
