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
// $Id: FillEpsilonPlot.cc,v 1.5 2013/03/18 22:48:42 lpernie Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include<iostream>
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

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"

#include "CalibCode/FillEpsilonPlot/interface/FillEpsilonPlot.h"
#include "CalibCode/CalibTools/interface/GlobalFunctions.h"
#include "CalibCode/CalibTools/interface/EcalRecHitCompare.h"
//@#include "CalibCode/CalibTools/interface/PreshowerCluster.h"
#include "CalibCode/CalibTools/interface/PreshowerTools.h"
#include "CalibCode/CalibTools/interface/GeometryService.h"
//Geom
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
//ES
#include "FWCore/Framework/interface/ESProducer.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
//HLT
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include <FWCore/Common/interface/TriggerNames.h>
#include <DataFormats/Common/interface/TriggerResults.h>


using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::max;

//JSON
//#include "CalibCode/FillEpsilonPlot/interface/JSON.h"

//MVA Stuff
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif

#ifdef MVA_REGRESSIO
#include "CalibCode/GBRTrain/interface/GBRApply.h"
#include "CalibCode/EgammaObjects/interface/GBRForest.h"
#include "Cintex/Cintex.h"
#include "TLorentzVector.h"
#endif

using namespace TMVA;

//Function
double max_array(double *A, int n);
double max(double x, double y);
int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3);

FillEpsilonPlot::FillEpsilonPlot(const edm::ParameterSet& iConfig)
{
    /// to be moved in parameters.py
    useMassInsteadOfEpsilon_ = 1;

    /// parameters from python
    EBRecHitCollectionTag_  = iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag");
    EERecHitCollectionTag_  = iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag");
    ESRecHitCollectionTag_  = iConfig.getUntrackedParameter<edm::InputTag>("ESRecHitCollectionTag");
    HLTResults_             = iConfig.getUntrackedParameter<bool>("HLTResults",false);
    l1TriggerTag_           = iConfig.getUntrackedParameter<edm::InputTag>("L1TriggerTag");
    triggerTag_             = iConfig.getUntrackedParameter<edm::InputTag>("triggerTag",edm::InputTag("TriggerResults"));
    outfilename_            = iConfig.getUntrackedParameter<std::string>("OutputFile");
    ebContainmentCorrections_       = iConfig.getUntrackedParameter<std::string>("EBContainmentCorrections");
    MVAEBContainmentCorrections_01_ = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_01");
    MVAEBContainmentCorrections_02_ = iConfig.getUntrackedParameter<std::string>("MVAEBContainmentCorrections_02");
    MVAEEContainmentCorrections_01_ = iConfig.getUntrackedParameter<std::string>("MVAEEContainmentCorrections_01");
    MVAEEContainmentCorrections_02_ = iConfig.getUntrackedParameter<std::string>("MVAEEContainmentCorrections_02");
    Endc_x_y_                       = iConfig.getUntrackedParameter<std::string>("Endc_x_y");
    ebPHIContainmentCorrections_    = iConfig.getUntrackedParameter<std::string>("EBPHIContainmentCorrections");
    eeContainmentCorrections_       = iConfig.getUntrackedParameter<std::string>("EEContainmentCorrections");
    useEBContainmentCorrections_    = iConfig.getUntrackedParameter<bool>("useEBContainmentCorrections");
    useEEContainmentCorrections_    = iConfig.getUntrackedParameter<bool>("useEEContainmentCorrections");
    ContCorr_EB_                    = iConfig.getUntrackedParameter<std::string>("ContCorr_EB");
    externalGeometry_               = iConfig.getUntrackedParameter<std::string>("ExternalGeometry");
    currentIteration_               = iConfig.getUntrackedParameter<int>("CurrentIteration");
    outputDir_                      = iConfig.getUntrackedParameter<std::string>("OutputDir");
    jsonFile_                       = iConfig.getUntrackedParameter<std::string>("json_file");
    calibMapPath_                   = iConfig.getUntrackedParameter<std::string>("calibMapPath");
    Barrel_orEndcap_                = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");
    pi0PtCut_[EcalBarrel]  = iConfig.getUntrackedParameter<double>("Pi0PtCutEB");
    pi0PtCut_[EcalEndcap]  = iConfig.getUntrackedParameter<double>("Pi0PtCutEE");
    gPtCut_[EcalBarrel]  = iConfig.getUntrackedParameter<double>("gPtCutEB");
    gPtCut_[EcalEndcap]  = iConfig.getUntrackedParameter<double>("gPtCutEE");
    pi0IsoCut_[EcalBarrel] = iConfig.getUntrackedParameter<double>("Pi0IsoCutEB");
    pi0IsoCut_[EcalEndcap] = iConfig.getUntrackedParameter<double>("Pi0IsoCutEE");
    nXtal_1_cut_[EcalEndcap] = iConfig.getUntrackedParameter<double>("nXtal_1_EE");
    nXtal_2_cut_[EcalEndcap] = iConfig.getUntrackedParameter<double>("nXtal_2_EE");
    nXtal_1_cut_[EcalBarrel] = iConfig.getUntrackedParameter<double>("nXtal_1_EB");
    nXtal_2_cut_[EcalBarrel] = iConfig.getUntrackedParameter<double>("nXtal_2_EB");
    S4S9_cut_[EcalBarrel] = iConfig.getUntrackedParameter<double>("S4S9_EB");
    S4S9_cut_[EcalEndcap] = iConfig.getUntrackedParameter<double>("S4S9_EE");
    SystOrNot_ = iConfig.getUntrackedParameter<double>("SystOrNot",0);

cout<<"Cus used: EB)"<<endl;
cout<<"Pt(pi0): "<<pi0PtCut_[EcalBarrel]<<", Pt(Clus): "<<gPtCut_[EcalBarrel]<<", Iso: "<<pi0IsoCut_[EcalBarrel]<<", Nxtal_1: "<<nXtal_1_cut_[EcalBarrel]<<", Nxtal_2: "<<nXtal_2_cut_[EcalBarrel]<<", S4S9: "<<S4S9_cut_[EcalBarrel]<<endl;
cout<<"Cus used: EE)"<<endl;
cout<<"Pt(pi0): "<<pi0PtCut_[EcalEndcap]<<", Pt(Clus): "<<gPtCut_[EcalEndcap]<<", Iso: "<<pi0IsoCut_[EcalEndcap]<<", Nxtal_1: "<<nXtal_1_cut_[EcalEndcap]<<", Nxtal_2: "<<nXtal_2_cut_[EcalEndcap]<<", S4S9: "<<S4S9_cut_[EcalEndcap]<<endl;
cout<<"The StatError option choose is: "<<SystOrNot_<<" [0= No error stat computation, 1 = yes only even events, 2 = yes only odd events]"<<endl;

    alcaL1TrigNames_ = iConfig.getUntrackedParameter< std::vector<std::string> >("AlcaL1TrigNames"); 

    useOnlyEEClusterMatchedWithES_ = iConfig.getUntrackedParameter<bool>("useOnlyEEClusterMatchedWithES"); 

    /// shower shape parameters
    PCparams_.param_LogWeighted_ = true;
    //PCparams_.param_T0_barl_     = 5.7;
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
    externalGeometryFile_ = TFile::Open(externalGeometry_.c_str());
    if(!externalGeometryFile_) cms::Exception("ExtGeom") << "External Geometry file (" << externalGeometry_ << ") not found" << endl;
    geom_ = ECALGeometry::getGeometry(externalGeometryFile_);
    GeometryService::setGeometryName(externalGeometry_);
    GeometryService::setGeometryPtr(geom_);
    
    /// containment corrections
#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)
    if(useEEContainmentCorrections_)
        containmentCorrections_.loadContainmentPointCorrectionsEE(eeContainmentCorrections_.c_str());
    if(useEBContainmentCorrections_){
        containmentCorrections_.loadContainmentCorrectionsEB(ebContainmentCorrections_.c_str());
        EBPHI_Cont_Corr_load( ebPHIContainmentCorrections_.c_str() );
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
    if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
    else if(currentIteration_ > 0)
    {
       char fileName[200];
       //sprintf(fileName,"%s/iter_%d/calibMap.root", outputDir_.c_str(), currentIteration_-1);
       cout << "FillEpsilonPlot:: loading calibraion map at " << calibMapPath_ << endl;
       sprintf(fileName,"%s", calibMapPath_.c_str());
       regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
    }

    /// epsilon histograms
    if(useMassInsteadOfEpsilon_)
    {
      if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","#pi^{0} Mass distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
      if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","#pi^{0} Mass distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
    }
    else
    {
      if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EB_h = initializeEpsilonHistograms("epsilon_EB_iR_","Epsilon distribution EB - iR ", regionalCalibration_->getCalibMap()->getNRegionsEB() );
      if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) )  epsilon_EE_h = initializeEpsilonHistograms("epsilon_EE_iR_","Epsilon distribution EE - iR ", regionalCalibration_->getCalibMap()->getNRegionsEE() );
    }

    allEpsilon_EB = new TH1F("allEpsilon_EB", "allEpsilon_EB",240,0.,0.5);
    allEpsilon_EE = new TH1F("allEpsilon_EE", "allEpsilon_EE",240,0.,0.5);
    entries_EEp   = new TH2F("entries_EEp","entries_EEp",101,-0.5,100.5,101,-0.5,100.5);
    entries_EEm   = new TH2F("entries_EEm","entries_EEm",101,-0.5,100.5,101,-0.5,100.5);

    pi0MassVsIetaEB = new TH2F("pi0MassVsIetaEB","#pi^{0} mass vs i#eta",85,0.5,85.5,120,0.,0.3);
    pi0MassVsIetaEB->GetXaxis()->SetTitle("i#eta");
    pi0MassVsIetaEB->GetYaxis()->SetTitle("#pi^{0} mass");
    pi0MassVsETEB = new TH2F("pi0MassVsETEB", "#pi^{0} mass vs E_{T}(pi^{0})",120,0.,20.,120,0.,0.3);
    pi0MassVsETEB->GetXaxis()->SetTitle("E_{T}(pi^{0})");
    pi0MassVsETEB->GetYaxis()->SetTitle("#pi^{0} mass");

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

    /// trigger histo
    triggerComposition = new TH1F("triggerComposition", "Trigger Composition",128,-0.5,127.5);
    areLabelsSet_ = false;

    //noCorrections_ = true;

#ifdef MVA_REGRESSIO
    EBweight_file_pi01 = TFile::Open( MVAEBContainmentCorrections_01_.c_str() );
    EBweight_file_pi02 = TFile::Open( MVAEBContainmentCorrections_02_.c_str() );
    forest_EB_pi01 = (GBRForest *)EBweight_file_pi01->Get("Correction");    
    forest_EB_pi02 = (GBRForest *)EBweight_file_pi02->Get("Correction");    
#endif
#ifdef MVA_REGRESSIO_EE
    EEweight_file_pi01 = TFile::Open( MVAEEContainmentCorrections_01_.c_str() );
    EEweight_file_pi02 = TFile::Open( MVAEEContainmentCorrections_02_.c_str() );
    forest_EE_pi01 = (GBRForest *)EEweight_file_pi01->Get("Correction");
    forest_EE_pi02 = (GBRForest *)EEweight_file_pi02->Get("Correction");
#endif
    //JSON
    //ev_TOT=0;
    //ev_JSON=0;
    //myjson = new JSON(jsonFile_.c_str());
    //hev_TOT  = new TH1F("hev_TOT", "Events before Json file selection",1, 0.5, 1.5);
    //hev_JSON = new TH1F("hev_JSON", "Events after Json file selection",1, 0.5, 1.5);
}

FillEpsilonPlot::~FillEpsilonPlot()
{
    externalGeometryFile_->Close();
    outfile_->Write();
    outfile_->Close();

    if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) deleteEpsilonPlot(epsilon_EB_h, regionalCalibration_->getCalibMap()->getNRegionsEB() );
    if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) deleteEpsilonPlot(epsilon_EE_h, regionalCalibration_->getCalibMap()->getNRegionsEE() );
    delete allEpsilon_EB;
    delete allEpsilon_EE;
    delete entries_EEp;
    delete entries_EEm;
    delete pi0MassVsIetaEB;
    delete pi0MassVsETEB;
    delete triggerComposition;

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

#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)
    delete EBPHI_ConCorr_p;
    delete EBPHI_ConCorr_m;
#endif
    //JSON
    //delete myjson;
    //delete hev_TOT;
    //delete hev_JSON;
#ifdef MVA_REGRESSIO
    delete forest_EB_pi01;
    delete forest_EB_pi02;
#endif
#ifdef MVA_REGRESSIO_EE
    delete forest_EE_pi01;
    delete forest_EE_pi02;
#endif
}


//
// member functions
//

// ------------ method called for each event  ------------
void
FillEpsilonPlot::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//cout<<"----------------------------Start Event--------------------------------------------------"<<endl;
   using namespace edm;
    //For Syst error SystOrNot_=1 or 2, for normal calib is 0
    if(SystOrNot_==1. && int(iEvent.id().event())%2!=0 ) return;
    else if(SystOrNot_==2. && int(iEvent.id().event())%2==0 ) return;

    //ES
    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    geometry = geoHandle.product();
    estopology_ = new EcalPreshowerTopology(geoHandle);
    esGeometry_ = (dynamic_cast<const EcalPreshowerGeometry*>( (CaloSubdetectorGeometry*) geometry->getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ));

   /// ------ EB CLUSTERS -------
   //if( !getTriggerResult(iEvent, iSetup)) return;

   //JSON FILE
   //ev_TOT++;
   //if ( !myjson->isGoodLS( iEvent.id().run() , iEvent.luminosityBlock() ) ) return;
   //ev_JSON++;
   //Vectors
   std::vector< CaloCluster > ebclusters; // contains the output clusters
   ebclusters.clear();
   vs4s9.clear();
   vs2s9.clear();
   vs2s9.clear();
#ifdef MVA_REGRESSIO_EE
   vs4s9EE.clear(); vs2s9EE.clear(); vs2s9EE.clear(); ESratio.clear();
#endif
   std::vector< CaloCluster > eseeclusters; eseeclusters.clear();
   std::vector< CaloCluster > eseeclusters_tot; eseeclusters_tot.clear();
   Ncristal_EB.clear(); Ncristal_EE.clear();

   bool EB_HLT=true, EE_HLT=true;
   if( HLTResults_ ){
       EB_HLT = GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
       EE_HLT = GetHLTResults(iEvent, "AlCa_EcalPi0EEonly.*");
   }

   if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) && EB_HLT ) fillEBClusters(ebclusters, iEvent);
   if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) && EE_HLT ) fillEEClusters(eseeclusters, eseeclusters_tot, iEvent);

   if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(ebclusters,EcalBarrel);
   if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(eseeclusters_tot,EcalEndcap);

   delete estopology_;
}


/*===============================================================*/
void FillEpsilonPlot::fillEBClusters(std::vector< CaloCluster > & ebclusters, const edm::Event& iEvent)
/*===============================================================*/
{
   edm::Handle< EBRecHitCollection > ebHandle;
   bool goodhandle = iEvent.getByLabel ( EBRecHitCollectionTag_, ebHandle);
   if(!goodhandle) cout << "problem in eventCont->getByLabel ( ebLabel, ebHandle)" << endl;

   std::vector<EcalRecHit> ebseeds;

   typedef std::set<EBDetId> XtalInUse;
   XtalInUse isUsed; // map of which xtals have been used

   int dc = 0;

   // sort by energy and find the seeds
   for(EBRecHitCollection::const_iterator itb= ebHandle->begin(); itb != ebHandle->end(); ++itb, ++dc) 
   {
      if(itb->energy() > 0.200)  ebseeds.push_back( *itb );
   }
   /// debug
   //cout << "ebseeds.size() = " << ebseeds.size() << " / " << dc << endl;

#ifdef SELECTION_TREE
   Fill_NSeeds_EB(ebseeds.size());
#endif

   sort(ebseeds.begin(), ebseeds.end(), ecalRecHitLess());
        int seed_c = 0;
      // loop over seeds and make clusters
      for (std::vector<EcalRecHit>::iterator itseed=ebseeds.begin(); itseed!=ebseeds.end(); itseed++, seed_c++) 
      {
         EBDetId seed_id( itseed->id() );

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
         for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++) 
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
           cout << "skipping cluster with negative energy " << simple_energy << endl; 
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

         // EB only for the time being
         //float etaSeed = geom->getPosition(seed_id).eta() ;
         //double T0(0.);
         //if( fabs(ctreta) < 1.479 )  T0 = param_T0_barl_;
         //else T0 = param_T0_endcES_;

         // Calculate shower depth
         float T0 = PCparams_.param_T0_barl_;
         float maxDepth = PCparams_.param_X0_ * ( T0 + log( posTotalEnergy ) );
         float maxToFront = geom_->getPosition(seed_id).mag(); // to front face
         double EnergyCristals[9] = {0.};

         // loop over xtals and compute energy and position
         for(unsigned int j=0; j<RecHitsInWindow.size();j++)
         {
            EBDetId det(RecHitsInWindow[j]->id());
            //cout << "id: " << det << " calib coeff: " << calibMap->coeff(RecHitsInWindow[j]->id()) << endl;

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
                                // NOTA BENE: sto usando le frazioni per salvare energia rechit
               isUsed.insert(RecHitsInWindow[j]->id());

            }

            // compute position
            if(en>0.) 
            {
               float weight = std::max( float(0.), PCparams_.param_W0_ + log(en/posTotalEnergy) );

               float depth = maxDepth + maxToFront - geom_->getPosition(det).mag() ;
               GlobalPoint posThis = geom_->getPosition(det,depth);

               xclu += weight*posThis.x(); 
               yclu += weight*posThis.y(); 
               zclu += weight*posThis.z(); 
               total_weight += weight;
            }

         } // loop over 3x3 rechits

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

          if(s4s9<S4S9_cut_[EcalBarrel]) continue;

#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)
         if(useEBContainmentCorrections_) 
         {
             e3x3 *=  containmentCorrections_.getContainmentCorrectionsEB(e3x3, seed_id.ieta() );
             e3x3 *=  EBPHI_Cont_Corr(e3x3/cosh(clusPos.eta()), seed_id.iphi()%20, seed_id.ieta() );
         }
#endif

         // compute pt of gamma and cut
         float ptClus = e3x3/cosh(clusPos.eta());
         
         if(ptClus<gPtCut_[EcalBarrel]) continue;

         // make calo clusters
#ifdef MVA_REGRESSIO
         vs4s9.push_back( s4s9 );
         vs1s9.push_back( itseed->energy()/e3x3 );
         double maxEne = max_array( EnergyCristals, 9 );
         for(int i=0; i<9; i++){ if( EnergyCristals[i]==maxEne ) EnergyCristals[i]=0.; }
         double maxEne2 = max_array( EnergyCristals, 9);
         vs2s9.push_back( (maxEne+maxEne2)/e3x3 );
#endif
         Ncristal_EB.push_back(RecHitsInWindow.size() );
         ebclusters.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_BARREL),
                                          enFracs, CaloCluster::undefined, seed_id ) );
      } //loop over seeds to make EB clusters

}

/*===============================================================*/
void FillEpsilonPlot::fillEEClusters(std::vector< CaloCluster > & eseeclusters, std::vector< CaloCluster > & eseeclusters_tot, const edm::Event& iEvent)
/*===============================================================*/
{
    edm::Handle< EERecHitCollection > eeHandle;                     
    bool goodhandle = iEvent.getByLabel ( EERecHitCollectionTag_, eeHandle);                     
    if(!goodhandle) throw cms::Exception("Handle") << "Bad EE Handle\n";

    edm::Handle< ESRecHitCollection > esHandle;                     
    goodhandle = iEvent.getByLabel ( ESRecHitCollectionTag_, esHandle);  
    if(!goodhandle) throw cms::Exception("Handle") << "Bad ES Handle\n";

    PreshowerTools esClusteringAlgo(geometry, estopology_, esHandle);

    std::vector<EcalRecHit> eeseeds;

    vector <double> eeclusterS4S9; eeclusterS4S9.clear();
#ifdef MVA_REGRESSIO_EE
    vector <double> eeclusterS1S9; eeclusterS1S9.clear();
    vector <double> eeclusterS2S9; eeclusterS2S9.clear();
#endif
    
    std::vector< CaloCluster > eeclusters; // contains the output eeclusters
    eeclusters.clear();

    int dc = 0;

    // sort by energy and find the eeseeds
    for(EERecHitCollection::const_iterator ite= eeHandle->begin(); ite != eeHandle->end(); ++ite, ++dc) {
       int idXtal= ite->id();
       float posTotalEnergy(ite->energy());
       float T0 = PCparams_.param_T0_barl_;
       float maxDepth = PCparams_.param_X0_ * ( T0 + log( posTotalEnergy ) );
       float maxToFront = geom_->getPosition(idXtal).mag(); // to front face
       float depth = maxDepth + maxToFront - geom_->getPosition(idXtal).mag() ;
       GlobalPoint posThis = geom_->getPosition(idXtal,depth);
       if(ite->energy()/cosh(posThis.eta()) > 0.5)      eeseeds.push_back( *ite );
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
       for (std::vector<DetId>::const_iterator det=clus_v.begin(); det!=clus_v.end(); det++) 
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
          cout << "skipping cluster with negative energy " << simple_energy << endl; 
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
       float maxToFront = geom_->getPosition(eeseed_id).mag(); // to front face
#ifdef MVA_REGRESSIO_EE
       double EnergyCristals[9] = {0.};
#endif
       // loop over xtals and compute energy and position
       for(unsigned int j=0; j<RecHitsInWindow.size();j++)
       {
          EEDetId det(RecHitsInWindow[j]->id());
          //cout << "id: " << det << " calib coeff: " << calibMap->coeff(RecHitsInWindow[j]->id()) << endl;

          int ix = det.ix();
          int iy = det.iy();
 
          // use calibration coeff for energy and position
          float en = RecHitsInWindow[j]->energy() * regionalCalibration_->getCalibMap()->coeff(RecHitsInWindow[j]->id());
          int dx = seed_ix-ix;
          int dy = seed_iy-iy;
#ifdef MVA_REGRESSIO_EE
          EnergyCristals[j] = en;
#endif
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
  
            float depth = maxDepth + maxToFront - geom_->getPosition(det).mag() ;
            GlobalPoint posThis = geom_->getPosition(det,depth);
   
            xclu += weight*posThis.x(); 
            yclu += weight*posThis.y(); 
            zclu += weight*posThis.z(); 
            total_weight += weight;
         }
      } // loop over 3x3 eerechits
      
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

      if(s4s9<S4S9_cut_[EcalEndcap]) continue; //original cut 0.95

      float ptClus = e3x3/cosh(clusPos.eta());

      if(ptClus<gPtCut_[EcalEndcap]) continue;  //original cut 0.6

      // make calo clusters

      //Filling the usedxtal vector
      for(unsigned int j=0; j<RecHitsInWindow.size();j++){
         EEXisUsed [RecHitsInWindow[j]->id()] = true;
      }

      Ncristal_EE.push_back( RecHitsInWindow.size() );
      eeclusters.push_back( CaloCluster( e3x3, clusPos, CaloID(CaloID::DET_ECAL_ENDCAP),
                                         enFracs, CaloCluster::undefined, eeseed_id ) );

      eeclusterS4S9.push_back(s4s9);
#ifdef MVA_REGRESSIO_EE
      eeclusterS1S9.push_back(eeitseed->energy()/e3x3);
      double maxEne = max_array( EnergyCristals, 9 );
      for(int i=0; i<9; i++){ if( EnergyCristals[i]==maxEne ) EnergyCristals[i]=0.; }
      double maxEne2 = max_array( EnergyCristals, 9);
      eeclusterS2S9.push_back( (maxEne+maxEne2)/e3x3 );
#endif

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

          PreshowerCluster preshowerclusterp1 = esClusteringAlgo.makeOnePreshowerCluster( PreshowerTools::clusterwindowsize_, &tmp1_conversion);
          PreshowerCluster preshowerclusterp2 = esClusteringAlgo.makeOnePreshowerCluster( PreshowerTools::clusterwindowsize_, &tmp2_conversion);


          double e1 = preshowerclusterp1.energy();
          double e2 = preshowerclusterp2.energy();
          // GeV to #MIPs
          e1 = e1 / PreshowerTools::mip_;
          e2 = e2 / PreshowerTools::mip_;
          double tempenergy = eeclus_iter->energy();

          //if(e1+e2 > 1.0e-10) 
          if(e1 > 2.0 || e2 > 2.0) /// cut @ 2 MIPs as suggested by Ming @ DPG/EGamma Joint Meeting 19.03.2012 
          {
             double deltaE = PreshowerTools::gamma_*(PreshowerTools::calib_planeX_*e1 + PreshowerTools::calib_planeY_*e2);

             tempenergy = deltaE + eeclus_iter->energy();
#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)
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
#ifdef MVA_REGRESSIO_EE
                  vs4s9EE.push_back( eeclusterS4S9[ind] );
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
#ifdef MVA_REGRESSIO_EE
    vs4s9EE.push_back( eeclusterS4S9[ind] );
    vs1s9EE.push_back( eeclusterS1S9[ind] );
    vs2s9EE.push_back( eeclusterS2S9[ind] );
    ESratio.push_back( (-1998.)/eeclus_iter->energy() );
#endif
}
    }//end of the matching loop

    Ncristal_EE.clear();
    Ncristal_EE = Nxtal_tot; 
    Nxtal.clear();
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

        //@@ h[jR] = new TH1F(name_c, title_c,120,-0.5,0.5);

        if(useMassInsteadOfEpsilon_)
        {
            h[jR] = new TH1F(name_c, title_c, 120,0.,0.5);
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



void FillEpsilonPlot::computeEpsilon(std::vector< CaloCluster > & clusters, int subDetId ) 
{

    if(subDetId!=EcalBarrel && subDetId != EcalEndcap) 
        throw cms::Exception("FillEpsilonPlot::computeEpsilon") << "Subdetector Id not recognized\n";

    int myNPi0 = 0;

    // loop over clusters to make Pi0
    size_t i=0;
    for(std::vector<CaloCluster>::const_iterator g1  = clusters.begin(); g1 != clusters.end(); ++g1, ++i) 
    {
      size_t j=i+1;
      for(std::vector<CaloCluster>::const_iterator g2 = g1+1; g2 != clusters.end(); ++g2, ++j ) {
          
         float Corr1 = 1., Corr2 = 1.;
#if !defined(NEW_CONTCORR) && defined(MVA_REGRESSIO)
     if( subDetId==EcalBarrel && (g1->seed().subdetId()==1) && (g2->seed().subdetId()==1) ){
   
       TLorentzVector G_Sort_1, G_Sort_2;
       int ind1 = i, ind2 = j;
       EBDetId  id_1(g1->seed()); int iEta1 = id_1.ieta(); int iPhi1 = id_1.iphi();
       EBDetId  id_2(g2->seed()); int iEta2 = id_2.ieta(); int iPhi2 = id_2.iphi();
       bool Inverted=false;
   
       if( g1->energy()/cosh(g1->eta()) > g2->energy()/cosh(g2->eta()) ){
         G_Sort_1.SetPtEtaPhiE( g1->energy()/cosh(g1->eta()) ,g1->eta(),g1->phi(),g1->energy() );
         G_Sort_2.SetPtEtaPhiE( g2->energy()/cosh(g2->eta()) ,g2->eta(),g2->phi(),g2->energy() );
       }
       else{
         G_Sort_1.SetPtEtaPhiE( g2->energy()/cosh(g2->eta()) ,g2->eta(),g2->phi(),g2->energy() );
         G_Sort_2.SetPtEtaPhiE( g1->energy()/cosh(g1->eta()) ,g1->eta(),g1->phi(),g1->energy() );
         iEta1=id_2.ieta(); iEta2 = id_1.ieta();
         iPhi1=id_2.iphi(); iPhi2 = id_1.iphi();
         ind1=j; ind2=i;
         Inverted=true;
       }

       float value_pi01[14];
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
       float Correct1 = forest_EB_pi01->GetResponse(value_pi01);
       //cout<<"Correction1: "<<Correct1<<" iEta: "<<iEta1<<" iPhi "<<iPhi1<<" E "<<G_Sort_1.E()<<endl;
       float value_pi02[14];
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
       float Correct2 = forest_EB_pi02->GetResponse(value_pi02);
       //cout<<"Correction2: "<<Correct2<<" iEta: "<<iEta2<<" iPhi "<<iPhi2<<" E "<<G_Sort_2.E()<<endl;

         if( !Inverted ){ Corr1 = Correct1; Corr2 = Correct2; }
         else           { Corr1 = Correct2; Corr2 = Correct1; }

         //Correction1_mva = Correct1; Correction2_mva = Correct2;
         //iEta1_mva = iEta1; iEta2_mva = iEta2; iPhi1_mva = iPhi1; iPhi2_mva = iPhi2; Pt1_mva = G_Sort_1.Pt(); Pt2_mva = G_Sort_2.Pt();

        // TLorentzVector mvag1P4; mvag1P4.SetPtEtaPhiE( Correct1*G_Sort_1.E()/cosh(G_Sort_1.Eta()), G_Sort_1.Eta(), G_Sort_1.Phi(), Correct1*G_Sort_1.E() );
        // TLorentzVector mvag2P4; mvag2P4.SetPtEtaPhiE( Correct2*G_Sort_2.E()/cosh(G_Sort_2.Eta()), G_Sort_2.Eta(), G_Sort_2.Phi(), Correct2*G_Sort_2.E() );

         //TLorentzVector mvaOrg1P4; mvaOrg1P4.SetPtEtaPhiE( G_Sort_1.E()/cosh(G_Sort_1.Eta()), G_Sort_1.Eta(), G_Sort_1.Phi(), G_Sort_1.E() );
         //TLorentzVector mvaOrg2P4; mvaOrg2P4.SetPtEtaPhiE( G_Sort_2.E()/cosh(G_Sort_2.Eta()), G_Sort_2.Eta(), G_Sort_2.Phi(), G_Sort_2.E() );

         //Mass_mva = (mvag1P4 + mvag2P4).M();
         //MassOr_mva = (mvaOrg1P4 + mvaOrg2P4).M();
         //TTree_JoshMva->Fill();   

     }
#endif
#if !defined(NEW_CONTCORR) && defined(MVA_REGRESSIO_EE)
     if( subDetId==EcalEndcap && (g1->seed().subdetId()==2) && (g2->seed().subdetId()==2) ){
  
       TLorentzVector G_Sort_1, G_Sort_2;
       int ind1 = i, ind2 = j;
       EEDetId  id_1(g1->seed()); int iX1 = id_1.ix(); int iY1 = id_1.iy(); int iZ1 = id_1.zside();
       EEDetId  id_2(g2->seed()); int iX2 = id_2.ix(); int iY2 = id_2.iy(); int iZ2 = id_2.zside();
       bool Inverted=false;

       if( g1->energy()/cosh(g1->eta()) > g2->energy()/cosh(g2->eta()) ){
         G_Sort_1.SetPtEtaPhiE( g1->energy()/cosh(g1->eta()) ,g1->eta(),g1->phi(),g1->energy() );
         G_Sort_2.SetPtEtaPhiE( g2->energy()/cosh(g2->eta()) ,g2->eta(),g2->phi(),g2->energy() );
       }
       else{
         G_Sort_1.SetPtEtaPhiE( g2->energy()/cosh(g2->eta()) ,g2->eta(),g2->phi(),g2->energy() );
         G_Sort_2.SetPtEtaPhiE( g1->energy()/cosh(g1->eta()) ,g1->eta(),g1->phi(),g1->energy() );
         iX1=id_2.ix(); iX2 = id_1.ix();
         iY1=id_2.iy(); iY2 = id_1.iy();
         iZ1=id_2.zside(); iZ2 = id_1.zside();
         ind1=j; ind2=i;
         Inverted=true;
       }
       int EtaRing_1=GetRing( iX1, iY1, VectRing, false), EtaRing_2=GetRing( iX2, iY2, VectRing, false);
       float value_pi01[10];
       value_pi01[0] = ( (G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta()) );
       value_pi01[1] = ( G_Sort_1.E()/((G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta())) );
       value_pi01[2] = ( G_Sort_1.Pt() );
       value_pi01[3] = ( Ncristal_EE[ind1] );
       value_pi01[4] = ( Ncristal_EE[ind2] );
       value_pi01[5] = ( vs4s9EE[ind1] );
       value_pi01[6] = ( vs1s9EE[ind1] );
       value_pi01[7] = ( vs2s9EE[ind1] );
       value_pi01[8] = ( ESratio[ind1] );
       value_pi01[9] = ( EtaRing_1 );
       float Correct1 = forest_EE_pi01->GetResponse(value_pi01);
       cout<<"Correction1: "<<Correct1<<" iX: "<<iX1<<" iY "<<iY1<<" Epi0 "<<(G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta())
       <<" ratio E "<< G_Sort_1.E()/((G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta()))<<" Pt "<<G_Sort_1.Pt()
       <<" xtal "<<Ncristal_EE[ind1]<<" vs4s9EE "<<vs4s9EE[ind1]<<" vs1s9EE "<<vs1s9EE[ind1]<<" vs2s9EE "<<vs2s9EE[ind1]
       <<" ESratio "<<ESratio[ind1]<<" EtaRing_1 "<<EtaRing_1<<endl;
 
       float value_pi02[10];
       value_pi02[0] = ( (G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta()) );
       value_pi02[1] = ( G_Sort_2.E()/((G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta())) );
       value_pi02[2] = ( G_Sort_2.Pt() );
       value_pi02[3] = ( Ncristal_EE[ind1] );
       value_pi02[4] = ( Ncristal_EE[ind2] );
       value_pi02[5] = ( vs4s9EE[ind2] );
       value_pi02[6] = ( vs1s9EE[ind2] );
       value_pi02[7] = ( vs2s9EE[ind2] );
       value_pi02[8] = ( ESratio[ind2] );
       value_pi02[9] = ( EtaRing_2 );
       float Correct2 = forest_EE_pi02->GetResponse(value_pi02);
       cout<<"Correction2: "<<Correct2<<" iX: "<<iX2<<" iY "<<iY2<<" Epi0 "<<(G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta())
       <<" ratio E "<< G_Sort_2.E()/((G_Sort_1+G_Sort_2).E()/cosh((G_Sort_1+G_Sort_2).Eta()))<<" Pt "<<G_Sort_2.Pt()
       <<" xtal "<<Ncristal_EE[ind2]<<" vs4s9EE "<<vs4s9EE[ind2]<<" vs1s9EE "<<vs1s9EE[ind2]<<" vs2s9EE "<<vs2s9EE[ind2]
       <<" ESratio "<<ESratio[ind2]<<" EtaRing_1 "<<EtaRing_2<<endl;

     //  if( !Inverted ){ Corr1 = Correct1; Corr2 = Correct2; }
     //  else           { Corr1 = Correct2; Corr2 = Correct1; }

       Correction1EE_mva = Correct1; Correction2EE_mva = Correct2;
       iX1_mva = iX1; iX2_mva = iX2; iY1_mva = iY1; iY1_mva = iY2; Pt1EE_mva = G_Sort_1.Pt(); Pt2EE_mva = G_Sort_2.Pt();
       EtaRing1_mva = EtaRing_1; EtaRing2_mva = EtaRing_2;

       TLorentzVector mvag1P4; mvag1P4.SetPtEtaPhiE( Correct1*G_Sort_1.E()/cosh(G_Sort_1.Eta()), G_Sort_1.Eta(), G_Sort_1.Phi(), Correct1*G_Sort_1.E() );
       TLorentzVector mvag2P4; mvag2P4.SetPtEtaPhiE( Correct2*G_Sort_2.E()/cosh(G_Sort_2.Eta()), G_Sort_2.Eta(), G_Sort_2.Phi(), Correct2*G_Sort_2.E() );

       TLorentzVector mvaOrg1P4; mvaOrg1P4.SetPtEtaPhiE( G_Sort_1.E()/cosh(G_Sort_1.Eta()), G_Sort_1.Eta(), G_Sort_1.Phi(), G_Sort_1.E() );
       TLorentzVector mvaOrg2P4; mvaOrg2P4.SetPtEtaPhiE( G_Sort_2.E()/cosh(G_Sort_2.Eta()), G_Sort_2.Eta(), G_Sort_2.Phi(), G_Sort_2.E() );

       MassEE_mva = (mvag1P4 + mvag2P4).M();
       MassEEOr_mva = (mvaOrg1P4 + mvaOrg2P4).M();
       TTree_JoshMva_EE->Fill();   
     }
#endif

           math::PtEtaPhiMLorentzVector g1P4( (Corr1*g1->energy())/cosh(g1->eta()), g1->eta(), g1->phi(), 0. );
           math::PtEtaPhiMLorentzVector g2P4( (Corr2*g2->energy())/cosh(g2->eta()), g2->eta(), g2->phi(), 0. );
    
           math::PtEtaPhiMLorentzVector pi0P4 = g1P4 + g2P4;

          //In case ES give same posizion for different clusters
          if( pi0P4.mass()<0.002 ) continue;
#ifdef SELECTION_TREE
          if( subDetId == EcalBarrel ){ 
          Fill_PtPi0_EB( pi0P4.Pt() );
          Fill_mpi0_EB( pi0P4.mass() );
          Fill_etapi0_EB( pi0P4.eta() );
          Fill_phipi0_EB( pi0P4.phi() );
          Fill_Epsilon_EB( 0.5 * ( pow(pi0P4.mass()/PI0MASS,2)  - 1. ) );
          Pi0Info_EB->Fill();
          }
          if( subDetId == EcalEndcap ){
          Fill_PtPi0_EE( pi0P4.Pt() );
          Fill_mpi0_EE( pi0P4.mass() );
          Fill_etapi0_EE( pi0P4.eta() );
          Fill_phipi0_EE( pi0P4.phi() );
          Fill_Epsilon_EE( 0.5 * ( pow(pi0P4.mass()/PI0MASS,2)  - 1. ) );
          Pi0Info_EE->Fill();
          }
#endif

          if(pi0P4.Pt() < pi0PtCut_[subDetId]) continue;
          if( g1P4.eta() == g2P4.eta() && g1P4.phi() == g2P4.phi() ) continue;

          float nextClu = 999., Drtmp = 999.;
          for(size_t ind=0; ind<clusters.size(); ++ind){
            const CaloCluster* Gtmp = &(clusters[ind]);
            double deltaR1 = GetDeltaR(Gtmp->eta(),g1P4.eta(),Gtmp->phi(),g1P4.phi());
            double deltaR2 = GetDeltaR(Gtmp->eta(),g2P4.eta(),Gtmp->phi(),g2P4.phi());
            if( ind!=i && ind!=j && (deltaR1<Drtmp || deltaR2<Drtmp ) ){
                nextClu = min(deltaR1,deltaR2);
                Drtmp = nextClu;
            }
          }
          if( nextClu<pi0IsoCut_[subDetId] ) continue;
          int Nxtal_EnergGamma=0;
          int Nxtal_EnergGamma2=0;
          if(subDetId==EcalEndcap){
            if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EE[i]; Nxtal_EnergGamma2 = Ncristal_EE[j]; }
            else                           {  Nxtal_EnergGamma = Ncristal_EE[j]; Nxtal_EnergGamma2 = Ncristal_EE[i]; }
          }
          else{
            if( g1->energy()>g2->energy() ){  Nxtal_EnergGamma = Ncristal_EB[i]; Nxtal_EnergGamma2 = Ncristal_EB[j]; }
            else                           {  Nxtal_EnergGamma = Ncristal_EB[j]; Nxtal_EnergGamma2 = Ncristal_EB[i]; }
          }

          if( Nxtal_EnergGamma < nXtal_1_cut_[subDetId] ) continue;
          if( Nxtal_EnergGamma2 < nXtal_2_cut_[subDetId] ) continue;

          /// debug
          if(subDetId==EcalBarrel)
          {
                pi0MassVsIetaEB->Fill( fabs(pi0P4.eta())/0.0174, pi0P4.mass());
                pi0MassVsETEB->Fill(pi0P4.Pt(), pi0P4.mass());
          }

#ifdef DEBUG
            cout << "---------------------------------------" << endl;
            cout << "---- pi0 mass: " << pi0P4.mass() << endl;
            cout << "---------------------------------------" << endl;
#endif
            // compute region weights
            RegionWeightVector w1 = regionalCalibration_->getWeights( &(*g1), subDetId ); // region weights W_j^k for clu1
            RegionWeightVector w2 = regionalCalibration_->getWeights( &(*g2), subDetId ); // region weights W_j^k for clu2

            // append w2 to w1
            w1.insert( w1.end(), w2.begin(), w2.end() );

            float r2 = pi0P4.mass()/PI0MASS;
                  r2 = r2*r2;
            //average <eps> for cand k
            float eps_k = 0.5 * ( r2 - 1. );

            myNPi0++;

            // compute quantities needed for <eps>_j in each region j
            for(RegionWeightVector::const_iterator it = w1.begin(); it != w1.end(); ++it) {
                const uint32_t& iR = (*it).iRegion;
                const float& w = (*it).value;

                if(subDetId==EcalBarrel){
                    epsilon_EB_h[iR]->Fill( useMassInsteadOfEpsilon_? pi0P4.mass() : eps_k, w );
                    allEpsilon_EB->Fill( pi0P4.mass(), w );
                }
                else 
                {
                    epsilon_EE_h[iR]->Fill( useMassInsteadOfEpsilon_? pi0P4.mass() : eps_k, w );
                    /// debug
                    //allEpsilon_EE->Fill( eps_k, w );
                    allEpsilon_EE->Fill( pi0P4.mass(), w );
                    std::vector<DetId> mioId(regionalCalibration_->allDetIdsInEERegion(iR));
                    for(unsigned int i=0; i<mioId.size(); ++i){
                    EEDetId tmp_id(mioId.at(i));
                      if( tmp_id.zside()==-1 ) entries_EEm->Fill( tmp_id.ix(), tmp_id.iy(), w );
                      else                     entries_EEp->Fill( tmp_id.ix(), tmp_id.iy(), w );
                      //cout<<tmp_id.ix()<<"  "<<tmp_id.iy()<<"  "<<tmp_id.zside()<<endl;
                    }
                }
            }
        } // loop over clusters (g2)
    } // loop over clusters to make pi0 
}


// ------------ method called once each job just before starting event loop  ------------
void 
FillEpsilonPlot::beginJob()
{
    // output file
    char fileName[200];
    sprintf(fileName,"%s/%s", outputDir_.c_str(), outfilename_.c_str());
    outfile_ = new TFile(fileName,"RECREATE");
    if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << string(fileName) << "\n";

    /// testing the EE eta ring
    TH2F eep("eep","EE+",102,0.5,101.5,102,-0.5,101.5);
    TH2F eem("eem","EE-",102,0.5,101.5,102,-0.5,101.5);
    cout << "beginJob:: about to call allDetIds" << endl;

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
 
#ifdef MVA_REGRESSIO_EE  
  ifstream file;
  file.open(Endc_x_y_.c_str(), ifstream::in);
  VectRing.clear();
  while ( !file.eof() ) {
    string Line;
    getline( file, Line);
    string value;
    stringstream MyLine(Line);

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
#endif

#ifdef MVA_REGRESSIO
   // TTree_JoshMva = new TTree("TTree_JoshMva","MVA corrections");
   // TTree_JoshMva->Branch("Correction1_mva", &Correction1_mva, "Correction1_mva/F");
   // TTree_JoshMva->Branch("Correction2_mva", &Correction2_mva, "Correction2_mva/F");
   // TTree_JoshMva->Branch("iEta1_mva", &iEta1_mva, "iEta1_mva/I");
   // TTree_JoshMva->Branch("iEta2_mva", &iEta2_mva, "iEta2_mva/I");
   // TTree_JoshMva->Branch("iPhi1_mva", &iPhi1_mva, "iPhi1_mva/I");
   // TTree_JoshMva->Branch("iPhi2_mva", &iPhi2_mva, "iPhi2_mva/I");
   // TTree_JoshMva->Branch("Pt1_mva", &Pt1_mva, "Pt1_mva/F");
   // TTree_JoshMva->Branch("Pt2_mva", &Pt2_mva, "Pt2_mva/F");
   // TTree_JoshMva->Branch("Mass_mva", &Mass_mva, "Mass_mva/F");
   // TTree_JoshMva->Branch("MassOr_mva", &MassOr_mva, "MassOr_mva/F");
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
  
   float diff = fabs(phi2 - phi1);
  
   while (diff >acos(-1)) diff -= 2*acos(-1);
   while (diff <= -acos(-1)) diff += 2*acos(-1);
  
   return diff; 
  
}

bool
FillEpsilonPlot::GetHLTResults(const edm::Event& iEvent, std::string s){

    edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
    iEvent.getByLabel(triggerTag_, hltTriggerResultHandle);

    edm::TriggerNames HLTNames;

    HLTNames = iEvent.triggerNames(*hltTriggerResultHandle);
    std::string tempnames;
    int hltCount = hltTriggerResultHandle->size();
    TRegexp reg(TString( s.c_str()) );

    for (int i = 0 ; i != hltCount; ++i) {
	  TString hltName_tstr(HLTNames.triggerName(i));
//cout<<"Trigger: "<<hltName_tstr<<endl;
	  std::string hltName_str(HLTNames.triggerName(i));
		if ( hltName_tstr.Contains(reg) ){
		//    cout<<hltTriggerResultHandle->accept(i)<<endl;
		    return hltTriggerResultHandle->accept(i);
            }
	  }
    return false;
} // HLT isValid

bool
FillEpsilonPlot::getTriggerByName( std::string s ) {
    std::map< std::string, int >::iterator currentTrigger;
    currentTrigger = l1TrigNames_.find(s);
    if(currentTrigger != l1TrigNames_.end())
	  return l1TrigBit_[currentTrigger->second];
    else 
	  std::cout << "Trigger Name not found" << std::endl;

    return false;
}

bool 
FillEpsilonPlot::getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // trigger Handles
    edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
    iEvent.getByLabel( l1TriggerTag_, gtReadoutRecord);

    const DecisionWord& gtDecisionWord = gtReadoutRecord->decisionWord();

    //int nL1bits = int( gtDecisionWord.size() );

    int thisBit =0;

    for (std::vector<bool>::const_iterator itBit = gtDecisionWord.begin(); itBit != gtDecisionWord.end(); ++itBit, ++thisBit) {
	  //l1TrigBit_[thisBit] = gtDecisionWord.at(thisBit);
	  if( gtDecisionWord.at(thisBit) )  
		triggerComposition->Fill(thisBit);
    }


    /// referece

    // edm::ESHandle<L1GtTriggerMenu> menuRcd;
    // iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    // const L1GtTriggerMenu* menu = menuRcd.product();

    //bool l1SingleEG2 = menu->gtAlgorithmResult("L1_SingleEG2", gtDecisionWord);

    //return gtDecisionWord.at(l1TrigNames_["L1_SingleEG20"]);
    return true;
}


// ------------ method called once each job just after ending the event loop  ------------
    void 
FillEpsilonPlot::endJob() 
{
    outfile_->cd();

#ifdef MVA_REGRESSIO
    //TTree_JoshMva->Write();
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
    //JSON
    //hev_TOT->SetBinContent(1,ev_TOT); hev_TOT->Write();
    //hev_JSON->SetBinContent(1,ev_JSON); hev_JSON->Write();

    allEpsilon_EB->Write();
    allEpsilon_EE->Write();
    entries_EEp->Write();
    entries_EEm->Write();
    pi0MassVsIetaEB->Write();
    pi0MassVsETEB->Write();
    triggerComposition->Write();
    if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EB_h, "Barrel" ,  regionalCalibration_->getCalibMap()->getNRegionsEB() );
    if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EE_h, "Endcap" ,  regionalCalibration_->getCalibMap()->getNRegionsEE() );

#ifdef MVA_REGRESSIO
    //delete TTree_JoshMva;
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
    else{                  cout<<"Cont. Correction too low... I'm using 1. Check if all is right please. (nBin = "<<nBin<<" )"<<endl;  return 1.;}
}

// ------------ method called when starting to processes a run  ------------
    void 
FillEpsilonPlot::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();


    std::map< std::string, int >::iterator currentTrigger;

    if(l1TrigNames_.size()>0) {
	  bool triggerChanged = false;
	  for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
		currentTrigger = l1TrigNames_.find((algo->second).algoName());
		if (currentTrigger == l1TrigNames_.end() || currentTrigger->second != (algo->second).algoBitNumber()) {
		    triggerChanged = true;
		    break;
		}
		// if(l1TrigNames_[(algo->second).algoName()] != (algo->second).algoBitNumber()) {
		//     cout << "Trigger numbering has changed" << endl;
		//     l1TrigNames_[(algo->second).algoName()] = (algo->second).algoBitNumber();
		// }
	  }
	  if(!triggerChanged) return;
	  cout << "beginRun:: Trigger names / ordering changed" << endl;
    }

    cout << "beginRun:: Filling trigger names" << endl;

    /// filling trigger map
    l1TrigNames_.clear();
    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
	  l1TrigNames_[(algo->second).algoName()] = (algo->second).algoBitNumber();

	  /// using same loop to set trigger histogram labels
	  if(!areLabelsSet_)
		triggerComposition->GetXaxis()->SetBinLabel((algo->second).algoBitNumber()+1,(algo->second).algoName().c_str());
    }

    if(!areLabelsSet_)
    {
	  areLabelsSet_ = true;
	  cout << "beginRun:: setting labels of triggerComposition histogram" << endl;
    }
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
    for( size_t i=0; i<VectRing.size(); i++){
	  if(  VectRing[i].iX != x || VectRing[i].iY != y ) index++;
	  if(  VectRing[i].iX == x && VectRing[i].iY == y ) break;
    }
    if(debug3) cout<<"Found : "<<VectRing[index].iX<<" and "<<VectRing[index].iY<<endl;
    return VectRing[index].Ring;
}


//define this as a plug-in
DEFINE_FWK_MODULE(FillEpsilonPlot);
