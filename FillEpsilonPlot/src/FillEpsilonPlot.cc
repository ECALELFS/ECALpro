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
// $Id: FillEpsilonPlot.cc,v 1.38 2013/01/30 15:55:29 lpernie Exp $
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
#include "CalibCode/FillEpsilonPlot/interface/JSON.h"

//MVA Stuff
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif

#ifdef MVA_REGRESSIO_EB
#include "CalibCode/GBRTrain/interface/GBRApply.h"
#include "CalibCode/EgammaObjects/interface/GBRForest.h"
#include "Cintex/Cintex.h"
#include "TLorentzVector.h"
#endif

using namespace TMVA;

//Function
double max_array(double *A, int n);
double max(double x, double y);

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
#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO_EB)
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

#ifdef MVA_REGRESSIO_EB
    EBweight_file_pi01 = TFile::Open( MVAEBContainmentCorrections_01_.c_str() );
    EBweight_file_pi02 = TFile::Open( MVAEBContainmentCorrections_02_.c_str() );
    forest_EB_pi01 = (GBRForest *)EBweight_file_pi01->Get("Correction");    
    forest_EB_pi02 = (GBRForest *)EBweight_file_pi02->Get("Correction");    
#endif
    //JSON
    ev_TOT=0;
    ev_JSON=0;
    myjson = new JSON(jsonFile_.c_str());
    hev_TOT  = new TH1F("hev_TOT", "Events before Json file selection",1, 0.5, 1.5);
    hev_JSON = new TH1F("hev_JSON", "Events after Json file selection",1, 0.5, 1.5);
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

#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO_EB)
    delete EBPHI_ConCorr_p;
    delete EBPHI_ConCorr_m;
#endif

#ifdef MVA_REGRESSIO_EB
    //delete TTree_JoshMva;
    //delete reader_;
#endif

    //JSON
    delete myjson;
    delete hev_TOT;
    delete hev_JSON;
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

   ///------- HANDLE ---------
    
    if(SystOrNot_==1. && int(iEvent.id().event())%2!=0 ) return;
    else if(SystOrNot_==2. && int(iEvent.id().event())%2==0 ) return;
    //cout<<"passed!"<<endl;
    //ES

    edm::ESHandle<CaloGeometry> geoHandle;
    iSetup.get<CaloGeometryRecord>().get(geoHandle);
    geometry = geoHandle.product();
    estopology_ = new EcalPreshowerTopology(geoHandle);
    const CaloSubdetectorGeometry *geometry_p = geometry->getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;
    esGeometry_ = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p));

   /// ------ EB CLUSTERS -------
   //if( !getTriggerResult(iEvent, iSetup)) return;

   //JSON FILE
   ev_TOT++;
   if ( !myjson->isGoodLS( iEvent.id().run() , iEvent.luminosityBlock() ) ) return;
   ev_JSON++;

   std::vector< CaloCluster > ebclusters; // contains the output clusters
   ebclusters.clear();
   vNxtal.clear();
   vs4s9.clear();
   vs2s9.clear();
   vs2s9.clear();

   Ncristal_EB.clear();
   bool EB_HLT=true, EE_HLT=true;
   if( HLTResults_ ){
       EB_HLT = GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
       EE_HLT = GetHLTResults(iEvent, "AlCa_EcalPi0EEonly.*");
   }
    
   if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) && EB_HLT ) fillEBClusters(ebclusters, iEvent);

   std::vector< CaloCluster > eseeclusters; 
   std::vector< CaloCluster > eseeclusters_tot; 
   eseeclusters.clear(); eseeclusters_tot.clear();
   Ncristal_EE.clear();  
   if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) && EE_HLT ) fillEEClusters(eseeclusters, eseeclusters_tot, iEvent);


   //EB
   if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(ebclusters,EcalBarrel);
   //EE -> choose if use normal cluster or those with ES position matching
   #ifndef NEWCUT_ESPOS
   if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(eseeclusters,EcalEndcap);
   #endif
   #ifdef NEWCUT_ESPOS
   if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) computeEpsilon(eseeclusters_tot,EcalEndcap);
   #endif

   delete estopology_;
   //delete geometry_p;
//   BeamSpot beamSpot;
//   edm::Handle<BeamSpot> beamSpotHandle;
//   eventCont->getByLabel("offlineBeamSpot", beamSpotHandle);
//   beamSpot = *beamSpotHandle;
//
//   double xPV = beamSpot.x0();
//   double yPV = beamSpot.y0();
//   double zPV = beamSpot.z0();

}


/*===============================================================*/
void FillEpsilonPlot::fillEBClusters(std::vector< CaloCluster > & ebclusters, const edm::Event& iEvent)
/*===============================================================*/
{
   edm::Handle< EBRecHitCollection > ebHandle;
   bool goodhandle = iEvent.getByLabel ( EBRecHitCollectionTag_, ebHandle);
   if(!goodhandle) cout << "problem in eventCont->getByLabel ( ebLabel, ebHandle)" << endl;

   std::vector<EcalRecHit> ebseeds;
   ebseeds.clear();

   typedef std::map< EBDetId, bool > XtalInUse;
   XtalInUse isUsed; // map of which xtals have been used

   std::vector<EBDetId> detIdEBRecHits;
   std::vector<EcalRecHit> EBRecHits;

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
         XtalInUse::const_iterator mapit = isUsed.find( seed_id );
         if( mapit != isUsed.end() ) continue; // seed already in use

         // find 3x3 matrix of xtals
         int clusEtaSize_(3), clusPhiSize_(3);
         std::vector<DetId> clus_v = ebtopology_->getWindow(seed_id,clusEtaSize_,clusPhiSize_);       
         // needed for position calculator
         std::vector<std::pair<DetId,float> > clus_used;

         /// debug
         //cout << "seed #" << seed_c << " clus_v.size():" << clus_v.size() << endl;

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
            XtalInUse::const_iterator mapit = isUsed.find( thisId );
            if( mapit != isUsed.end() ) continue; // xtal already used

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

         float s4s9_tmp[4];
         for(int i=0;i<4;i++){ 
            s4s9_tmp[i]= 0;
         }

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
               isUsed[ RecHitsInWindow[j]->id() ] = true;

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

#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO_EB)
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
#ifdef MVA_REGRESSIO_EB
         vNxtal.push_back( RecHitsInWindow.size() );
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
    eeseeds.clear();

    vector <double> eeclusterS4S9; //s4s9 stored for each 3x3 eeclusters
    vector <bool> eeusedcluster;  //used for debugging session, allows to check if the cluster has been used inside the printing procedure or not
    eeusedcluster.clear();
    
    std::vector< CaloCluster > eeclusters; // contains the output eeclusters
    //std::vector< CaloCluster > eseeclusters;
    eeclusters.clear();
    //eseeclusters.clear();     

    int dc = 0;

    // sort by energy and find the eeseeds
    for(EERecHitCollection::const_iterator ite= eeHandle->begin(); ite != eeHandle->end(); ++ite, ++dc) {
       #ifndef NEWCUT_ESPOS
       if(ite->energy() > 0.800)  eeseeds.push_back( *ite ); //Original cut 0.800
       #endif
       //----or cut in ET-----
       #ifdef NEWCUT_ESPOS 
       int idXtal= ite->id();
       float posTotalEnergy(ite->energy());
       float T0 = PCparams_.param_T0_barl_;
       float maxDepth = PCparams_.param_X0_ * ( T0 + log( posTotalEnergy ) );
       float maxToFront = geom_->getPosition(idXtal).mag(); // to front face
       float depth = maxDepth + maxToFront - geom_->getPosition(idXtal).mag() ;
       GlobalPoint posThis = geom_->getPosition(idXtal,depth);
       if(ite->energy()/cosh(posThis.eta()) > 0.5)      eeseeds.push_back( *ite );
       #endif
     } // loop over xtals

    /// debug
    // cout << "eeseeds.size() = " << eeseeds.size() << " / " << dc << endl;

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

      eeusedcluster.push_back(false);
    } //loop over seeds to make eeclusters

    //cout << "noES.size() = " << eeclusters.size() << endl;
    // for(int i=0; i<int(eeclusters.size()); ++i)
    // {
    //     cout << "---clus NOES  #" << i << " energy: " << eeclusters.at(i).energy() << endl;
    // }


    //cout << "fillEEClusters :: eeclusters.size() = " << eeclusters.size() << endl;

    /************************** ENDCAP-PRESHOWER MATCHING ************************/

    //loop over eecluster to find matches with preshower
    int ind=0;
    std::vector<int> Nxtal; Nxtal.clear();
    std::vector<int> Nxtal_tot; Nxtal_tot.clear();

    for( std::vector<CaloCluster>::const_iterator eeclus_iter  = eeclusters.begin(); eeclus_iter != eeclusters.end(); ++eeclus_iter, ++ind)
    {
//@
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
#if defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO_EB)
             if(useEEContainmentCorrections_) tempenergy *= containmentCorrections_.getContainmentPointCorrectionsEE( tempenergy , (eeclus_iter->position()).eta() );
#endif

             eseeclusters.push_back( CaloCluster( tempenergy, eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
             Nxtal.push_back(Ncristal_EE[ind]);
             //The perfect cluster
             //double DZ2 = (preshowerclusterp2.get_z()-preshowerclusterp1.get_z() )/2.;
             //GlobalPoint posClu(preshowerclusterp1.get_x()*(1.+DZ2/preshowerclusterp1.get_z() ),preshowerclusterp2.get_y()*(1.-DZ2/preshowerclusterp2.get_z()),(preshowerclusterp1.get_z()+preshowerclusterp2.get_z() )/2. );
             double DZ2 = (preshowerclusterp2.z()-preshowerclusterp1.z() )/2.;
             GlobalPoint posClu(preshowerclusterp1.x()*(1.+DZ2/preshowerclusterp1.z() ),preshowerclusterp2.y()*(1.-DZ2/preshowerclusterp2.z()),(preshowerclusterp1.z()+preshowerclusterp2.z() )/2. );

                if( fabs(preshowerclusterp1.z())>30  && fabs(preshowerclusterp2.z())>30){

                  math::XYZPoint posit(posClu.x(),posClu.y(),posClu.z());
                  eseeclusters_tot.push_back( CaloCluster( tempenergy, posit, CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
                  Nxtal_tot.push_back(Ncristal_EE[ind]);
                }
           }
        }
}

else{
    eseeclusters_tot.push_back( CaloCluster( eeclus_iter->energy(), eeclus_iter->position(), CaloID(CaloID::DET_ECAL_ENDCAP),  eeclus_iter->hitsAndFractions(), CaloCluster::undefined, eeclus_iter->seed() ) );
    Nxtal_tot.push_back(Ncristal_EE[ind]);
}
    }//end of the matching loop

    Ncristal_EE.clear();
    #ifndef NEWCUT_ESPOS
    Ncristal_EE = Nxtal;
    #endif
    #ifdef NEWCUT_ESPOS
    Ncristal_EE = Nxtal_tot; 
    #endif
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
#if !defined(NEW_CONTCORR) && defined(MVA_REGRESSIO_EB)
     //JOSHMVA
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
         vector<float> value_pi01;
         value_pi01.push_back( G_Sort_1.E()/G_Sort_2.E() );
         value_pi01.push_back( G_Sort_1.Pt() );
         value_pi01.push_back( vNxtal[ind1] );
         value_pi01.push_back( vNxtal[ind2] );
         value_pi01.push_back( vs4s9[ind1] );
         value_pi01.push_back( vs1s9[ind1] );
         value_pi01.push_back( vs2s9[ind1] );
         value_pi01.push_back( iEta1 );
         value_pi01.push_back( iPhi1 );
         value_pi01.push_back( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
         value_pi01.push_back( iEta1%5 );
         value_pi01.push_back( iPhi1%2 );
         value_pi01.push_back( (TMath::Abs(iEta1)<=25)*(iEta1%25) + (TMath::Abs(iEta1)>25)*((iEta1-25*TMath::Abs(iEta1)/iEta1)%20) );
         value_pi01.push_back( iPhi1%20 );
         float Correct1 = gbrapply->ApplyAsFriend_calib( forest_EB_pi01, value_pi01);
         //cout<<"Correction1: "<<Correct1<<" iEta: "<<iEta1<<" iPhi "<<iPhi1<<" E "<<G_Sort_1.E()<<endl;
   
         vector<float> value_pi02;
         value_pi02.push_back( G_Sort_1.E()/G_Sort_2.E() );
         value_pi02.push_back( G_Sort_2.Pt() );
         value_pi02.push_back( vNxtal[ind1] );
         value_pi02.push_back( vNxtal[ind2] );
         value_pi02.push_back( vs4s9[ind2] );
         value_pi02.push_back( vs1s9[ind2] );
         value_pi02.push_back( vs2s9[ind2] );
         value_pi02.push_back( iEta2 );
         value_pi02.push_back( iPhi2 );
         value_pi02.push_back( sqrt(pow((iEta1-iEta2),2)+pow((iPhi1-iPhi2),2)));
         value_pi02.push_back( iEta2%5 );
         value_pi02.push_back( iPhi2%2 );
         value_pi02.push_back( (TMath::Abs(iEta2)<=25)*(iEta2%25) + (TMath::Abs(iEta2)>25)*((iEta2-25*TMath::Abs(iEta2)/iEta2)%20) );
         value_pi02.push_back( iPhi2%20 );
         float Correct2 = gbrapply->ApplyAsFriend_calib( forest_EB_pi02, value_pi02);
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

           // ptPi0 cut: default alca > 1.6
           if(pi0P4.Pt() < pi0PtCut_[subDetId]) continue;
           if( g1P4.eta() == g2P4.eta() && g1P4.phi() == g2P4.phi() ) continue;
#ifndef NEWCUT_ESPOS  
          float isolation = 0;
          for( std::vector<CaloCluster>::const_iterator g3  = clusters.begin(); g3 != clusters.end(); ++g3){
             if( g3->seed() == g1->seed() || g3->seed() == g2->seed()) continue;
             math::PtEtaPhiMLorentzVector g3P4( g3->energy()/cosh(g3->eta()), g3->eta(), g3->phi(), 0. );
             float drcl = GetDeltaR(pi0P4.eta(),g3P4.eta(),pi0P4.phi(),g3P4.phi()); 
             float dretacl = fabs(pi0P4.eta() - g3P4.eta());
             if( drcl > 0.2 ||  dretacl > 0.05 ) continue; 
             isolation += g3P4.Pt();
          }//End of loop over gammas to measure isolation

            if(isolation/pi0P4.Pt() > pi0IsoCut_[subDetId]) continue;  //Pi0 isolation cut
#endif

#ifdef NEWCUT_ESPOS
          //if(subDetId==EcalEndcap){

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
            //if( nextClu<0.2 ) continue;
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

            //###if( Nxtal_EnergGamma >= 8 ) continue;
            if( Nxtal_EnergGamma < nXtal_1_cut_[subDetId] ) continue;
            if( Nxtal_EnergGamma2 < nXtal_2_cut_[subDetId] ) continue;
          //}
#endif

          /// debug
          if(subDetId==EcalBarrel)
          {
                pi0MassVsIetaEB->Fill( fabs(pi0P4.eta())/0.0174, pi0P4.mass());
                pi0MassVsETEB->Fill(pi0P4.Pt(), pi0P4.mass());
          }
    
          // use candidates in the signal mass window just for the Barrel
	// REMOVED FOR HAVING SIDEBANDS @@
	// if(subDetId==EcalBarrel && (pi0P4.mass()<0.095 || pi0P4.mass()>0.175) ) continue;

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
    
#ifdef MVA_REGRESSIO_EB
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

#ifdef SELECTION_TREE
    CutVariables_EB->Write();
    CutVariables_EE->Write();
    Pi0Info_EB->Write();
    Pi0Info_EE->Write();
#endif
    //JSON
    hev_TOT->SetBinContent(1,ev_TOT); hev_TOT->Write();
    hev_JSON->SetBinContent(1,ev_JSON); hev_JSON->Write();

    allEpsilon_EB->Write();
    allEpsilon_EE->Write();
    entries_EEp->Write();
    entries_EEm->Write();
    pi0MassVsIetaEB->Write();
    pi0MassVsETEB->Write();
    triggerComposition->Write();
    if( (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EB_h, "Barrel" ,  regionalCalibration_->getCalibMap()->getNRegionsEB() );
    if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) writeEpsilonPlot(epsilon_EE_h, "Endcap" ,  regionalCalibration_->getCalibMap()->getNRegionsEE() );

#ifdef MVA_REGRESSIO_EB
    //TTree_JoshMva->Write();
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

//define this as a plug-in
DEFINE_FWK_MODULE(FillEpsilonPlot);
