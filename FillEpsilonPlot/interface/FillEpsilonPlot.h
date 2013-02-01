#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "CalibCode/CalibTools/interface/PosCalcParams.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"
#include "CalibCode/CalibTools/interface/EcalEnerCorr.h"
#include "CalibCode/CalibTools/interface/EcalCalibTypes.h"
#include "CalibCode/CalibTools/interface/EcalRegionalCalibration.h"
#include "CalibCode/CalibTools/interface/EcalPreshowerHardcodedTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//JSON
#include "CalibCode/FillEpsilonPlot/interface/JSON.h"

//#define SELECTION_TREE
#define NEWCUT_ESPOS
//#define RESIDUALCORR
//#define NEW_CONTCORR
#define MVA_REGRESSIO_EB
//#define DEBUG_HIGHETA
//To see why at high Eta you have a hot eta-ring

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
#endif

enum calibGranularity{ xtal, tt, etaring };
//enum subdet{ thisIsEE, thisIsEB }; 

using namespace reco;

class FillEpsilonPlot : public edm::EDAnalyzer {
   public:
      explicit FillEpsilonPlot(const edm::ParameterSet&);
      ~FillEpsilonPlot();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ---------- user defined ------------------------
      void fillEBClusters(std::vector< CaloCluster > & ebclusters, const edm::Event& iEvent);
      void fillEEClusters(std::vector< CaloCluster > & eseeclusters,std::vector< CaloCluster > & eseeclusters_tot, const edm::Event& iEvent);
      void computeEpsilon(std::vector< CaloCluster > & clusters, int subDetId);
      float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
      float DeltaPhi(float phi1, float phi2);
      double min( double a, double b);

      TH1F** initializeEpsilonHistograms(const char *name, const char *title, int size );
      void  deleteEpsilonPlot(TH1F **h, int size);
      void  writeEpsilonPlot(TH1F **h, const char *folder, int size);
      bool getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool getTriggerByName( std::string s );
      bool GetHLTResults(const edm::Event& iEvent, std::string s);

      float EB_Cont_Corr( float PT, int iEta);
      float EB_Resid_Corr(int iEta, int iPhi);
      float EBPHI_Cont_Corr(float PT, int giPhi, int ieta);
      void  EB_Cont_Corr_load(std::string FileName );
      void  EB_Resid_Corr_load( bool isETA );
      void  EBPHI_Cont_Corr_load(std::string FileName );
      TH1F * EB_ConCorr_1;
      TH1F * EB_ConCorr_2;
      TH1F * EB_ConCorr_3;
      TH1F * EB_ConCorr_4;
      TH1F * EB_ConCorr_5;
      TH1F * EB_ConCorr_6;
      TH1F * EB_ConCorr_7;
      TH1F * EBPHI_ConCorr_p;
      TH1F * EBPHI_ConCorr_m;
      //Residual
      TH1F * EB_Residual_Eta;
      TH1F * EB_Residual_Phi;

      // ----------member data ---------------------------
 
      const EcalPreshowerGeometry *esGeometry_;     
      const CaloGeometry* geometry;

      std::string outfilename_;
      std::string externalGeometry_;
      std::string calibMapPath_; 
      std::string jsonFile_; 
      std::string ebContainmentCorrections_;
      std::string MVAEBContainmentCorrections_01_;
      std::string MVAEBContainmentCorrections_02_;
      std::string ebPHIContainmentCorrections_;
      std::string eeContainmentCorrections_;
      std::string Barrel_orEndcap_;
      bool useEBContainmentCorrections_;
      bool useEEContainmentCorrections_;
      bool useOnlyEEClusterMatchedWithES_;
      bool Is2012_;

      edm::InputTag EBRecHitCollectionTag_;
      edm::InputTag EERecHitCollectionTag_;
      edm::InputTag ESRecHitCollectionTag_;
      edm::InputTag l1TriggerTag_;
      edm::InputTag triggerTag_;

      PosCalcParams PCparams_;
      //const double preshowerStartEta_ =  1.653;

      ECALGeometry* geom_;
      CaloTopology *ebtopology_;
      CaloTopology *eetopology_;
      CaloSubdetectorTopology *estopology_;
      
      EcalEnerCorr containmentCorrections_;
      //bool noCorrections_;

      std::string calibTypeString_;
      calibGranularity calibTypeNumber_;

      // selection criteria
      double gPtCut_[3];
      double pi0PtCut_[3];
      double pi0IsoCut_[3];
      double nXtal_1_cut_[3];
      double nXtal_2_cut_[3];
      double S4S9_cut_[3];
      double SystOrNot_;

      /// all the three options have to be instantiated to allow the
      //choice at runtime
      EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib;
      EcalRegionalCalibration<EcalCalibType::EtaRing> etaCalib;
      EcalRegionalCalibration<EcalCalibType::TrigTower> TTCalib;

      EcalRegionalCalibrationBase *regionalCalibration_;

      int currentIteration_;
      string outputDir_;

      TFile *outfile_;
      TFile *externalGeometryFile_;

      std::vector<int> Ncristal_EE;
      std::vector<int> Ncristal_EB;

      TH1F **epsilon_EB_h;  // epsilon distribution by region
      TH1F **epsilon_EE_h;  // epsilon distribution in EE
      TH1F *allEpsilon_EE; 
      TH1F *allEpsilon_EB;
      TH2F *entries_EEp;
      TH2F *entries_EEm;
      TH2F *pi0MassVsIetaEB;
      TH2F *pi0MassVsETEB;
      bool useMassInsteadOfEpsilon_;

#ifdef SELECTION_TREE
      Float_t NSeeds_EB, Xclus_EB, Yclus_EB, Zclus_EB, e3x3_EB, S4S9_EB, PTClus_EB;//@@@@
      void Fill_NSeeds_EB(float nSeed){ NSeeds_EB=nSeed; };
      void Fill_xClus_EB(float x){ Xclus_EB=x; };
      void Fill_yClus_EB(float y){ Xclus_EB=y; };
      void Fill_zClus_EB(float z){ Xclus_EB=z; };
      void Fill_e3x3_EB(float E3x3){ e3x3_EB=E3x3; };
      void Fill_S4S9_EB(float s4s9){ S4S9_EB=s4s9; };
      void Fill_PtClus_EB(float clus){ PTClus_EB=clus; };
      TTree *CutVariables_EB; 
      Float_t NSeeds_EE, Xclus_EE, Yclus_EE, Zclus_EE, e3x3_EE, S4S9_EE, PTClus_EE;
      void Fill_NSeeds_EE(float nSeed){ NSeeds_EE=nSeed; };
      void Fill_xClus_EE(float x){ Xclus_EE=x; };
      void Fill_yClus_EE(float y){ Xclus_EE=y; };
      void Fill_zClus_EE(float z){ Xclus_EE=z; };
      void Fill_e3x3_EE(float E3x3){ e3x3_EE=E3x3; }; 
      void Fill_S4S9_EE(float s4s9){ S4S9_EE=s4s9; };
      void Fill_PtClus_EE(float clus){ PTClus_EE=clus; };
      TTree *CutVariables_EE;

      Float_t PtPi0_EB, mpi0_EB, Etapi0_EB, Phipi0_EB, Epsilon_EB;
      void Fill_PtPi0_EB(float pt){ PtPi0_EB=pt; };
      void Fill_mpi0_EB(float m){ mpi0_EB=m; };
      void Fill_etapi0_EB( float eta){ Etapi0_EB =eta; };
      void Fill_phipi0_EB( float phi){ Phipi0_EB =phi; };
      void Fill_Epsilon_EB(float eps ){ Epsilon_EB=eps; };
      TTree *Pi0Info_EB;
      Float_t PtPi0_EE, mpi0_EE,Etapi0_EE, Phipi0_EE, Epsilon_EE;
      void Fill_PtPi0_EE(float pt){ PtPi0_EE=pt; };
      void Fill_mpi0_EE(float m){ mpi0_EE=m; };
      void Fill_etapi0_EE( float eta){ Etapi0_EE =eta; };
      void Fill_phipi0_EE( float phi){ Phipi0_EE =phi; };
      void Fill_Epsilon_EE(float eps ){ Epsilon_EE=eps; };
      TTree *Pi0Info_EE;
#endif

      std::string ContCorr_EB_;
      TH1F *triggerComposition;
      bool areLabelsSet_;

      std::vector<std::string> alcaL1TrigNames_;
      //std::vector<std::string> l1TrigNames_(128,"");
      std::map< std::string, int > l1TrigNames_;
      bool l1TrigBit_[128];
#ifdef MVA_REGRESSIO_EB
      Float_t   MVA_Pt_1, MVA_S1S9_1, MVA_S2S9_1, MVA_S4S9_1;
      //Int_t   MVA_Nxtal_1, MVA_Eta_1, MVA_Phi_1, MVA_Eta_1on5, MVA_Phi_1on2, MVA_Eta_1on2520, MVA_Phi_1on20;
      Float_t   MVA_Nxtal_1, MVA_Eta_1, MVA_Phi_1, MVA_Eta_1on5, MVA_Phi_1on2, MVA_Eta_1on2520, MVA_Phi_1on20;
      Float_t MVA_E3x3_1;
      Float_t MVA_E3x3MC_1;
      //TMVA::Reader *reader_;
      vector<float> vNxtal;
      vector<float> vs4s9;
      vector<float> vs1s9;
      vector<float> vs2s9;
      TFile *EBweight_file_pi01;
      TFile *EBweight_file_pi02;
      const GBRForest *forest_EB_pi01;
      const GBRForest *forest_EB_pi02;
      GBRApply *gbrapply;
      //TTree *TTree_JoshMva;
      //Float_t Correction1_mva, Correction2_mva, Pt1_mva, Pt2_mva, Mass_mva, MassOr_mva;
      //Int_t   iEta1_mva, iPhi1_mva, iEta2_mva, iPhi2_mva;
#endif   
#ifdef DEBUG_HIGHETA
      static const int nMaxMC = 350;
      TTree *Tree_HighEta;
      Int_t nPi0HighEta;
      Float_t ESpos_eta_1[nMaxMC];
      Float_t ESpos_phi_1[nMaxMC];
      Float_t ESpos_eta_2[nMaxMC];
      Float_t ESpos_phi_2[nMaxMC];
      Float_t ESpos_mass[nMaxMC];
      TTree *Tree_HighEta_clus;
      Int_t nPi0HighEta_clus;
      Float_t ESpos_eta_1_clus[nMaxMC];
      Float_t ESpos_phi_1_clus[nMaxMC];
      Float_t ESpos_eta_2_clus[nMaxMC];
      Float_t ESpos_phi_2_clus[nMaxMC];
#endif
      //JSON
      int ev_TOT;
      int ev_JSON;
      JSON* myjson;
      TH1F *hev_TOT;
      TH1F *hev_JSON;    
};
