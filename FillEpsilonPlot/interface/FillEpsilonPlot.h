#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLorentzVector.h"

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
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CalibCode/FillEpsilonPlot/interface/JSON.h"
// to get L1 info
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h" // included to get L1 info
//L1                                                                                                                                         
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#define NPI0MAX 30000
#define NL1SEED GlobalAlgBlk::maxPhysicsTriggers  // was 128
//#define SELECTION_TREE
//#define NEW_CONTCORR    // to use Yong's parametric CC, act on both EE and EB
#define MVA_REGRESSIO     // to use regression in EB
//#define MVA_REGRESSIO_Tree  // when using regression (defined MVA_REGRESSIO), decide to store some variables in a tree. This is for EB
//#define MVA_REGRESSIO_EE    // should be as MVA_REGRESSIO but actually it also act as MVA_REGRESSIO_Tree for EE (define it to use regression in EE)
//#define MVA_REGRESSIO_EE_Tree  // not used anywere apparently

// developing new feature to have Yong's parametric containment corrections in EE and MVA regression in EB
//
// We would use regression in 2012 (or 2016) for EB and parametric containment corrections in EE
// To do it, you can uncomment the following directive and also uncomment MVA_REGRESSIO while keeping NEW_CONTCORR commented
// the temporary implementation of this solution is done so that REGRESS_AND_PARAM_CONTCORR substitutes MVA_REGRESSIO and NEW_CONTCORR, but it is not harmful
// to keep MVA_REGRESSIO uncommented

//#define REGRESS_AND_PARAM_CONTCORR

// then in parameters.py you'll have 
//    if ContainmentCorrection == 'mixed':
//       useEBContainmentCorrections = 'False'  // no parametric CC in EB
//       useEEContainmentCorrections = 'True'   // parametric CC in EB
//       useMVAContainmentCorrections = True    //  regression (eventually used only in EB)
//       new_pi0ContainmentCorrections = False  // use new 2016 regression: True to use it, False to use old one (useMVAContainmentCorrections must be true anyway)


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

enum calibGranularity{ xtal, tt, etaring };
//enum subdet{ thisIsEE, thisIsEB }; 

struct iXiYtoRing {
    int iX;
    int iY;
    int sign;
    int Ring;
    iXiYtoRing() : iX(0), iY(0), sign(0), Ring(-1) { }
} GiveRing;

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
      void fillEBClusters(std::vector< CaloCluster > & ebclusters, const edm::Event& iEvent, const EcalChannelStatus &channelStatus);
      void fillEEClusters(std::vector< CaloCluster > & eseeclusters,std::vector< CaloCluster > & eseeclusters_tot, const edm::Event& iEvent, const EcalChannelStatus &channelStatus);
      std::vector< CaloCluster > MCTruthAssociate(std::vector< CaloCluster > & clusters, double deltaR, bool isEB);
      std::vector< CaloCluster > MCTruthAssociateMultiPi0(std::vector< CaloCluster > & clusters, double deltaR, bool isEB);
      void computeEpsilon(std::vector< CaloCluster > & clusters, int subDetId);
      bool checkStatusOfEcalRecHit(const EcalChannelStatus &channelStatus,const EcalRecHit &rh);
      bool isInDeadMap( bool isEB, const EcalRecHit &rh );
      float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
      float DeltaPhi(float phi1, float phi2);
      double min( double a, double b);

      TH1F** initializeEpsilonHistograms(const char *name, const char *title, int size );
      void deleteEpsilonPlot(TH1F **h, int size);
      void writeEpsilonPlot(TH1F **h, const char *folder, int size);
      bool getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool getTriggerByName( std::string s );
      bool GetHLTResults(const edm::Event& iEvent, std::string s);

      float EBPHI_Cont_Corr(float PT, int giPhi, int ieta);
      void  EBPHI_Cont_Corr_load(std::string FileName );
      TFile* DeadMap;
      TH2F * EBMap_DeadXtal;
      TH2F * EEmMap_DeadXtal;
      TH2F * EEpMap_DeadXtal;
      TH1F * EBPHI_ConCorr_p;
      TH1F * EBPHI_ConCorr_m;
#if (defined(NEW_CONTCORR) && !defined(MVA_REGRESSIO)) || defined(REGRESS_AND_PARAM_CONTCORR)
      EcalEnerCorr containmentCorrections_;
#endif
      // ----------member data ---------------------------
      edm::Handle< EBRecHitCollection > ebHandle;
      edm::Handle< EBRecHitCollection > eeHandle;
      edm::Handle< ESRecHitCollection > esHandle;
      // edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esHandle;

      const EcalPreshowerGeometry *esGeometry_;     
      const CaloGeometry* geometry;
      bool GeometryFromFile_;

      std::string outfilename_;
      std::string externalGeometry_;
      std::string calibMapPath_; 
      std::string jsonFile_; 
      std::string ebContainmentCorrections_;
      std::string MVAEBContainmentCorrections_01_;
      std::string MVAEBContainmentCorrections_02_;
      std::string MVAEEContainmentCorrections_01_;
      std::string MVAEEContainmentCorrections_02_;
      std::string MVAEBContainmentCorrections_eta01_;
      std::string MVAEBContainmentCorrections_eta02_;
      std::string Endc_x_y_;
      bool        EtaRingCalibEB_;
      bool        SMCalibEB_;
      bool        EtaRingCalibEE_;
      bool        SMCalibEE_;
      std::string CalibMapEtaRing_;
      std::string ebPHIContainmentCorrections_;
      std::string eeContainmentCorrections_;
      std::string Barrel_orEndcap_;
      bool useEBContainmentCorrections_;
      bool useEEContainmentCorrections_;
      bool useOnlyEEClusterMatchedWithES_;
      bool HLTResults_;
      std::string HLTResultsNameEB_;
      std::string HLTResultsNameEE_;
      bool RemoveDead_Flag_;
      TString RemoveDead_Map_;
      TString L1_Bit_Sele_;
      //float L1BitCollection_[NL1SEED];

      bool Are_pi0_;
      bool useMVAContainmentCorrections_;
      bool new_pi0ContainmentCorrections_;

      bool L1TriggerInfo_;
      edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_;
      edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_;
      edm::EDGetTokenT<ESRecHitCollection> ESRecHitCollectionToken_;
      edm::InputTag l1TriggerTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
      edm::EDGetTokenT<GlobalAlgBlkBxCollection> L1GTobjmapToken_;
      edm::InputTag l1InputTag_;
      //std::map<string,int> L1_nameAndNumb;
      edm::EDGetTokenT<GenParticleCollection> GenPartCollectionToken_;

      edm::EDGetTokenT<edm::SimTrackContainer>  g4_simTk_Token_;
      edm::EDGetTokenT<edm::SimVertexContainer> g4_simVtx_Token_;      

      PosCalcParams PCparams_;
      //const double preshowerStartEta_ =  1.653;

      ECALGeometry* geom_;
      CaloTopology *ebtopology_;
      CaloTopology *eetopology_;
      CaloSubdetectorTopology *estopology_;
      
      std::string calibTypeString_;
      calibGranularity calibTypeNumber_;

      // Seeds
      double EB_Seed_E_;
      bool useEE_EtSeed_;
      double EE_Seed_E_;
      double EE_Seed_Et_;
      // selection criteria
      double gPtCut_low_[3];
      double gPtCut_high_[3];
      double pi0PtCut_low_[3];
      double pi0PtCut_high_[3];

      double pi0IsoCut_low_[3];
      double pi0IsoCut_high_[3];
      bool   CutOnHLTIso_;
      double pi0HLTIsoCut_low_[3];
      double pi0HLTIsoCut_high_[3];

      double nXtal_1_cut_low_[3];
      double nXtal_1_cut_high_[3];
      double nXtal_2_cut_low_[3];
      double nXtal_2_cut_high_[3];
      double S4S9_cut_low_[3];
      double S4S9_cut_high_[3];
      double SystOrNot_;
      bool isMC_;
      bool MC_Assoc_;
      double MC_Assoc_DeltaR;
      math::XYZPoint Gamma1MC;
      math::XYZPoint Gamma2MC;
      bool isCRAB_;
      bool MakeNtuple4optimization_;
      bool isDebug_; 
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

      std::vector<int> Ncristal_EE, Ncristal_EE_used;
      std::vector<int> Ncristal_EB, Ncristal_EB_used;

      TH1F *EventFlow_EB;
      TH1F *EventFlow_EE;
      TH1F *EventFlow_EB_debug;
      TH1F *EventFlow_EE_debug;
      TH1F **epsilon_EB_h;  // epsilon distribution by region
      TH1F **epsilon_EE_h;  // epsilon distribution in EE
      TH1F *allEpsilon_EE; 
      TH1F *allEpsilon_EEnw; 
      TH1F *allEpsilon_EB;
      TH1F *allEpsilon_EBnw;
      TH2F *entries_EEp;
      TH2F *entries_EEm;
      TH2F *entries_EB;
      TH2F *Occupancy_EEp;
      TH2F *Occupancy_EEm;
      TH2F *Occupancy_EB;
      TH2F *pi0MassVsIetaEB;
      TH2F *pi0MassVsETEB;
      bool useMassInsteadOfEpsilon_;

#ifdef SELECTION_TREE
      Float_t NSeeds_EB, Xclus_EB, Yclus_EB, Zclus_EB, e3x3_EB, S4S9_EB, PTClus_EB;
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
      //adding these variables
      Float_t PtGamma1_EB, PtGamma2_EB, EtaGamma1_EB, EtaGamma2_EB, NxtalGamma1_EB, NxtalGamma2_EB, S4S9Gamma1_EB, S4S9Gamma2_EB;
      void Fill_PtPi0_EB(float pt){ PtPi0_EB=pt; };
      void Fill_mpi0_EB(float m){ mpi0_EB=m; };
      void Fill_etapi0_EB( float eta){ Etapi0_EB =eta; };
      void Fill_phipi0_EB( float phi){ Phipi0_EB =phi; };
      //adding these method to fill variables
      void Fill_PtGamma_EB(float pt1, float pt2){ PtGamma1_EB = pt1; PtGamma2_EB = pt2;};
      void Fill_EtaGamma_EB(float eta1, float eta2){ EtaGamma1_EB = eta1; EtaGamma2_EB = eta2;};
      void Fill_NcrystalUsedGamma_EB(float Nxtal1, float Nxtal2) { NxtalGamma1_EB = Nxtal1 ; NxtalGamma2_EB = Nxtal2;};
      void Fill_S4S9Gamma_EB(float s4s9g1, float s4s9g2) { S4S9Gamma1_EB = s4s9g1; S4S9Gamma2_EB = s4s9g2;};
      //
      void Fill_Epsilon_EB(float eps ){ Epsilon_EB=eps; };
      TTree *Pi0Info_EB;

      Float_t PtPi0_EE, mpi0_EE,Etapi0_EE, Phipi0_EE, Epsilon_EE;
      //adding these variables
      Float_t PtGamma1_EE, PtGamma2_EE, EtaGamma1_EE, EtaGamma2_EE, NxtalGamma1_EE, NxtalGamma2_EE, S4S9Gamma1_EE, S4S9Gamma2_EE;
      void Fill_PtPi0_EE(float pt){ PtPi0_EE=pt; };
      void Fill_mpi0_EE(float m){ mpi0_EE=m; };
      void Fill_etapi0_EE( float eta){ Etapi0_EE =eta; };
      void Fill_phipi0_EE( float phi){ Phipi0_EE =phi; };
      //adding these methods to fill variables
      void Fill_PtGamma_EE(float pt1, float pt2){ PtGamma1_EE = pt1; PtGamma2_EE = pt2;};
      void Fill_EtaGamma_EE(float eta1, float eta2){ EtaGamma1_EE = eta1; EtaGamma2_EE = eta2;};
      void Fill_NcrystalUsedGamma_EE(float Nxtal1, float Nxtal2) { NxtalGamma1_EE = Nxtal1 ; NxtalGamma2_EE = Nxtal2;};
      void Fill_S4S9Gamma_EE(float s4s9g1, float s4s9g2) { S4S9Gamma1_EE = s4s9g1; S4S9Gamma2_EE = s4s9g2;};
      //
      void Fill_Epsilon_EE(float eps ){ Epsilon_EE=eps; };
      TTree *Pi0Info_EE;
#endif
      TTree*  Tree_Optim;
      Int_t   nPi0;
      //Int_t   Op_L1Seed[NL1SEED];
      Int_t   Op_NPi0_rec;
      Int_t   Op_Pi0recIsEB[NPI0MAX];
      Float_t Op_IsoPi0_rec[NPI0MAX];
      Float_t Op_HLTIsoPi0_rec[NPI0MAX];
      Int_t   Op_n1CrisPi0_rec[NPI0MAX];
      Int_t   Op_n2CrisPi0_rec[NPI0MAX];
      Float_t Op_mPi0_rec[NPI0MAX];
      Float_t Op_enG1_rec[NPI0MAX];
      Float_t Op_enG2_rec[NPI0MAX];
      Float_t Op_etaPi0_rec[NPI0MAX];
      Float_t Op_ptPi0_rec[NPI0MAX];
      Float_t Op_DeltaRG1G2[NPI0MAX];
      Float_t Op_Es_e1_1[NPI0MAX];
      Float_t Op_Es_e1_2[NPI0MAX];
      Float_t Op_Es_e2_1[NPI0MAX];
      Float_t Op_Es_e2_2[NPI0MAX];
      Float_t Op_S4S9_1[NPI0MAX];
      Float_t Op_S4S9_2[NPI0MAX];
      Float_t Op_S1S9_1[NPI0MAX];
      Float_t Op_S1S9_2[NPI0MAX];
      Float_t Op_S2S9_1[NPI0MAX];
      Float_t Op_S2S9_2[NPI0MAX];
      Float_t Op_Eta_1[NPI0MAX];
      Float_t Op_Eta_2[NPI0MAX];
      Float_t Op_Phi_1[NPI0MAX];
      Float_t Op_Phi_2[NPI0MAX];
      Float_t Op_Time_1[NPI0MAX];
      Float_t Op_Time_2[NPI0MAX];
      Float_t Op_DeltaR_1[NPI0MAX];
      Float_t Op_DeltaR_2[NPI0MAX];
      Float_t Op_enG1_nocor[NPI0MAX];
      Float_t Op_enG2_nocor[NPI0MAX];
      Float_t Op_ptPi0_nocor[NPI0MAX];
      Float_t Op_mPi0_nocor[NPI0MAX];
      Float_t Op_enG1_true[NPI0MAX];
      Float_t Op_enG2_true[NPI0MAX];
      Int_t Op_Nxtal_1[NPI0MAX];
      Int_t Op_Nxtal_2[NPI0MAX];
      Int_t Op_iEtaiX_1[NPI0MAX];
      Int_t Op_iEtaiX_2[NPI0MAX];
      Int_t Op_iPhiiY_1[NPI0MAX];
      Int_t Op_iPhiiY_2[NPI0MAX];
      Int_t Op_iEta_1on5[NPI0MAX];
      Int_t Op_iEta_2on5[NPI0MAX];
      Int_t Op_iPhi_1on2[NPI0MAX];
      Int_t Op_iPhi_2on2[NPI0MAX];
      Int_t Op_iEta_1on2520[NPI0MAX];
      Int_t Op_iEta_2on2520[NPI0MAX];
      Int_t Op_iPhi_1on20[NPI0MAX];
      Int_t Op_iPhi_2on20[NPI0MAX];

      vector<float> Es_1;
      vector<float> Es_2;

      std::string ContCorr_EB_;
      TH1F *triggerComposition;      
      TH1F *triggerComposition_EB; // require that HLT in EB fired
      TH1F *triggerComposition_EE; // require that HLT in EE fired
      bool areLabelsSet_;

      std::map< std::string, int > l1TrigNames_;
      bool l1TrigBit_[128];
      vector<float> vs4s9;
      vector<float> vs1s9;
      vector<float> vs2s9;
      TFile *EBweight_file_1;
      TFile *EBweight_file_2;
      const GBRForest *forest_EB_1;
      const GBRForestD *forestD_EB_1;
      const GBRForest *forest_EB_2;
      const GBRForestD *forestD_EB_2;
      GBRApply *gbrapply;
#if defined(MVA_REGRESSIO_Tree) && defined(MVA_REGRESSIO)
      TTree *TTree_JoshMva;
      Float_t Correction1_mva, Correction2_mva, Pt1_mva, Pt2_mva, Mass_mva, MassOr_mva, pi0Eta;
      Int_t   iEta1_mva, iPhi1_mva, iEta2_mva, iPhi2_mva, iSM1_mva, iSM2_mva;
#endif
      vector<iXiYtoRing> VectRing;
      std::map<int,vector<int>> ListEtaFix_xtalEB;
      std::map<int,vector<int>> ListSMFix_xtalEB;
      std::map<int,vector<int>> ListEtaFix_xtalEEm;
      std::map<int,vector<int>> ListEtaFix_xtalEEp;
      std::map<int,vector<int>> ListQuadFix_xtalEEm;
      std::map<int,vector<int>> ListQuadFix_xtalEEp;
      std::map<int,vector<int>> List_IR_EtaPhi;
      std::map<int,vector<int>> List_IR_XYZ;
      vector<float> vs4s9EE;
      vector<float> vSeedTime;
      vector<float> vSeedTimeEE;
      vector<float> vs1s9EE;
      vector<float> vs2s9EE;
      vector<float> ESratio;
#ifdef MVA_REGRESSIO_EE
      TFile *EEweight_file_pi01;
      TFile *EEweight_file_pi02;
      const GBRForest *forest_EE_pi01;
      const GBRForestD *forestD_EE_pi01;
      const GBRForest *forest_EE_pi02;
      const GBRForestD *forestD_EE_pi02;

      TTree *TTree_JoshMva_EE;
      Float_t Correction1EE_mva, Correction2EE_mva, Pt1EE_mva, Pt2EE_mva, MassEE_mva, MassEEOr_mva;
      Int_t   iX1_mva, iY1_mva, iX2_mva, iY2_mva, EtaRing1_mva, EtaRing2_mva;
#endif
      //JSON
      std::string JSONfile_;
      JSON* myjson;
      int Num_Fail_Sel;
      int Num_Fail_tot;
      TH1F *Selec_Efficiency;
      //Preselection
      //int Num_Fail_Presel;
      //bool FailPreselEB;
      //bool FailPreselEE;
      //std::map<int,bool>  PassPreselection;
      
      //Containment correction
      /*constexpr*/ double meanlimlow  = 0.2;
      /*constexpri*/ double meanlimhigh = 2.0;
      /*constexpr*/ double meanoffset  = meanlimlow + 0.5*(meanlimhigh-meanlimlow);
      /*constexpr*/ double meanscale   = 0.5*(meanlimhigh-meanlimlow);

      // for L1
      short *l1flag;
      TString* algoBitToName;
      std::string L1SeedsPi0Stream_;
      int nL1SeedsPi0Stream_; // number of seeds used by the stream (given L1SeedsPi0Stream_, it is the number of " OR " +1, e.g. "seed1 OR seed2 OR seed3" has 3 seeds
      int *seedIsInStream;

      // store for each event if AlCa_EcalPi0(Eta)EB(EE)only_v* fired
      bool EB_HLT, EE_HLT;

      // event info
      ULong64_t myEvent;
      int myLumiBlock;
      int myRun;
      int myBunchCrossing;

};
