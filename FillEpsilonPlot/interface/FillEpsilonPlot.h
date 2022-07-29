#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "CalibCode/CalibTools/interface/PosCalcParams.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"
#include "CalibCode/CalibTools/interface/EcalEnerCorr.h"
#include "CalibCode/CalibTools/interface/EndcapTools.h"
#include "CalibCode/CalibTools/interface/EcalCalibTypes.h"
#include "CalibCode/CalibTools/interface/EcalRegionalCalibration.h"
#include "CalibCode/CalibTools/interface/EcalPreshowerHardcodedTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
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

//#include "CalibCode/GBRTrain/interface/GBRApply.h"
//#include "CalibCode/EgammaObjects/interface/GBRForest.h"
//#include "CondFormats/EgammaObjects/interface/GBRForestD.h"

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

class FillEpsilonPlot : public edm::one::EDAnalyzer<> {
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
    //std::vector< CaloCluster > MCTruthAssociate(std::vector< CaloCluster > & clusters, double deltaR, bool isEB);
    std::vector< CaloCluster > MCTruthAssociateMultiPi0(std::vector< CaloCluster > & clusters, int& retNumberUnmergedGen, int& retNumberMatchedGen, std::vector<TLorentzVector*>& retClusters_matchedGenPhotonEnergy, const double deltaR, const bool isEB);
    // void computePairProperties(std::vector<CaloCluster>::const_iterator g1, std::vector<CaloCluster>::const_iterator g2, math::XYZVector &tmp_photon1, math::XYZVector &tmp_photon2, float &m_pair, float &pt_pair, float &eta_pair, float &phi_pair);
    void computePairProperties(const CaloCluster* g1, const CaloCluster* g2, math::XYZVector &tmp_photon1, math::XYZVector &tmp_photon2, float &m_pair, float &pt_pair, float &eta_pair, float &phi_pair);
    void computeEpsilon(std::vector< CaloCluster > & clusters, std::vector<TLorentzVector*>& clusters_matchedGenPhoton, int subDetId);
    void computeEoverEtrue(std::vector< CaloCluster > & clusters, std::vector<TLorentzVector*>& clusters_matchedGenPhoton, int subDetId);
    bool checkStatusOfEcalRecHit(const EcalChannelStatus &channelStatus,const EcalRecHit &rh);
    bool isInDeadMap( bool isEB, const EcalRecHit &rh );
    float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
    float DeltaPhi(float phi1, float phi2);
    double min( double a, double b);
    int getNumberOverlappingCrystals(std::vector<CaloCluster>::const_iterator g1, std::vector<CaloCluster>::const_iterator g2, const bool isEB);


    TH2F* initializeEpsilonHistograms2D(const char *name, const char *title, int size );
    void deleteEpsilonPlot2D(TH2F *h);
    void writeEpsilonPlot2D(TH2F *h);

    bool getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    //bool getTriggerByName( std::string s ); not used anymore
    bool GetHLTResults(const edm::Event& iEvent, std::string s);
     
    TFile* DeadMap;
    TH2F * EBMap_DeadXtal;
    TH2F * EEmMap_DeadXtal;
    TH2F * EEpMap_DeadXtal;

    // for containment corrections based on E/Etrue in MC
    TH2F* hCC_EoverEtrue_g1 = nullptr;
    TH2F* hCC_EoverEtrue_g2 = nullptr;
    void loadEoverEtrueContainmentCorrections(const string& fileName);
    CaloCluster getClusterAfterContainmentCorrections(std::vector<CaloCluster>::const_iterator, const bool isSecondPhoton, const bool isEB);

    // ----------member data ---------------------------
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geoToken_;
    edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> chStatusToken_;    
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
    std::string fileEoverEtrueContainmentCorrections_;
    double scalingEoverEtrueCC_g1_;
    double scalingEoverEtrueCC_g2_;
    std::string Endc_x_y_;
    bool        EtaRingCalibEB_;
    bool        SMCalibEB_;
    bool        EtaRingCalibEE_;
    bool        SMCalibEE_;
    std::string CalibMapEtaRing_;
    std::string Barrel_orEndcap_;
    bool useContainmentCorrectionsFromEoverEtrue_;
    bool useOnlyEEClusterMatchedWithES_;
    bool HLTResults_;
    std::string HLTResultsNameEB_;
    std::string HLTResultsNameEE_;
    bool RemoveDead_Flag_;
    TString RemoveDead_Map_;
    bool RemoveSeedsCloseToDeadXtal_;
    TString L1_Bit_Sele_;
    //float L1BitCollection_[NL1SEED];

    bool Are_pi0_;

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
    edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileupSummaryToken_;
    edm::Handle<reco::GenParticleCollection> genParticles;

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
    std::vector<double> gPtCut_;
    double pi0PtCut_low_[3];
    double pi0PtCut_high_[3];
    std::vector<double> pi0PtCut_;

    double pi0IsoCut_low_[3];
    double pi0IsoCut_high_[3];
    std::vector<double> pi0IsoCut_;
    bool   CutOnHLTIso_;
    double pi0HLTIsoCut_low_[3];
    double pi0HLTIsoCut_high_[3];
    std::vector<double> pi0HLTIsoCut_;

    int nXtal_1_cut_low_[3];
    int nXtal_1_cut_high_[3];
    int nXtal_2_cut_low_[3];
    int nXtal_2_cut_high_[3];
    std::vector<int> nXtal_1_cut_;
    std::vector<int> nXtal_2_cut_;
    double S4S9_cut_low_[3];
    double S4S9_cut_high_[3];
    std::vector<double> S4S9_cut_;
    int SystOrNot_;

    // MC stuff
    bool isMC_;
    bool MC_Assoc_;
    double MC_Assoc_DeltaR;
    bool isEoverEtrue_;  // to run E/Etrue flow (with MC only)
    math::XYZPoint Gamma1MC;
    math::XYZPoint Gamma2MC;
    // for MC truth with more pi0
    // start with TLorentzVector, then I will optimize
    /* vector<math::XYZPoint> vecGamma1MC; */
    /* vector<math::XYZPoint> vecGamma2MC; */
    /* vector<TLorentzVector> vecGamma1MC; */
    /* vector<TLorentzVector> vecGamma2MC; */
    vector<TLorentzVector> vecGamma1MC_EB;
    vector<TLorentzVector> vecGamma2MC_EB;
    vector<TLorentzVector> vecGamma1MC_EE;
    vector<TLorentzVector> vecGamma2MC_EE;
    std::vector< TLorentzVector* > ebclusters_matchedGenPhoton;  // will store the gen photon corresponding to a given reco cluster (ordered pairs with seed energy)   
    std::vector< TLorentzVector* > eeclusters_matchedGenPhoton;  // will store the gen photon corresponding to a given reco cluster (ordered pairs with seed energy)
    TH1F* h_numberUnmergedGenPhotonPairs_EB; // fraction of gen photon pairs that are not merged (i.e. the photons are separated by a DR defined in .cc)
    TH1F* h_numberMatchedGenPhotonPairs_EB;  // fraction of gen photon pairs that are succesfully matched to reco clusters
    TH1F* h_numberUnmergedGenPhotonPairs_EE; 
    TH1F* h_numberMatchedGenPhotonPairs_EE;  
    TH1F* h_numberUnmergedGenPhotonPairs; // absolute number without separating EB and EE 
    TH1F* h_numberMatchedGenPhotonPairs;  
    TH1F* g1RecoGenDR_EB;
    TH1F* g2RecoGenDR_EB;
    TH1F* diff_g2Recog1GenDR_g2RecoGenDR_EB;
    // for E/Etrue with MC
    TH1F **EoverEtrue_g1_EB_h;  
    TH1F **EoverEtrue_g1_EE_h;  
    TH1F **EoverEtrue_g2_EB_h;  
    TH1F **EoverEtrue_g2_EE_h;  
    TH2F *EoverEtrue_g1_EB_h2D;  
    TH2F *EoverEtrue_g2_EB_h2D;  
    TH2F *EoverEtrue_g1_EE_h2D;  
    TH2F *EoverEtrue_g2_EE_h2D;  
    /* TH1F *allEoverEtrue_g1_EE;  */
    /* TH1F *allEoverEtrue_g1_EEnw;  */
    /* TH1F *allEoverEtrue_g1_EB; */
    /* TH1F *allEoverEtrue_g1_EBnw; */
    /* TH1F *allEoverEtrue_g2_EE;  */
    /* TH1F *allEoverEtrue_g2_EEnw;  */
    /* TH1F *allEoverEtrue_g2_EB; */
    /* TH1F *allEoverEtrue_g2_EBnw; */

    // Some kinematic variables (use option in parameters.py to choose whether to fill and save them)
    bool fillKinematicVariables_;
    int whichRegionEcalStreamPi0; // will be used to say in which region we are based on eta of pi0
    TH2F* seedEnergyInCluster;
    TH2F* pi0pt_afterCuts;  // 5 regions (2 in EB and 3 in EE, last 2 in EE could be merged)
    TH2F* g1pt_afterCuts;
    TH2F* g2pt_afterCuts;
    TH2F* g1Nxtal_afterCuts;
    TH2F* g2Nxtal_afterCuts;
    TH2F* pi0PhotonsNoverlappingXtals_afterCuts;
    TH2F* g1g2DR_afterCuts;
    std::vector<TH2F*> pi0MassVsPU;  // BX 0
    //std::vector<TH2F*> pi0MassVsPU_BXm1;
    //std::vector<TH2F*> pi0MassVsPU_BXm2;
    //std::vector<TH2F*> pi0MassVsPU_BXp1;

    /////////
    bool isCRAB_;
    bool MakeNtuple4optimization_;
    bool isDebug_;
    /// all the three options have to be instantiated to allow the
    //choice at runtime
    EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib;
    EcalRegionalCalibration<EcalCalibType::EtaRing> etaCalib;
    EcalRegionalCalibration<EcalCalibType::TrigTower> TTCalib;
    EcalRegionalCalibrationBase *regionalCalibration_;  // use it for pi0 mass or first photon with E/overEtrue

    // for second photon with E/Etrue (MC only)
    // I create them with "regionalCalibration_g2_ = new EcalRegionalCalibration<EcalCalibType::Xtal>()" directly in the source 
    /* EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib_g2; */
    /* EcalRegionalCalibration<EcalCalibType::EtaRing> etaCalib_g2; */
    /* EcalRegionalCalibration<EcalCalibType::TrigTower> TTCalib_g2; */
    EcalRegionalCalibrationBase *regionalCalibration_g2_; // use it for second gen photon with E/Etrue

    int currentIteration_;
    string outputDir_;

    TFile *outfile_;
    TFile *externalGeometryFile_;

    std::vector<int> Ncristal_EE, Ncristal_EE_used;
    std::vector<int> Ncristal_EB, Ncristal_EB_used;

    TH1F **epsilon_EB_h;  // epsilon distribution by region
    TH1F **epsilon_EE_h;  // epsilon distribution in EE
    TH2F *epsilon_EB_h2D;  // epsilon distribution by region
    TH2F *epsilon_EE_h2D;  // epsilon distribution by region
    TH2F *pi0MassVsIetaEB;
    TH2F *pi0MassVsETEB;
    TH2F *photonDeltaRVsIetaEB;

    ///SJ
    TH1F *pi0_mass_EB;
    TH1F *pi0_mass_EEP;
    TH1F *pi0_mass_EEM;

    bool useMassInsteadOfEpsilon_;

    TTree*  Tree_Optim;
    Int_t   nPi0;
    //Int_t   Op_L1Seed[NL1SEED];
    Int_t   Op_NPi0;
    std::vector<Int_t> Op_Pi0recIsEB;
    std::vector<Float_t> Op_ClusIsoPi0;
    std::vector<Float_t> Op_HLTIsoPi0;
    std::vector<Int_t> Op_nCrisG1;
    std::vector<Int_t> Op_nCrisG2;
    std::vector<Float_t> Op_enG1_cor;
    std::vector<Float_t> Op_enG2_cor;
    std::vector<Float_t> Op_etaG1_cor;
    std::vector<Float_t> Op_etaG2_cor;
    std::vector<Float_t> Op_phiG1_cor;
    std::vector<Float_t> Op_phiG2_cor;
    std::vector<Float_t> Op_mPi0_cor;
    std::vector<Float_t> Op_etaPi0_cor;
    std::vector<Float_t> Op_ptPi0_cor;
    std::vector<Float_t> Op_phiPi0_cor;
    std::vector<Float_t> Op_DeltaRG1G2;
    std::vector<Float_t> Op_Es_e1_1;
    std::vector<Float_t> Op_Es_e1_2;
    std::vector<Float_t> Op_Es_e2_1;
    std::vector<Float_t> Op_Es_e2_2;
    std::vector<Float_t> Op_S4S9_1;
    std::vector<Float_t> Op_S4S9_2;
    std::vector<Float_t> Op_S1S9_1;
    std::vector<Float_t> Op_S1S9_2;
    std::vector<Float_t> Op_S2S9_1;
    std::vector<Float_t> Op_S2S9_2;
    std::vector<Float_t> Op_Time_1;
    std::vector<Float_t> Op_Time_2;
    std::vector<Float_t> Op_DeltaR_1;
    std::vector<Float_t> Op_DeltaR_2;
    std::vector<Float_t> Op_enG1_nocor;
    std::vector<Float_t> Op_enG2_nocor;
    std::vector<Float_t> Op_etaG1_nocor;
    std::vector<Float_t> Op_etaG2_nocor;
    std::vector<Float_t> Op_phiG1_nocor;
    std::vector<Float_t> Op_phiG2_nocor;
    std::vector<Float_t> Op_ptPi0_nocor;
    std::vector<Float_t> Op_etaPi0_nocor;
    std::vector<Float_t> Op_phiPi0_nocor;
    std::vector<Float_t> Op_mPi0_nocor;
    std::vector<Float_t> Op_enG1_true;
    std::vector<Float_t> Op_enG2_true;
    std::vector<Int_t> Op_iEtaiX_1;
    std::vector<Int_t> Op_iEtaiX_2;
    std::vector<Int_t> Op_iPhiiY_1;
    std::vector<Int_t> Op_iPhiiY_2;
    std::vector<Int_t> Op_iEta_1on5;
    std::vector<Int_t> Op_iEta_2on5;
    std::vector<Int_t> Op_iPhi_1on2;
    std::vector<Int_t> Op_iPhi_2on2;
    std::vector<Int_t> Op_iEta_1on2520;
    std::vector<Int_t> Op_iEta_2on2520;
    std::vector<Int_t> Op_iPhi_1on20;
    std::vector<Int_t> Op_iPhi_2on20;
    // Optmization tree's variables


    vector<float> Es_1;
    vector<float> Es_2;

    std::string ContCorr_EB_;
    TH1F *triggerComposition;      
    TH1F *triggerComposition_EB; // require that HLT in EB fired
    TH1F *triggerComposition_EE; // require that HLT in EE fired
    bool areLabelsSet_;

    /* std::map< std::string, int > l1TrigNames_; */
    /* bool l1TrigBit_[128]; */
    vector<float> vs4s9;
    vector<float> vs1s9;
    vector<float> vs2s9;
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

    // PU info for MC
    Float_t nPUtrue_;
    std::map<Int_t,Int_t> nPUobs_;
    Int_t nPUobs_BX0_;

    //variables for monitoring tree
    TTree*  tree_mon;
    int event_year, event_month, event_day, isPi0EB;
    double event_time, pi0_mass, pho1_eta, pho2_eta;

};
