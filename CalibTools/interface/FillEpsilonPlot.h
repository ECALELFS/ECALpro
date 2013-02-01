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

#include "CalibCode/FillEpsilonPlot/interface/PosCalcParams.h"
#include "CalibCode/FillEpsilonPlot/interface/ECALGeometry.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalEnerCorr.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalCalibTypes.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalRegionalCalibration.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalPreshowerHardcodedGeometry.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalPreshowerHardcodedTopology.h"

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
      void fillEEClusters(std::vector< CaloCluster > & eseeclusters, const edm::Event& iEvent);
      void computeEpsilon(std::vector< CaloCluster > & clusters, int subDetId);
      float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
      float DeltaPhi(float phi1, float phi2);

      TH1F** initializeEpsilonHistograms(const char *name, const char *title, int size );
      void  deleteEpsilonPlot(TH1F **h, int size);
      void  writeEpsilonPlot(TH1F **h, const char *folder, int size);
      bool getTriggerResult(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool getTriggerByName( std::string s );
      // ----------member data ---------------------------
      
      std::string outfilename_;
      std::string ebContainmentCorrections_;
      std::string eeContainmentCorrections_;
      std::string externalGeometry_;

      edm::InputTag EBRecHitCollectionTag_;
      edm::InputTag EERecHitCollectionTag_;
      edm::InputTag ESRecHitCollectionTag_;
      edm::InputTag l1TriggerTag_,;

      PosCalcParams PCparams_;

      ECALGeometry* geom_;
      CaloTopology *ebtopology_;
      CaloTopology *eetopology_;
      CaloSubdetectorTopology *estopology_;
      EcalPreshowerHardcodedGeometry* hardcodedPreshowerGeometry_;
      
      EcalEnerCorr energycorr;
      bool noCorrections_;

      std::string calibTypeString_;
      calibGranularity calibTypeNumber_;

      // selection criteria
      double pi0PtCut_[3];
      double pi0IsoCut_[3];

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

      TH1F **epsilon_EB_h;  // epsilon distribution by region
      TH1F **epsilon_EE_h;  // epsilon distribution in EE
      TH1F *allEpsilon_EE; // debug
      TH2F *entries_EEp;
      TH2F *entries_EEm;
      TH2F *pi0MassVsIetaEB;
      TH2F *pi0MassVsETEB;

      TH1F *triggerComposition;
      bool areLabelsSet_;

      std::vector<std::string> alcaL1TrigNames_;
      //std::vector<std::string> l1TrigNames_(128,"");
      std::map< std::string, int > l1TrigNames_;
      bool l1TrigBit_[128];
};
