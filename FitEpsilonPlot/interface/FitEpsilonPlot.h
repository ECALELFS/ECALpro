#include <memory>

#include "RooRealVar.h"
#include "RooFitResult.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibCode/CalibTools/interface/EcalRegionalCalibration.h"
#include "CalibCode/CalibTools/interface/EcalCalibTypes.h"


enum calibGranularity{ xtal, tt, etaring };

struct Pi0FitResult {
   RooFitResult* res;
   float S;     // signal in 3 sigma region
   float Serr;
   float B;     // bkg in 3 sigma region
   float Berr; 
   float  SoB;  // S/B
   float  SoBerr;
   float chi2;
   int dof;
   float probchi2;
};

class FitEpsilonPlot : public edm::EDAnalyzer {
   public:
      enum FitMode{ Eta=0, Pt, GausPol3, GausEndpoint, Pi0EB, Pi0EE, EtaEB };
      explicit FitEpsilonPlot(const edm::ParameterSet&);
      ~FitEpsilonPlot();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void loadEpsilonPlot(char *filename);
      void saveCoefficients();
      void IterativeFit(TH1F* h, TF1 & ffit); 
      void deleteEpsilonPlot(TH1F **h, int size);

      Pi0FitResult FitMassPeakRooFit(TH1F* h,double xlo, double xhi, uint32_t HistoIndex, int ngaus=1, FitMode mode=Pi0EB, int niter=0, bool isNot_2010_=true);

      // ----------member data ---------------------------

      EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib;
      EcalRegionalCalibration<EcalCalibType::EtaRing> etaCalib;
      EcalRegionalCalibration<EcalCalibType::TrigTower> TTCalib;

      EcalRegionalCalibrationBase *regionalCalibration_;

      int currentIteration_;
      std::string outputDir_;
      std::string outfilename_;
      std::string calibTypeString_;
      std::string epsilonPlotFileName_;
      std::string calibMapPath_; 
      std::string Barrel_orEndcap_; 

      std::string EEoEB_; 
      bool isNot_2010_; 
      bool Are_pi0_; 
      bool StoreForTest_; 
      int inRangeFit_; 
      int finRangeFit_; 

      calibGranularity calibTypeNumber_;

      TH1F **epsilon_EB_h;  // epsilon distribution by region
      TH1F **epsilon_EE_h;  // epsilon distribution in EE

      TFile *inputEpsilonFile_;
      TFile *outfile_;
      TFile *outfileTEST_;

      bool useMassInsteadOfEpsilon_;

      std::map<int,float> EBmap_Signal;//#
      std::map<int,float> EBmap_Backgr;
      std::map<int,float> EBmap_Chisqu;
      std::map<int,float> EBmap_ndof;
      std::map<int,float> EBmap_mean;
      std::map<int,float> EBmap_mean_err;
      std::map<int,float> EBmap_sigma;
      std::map<int,float> EBmap_Snorm;
      std::map<int,float> EBmap_b0;
      std::map<int,float> EBmap_b1;
      std::map<int,float> EBmap_b2;
      std::map<int,float> EBmap_b3;
      std::map<int,float> EBmap_Bnorm;

      std::map<int,float> EEmap_Signal;
      std::map<int,float> EEmap_Backgr;
      std::map<int,float> EEmap_Chisqu;
      std::map<int,float> EEmap_ndof;
      std::map<int,float> EEmap_mean;
      std::map<int,float> EEmap_mean_err;
      std::map<int,float> EEmap_sigma;
      std::map<int,float> EEmap_Snorm;
      std::map<int,float> EEmap_b0;
      std::map<int,float> EEmap_b1;
      std::map<int,float> EEmap_b2;
      std::map<int,float> EEmap_b3;
      std::map<int,float> EEmap_Bnorm;





};



