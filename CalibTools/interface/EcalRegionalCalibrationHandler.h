#include "CalibCode/FillEpsilonPlot/interface/EcalCalibTypes.h"
#include "CalibCode/FillEpsilonPlot/interface/EcalRegionalCalibration.h"

enum calibGranularity{ xtal, tt, etaring };

class EcalRegionalCalibrationHandler{
    public:
        EcalRegionalCalibrationHandler();
        ~EcalRegionalCalibrationHandler();
        void setCa
    private:
        EcalRegionalCalibration<EcalCalibType::Xtal> *xtalCalib_;
        EcalRegionalCalibration<EcalCalibType::EtaRing> *etaCalib_;
        EcalRegionalCalibration<EcalCalibType::TrigTower> *TTCalib_;
        calibGranularity calibType_;
        bool isSet_;
};

EcalRegionalCalibrationHandler::EcalRegionalCalibrationHandler() {
    xtalCalib_ = new  EcalRegionalCalibration<EcalCalibType::Xtal>;
    etaCalib_ = new EcalRegionalCalibration<EcalCalibType::EtaRing>;
    TTCalib_ = new EcalRegionalCalibration<EcalCalibType::TrigTower>;
    isSet_ = false;
    calibGranularity = xtal;
}

EcalRegionalCalibrationHandler::~EcalRegionalCalibrationHandler() {
    delete xtalCalib_;
    delete etaCalib_;
    delete TTCalib_;
}
