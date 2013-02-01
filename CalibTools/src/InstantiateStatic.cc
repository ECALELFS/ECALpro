#include <string>

#include "CalibCode/CalibTools/interface/EndcapTools.h"
#include "CalibCode/CalibTools/interface/GeometryService.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"
// #include "CalibCode/FillEpsilonPlot/interface/EcalCalibMap.h"

using std::string;

bool GeometryService::isNameSet_ = false;
bool GeometryService::isPtrSet_ = false;
string GeometryService::geometryFileName_;
ECALGeometry* GeometryService::geometryPtr_ = 0;

bool EndcapTools::isInitializedFromGeometry_ = false;
int EndcapTools::endcapRingIndex_[EEDetId::IX_MAX][EEDetId::IY_MAX]; 
ECALGeometry* EndcapTools::caloGeometry_ = 0;
TFile* EndcapTools::externalGeometryFile_ = 0;


// template<typename Type> float EcalCalibMap<Type>::mapEB[Type::nRegions];
// template<typename Type> float EcalCalibMap<Type>::mapEE[Type::nRegionsEE];
// template<typename Type> float EcalCalibMap<Type>::bad_coeff;
