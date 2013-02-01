#ifndef _ECALGeometry_
#define _ECALGeometry_


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <map>

class TFile;


class ECALGeometry {
public:
  static ECALGeometry* getGeometry(TFile* f=0);

  GlobalPoint getPosition(DetId id, float depth=0.) const;
  GlobalVector   getAxis(DetId id) const;

  const std::map<DetId,GlobalPoint> & getPositionMap() const { return posMap;}
  const std::map<DetId,GlobalVector> & getAxisMap() const { return axisMap;}

  void print() const;



private:
  ECALGeometry() { }
  ECALGeometry(TFile* f);
 
  std::map<DetId,GlobalPoint> posMap;
  std::map<DetId,GlobalVector> axisMap;
  static ECALGeometry* instance;

};
#endif
