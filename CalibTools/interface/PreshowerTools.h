#ifndef PreshowerTools_H
#define PreshowerTools_H

#include <set>
#include <map>
#include <vector>
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Analysis/Modules/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"

//using namespace reco;

typedef std::map<DetId, EcalRecHit> RecHitsMap;
typedef std::set<DetId> HitsID;

using namespace reco;

class PreshowerTools{
    public:
      //PreshowerTools(CaloSubdetectorGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle); 
      PreshowerTools(const CaloGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle); 
      PreshowerCluster makeOnePreshowerCluster(int stripwindow, ESDetId *strip);

      /// preshower calibration constants
      static const double mip_;
      static const double gamma_;
      static const double calib_planeX_;
      static const double calib_planeY_;
      static const int clusterwindowsize_;


    private:
      const CaloGeometry* geom_;
      CaloSubdetectorTopology *estopology_;

      std::vector<ESDetId> esroad_2d;

      // The set of used DetID's for a given event used in the preshower clustering algo:
      HitsID  used_strips;

      // creates the map of rechits used in the preshower clustering algo
      RecHitsMap  rechits_map;

      void findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane);
      bool goodStrip(RecHitsMap::iterator candidate_it);
};

#endif
