#ifndef EcalRegionalCalibration_H
#define EcalRegionalCalibration_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <string>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "CalibCode/CalibTools/interface/EcalCalibTypes.h"
#include "CalibCode/CalibTools/interface/EcalCalibMap.h"

#define PI0MASS 0.135

//using namespace std;

class RegionWeight {
  public:
    uint32_t   iRegion;
    float value;
};

typedef std::vector<RegionWeight> RegionWeightVector;

typedef std::pair<DetId, float> EnergyFraction;
typedef std::vector< EnergyFraction > EnergyFractionVector;

/// interface to avoid dealing with templated objects
class EcalRegionalCalibrationBase{
    public:
        virtual RegionWeightVector getWeights(const reco::CaloCluster* clus, int subDetId ) const =0;
        //EcalCalibMap<Type>* getCalibMap() =0;
        virtual  EcalCalibMapBase* getCalibMap() =0;
        virtual std::string printType() =0;
        virtual std::vector<DetId> allDetIdsInEBRegion(uint32_t iR) =0;
        virtual std::vector<DetId> allDetIdsInEERegion(uint32_t iR) =0;
};

template<class Type> class EcalRegionalCalibration : public EcalRegionalCalibrationBase {
    public:
        EcalRegionalCalibration();
        ~EcalRegionalCalibration(){}

        RegionWeightVector getWeights(const reco::CaloCluster* clus, int subDetId ) const;

        //EcalCalibMap<Type>* getCalibMap() { 
        EcalCalibMapBase* getCalibMap() { 
            EcalCalibMapBase *basePtr = calibMap; 
            return basePtr;
        }
        std::string printType() { return std::string(Type::printType()); }

        //RegionWeightVector getWeightsEE(const reco::CaloCluster* clus) const;
 
        std::vector<DetId> allDetIdsInEBRegion(uint32_t iR) {
            return Type::allDetIdsInRegion(iR);
        }

        std::vector<DetId> allDetIdsInEERegion(uint32_t iR) {
            return Type::allDetIdsInEERegion(iR);
        }

    private:
        EcalCalibMap<Type>*   calibMap;

};


//================================================================================================
template<class Type>
EcalRegionalCalibration<Type>::EcalRegionalCalibration()
//================================================================================================
{
    //std::cout << "EcalRegionalCalibration<Type>::EcalRegionalCalibration" << std::endl;
    calibMap = new EcalCalibMap<Type>();
}



//================================================================================================
template<class Type>
RegionWeightVector EcalRegionalCalibration<Type>::getWeights(const reco::CaloCluster* clus, int subDetId ) const {             
//================================================================================================
    RegionWeightVector weights;
    
    bool isEB = true;

    if( subDetId == EcalBarrel ) { isEB = true; }
    else if( subDetId == EcalEndcap ) { isEB = false; }
    else throw cms::Exception("EcalRegionalCalibration::getWeights") << "Subdetector Id not recognized\n";

    // map to compute energy&weight for each region
    std::map<uint32_t,float> regionWeightMap;
    const EnergyFractionVector& enHits = clus->hitsAndFractions();

    #ifdef DEBUG
    std::cout << "   -- cluster energy: " << clus->energy() << std::endl;
    #endif

    for(EnergyFractionVector::const_iterator it =  enHits.begin(); it != enHits.end(); ++it) {
        uint32_t iR = isEB ? Type::iRegion( it->first ) : Type::iRegionEE( it->first );
        float energy = it->second;
        if(energy==0.) continue;

        #ifdef DEBUG
        std::cout << "id: " << isEB ? EBDetId(it->first) : EEDetId(it->first);
        std::cout << " iR: " << iR << " E: " << energy
             << std::endl;
        #endif

        regionWeightMap[iR] += energy;
    } 

    for(std::map<uint32_t,float>::const_iterator it2 = regionWeightMap.begin(); it2 != regionWeightMap.end(); ++it2) 
    {
        RegionWeight w; 
        w.iRegion = it2->first;
        // only positive weights
        w.value   = (it2->second<0.) ? 0. : it2->second/clus->energy();
        // no w>1. weight!
        w.value   =  (w.value>1.) ? 1. : w.value;
        weights.push_back( w );

        #ifdef DEBUG
        std::cout << "iR: " << w.iRegion << " energy: " << it2->second
             << " weight: " << w.value
             << std::endl;
        #endif
    }
    return weights;
}



#endif
