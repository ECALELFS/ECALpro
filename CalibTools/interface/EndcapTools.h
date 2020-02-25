#ifndef EndcapTools_h
#define EndcapTools_h

#include <vector>
#include <iostream>

#include "TFile.h"

#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"

class EndcapTools
{
    public:
        EndcapTools() { std::cout << "EndcapTools::EndcapTools" << std::endl;}
        ~EndcapTools(){}
            
        static std::vector<DetId> getDetIdsInRing(int aRingIndex);
        static int getRingIndex(DetId aDetId); 

        static void freeMemory();

        //void setCaloGeometry(const ECALGeometry* geometry) { caloGeometry_ = geometry; };
        //void setCaloGeometry(string & geometryName);

        static const int N_RING_ENDCAP = 78;
        static const int N_RING_ENDCAP_SIDE = 39;

    private:
        static void initializeFromGeometry(); 
  
        static bool isInitializedFromGeometry_;
        static int endcapRingIndex_[EEDetId::IX_MAX][EEDetId::IY_MAX]; 
        static ECALGeometry* caloGeometry_;
        static TFile *externalGeometryFile_;
};

#endif
