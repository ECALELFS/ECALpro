

//Includes need to read from geometry
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

//#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TFile.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "CalibCode/CalibTools/interface/EndcapTools.h"
#include "CalibCode/CalibTools/interface/GeometryService.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"

// /*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
// EndcapTools::EndcapTools()
// /*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
// {
//     isInitializedFromGeometry_ = false;
//     caloGeometry_ = 0;
//     //endcapRingIndex_[EEDetId::IX_MAX][EEDetId::IY_MAX];
// }


/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
int EndcapTools::getRingIndex(DetId id) 
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
{
    if (id.subdetId() != EcalEndcap)
        throw cms::Exception("EndcapTools") << "Current DetId does not belong to Ecal Endcap \n";

    if (!isInitializedFromGeometry_)
        initializeFromGeometry();

    EEDetId eid(id);
    int endcapRingIndex = endcapRingIndex_[eid.ix()-1][eid.iy()-1];
    if (eid.zside() == 1) endcapRingIndex += N_RING_ENDCAP_SIDE;
    return endcapRingIndex;
}



/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
std::vector<DetId> EndcapTools::getDetIdsInRing(int etaIndex) 
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
{
    std::vector<DetId> ringIds;

    if (etaIndex < 0 || etaIndex > N_RING_ENDCAP)
        throw cms::Exception("EndcapTools") << "Initializing without geometry handle" ;
  
    if (!isInitializedFromGeometry_)
	    initializeFromGeometry();
      
    int zside= (etaIndex < N_RING_ENDCAP_SIDE ) ? -1 : 1;
    int eeEtaIndex = (etaIndex)%N_RING_ENDCAP_SIDE; 

    for (int ix=0;ix<EEDetId::IX_MAX;++ix)
	    for (int iy=0;iy<EEDetId::IY_MAX;++iy)
	        if (endcapRingIndex_[ix][iy] == eeEtaIndex)
	            ringIds.push_back(EEDetId(ix+1,iy+1,zside));
  
    return ringIds;
} 


/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
void EndcapTools::initializeFromGeometry()
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/
{
    std::cout << "EndcapTools::initializeFromGeometry() " << std::endl;
    // TFile* externalGeometryFile_ = TFile::Open((GeometryService::getGeometryName()).c_str());
    // if(!externalGeometryFile_) cms::Exception("ExtGeom") << "External Geometry file (" << GeometryService::getGeometryName() << ") not found\n";
    // caloGeometry_ = ECALGeometry::getGeometry(externalGeometryFile_);
    caloGeometry_ = GeometryService::getGeometryPtr();
    std::cout << "EndcapTools::caloGeometry_ = " << caloGeometry_ << std::endl;

    if (!caloGeometry_)
        throw cms::Exception("EndcapTools") << "Initializing without geometry handle" ;

    float m_cellPosEta[EEDetId::IX_MAX][EEDetId::IY_MAX];
    for (int ix=0; ix<EEDetId::IX_MAX; ++ix) 
        for (int iy=0; iy<EEDetId::IY_MAX; ++iy) 
        {
            m_cellPosEta[ix][iy] = -1.;
	        endcapRingIndex_[ix][iy]=-9;
        }
  
  
    //const CaloSubdetectorGeometry *endcapGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

    //if (!endcapGeometry)
        //throw cms::Exception("EndcapTools") << "Ecal Endcap geometry not found" ;

    /// all the EE det ids
    //const std::vector<DetId>& m_endcapCells= caloGeometry_->getValidDetIds(DetId::Ecal, EcalEndcap);
    const std::map<DetId,GlobalPoint> & m_allCells = caloGeometry_->getPositionMap();


    /// find eta value for all the cells in EE+
    for (std::map<DetId,GlobalPoint>::const_iterator cellIt = m_allCells.begin();
         cellIt!=m_allCells.end();
         ++cellIt)
    {
        if( cellIt->first.subdetId() != EcalEndcap ) continue;
        EEDetId ee(cellIt->first);
        if (ee.zside() == -1) continue; //Just using +side to fill absEta x,y map
        // const CaloCellGeometry *cellGeometry = endcapGeometry->getGeometry(*endcapIt) ;
        int ics=ee.ix() - 1 ;
        int ips=ee.iy() - 1 ;
        m_cellPosEta[ics][ips] = fabs(cellIt->second.eta());

        //std::cout<<"EE Xtal, |eta| is "<<fabs(cellGeometry->getPosition().eta())<<std::endl;
    }

    float eta_ring[N_RING_ENDCAP_SIDE];

    /// span half of the y axis (x value is fixed)
    for (int ring=0; ring<N_RING_ENDCAP_SIDE; ++ring)
        eta_ring[ring]=m_cellPosEta[ring][50];

    double etaBoundary[N_RING_ENDCAP_SIDE + 1];
    etaBoundary[0]=1.47;
    etaBoundary[N_RING_ENDCAP/2]=4.0;

    for (int ring=1; ring<N_RING_ENDCAP_SIDE; ++ring)
        etaBoundary[ring]=(eta_ring[ring]+eta_ring[ring-1])/2.;
  
    for (int ring=0; ring<N_RING_ENDCAP_SIDE; ring++){
    // std::cout<<"***********************EE ring: "<<ring<<" eta "<<(etaBoundary[ring] + etaBoundary[ring+1])/2.<<std::endl;
        for (int ix=0; ix<EEDetId::IX_MAX; ix++)
            for (int iy=0; iy<EEDetId::IY_MAX; iy++)
	            if (m_cellPosEta[ix][iy]>etaBoundary[ring] && m_cellPosEta[ix][iy]<etaBoundary[ring+1]) {
	                endcapRingIndex_[ix][iy]=ring;
	                //std::cout<<"endcapRing_["<<ix+1<<"]["<<iy+1<<"] = "<<ring<<";"<<std::endl;  
	            }
    }

    // const std::vector<DetId>& m_barrelCells= caloGeometry_->getValidDetIds(DetId::Ecal, EcalBarrel);

    // for (std::vector<DetId>::const_iterator barrelIt = m_barrelCells.begin();
    //      barrelIt!=m_barrelCells.end();
    //      ++barrelIt)
    //   {
    //     EBDetId eb(*barrelIt);
    //   }


    // //EB

  isInitializedFromGeometry_ = true;

}



void EndcapTools::freeMemory()
{
    //if(caloGeometry_) delete caloGeometry_;
    if(externalGeometryFile_)
    {
        externalGeometryFile_->Close();
        externalGeometryFile_->Delete();
    }
    isInitializedFromGeometry_ = false;

}
