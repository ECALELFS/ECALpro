// $Id: DumpCaloGeometry.cc,v 1.2 2012/03/16 12:22:09 mgrassi Exp $
// Shahram Rahatlou, Sapienza Universita` di Roma & INFN
// May 2010
//
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CalibCode/DumpCaloGeometry/src/DumpCaloGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Geometry
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

// ROOT stuff
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

#include <iostream>
#include <cmath>


using namespace std;
using namespace edm;
using namespace reco;


// ******************************************
// constructors
// *****************************************

DumpCaloGeometry::DumpCaloGeometry( const edm::ParameterSet& ps ) 
{
 
  /// Ttree
  m_outfilename = ps.getUntrackedParameter<string>("OutputFile","caloGeometry.root");


}


DumpCaloGeometry::~DumpCaloGeometry()
{
}


//--------------------------------------------------------
void DumpCaloGeometry::beginRun(const edm::Run& r, const EventSetup& context) {


  // TTree: Pi0 EB
  

  // create and cd into new folder


  m_file = TFile::Open(m_outfilename.c_str(),"RECREATE");
  m_tree = new TTree("Geometry","ECAL Geometry Tree");

  // GENERAL block branches
  m_tree->Branch("id",&id,"id/i");
  m_tree->Branch("xXtal",&xtalPos[0],"xXtal/F");
  m_tree->Branch("yXtal",&xtalPos[1],"yXtal/F");
  m_tree->Branch("zXtal",&xtalPos[2],"zXtal/F");
  m_tree->Branch("xAxisXtal",&xtalAxis[0],"xAxisXtal/F");
  m_tree->Branch("yAxisXtal",&xtalAxis[1],"yAxisXtal/F");
  m_tree->Branch("zAxisXtal",&xtalAxis[2],"zAxisXtal/F");


}

//--------------------------------------------------------
//--------------------------------------------------------
void DumpCaloGeometry::beginLuminosityBlock(const LuminosityBlock& lumiSeg, 
     const EventSetup& context) {
  
}

//-------------------------------------------------------------

void DumpCaloGeometry::analyze(const Event& iEvent, 
			       const EventSetup& iSetup ){  

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);     
  const CaloSubdetectorGeometry *geoEB = geoHandle->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  const CaloSubdetectorGeometry *geoEE = geoHandle->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
  const CaloSubdetectorGeometry *geoES = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);

  int n=0;
  std::vector<DetId> eb_ids = geoEB->getValidDetIds(DetId::Ecal,EcalBarrel);
  std::vector<DetId> ee_ids = geoEE->getValidDetIds(DetId::Ecal,EcalEndcap);
  std::vector<DetId> es_ids = geoES->getValidDetIds(DetId::Ecal,EcalPreshower);

  for (std::vector<DetId>::iterator i=eb_ids.begin(); i!=eb_ids.end(); i++) {
    n++;
    // since at least CMSSW_10_1_1, the object geoEB->getGeometry(*i) return an std::shared_ptr<const CaloCellGeometry>
    // probably the type was changed wrt to release 94X, because now code does not compile, saying cannot convert it to const CaloCellGeometry*
    // the solution is to use get() method of std::shared_ptr
    const CaloCellGeometry* cell = geoEB->getGeometry(*i).get();

    id = i->rawId();
    xtalPos[0] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).x();
    xtalPos[1] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).y();
    xtalPos[2] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).z();

    xtalAxis[0] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().x();
    xtalAxis[1] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().y();
    xtalAxis[2] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().z();

    m_tree->Fill();
  }

  for (std::vector<DetId>::iterator i=ee_ids.begin(); i!=ee_ids.end(); i++) {
    n++;
    const CaloCellGeometry* cell= geoEE->getGeometry(*i).get();

    id = i->rawId();
    xtalPos[0] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).x();
    xtalPos[1] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).y();
    xtalPos[2] = dynamic_cast<const TruncatedPyramid*>(cell)->getPosition(0.).z();

    xtalAxis[0] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().x();
    xtalAxis[1] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().y();
    xtalAxis[2] = dynamic_cast<const TruncatedPyramid*>(cell)->axis().z();

    m_tree->Fill();
  }


  for (std::vector<DetId>::iterator i=es_ids.begin(); i!=es_ids.end(); i++) {
    n++;
    const CaloCellGeometry* cell=geoES->getGeometry(*i).get();
    GlobalPoint position = cell->getPosition();

    id = i->rawId();
    xtalPos[0] = position.x();
    xtalPos[1] = position.y();
    xtalPos[2] = position.z();

    xtalAxis[0] = -999.;
    xtalAxis[1] = -999.;
    xtalAxis[2] = -999.;

    m_tree->Fill();
  }


  /*** temporary ***/
/*
  cout << "index\tieta\tiphi" << endl;
  for(int i=0; i < (EBDetId::MAX_SM*EBDetId::kCrystalsInPhi*EBDetId::kCrystalsInEta); ++i)
  {
      EBDetId myid(EBDetId::detIdFromDenseIndex( i ));
      cout << i << "\t" << myid.ieta() << "\t" << myid.iphi() << endl;
  }
*/

}


//--------------------------------------------------------
void DumpCaloGeometry::endLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {
}
//--------------------------------------------------------

//--------------------------------------------------------
void DumpCaloGeometry::endRun(const Run& r, const EventSetup& context){


  m_file->Write();
  m_file->Close();
  m_file->Delete();



}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DumpCaloGeometry);
