#include "CalibCode/FillEpsilonPlot/plugins/CleanedDigiCollectionProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"



#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "FWCore/Utilities/interface/transform.h"

#include <iostream>

CleanedDigiCollectionProducer::CleanedDigiCollectionProducer(const edm::ParameterSet& iConfig) 
{

  ebDigisToken_ = consumes<EBDigiCollection>(iConfig.getParameter< edm::InputTag > ("ebDigis"));
  eeDigisToken_ = consumes<EEDigiCollection>(iConfig.getParameter< edm::InputTag > ("eeDigis"));
	  
  cleanedEBDigiCollection_ = iConfig.getParameter<std::string>("cleanedEBDigiCollection");
  cleanedEEDigiCollection_ = iConfig.getParameter<std::string>("cleanedEEDigiCollection");
  filter_ = iConfig.getParameter<bool>("filter");
  
   //register your products
  produces< EBDigiCollection > (cleanedEBDigiCollection_) ;
  produces< EEDigiCollection > (cleanedEEDigiCollection_) ;
  
}


CleanedDigiCollectionProducer::~CleanedDigiCollectionProducer()
{}


// ------------ method called to produce the data  ------------
bool
CleanedDigiCollectionProducer::filter (edm::Event& iEvent, 
                                const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<EBDigiCollection> ebDigisHandle;
   Handle<EEDigiCollection> eeDigisHandle;
   iEvent.getByToken(ebDigisToken_,ebDigisHandle);
   iEvent.getByToken(eeDigisToken_,eeDigisHandle);
   if( !ebDigisHandle.isValid() || !eeDigisHandle.isValid() ) 
     {
       edm::LogError("CleanedDigiCollectionProducer") << "Barrel or Endcap digi collection not found";
       return false;
     }
   
   //Create empty output collections
   std::unique_ptr< EBDigiCollection > cleanedEBDigiCollection (new EBDigiCollection()) ;
   std::unique_ptr< EEDigiCollection > cleanedEEDigiCollection (new EEDigiCollection()) ;
   
   for(EBDigiCollection::const_iterator itdg=ebDigisHandle->begin(); itdg!=ebDigisHandle->end(); ++itdg) {
     DetId detid(itdg->id());
     if(detid.subdetId()==EcalBarrel)
       cleanedEBDigiCollection->push_back(itdg->id(),itdg->begin());
   }
   for(EEDigiCollection::const_iterator itdg=eeDigisHandle->begin(); itdg!=eeDigisHandle->end(); ++itdg) {
     DetId detid(itdg->id());
     if(detid.subdetId()==EcalEndcap)
       cleanedEEDigiCollection->push_back(itdg->id(),itdg->begin());
   }

   // Populate the new collections
   bool pass = (cleanedEBDigiCollection->size() > 0 && cleanedEEDigiCollection->size() > 0);
   iEvent.put( std::move(cleanedEBDigiCollection),cleanedEBDigiCollection_ );
   iEvent.put( std::move(cleanedEEDigiCollection),cleanedEEDigiCollection_ );
   
   if(filter_) return pass;
   else return true;
}

DEFINE_FWK_MODULE( CleanedDigiCollectionProducer );
