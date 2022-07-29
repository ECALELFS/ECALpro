// -*- C++ -*-
//
// Package:    ECALlite/DummyRechitDigisPi0
// Class:      DummyRechitDigisPi0
// 
/**\class DummyRechitDigisPi0 DummyRechitDigisPi0.cc ECALlite/DummyRechitDigisPi0/plugins/DummyRechitDigisPi0.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Joshua Robert Hardenbrook
//         Created:  Mon, 01 Jun 2015 06:58:24 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

//
// class declaration
//

class DummyRechitDigisPi0 : public edm::EDProducer {
public:
  explicit DummyRechitDigisPi0(const edm::ParameterSet&);
  ~DummyRechitDigisPi0();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
      
  // ----------member data ---------------------------
  edm::InputTag tag_barrelHitProducer_;
  edm::InputTag tag_endcapHitProducer_;
  const std::string barrelRecHitCollection_;
  const std::string endcapRecHitCollection_;
  edm::InputTag tag_barrelDigiProducer_;
  edm::InputTag tag_endcapDigiProducer_;
  const std::string barrelDigiCollection_;
  const std::string endcapDigiCollection_;
  const bool doDigi_;
  edm::EDGetTokenT<EcalRecHitCollection>  barrelRecHit_pi0_token_;
  edm::EDGetTokenT<EcalRecHitCollection>  endcapRecHit_pi0_token_;
  edm::EDGetTokenT<EBDigiCollection>  barrelDigis_pi0_token_;
  edm::EDGetTokenT<EEDigiCollection>  endcapDigis_pi0_token_;

};

DummyRechitDigisPi0::DummyRechitDigisPi0(const edm::ParameterSet& iConfig):
  tag_barrelHitProducer_       (iConfig.getParameter< edm::InputTag > ("barrelHitProducer")),
  tag_endcapHitProducer_       (iConfig.getParameter< edm::InputTag > ("endcapHitProducer")),
  barrelRecHitCollection_  (iConfig.getUntrackedParameter<std::string>("barrelRecHitCollection")),
  endcapRecHitCollection_  (iConfig.getUntrackedParameter<std::string>("endcapRecHitCollection")),
  tag_barrelDigiProducer_      (iConfig.getParameter< edm::InputTag > ("barrelDigis")),
  tag_endcapDigiProducer_      (iConfig.getParameter< edm::InputTag > ("endcapDigis")),
  barrelDigiCollection_ (iConfig.getUntrackedParameter<std::string>("barrelDigiCollection")),
  endcapDigiCollection_ (iConfig.getUntrackedParameter<std::string>("endcapDigiCollection")),
  doDigi_ (iConfig.getUntrackedParameter<bool>("doDigi"))
{
  if(doDigi_) { 
    produces< EBDigiCollection >(barrelDigiCollection_);
    produces< EEDigiCollection >(endcapDigiCollection_);
    barrelDigis_pi0_token_ = consumes<EBDigiCollection> (tag_barrelDigiProducer_);
    endcapDigis_pi0_token_ = consumes<EEDigiCollection> (tag_endcapDigiProducer_);
  }
  else {
    produces< EcalRecHitCollection >(barrelRecHitCollection_);
    produces< EcalRecHitCollection >(endcapRecHitCollection_);
    barrelRecHit_pi0_token_ = consumes<EcalRecHitCollection> (tag_barrelHitProducer_);
    endcapRecHit_pi0_token_ = consumes<EcalRecHitCollection> (tag_endcapHitProducer_);
   }
}

DummyRechitDigisPi0::~DummyRechitDigisPi0(){ }

void DummyRechitDigisPi0::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  //std::cout << "\n-----------New Event ----------------\n " << std::endl;
   using namespace edm;   
   // build an empty collection

   // fake rechits
   // handle to try to fill
   edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
   // dummy collection to Put()
   std::unique_ptr< EcalRecHitCollection > rechits_temp( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_temp2( new EcalRecHitCollection);

   //   EcalRecHit::EcalRecHit(const DetId& id, float energy, float time, uint32_t flags, uint32_t flagBits)
   // add one rechit
   EcalRecHit zero_rechit(0,0,0,0,0);
   EcalRecHitCollection zero_collection;
   zero_collection.push_back(zero_rechit);
   *rechits_temp  = zero_collection;
   *rechits_temp2 = zero_collection;

   std::unique_ptr< EcalRecHitCollection > rechits_eb( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > rechits_ee( new EcalRecHitCollection);

   // fake digis
   // handle to try to fill
   edm::Handle<EBDigiCollection> digisEBHandle;
   edm::Handle<EEDigiCollection> digisEEHandle;
   // dummy collection to Put()
   std::unique_ptr<EBDigiCollection> outputEBDigiCollection( new EBDigiCollection );
   std::unique_ptr<EEDigiCollection> outputEEDigiCollection( new EEDigiCollection );

   std::unique_ptr< EcalRecHitCollection > outputEBRecHitCollection( new EcalRecHitCollection);
   std::unique_ptr< EcalRecHitCollection > outputEERecHitCollection( new EcalRecHitCollection);

   //Digi zero_digi;
   EBDigiCollection ebfakecol;
   EEDigiCollection eefakecol;
   //ebfakecol.push_back(zerodigi);
   //eefakecol.push_back(zerodigi);

   // fake empty collections
   std::unique_ptr<EBDigiCollection> fakeEBDigiCollection( new EBDigiCollection );
   *fakeEBDigiCollection = ebfakecol;
   std::unique_ptr<EEDigiCollection> fakeEEDigiCollection( new EEDigiCollection) ;
   *fakeEEDigiCollection = eefakecol;


   if(doDigi_) { 

     iEvent.getByToken(barrelDigis_pi0_token_, digisEBHandle);
     if (digisEBHandle.isValid()) {
       *outputEBDigiCollection = *(digisEBHandle.product());
       // insert the EB collection
       iEvent.put(std::move(outputEBDigiCollection), barrelDigiCollection_);
     } else {
       iEvent.put(std::move(fakeEBDigiCollection), barrelDigiCollection_);

     }

     iEvent.getByToken(endcapDigis_pi0_token_, digisEEHandle);
     if (digisEEHandle.isValid()) {
       *outputEEDigiCollection = *(digisEEHandle.product());
       iEvent.put(std::move(outputEEDigiCollection), endcapDigiCollection_);
     } else {
       iEvent.put(std::move(fakeEEDigiCollection), endcapDigiCollection_);
     }
     
   } // end dummy digis
   
   // Build fake rechits collections
   else {    

     iEvent.getByToken(barrelRecHit_pi0_token_, barrelRecHitsHandle);
     // if you dont find the barrel rechits youre looking for, put in a fake one
     if (barrelRecHitsHandle.isValid()) {
       *rechits_eb = *(barrelRecHitsHandle.product());
       iEvent.put(std::move(rechits_eb), barrelRecHitCollection_);
     } else {
       iEvent.put(std::move(rechits_temp), barrelRecHitCollection_);
     }
     
     iEvent.getByToken(endcapRecHit_pi0_token_, endcapRecHitsHandle);
     // if you dont find the endcap rechits youre looking for, put in a fake one
     if (endcapRecHitsHandle.isValid()) {
       *rechits_ee = *(endcapRecHitsHandle.product());
       iEvent.put(std::move(rechits_ee), endcapRecHitCollection_);
     } else {
       iEvent.put(std::move(rechits_temp2), endcapRecHitCollection_);
     }

   }
   //std::cout << "-----------End Dummy Rechits ---------- " << std::endl;
}

void DummyRechitDigisPi0::beginJob(){ }
void DummyRechitDigisPi0::endJob(){ }
void DummyRechitDigisPi0::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DummyRechitDigisPi0);
