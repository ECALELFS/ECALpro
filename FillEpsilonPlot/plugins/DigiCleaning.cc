// -*- C++ -*-
//
// Package:    ECALpro/DigiCleaning
// Class:      DigiCleaning
// 
/**\class DigiCleaning DigiCleaning.cc ECALpro/DigiCleaning/plugins/DigiCleaning.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shervin Nourbakhsh
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

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include <cassert>
#define DEBUG

//
// class declaration
//

class DigiCleaning : public edm::EDProducer {
public:
  explicit DigiCleaning(const edm::ParameterSet&);
  ~DigiCleaning();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
      
  // ----------member data ---------------------------
  edm::InputTag tag_barrelDigiProducer_;
  edm::InputTag tag_endcapDigiProducer_;
  const std::string barrelDigiCollection_;
  const std::string endcapDigiCollection_;
};

DigiCleaning::DigiCleaning(const edm::ParameterSet& iConfig):
  tag_barrelDigiProducer_      (iConfig.getParameter< edm::InputTag > ("barrelDigis")),
  tag_endcapDigiProducer_      (iConfig.getParameter< edm::InputTag > ("endcapDigis")),
  barrelDigiCollection_ (iConfig.getParameter<std::string>("cleanedEBDigiCollection")),
  endcapDigiCollection_ (iConfig.getParameter<std::string>("cleanedEEDigiCollection"))
{


    produces< EBDigiCollection >(barrelDigiCollection_);
    produces< EEDigiCollection >(endcapDigiCollection_);

}

DigiCleaning::~DigiCleaning(){ }

void DigiCleaning::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  //std::cout << "\n-----------New Event ----------------\n " << std::endl;
   using namespace edm;   

   Handle<EBDigiCollection> digisEBHandle;
   Handle<EEDigiCollection> digisEEHandle;

   std::auto_ptr<EBDigiCollection> outputEBDigiCollection( new EBDigiCollection );
   std::auto_ptr<EEDigiCollection> outputEEDigiCollection( new EEDigiCollection );
   
//     bool foundEBDigi = true;
//     bool foundEEDigi = true;

   iEvent.getByLabel(tag_barrelDigiProducer_, digisEBHandle);
   iEvent.getByLabel(tag_endcapDigiProducer_, digisEEHandle);

   if(digisEBHandle.isValid()){
	   unsigned long long int size= digisEBHandle->size();
	   
	   for(unsigned long long int i=0; i < size; ++i
			   //EBDigiCollection::const_iterator digi_itr = digisEBHandle->begin(); 
			   //digi_itr != digisEBHandle->end(); 
			   //++digi_itr
		   ){
		   DataFrameContainer::const_IterPair itr = digisEBHandle->pair(i);
		   if(DetId(*(itr.first)).subdetId()!=EcalBarrel) continue; // if invalid is skipped
		   outputEBDigiCollection->push_back(*(itr.first),  &(*itr.second));
	   }
   }

   if(digisEEHandle.isValid()){
	   unsigned long long int size= digisEEHandle->size();
	   
	   for(unsigned long long int i=0; i < size; ++i){
		   DataFrameContainer::const_IterPair itr = digisEEHandle->pair(i);
		   if(DetId(*(itr.first)).subdetId()!=EcalEndcap){
#ifdef DEBUG
			   edm::LogWarning("Invalid EE digi")<<  "\t SubdetId = " << DetId(*(itr.first)).subdetId();
#endif
			   continue; // if invalid is skipped
		   }
		   outputEEDigiCollection->push_back(*(itr.first),  &(*itr.second));
	   }
   }


   iEvent.put(outputEBDigiCollection, barrelDigiCollection_);
   iEvent.put(outputEEDigiCollection, endcapDigiCollection_);

}

void DigiCleaning::beginJob(){ }
void DigiCleaning::endJob(){ }
void DigiCleaning::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  //desc.add<std::string>("cleanedEBDigiCollection");
  descriptions.addDefault(desc);
  //descriptions.add("digiCleaning", desc);
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(DigiCleaning);
