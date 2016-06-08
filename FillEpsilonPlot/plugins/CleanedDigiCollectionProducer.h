#ifndef _CLEANEDDIGICOLLECTIONPRODUCER_H
#define _CLEANEDDIGICOLLECTIONPRODUCER_H

// -*- C++ -*-
//
// Package:    CleanedDigiCollectionProducer
// Class:      CleanedDigiCollectionProducer
// 




// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

class CleanedDigiCollectionProducer : public edm::stream::EDFilter<> {
   public:
      //! ctor
      explicit CleanedDigiCollectionProducer(const edm::ParameterSet&);
      ~CleanedDigiCollectionProducer();
      //! producer
      virtual bool filter(edm::Event &, const edm::EventSetup&);

   private:
      // ----------member data ---------------------------
      edm::EDGetTokenT<EBDigiCollection>     ebDigisToken_;
      edm::EDGetTokenT<EEDigiCollection>     eeDigisToken_;
      std::string cleanedEBDigiCollection_;
      std::string cleanedEEDigiCollection_;
      bool filter_;
  
};

#endif
