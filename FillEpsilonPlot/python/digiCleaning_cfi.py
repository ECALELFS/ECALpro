import FWCore.ParameterSet.Config as cms

digiCleaning = cms.EDProducer("DigiCleaning",
                              cleanedEBDigiCollection = cms.string("cleanedEBDigiCollection"),
                              cleanedEEDigiCollection = cms.string("cleanedEEDigiCollection"),
                              barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","HLT"),
                              endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","HLT"),
)
