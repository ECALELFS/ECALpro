import FWCore.ParameterSet.Config as cms

cleanedEcalDigis = cms.EDFilter("CleanedDigiCollectionProducer",
    ebDigis = cms.InputTag("ecalDigis","ebDigis"),
    eeDigis = cms.InputTag("ecalDigis","eeDigis"),
    cleanedEBDigiCollection = cms.string('ebCleanedDigis'),
    cleanedEEDigiCollection = cms.string('eeCleanedDigis'),
    filter = cms.bool( True )
)
