
useHLTFilter = True

correctHits = True

import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

process = cms.Process("analyzerFillEpsilon")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


### Recalibration Module to apply laser corrections on the fly
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V22::All'
process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                tag = cms.string("EcalIntercalibConstants_Bon_V20110616_weightedAverage_2010V3_ratio"),
                connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
        ),
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                tag = cms.string("EcalLaserAPDPNRatios_data_20111122_158851_180363"),
                connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                #connect = cms.untracked.string("oracle://cms_orcoff_prep/CMS_COND_ECAL")
        ),
        cms.PSet(record = cms.string("EcalLaserAlphasRcd"),
                tag = cms.string("EcalLaserAlphas_lto420-620_progr_data_20111122"),
                connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
                #connect = cms.untracked.string("oracle://cms_orcoff_prep/CMS_COND_ECAL")
        )
)


### Running on AlcaRAW requires filtering AlcaPi0 events from AlcaEta events
if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(False),
        doIntercalib = cms.bool(False),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEB","HLT"),
        EERecHitCollection = cms.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEE","HLT"),
        EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
        )

if useHLTFilter:
     import copy
     from HLTrigger.HLTfilters.hltHighLevel_cfi import *
     process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
     process.AlcaP0Filter.throw = cms.bool(False)
     process.AlcaP0Filter.HLTPaths = ["AlCa_EcalPi0_*"]


### Paramters of the FillEpsilonPlot module
process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')

process.analyzerFillEpsilon.OutputFile = cms.untracked.string('EcalNtp.root')
process.analyzerFillEpsilon.EBContainmentCorrections = cms.untracked.string('totNewPi0TupleMB_fillingTot.fittedcorrectionsEB.root')
process.analyzerFillEpsilon.EEContainmentCorrections = cms.untracked.string('totNewPi0TupleMB_fillingTot.fittedcorrectionsEE.root')
process.analyzerFillEpsilon.EEContainmentCorrections = cms.untracked.string('totNewPi0TupleMB_fillingTot.fittedcorrectionsEE.root')
process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('caloGeometry.root')
process.analyzerFillEpsilon.L1TriggerTag = cms.untracked.InputTag("hltGtDigis")
process.analyzerFillEpsilon.CalibType = cms.untracked.string('xtal')
process.analyzerFillEpsilon.CurrentIteration = cms.untracked.int32(0) # <--- iteration number
process.analyzerFillEpsilon.OutputDir = cms.untracked.string('/afs/cern.ch/user/m/mgrassi/scratch0/calib/428-parallelCalib/src/CalibCode/FillEpsilonPlot/test')
process.analyzerFillEpsilon.Pi0PtCutEB = cms.untracked.double(1.7)
process.analyzerFillEpsilon.Pi0PtCutEE = cms.untracked.double(2.5)
process.analyzerFillEpsilon.Pi0IsoCutEB = cms.untracked.double(0.5)
process.analyzerFillEpsilon.Pi0IsoCutEE = cms.untracked.double(0.5)
process.analyzerFillEpsilon.AlcaL1TrigNames = cms.untracked.vstring("L1_SingleIsoEG5","L1_SingleIsoEG8","L1_SingleIsoEG10","L1_SingleIsoEG12","L1_SingleIsoEG15","L1_SingleEG2","L1_SingleEG5","L1_SingleEG8","L1_SingleEG10","L1_SingleEG12","L1_SingleEG15","L1_SingleEG20","L1_SingleJet6U","L1_SingleJet10U","L1_SingleJet20U","L1_SingleJet30U","L1_SingleJet40U","L1_SingleJet50U","L1_DoubleJet30U","L1_DoubleEG5","L1_DoubleEG2")

### choosing proper input tag (recalibration module changes the collection names)
if correctHits:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEB")
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEE")
else:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEB", "HLT")
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEE", "HLT")

process.analyzerFillEpsilon.ESRecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsES", "HLT")

# input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/data/Run2010A/AlCaP0/ALCARECO/v4/000/144/114/A63B55F1-3CB4-DF11-9D8F-0030487A1884.root'
        # 'file:/cmshome/mgrassi/calib/399-shacalib/src/data/A63B55F1-3CB4-DF11-9D8F-0030487A1884.root'
        '/store/data/Run2011A/AlCaP0/RAW/v1/000/173/692/BEDBB812-91CC-E011-83AF-BCAEC53296F2.root'
    )
)

process.p = cms.Path()

if useHLTFilter:
    process.p *= process.AlcaP0Filter

if correctHits:
    print "ADDING RECALIB RECHIT MODULE WITH PARAMETERS"
    print "ENERGY SCALE "+str(process.ecalPi0ReCorrected.doEnergyScale)
    print "INTERCALIBRATION "+str(process.ecalPi0ReCorrected.doIntercalib)
    print "LASER "+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected

process.p *= process.analyzerFillEpsilon


