useHLTFilter = True
correctHits = False

import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi
import os, sys, imp, re
CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("analyzerFillEpsilon")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'
#DUMMY RECHIT
process.dummyHits = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB' ,'HLT'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE' ,'HLT'),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","HLT"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","HLT"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))
process.load('CalibCode.FillEpsilonPlot.digiCleaning_cfi')

#RAW to DIGI'
#https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/RecoLocalCalo/EcalRecProducers/test/testMultipleEcalRecoLocal_cfg.py
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.raw2digi_step = cms.Sequence(process.RawToDigi)
#DIGI to UNCALIB
process.load('Configuration.StandardSequences.Reconstruction_cff')
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
process.ecalMultiFitUncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
#process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','analyzerFillEpsilon')
#process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','analyzerFillEpsilon')
process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('digiCleaning','cleanedEBDigiCollection')
process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('digiCleaning','cleanedEEDigiCollection')
#UNCALIB to CALIB
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *
process.ecalDetIdToBeRecovered =  RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi.ecalDetIdToBeRecovered.clone()
process.ecalRecHit.killDeadChannels = cms.bool( False )
process.ecalRecHit.recoverEBVFE = cms.bool( False )
process.ecalRecHit.recoverEEVFE = cms.bool( False )
process.ecalRecHit.recoverEBFE = cms.bool( False )
process.ecalRecHit.recoverEEFE = cms.bool( False )
process.ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )
process.ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )
process.ecalLocalRecoSequence = cms.Sequence(ecalRecHit)
### Recalibration Module to apply laser corrections on the fly
if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(False),
        doIntercalib = cms.bool(False),
        doLaserCorrections = cms.bool(False),
        EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon"),
        EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon"),
        EBRecalibRecHitCollection = cms.string("pi0EcalRecHitsEB"),
        EERecalibRecHitCollection = cms.string("pi0EcalRecHitsEE")
    )

### Running on AlcaRAW requires filtering AlcaPi0 events from AlcaEta events
if useHLTFilter:
    import copy
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
    process.AlcaP0Filter.throw = cms.bool(False)
    process.AlcaP0Filter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
    process.AlcaP0Filter.HLTPaths = ["AlCa_EcalPi0E*"]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr = cms.untracked.PSet(
#        threshold  = cms.untracked.string('WARNING'),
#        ERROR      = cms.untracked.PSet (
#                                         limit = cms.untracked.int32(1)
#        )
#)

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   #SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')
)
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/data/Run2015B/AlCaP0/RAW/v1/000/251/548/00000/36828250-1A28-E511-B1A9-02163E013515.root',
#        'root://eoscms//eos/cms/store/data/Run2015B/AlCaP0/RAW/v1/000/251/548/00000/A85DBC95-1628-E511-8956-02163E0136C5.root',
#        'root://eoscms//eos/cms/store/data/Run2015B/AlCaP0/RAW/v1/000/251/559/00000/066C6BFD-E52A-E511-862A-02163E014181.root'
    )
)

process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')
process.analyzerFillEpsilon.OutputDir = cms.untracked.string('/tmp/')
process.analyzerFillEpsilon.OutputFile = cms.untracked.string('2015B_TEcalNtp_2.root')
process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('CalibCode/FillEpsilonPlot/data/caloGeometry.root')
process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms/store/caf/user/cmackay/ALL_2015B_Test2/iter_-1/2015B_TcalibMap.root')
process.analyzerFillEpsilon.useEBContainmentCorrections = cms.untracked.bool(True)
process.analyzerFillEpsilon.useEEContainmentCorrections = cms.untracked.bool(False)
process.analyzerFillEpsilon.EBContainmentCorrections = cms.untracked.string('CalibCode/FillEpsilonPlot/data/totNewPi0TupleMB_fillingTot.fittedcorrectionsEB.root')
process.analyzerFillEpsilon.MVAEBContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_pi01_Mediumtrain.root')
process.analyzerFillEpsilon.MVAEBContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_pi02_Mediumtrain.root')
process.analyzerFillEpsilon.MVAEEContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_pi01_Mediumtrain_EE.root')
process.analyzerFillEpsilon.MVAEEContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_pi02_Mediumtrain_EE.root')
process.analyzerFillEpsilon.MVAEBContainmentCorrections_eta01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_eta1_Mediumtrain.root')
process.analyzerFillEpsilon.MVAEBContainmentCorrections_eta02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/JOSH_MVA_eta2_Mediumtrain.root')
process.analyzerFillEpsilon.Endc_x_y                        = cms.untracked.string('CalibCode/FillEpsilonPlot/data/Endc_x_y_ring.txt')
process.analyzerFillEpsilon.EBPHIContainmentCorrections = cms.untracked.string('CalibCode/FillEpsilonPlot/data/correctionsEB_PHI.root')
process.analyzerFillEpsilon.EEContainmentCorrections    = cms.untracked.string('CalibCode/FillEpsilonPlot/data/totNewPi0TupleMB_fillingTot.fittedcorrectionsEE.root')
process.analyzerFillEpsilon.ContCorr_EB                 = cms.untracked.string('CalibCode/FillEpsilonPlot/data/correctionsEB.root')
process.analyzerFillEpsilon.HLTResults                  = cms.untracked.bool(True)
process.analyzerFillEpsilon.HLTResultsNameEB            = cms.untracked.string('AlCa_EcalPi0EB')
process.analyzerFillEpsilon.HLTResultsNameEE            = cms.untracked.string('AlCa_EcalPi0EE')
process.analyzerFillEpsilon.RemoveDead_Flag             = cms.untracked.bool(True)
process.analyzerFillEpsilon.RemoveDead_Map              = cms.untracked.string('')
process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(True)
process.analyzerFillEpsilon.useOnlyEEClusterMatchedWithES = cms.untracked.bool(True)

### choosing proper input tag (recalibration module changes the collection names)
if correctHits:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')
else:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")
process.analyzerFillEpsilon.ESRecHitCollectionTag = cms.untracked.InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES','HLT')
process.analyzerFillEpsilon.L1TriggerTag = cms.untracked.InputTag('simGtDigis','','HLT')
process.analyzerFillEpsilon.triggerTag   = cms.untracked.InputTag("TriggerResults","","HLT")
process.analyzerFillEpsilon.hltL1GtObjectMap   = cms.untracked.InputTag("hltL1GtObjectMap","","HLT")
process.analyzerFillEpsilon.CalibType    = cms.untracked.string('xtal')
process.analyzerFillEpsilon.CurrentIteration = cms.untracked.int32(0)
process.analyzerFillEpsilon.EB_Seed_E = cms.untracked.double(0.5)
process.analyzerFillEpsilon.useEE_EtSeed = cms.untracked.bool(False)
process.analyzerFillEpsilon.EE_Seed_E = cms.untracked.double(1.5)
process.analyzerFillEpsilon.EE_Seed_Et = cms.untracked.double(0.5)
process.analyzerFillEpsilon.Pi0PtCutEB_low = cms.untracked.double(1.8)
process.analyzerFillEpsilon.Pi0PtCutEB_high = cms.untracked.double(2.6)
process.analyzerFillEpsilon.Pi0PtCutEE_low = cms.untracked.double(3.6)
process.analyzerFillEpsilon.Pi0PtCutEE_high = cms.untracked.double(3.6)
process.analyzerFillEpsilon.gPtCutEB_low = cms.untracked.double(0.6)
process.analyzerFillEpsilon.gPtCutEB_high = cms.untracked.double(0.6)
process.analyzerFillEpsilon.gPtCutEE_low = cms.untracked.double(1.)
process.analyzerFillEpsilon.gPtCutEE_high = cms.untracked.double(1.)
process.analyzerFillEpsilon.Pi0IsoCutEB_low = cms.untracked.double(0.2)
process.analyzerFillEpsilon.Pi0IsoCutEB_high = cms.untracked.double(0.05)
process.analyzerFillEpsilon.Pi0IsoCutEE_low = cms.untracked.double(0.3)
process.analyzerFillEpsilon.Pi0IsoCutEE_high = cms.untracked.double(0.3)
process.analyzerFillEpsilon.CutOnHLTIso = cms.untracked.bool(False)
process.analyzerFillEpsilon.Pi0HLTIsoCutEB_low = cms.untracked.double(999)
process.analyzerFillEpsilon.Pi0HLTIsoCutEB_high = cms.untracked.double(999)
process.analyzerFillEpsilon.Pi0HLTIsoCutEE_low = cms.untracked.double(999)
process.analyzerFillEpsilon.Pi0HLTIsoCutEE_high = cms.untracked.double(999)
process.analyzerFillEpsilon.nXtal_1_EB_low = cms.untracked.double(4)
process.analyzerFillEpsilon.nXtal_1_EB_high = cms.untracked.double(4)
process.analyzerFillEpsilon.nXtal_2_EB_low = cms.untracked.double(5)
process.analyzerFillEpsilon.nXtal_2_EB_high = cms.untracked.double(5)
process.analyzerFillEpsilon.nXtal_1_EE_low = cms.untracked.double(4)
process.analyzerFillEpsilon.nXtal_1_EE_high = cms.untracked.double(4)
process.analyzerFillEpsilon.nXtal_2_EE_low = cms.untracked.double(5)
process.analyzerFillEpsilon.nXtal_2_EE_high = cms.untracked.double(5)
process.analyzerFillEpsilon.S4S9_EB_low = cms.untracked.double(0.6)
process.analyzerFillEpsilon.S4S9_EB_high = cms.untracked.double(0.75)
process.analyzerFillEpsilon.S4S9_EE_low = cms.untracked.double(0.8)
process.analyzerFillEpsilon.S4S9_EE_high = cms.untracked.double(0.8)
process.analyzerFillEpsilon.Barrel_orEndcap = cms.untracked.string('ALL_PLEASE')
process.analyzerFillEpsilon.JSONfile = cms.untracked.string('CalibCode/FillEpsilonPlot/data/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_v2.txt')
process.p = cms.Path(process.bunchSpacingProducer)
if useHLTFilter:
    process.p *= process.AlcaP0Filter
if correctHits:
    print 'ADDING RECALIB RECHIT MODULE WITH PARAMETERS'
    print 'ENERGY SCALE '+str(process.ecalPi0ReCorrected.doEnergyScale)
    print 'INTERCALIBRATION '+str(process.ecalPi0ReCorrected.doIntercalib)
    print 'LASER '+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected
process.p *= process.dummyHits
process.p *= process.digiCleaning
process.p *= process.ecalMultiFitUncalibRecHit
process.p *= process.ecalLocalRecoSequence
process.p *= process.analyzerFillEpsilon


process.outputALCARAW = cms.OutputModule("PoolOutputModule",
    # SelectEvents = cms.untracked.PSet(
    #     SelectEvents = cms.vstring('pathALCARECOEcalUncalSingleElectron')
    # ),
    # dataset = cms.untracked.PSet(
    #     dataTier = cms.untracked.string('ALCARECO'),
    #     filterName = cms.untracked.string('')
    # ),
                                         fileName = cms.untracked.string('alcaraw.root'),
    # maxSize = cms.untracked.int32(5120000),
                                         outputCommands = cms.untracked.vstring('keep *'),
                                         #drop *', 
    #     'keep uint_bunchSpacingProducer_*_*', 
    #     'keep *_pfMet_*_*', 
)

process.ALCARECOoutput_step = cms.EndPath(process.outputALCARAW)

