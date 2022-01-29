useHLTFilter = True
correctHits = False

import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi
import os, sys, imp, re
CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("analyzerFillEpsilon")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag  

from Configuration.AlCa.GlobalTag import GlobalTag


process.GlobalTag.globaltag = '105X_dataRun2_v8'
#DUMMY RECHIT
process.dummyHits = cms.EDProducer('DummyRechitDigis',
                                     doDigi = cms.untracked.bool(True),
                                     # rechits
                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB'),
                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE'),
                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),
                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),
                                     # digis
                                     barrelDigis            = cms.InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis"),
                                     endcapDigis            = cms.InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis"),
                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),
                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))

#RAW to DIGI'
#https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/RecoLocalCalo/EcalRecProducers/test/testMultipleEcalRecoLocal_cfg.py
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.raw2digi_step = cms.Sequence(process.RawToDigi)
#DIGI to UNCALIB
process.load('Configuration.StandardSequences.Reconstruction_cff')
import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi
process.ecalMultiFitUncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()
process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','analyzerFillEpsilon')
process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','analyzerFillEpsilon')
process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = False #added this line to make code run
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
process.GlobalTag.toGet = cms.VPSet(
     cms.PSet(record = cms.string('EcalLaserAPDPNRatiosRcd'),
              tag = cms.string('EcalLaserAPDPNRatios_rereco2018_v3'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalPFRecHitThresholdsRcd'),
              tag = cms.string('EcalPFRecHitThresholds_UL_2018_2e3sig'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalPulseShapesRcd'),
              tag = cms.string('EcalPulseShapes_UltraLegacy2018_calib'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalPedestalsRcd'),
              tag = cms.string('EcalPedestals_timestamp_2018_18January2019_collisions_blue_laser'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalLaserAlphasRcd'),
              tag = cms.string('EcalLaserAlphas_EB152-150_EEoptimized18'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalIntercalibConstantsRcd'),
              tag = cms.string('EcalIntercalibConstants_Run2018ABCD_run297056_eopPNEB_v1'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
     cms.PSet(record = cms.string('EcalChannelStatusRcd'),
              tag = cms.string('EcalChannelStatus_v13_offline'),
              connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
     ),
)

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
)
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        #'root://cms-xrd-global.cern.ch//store/data/Run2018A/AlCaP0/RAW/v1/000/315/257/00000/84307A29-9C49-E811-8737-FA163E51C596.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018A/AlCaP0/RAW/v1/000/315/257/00000/72EB1782-9849-E811-BFF6-FA163E9393A6.root'
        'root://cms-xrd-global.cern.ch//store/data/Run2018D/AlCaP0/RAW/v1/000/321/396/00000/248ADA80-FBA1-E811-8742-FA163E7B2F96.root'
    ),
    skipBadFiles = cms.untracked.bool(True)
)
import FWCore.PythonUtilities.LumiList as LumiList
json_file = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
process.source.lumisToProcess = LumiList.LumiList(filename = json_file).getVLuminosityBlockRange()

process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')
process.analyzerFillEpsilon.OutputDir = cms.untracked.string('./')
process.analyzerFillEpsilon.OutputFile = cms.untracked.string('AlCaP0_2018_testCFG_EcalNtp_0.root')
process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('CalibCode/FillEpsilonPlot/data/caloGeometry.root')
process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/AlCaP0_2018_testCFG/iter_-1/AlCaP0_2018_testCFG_calibMap.root')
process.analyzerFillEpsilon.Endc_x_y                        = cms.untracked.string('CalibCode/FillEpsilonPlot/data/Endc_x_y_ring.txt')
process.analyzerFillEpsilon.HLTResults                  = cms.untracked.bool(True)
process.analyzerFillEpsilon.HLTResultsNameEB            = cms.untracked.string('AlCa_EcalPi0EB')
process.analyzerFillEpsilon.HLTResultsNameEE            = cms.untracked.string('AlCa_EcalPi0EE')
process.analyzerFillEpsilon.RemoveDead_Flag             = cms.untracked.bool(True)
process.analyzerFillEpsilon.RemoveDead_Map              = cms.untracked.string('')
process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(True)
process.analyzerFillEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( True )
process.analyzerFillEpsilon.scalingEoverEtrueCC_g1 = cms.untracked.double(1.0)
process.analyzerFillEpsilon.scalingEoverEtrueCC_g2 = cms.untracked.double(1.0)
process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/pi0CC_2018_EoverEtrue_foldSM_nFit10_onlyEB_fixGamma2EoverEtrue/iter_0/pi0CC_2018_EoverEtrue_foldSM_nFit10_onlyEB_fixGamma2EoverEtrue_calibMap.root")
process.analyzerFillEpsilon.useOnlyEEClusterMatchedWithES = cms.untracked.bool(True)

### choosing proper input tag (recalibration module changes the collection names)
if correctHits:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')
else:
    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")
    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")
process.analyzerFillEpsilon.ESRecHitCollectionTag = cms.untracked.InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES')
process.analyzerFillEpsilon.triggerTag   = cms.untracked.InputTag("TriggerResults","","HLT")
process.analyzerFillEpsilon.L1GTobjmapTag   = cms.untracked.InputTag("hltGtStage2Digis")
process.analyzerFillEpsilon.CalibType    = cms.untracked.string('xtal')
process.analyzerFillEpsilon.CurrentIteration = cms.untracked.int32(0)
process.analyzerFillEpsilon.EB_Seed_E = cms.untracked.double(0.5)
process.analyzerFillEpsilon.useEE_EtSeed = cms.untracked.bool(False)
process.analyzerFillEpsilon.EE_Seed_E = cms.untracked.double(1.0)
process.analyzerFillEpsilon.EE_Seed_Et = cms.untracked.double(0.0)
process.analyzerFillEpsilon.Pi0PtCutEB_low = cms.untracked.double(2.0)
process.analyzerFillEpsilon.Pi0PtCutEB_high = cms.untracked.double(1.75)
process.analyzerFillEpsilon.Pi0PtCutEE_low = cms.untracked.double(3.75)
process.analyzerFillEpsilon.Pi0PtCutEE_high = cms.untracked.double(2.0)
process.analyzerFillEpsilon.gPtCutEB_low = cms.untracked.double(0.65)
process.analyzerFillEpsilon.gPtCutEB_high = cms.untracked.double(0.65)
process.analyzerFillEpsilon.gPtCutEE_low = cms.untracked.double(1.1)
process.analyzerFillEpsilon.gPtCutEE_high = cms.untracked.double(0.95)
process.analyzerFillEpsilon.Pi0IsoCutEB_low = cms.untracked.double(0.2)
process.analyzerFillEpsilon.Pi0IsoCutEB_high = cms.untracked.double(0.2)
process.analyzerFillEpsilon.Pi0IsoCutEE_low = cms.untracked.double(0.2)
process.analyzerFillEpsilon.Pi0IsoCutEE_high = cms.untracked.double(0.2)
process.analyzerFillEpsilon.CutOnHLTIso = cms.untracked.bool(True)
process.analyzerFillEpsilon.Pi0HLTIsoCutEB_low = cms.untracked.double(0.5)
process.analyzerFillEpsilon.Pi0HLTIsoCutEB_high = cms.untracked.double(0.5)
process.analyzerFillEpsilon.Pi0HLTIsoCutEE_low = cms.untracked.double(0.5)
process.analyzerFillEpsilon.Pi0HLTIsoCutEE_high = cms.untracked.double(0.5)
process.analyzerFillEpsilon.nXtal_1_EB_low = cms.untracked.int32(7)
process.analyzerFillEpsilon.nXtal_1_EB_high = cms.untracked.int32(7)
process.analyzerFillEpsilon.nXtal_2_EB_low = cms.untracked.int32(7)
process.analyzerFillEpsilon.nXtal_2_EB_high = cms.untracked.int32(7)
process.analyzerFillEpsilon.nXtal_1_EE_low = cms.untracked.int32(6)
process.analyzerFillEpsilon.nXtal_1_EE_high = cms.untracked.int32(6)
process.analyzerFillEpsilon.nXtal_2_EE_low = cms.untracked.int32(6)
process.analyzerFillEpsilon.nXtal_2_EE_high = cms.untracked.int32(6)
process.analyzerFillEpsilon.S4S9_EB_low = cms.untracked.double(0.88)
process.analyzerFillEpsilon.S4S9_EB_high = cms.untracked.double(0.9)
process.analyzerFillEpsilon.S4S9_EE_low = cms.untracked.double(0.85)
process.analyzerFillEpsilon.S4S9_EE_high = cms.untracked.double(0.92)
process.analyzerFillEpsilon.Barrel_orEndcap = cms.untracked.string('ALL_PLEASE')
process.analyzerFillEpsilon.useMassInsteadOfEpsilon = cms.untracked.bool(True)
process.analyzerFillEpsilon.fillKinematicVariables = cms.untracked.bool(True)
process.analyzerFillEpsilon.L1TriggerInfo = cms.untracked.bool(False)
process.analyzerFillEpsilon.L1SeedsPi0Stream = cms.untracked.string("")
process.analyzerFillEpsilon.nL1SeedsPi0Stream = cms.untracked.int32(0)
process.p = cms.Path()
if useHLTFilter:
    process.p *= process.AlcaP0Filter
if correctHits:
    print 'ADDING RECALIB RECHIT MODULE WITH PARAMETERS'
    print 'ENERGY SCALE '+str(process.ecalPi0ReCorrected.doEnergyScale)
    print 'INTERCALIBRATION '+str(process.ecalPi0ReCorrected.doIntercalib)
    print 'LASER '+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected
process.p *= process.dummyHits
process.p *= process.ecalMultiFitUncalibRecHit
process.p *= process.ecalLocalRecoSequence
process.p *= process.analyzerFillEpsilon
process.endp = cms.EndPath()
