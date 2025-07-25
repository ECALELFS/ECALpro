useHLTFilter = True
correctHits = False

import FWCore.ParameterSet.Config as cms
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi
import os, sys, imp, re

from FWCore.ParameterSet.VarParsing import VarParsing

CMSSW_VERSION=os.getenv("CMSSW_VERSION")
process = cms.Process("analyzerFillEpsilon")
process.load("FWCore.MessageService.MessageLogger_cfi")

# Cmd line options
#options = VarParsing('analysis')
options = VarParsing('standard')
options.register('inputFiles', 
                 [''], 
                 VarParsing.multiplicity.list, 
                 VarParsing.varType.string, 
                 "Input file names")

options.register ('jsonFile', 
                  '', 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.string, 
                  "Input JSON")

options.register ('globalTag', 
                  '', 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.string, 
                  "Input global tag")

options.register ('outputDir', 
                  './output', 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.string, 
                  "Output directory")

options.register ('outputFile', 
                  '', 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.string, 
                  "Output file")

options.register ('conditionsFileName', 
                  '', 
                  VarParsing.multiplicity.singleton, 
                  VarParsing.varType.string, 
                  "Conditions file name")


options.parseArguments()



process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if options.globalTag != '':
    process.GlobalTag.globaltag = options.globalTag
else:
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_prompt', '')

###config with the conditions we want (REVISE LATER, for now run on prompt GT)
#process.load("CalibCode.FillEpsilonPlot.GTconditions_cff")
#process.load("CalibCode.FillEpsilonPlot.%s"%options.conditionsFileName) 
# process.GlobalTag.toGet = process.GTconditions

#DUMMY RECHIT
process.dummyHits = cms.EDProducer('DummyRechitDigisPi0',
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
process.ecalRecHit.cpu.killDeadChannels = cms.bool( False )
process.ecalRecHit.cpu.recoverEBVFE = cms.bool( False )
process.ecalRecHit.cpu.recoverEEVFE = cms.bool( False )
process.ecalRecHit.cpu.recoverEBFE = cms.bool( False )
process.ecalRecHit.cpu.recoverEEFE = cms.bool( False )
process.ecalRecHit.cpu.recoverEEIsolatedChannels = cms.bool( False )
process.ecalRecHit.cpu.recoverEBIsolatedChannels = cms.bool( False )
process.ecalLocalRecoSequence = cms.Sequence(process.ecalRecHit)

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
    process.AlcaP0Filter.HLTPaths = ["AlCa_EcalPi0E*", "AlCa_HIEcalPi0E*"]

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
)

process.source = cms.Source('PoolSource',
                            
                        fileNames = cms.untracked.vstring(
                            
                        #'root://cms-xrd-global.cern.ch//store/data/Run2018D/AlCaP0/RAW/v1/000/321/396/00000/248ADA80-FBA1-E811-8742-FA163E7B2F96.root'
                        ),
                            skipBadFiles = cms.untracked.bool(False)
)


process.source.fileNames = options.inputFiles 
###SJ - use the present available file to store the monitoring tree
#process.TFileService = cms.Service("TFileService", fileName = cms.string('monitoringTree.root'))


import FWCore.PythonUtilities.LumiList as LumiList
#json_file = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
json_file = options.jsonFile

process.source.lumisToProcess = LumiList.LumiList(filename = json_file).getVLuminosityBlockRange()

process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')
#process.analyzerFillEpsilon.OutputDir = cms.untracked.string('./')
#process.analyzerFillEpsilon.OutputFile = cms.untracked.string('AlCaP0_2018_testCFG_EcalNtp_0.root')
process.analyzerFillEpsilon.OutputDir = cms.untracked.string(options.outputDir)
process.analyzerFillEpsilon.OutputFile = cms.untracked.string(options.outputFile)

process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('CalibCode/FillEpsilonPlot/data/caloGeometry.root')
process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/AlCaP0_2018_testCFG/iter_-1/AlCaP0_2018_testCFG_calibMap.root') ####from M.Cipriani: actually, this is a file that doesn't really exist: when the very first iteration is run, there is no initial calibration map (the "iter_-1" somehow signals it to the code), while for the following iteration the map created by the previous step and stored on eos will be used

process.analyzerFillEpsilon.Endc_x_y                        = cms.untracked.string('CalibCode/FillEpsilonPlot/data/Endc_x_y_ring.txt')
process.analyzerFillEpsilon.HLTResults                  = cms.untracked.bool(True)
process.analyzerFillEpsilon.HLTResultsNameEB            = cms.untracked.string('EcalPi0EB') # the HLT path should contain this pattern. The Alca_ prefix is omitted to allow paths like Alca_HIEcalPi0EB as well
process.analyzerFillEpsilon.HLTResultsNameEE            = cms.untracked.string('EcalPi0EE') # the HLT path should contain this pattern. The Alca_ prefix is omitted to allow paths like Alca_HIEcalPi0EE as well
process.analyzerFillEpsilon.RemoveDead_Flag             = cms.untracked.bool(True)
process.analyzerFillEpsilon.RemoveDead_Map              = cms.untracked.string('')
process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(True)
#process.analyzerFillEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( True )
process.analyzerFillEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( False )
process.analyzerFillEpsilon.scalingEoverEtrueCC_g1 = cms.untracked.double(1.0)
process.analyzerFillEpsilon.scalingEoverEtrueCC_g2 = cms.untracked.double(1.0)
#process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/pi0CC_2018_EoverEtrue_foldSM_nFit10_onlyEB_fixGamma2EoverEtrue/iter_0/pi0CC_2018_EoverEtrue_foldSM_nFit10_onlyEB_fixGamma2EoverEtrue_calibMap.root")
process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string("")
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
    print('ADDING RECALIB RECHIT MODULE WITH PARAMETERS')
    print('ENERGY SCALE %s'%(str(process.ecalPi0ReCorrected.doEnergyScale)))
    print('INTERCALIBRATION %s'%(str(process.ecalPi0ReCorrected.doIntercalib)))
    print('LASER %s'%(str(process.ecalPi0ReCorrected.doLaserCorrections)))
    process.p *= process.ecalPi0ReCorrected
process.p *= process.dummyHits
process.p *= process.ecalMultiFitUncalibRecHit
process.p *= process.ecalLocalRecoSequence
process.p *= process.analyzerFillEpsilon


# process.out = cms.OutputModule("PoolOutputModule",
#                                fileName = cms.untracked.string('myOutputFile.root'),
#                                SelectEvents = cms.untracked.PSet( 
#                                    SelectEvents = cms.vstring("p") 
#                                ),
#                                outputCommands = cms.untracked.vstring('keep *')
# )


#process.endp = cms.EndPath(process.out)


# process.s = cms.Schedule(process.p, process.endp)

#print(process.dumpPython())
