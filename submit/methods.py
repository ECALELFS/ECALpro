from parameters import *
####from parameters_NEWESTCRAB import *

def printFillCfg1( outputfile ):
    outputfile.write("useHLTFilter = " + useHLTFilter + "\n")
    outputfile.write("correctHits = " + correctHits + "\n\n")
    outputfile.write('import FWCore.ParameterSet.Config as cms\n')
    outputfile.write('import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi\n')
    outputfile.write("import os, sys, imp, re\n")
    outputfile.write('CMSSW_VERSION=os.getenv("CMSSW_VERSION")\n')
    outputfile.write('process = cms.Process("analyzerFillEpsilon")\n')
    outputfile.write('process.load("FWCore.MessageService.MessageLogger_cfi")\n\n')
    outputfile.write('process.load("Configuration.Geometry.GeometryIdeal_cff")\n')

    if(globaltag_New):
       outputfile.write('process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")\n')
    else:
       outputfile.write('process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")\n')
    outputfile.write("process.GlobalTag.globaltag = '" + globaltag + "'\n")
    #From DIGI
    if (FROMDIGI):
        outputfile.write("#DUMMY RECHIT\n")
        outputfile.write("process.dummyHits = cms.EDProducer('DummyRechitDigis',\n")
        outputfile.write("                                     doDigi = cms.untracked.bool(True),\n")
        outputfile.write("                                     # rechits\n")                                                                                                
        outputfile.write("                                     barrelHitProducer      = cms.InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB'),\n")
        outputfile.write("                                     endcapHitProducer      = cms.InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE'),\n")
        outputfile.write("                                     barrelRecHitCollection = cms.untracked.string('dummyBarrelRechits'),\n")
        outputfile.write("                                     endcapRecHitCollection = cms.untracked.string('dummyEndcapRechits'),\n")
        outputfile.write("                                     # digis\n")                                                                                                                               
        outputfile.write("                                     barrelDigis            = cms." + EBdigi + ",\n")
        outputfile.write("                                     endcapDigis            = cms." + EEdigi + ",\n")
        outputfile.write("                                     barrelDigiCollection   = cms.untracked.string('dummyBarrelDigis'),\n")
        outputfile.write("                                     endcapDigiCollection   = cms.untracked.string('dummyEndcapDigis'))\n")
        outputfile.write("\n")
        if(FixGhostDigis):
            outputfile.write("# GHOST DIGIS CLEANER (FOR 2015 STREAM DATA)\n")
            outputfile.write("process.load(\"CalibCode.FillEpsilonPlot.cleanedDigiCollectionProducer_cfi\")\n")
            outputfile.write("process.cleanedEcalDigis.ebDigis = cms.InputTag('dummyHits','dummyBarrelDigis')\n")
            outputfile.write("process.cleanedEcalDigis.eeDigis = cms.InputTag('dummyHits','dummyEndcapDigis')\n")
            outputfile.write("\n")
        outputfile.write("#RAW to DIGI'\n")
        outputfile.write("#https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/RecoLocalCalo/EcalRecProducers/test/testMultipleEcalRecoLocal_cfg.py\n")
        outputfile.write("#process.load('Configuration.StandardSequences.RawToDigi_cff')\n")
        outputfile.write("#process.raw2digi_step = cms.Sequence(process.RawToDigi)\n")
        outputfile.write("#DIGI to UNCALIB\n")
        outputfile.write("process.load('Configuration.StandardSequences.Reconstruction_cff')\n")
        if(MULTIFIT):
           outputfile.write("import RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi\n")
           outputfile.write("process.ecalMultiFitUncalibRecHit =  RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHit_cfi.ecalMultiFitUncalibRecHit.clone()\n")
           if( DigiCustomization ):
               outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = cms.bool( False ) # To read the conditions from the header\n") 
           if( is50ns and DigiCustomization ):
               outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.activeBXs = cms.vint32(-4,-2,0,2,4) #Are 10 (-5-5). For 50ns is (-4,-2,0,2,4) #No .algoPSet. in old releases\n")
           if( not is50ns and DigiCustomization ):
               outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4) #Are 10 (-5-5). For 50ns is (-4,-2,0,2,4) #No .algoPSet. in old releases\n")
           if(FixGhostDigis):
               outputfile.write("process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('cleanedEcalDigis','ebCleanedDigis')\n")
               outputfile.write("process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('cleanedEcalDigis','eeCleanedDigis')\n")
           else:
               outputfile.write("process.ecalMultiFitUncalibRecHit.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','analyzerFillEpsilon')\n")
               outputfile.write("process.ecalMultiFitUncalibRecHit.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','analyzerFillEpsilon')\n")    
    
           outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = False #added this line to make code run\n") #can enable setting --> DigiCustomization = True <-- in parameters.py, but this also set --> outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.activeBXs = cms.vint32(-5,-4,-3,-2,-1,0,1,2,3,4) #Are 10 (-5-5). For 50ns is (-4,-2,0,2,4) \
#No .algoPSet. in old releases\n")  <-- line above , so I prefer to add it here
        if(WEIGHTS):
           outputfile.write("import RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi\n")
           outputfile.write("process.load('RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi')\n")
           outputfile.write("process.ecalweight =  RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi.ecalGlobalUncalibRecHit.clone()\n")
           if(FixGhostDigis):
               outputfile.write("process.ecalweight.EBdigiCollection = cms.InputTag('cleanedEcalDigis','ebCleanedDigis')\n")
               outputfile.write("process.ecalweight.EEdigiCollection = cms.InputTag('cleanedEcalDigis','eeCleanedDigis')\n")
           else:
               outputfile.write("process.ecalweight.EBdigiCollection = cms.InputTag('dummyHits','dummyBarrelDigis','analyzerFillEpsilon')\n")
               outputfile.write("process.ecalweight.EEdigiCollection = cms.InputTag('dummyHits','dummyEndcapDigis','analyzerFillEpsilon')\n")               
        outputfile.write("#UNCALIB to CALIB\n")
        outputfile.write("from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import *\n")
        outputfile.write("process.ecalDetIdToBeRecovered =  RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi.ecalDetIdToBeRecovered.clone()\n")
        outputfile.write("process.ecalRecHit.killDeadChannels = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEBVFE = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEEVFE = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEBFE = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEEFE = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEEIsolatedChannels = cms.bool( False )\n")
        outputfile.write("process.ecalRecHit.recoverEBIsolatedChannels = cms.bool( False )\n")
        if(WEIGHTS):
           outputfile.write("process.ecalRecHit.EEuncalibRecHitCollection =  cms.InputTag('ecalweight','EcalUncalibRecHitsEE')\n")
           outputfile.write("process.ecalRecHit.EBuncalibRecHitCollection =  cms.InputTag('ecalweight','EcalUncalibRecHitsEB')\n")
        outputfile.write("process.ecalLocalRecoSequence = cms.Sequence(ecalRecHit)\n")

    if (overWriteGlobalTag):        
        outputfile.write("process.GlobalTag.toGet = cms.VPSet(\n")
        if not(laserTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + laserTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + laserTag + "'),\n")
            outputfile.write("              connect = cms.string('" + laserDB + "')\n")
            outputfile.write('     ),\n')
        if not(alphaTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + alphaTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + alphaTag + "'),\n")
            outputfile.write("              connect = cms.string('" + alphaDB + "')\n")
            outputfile.write('     ),\n')
        if not(GeVTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + GeVTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + GeVTag + "'),\n")
            outputfile.write("              connect = cms.string('" + GeVDB + "')\n")
            outputfile.write('     ),\n')
        if not(pulseShapeTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + pulseShapeTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + pulseShapeTag + "'),\n")
            outputfile.write("              connect = cms.string('" + pulseShapeDB + "')\n")
            outputfile.write('     ),\n')
        if not(pedestalTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + pedestalTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + pedestalTag + "'),\n")
            outputfile.write("              connect = cms.string('" + pedestalDB + "')\n")
            outputfile.write('     ),\n')
        if not(laserAlphaTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + laserAlphaTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + laserAlphaTag + "'),\n")
            outputfile.write("              connect = cms.string('" + laserAlphaDB + "')\n")
            outputfile.write('     ),\n')
        if not(intercalibTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + intercalibTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + intercalibTag + "'),\n")
            outputfile.write("              connect = cms.string('" + intercalibDB + "')\n")
            outputfile.write('     ),\n')
        if not(ESEEIntercalibTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + ESEEIntercalibTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + ESEEIntercalibTag + "'),\n")
            outputfile.write("              connect = cms.string('" + ESEEIntercalibDB + "')\n")
            outputfile.write('     ),\n')
        if not(intercalibTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + intercalibTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + intercalibTag + "'),\n")
            outputfile.write("              connect = cms.string('" + intercalibDB + "')\n")
            outputfile.write('     ),\n')
        if not(linearCorrectionsTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + linearCorrectionsTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + linearCorrectionsTag + "'),\n")
            outputfile.write("              connect = cms.string('" + linearCorrectionsDB + "')\n")
            outputfile.write('     ),\n')
        outputfile.write(')\n\n')

    outputfile.write('### Recalibration Module to apply laser corrections on the fly\n')
    outputfile.write('if correctHits:\n')
    outputfile.write('    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(\n')
    outputfile.write('        doEnergyScale = cms.bool(' + doEnenerScale + '),\n')
    outputfile.write('        doIntercalib = cms.bool(' + doIC + '),\n')
    outputfile.write('        doLaserCorrections = cms.bool(' + doLaserCorr + '),\n')
    outputfile.write("        EBRecHitCollection = cms." + ebInputTag +",\n")
    outputfile.write("        EERecHitCollection = cms." + eeInputTag +",\n")
    outputfile.write('        EBRecalibRecHitCollection = cms.string("pi0EcalRecHitsEB"),\n')
    outputfile.write('        EERecalibRecHitCollection = cms.string("pi0EcalRecHitsEE")\n')
    outputfile.write('    )\n\n')

    outputfile.write('### Running on AlcaRAW requires filtering AlcaPi0 events from AlcaEta events\n')
    outputfile.write('if useHLTFilter:\n')
    outputfile.write('    import copy\n')
    outputfile.write('    from HLTrigger.HLTfilters.hltHighLevel_cfi import *\n')
    outputfile.write('    process.AlcaP0Filter = copy.deepcopy(hltHighLevel)\n')
    outputfile.write('    process.AlcaP0Filter.throw = cms.bool(False)\n')
    if( triggerTag!='' ):
        outputfile.write('    process.AlcaP0Filter.TriggerResultsTag = cms.' + triggerTag + '\n')
    outputfile.write('    process.AlcaP0Filter.HLTPaths = ["' + HLTPaths + '"]\n\n')

    outputfile.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(" + nEventsPerJob +") )\n")
    outputfile.write("process.MessageLogger.cerr.FwkReport.reportEvery = 100000\n")
    # outputfile.write("process.MessageLogger.cerr = cms.untracked.PSet(\n")
    # outputfile.write("        threshold  = cms.untracked.string('WARNING'),\n")
    # outputfile.write("        ERROR      = cms.untracked.PSet (\n")
    # outputfile.write("                                         limit = cms.untracked.int32(1)\n")
    # outputfile.write("        )\n")
    # outputfile.write(")\n")
    outputfile.write("process.options = cms.untracked.PSet(\n")
    outputfile.write("   wantSummary = cms.untracked.bool(True),\n")
    outputfile.write("   SkipEvent = cms.untracked.vstring('ProductNotFound','CrystalIDError')\n")
    outputfile.write(")\n")
    outputfile.write("process.source = cms.Source('PoolSource',\n")
    #outputfile.write("                            inputCommands = cms.untracked.vstring( #type_Module_instance_process\n")
    #outputfile.write("                                'drop *',\n")
    #outputfile.write("                                'keep E*DigiCollection_hltAlCa*RechitsToDigis_*_*',\n")
    #outputfile.write("                                'keep *EcalRecHit*_hltAlCa*RecHitsFilterEEonlyRegional_*_*',\n")
    #outputfile.write("                                'keep L1GlobalTriggerReadoutRecord_*_*_*',\n")
    #outputfile.write("                                'keep *TriggerResults_*_*_*'\n")
    #outputfile.write("                            ),\n")
    outputfile.write("    fileNames = cms.untracked.vstring(\n")

def printFillCfg2( outputfile, pwd , iteration, outputDir, ijob ):
    outputfile.write("    )\n")
    outputfile.write(")\n")
#    outputfile.write("\n")
#    if(len(json_file)>0):
#       outputfile.write('if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):\n')
#       outputfile.write("   import FWCore.PythonUtilities.LumiList as LumiList\n")
#       if (isCRAB):
#           outputfile.write("   process.source.lumisToProcess = LumiList.LumiList(filename = 'CalibCode/FillEpsilonPlot/data/" + json_file + "').getVLuminosityBlockRange()\n")
#       else:
#           outputfile.write("   process.source.lumisToProcess = LumiList.LumiList(filename = '" + pwd + "/../../CalibCode/FillEpsilonPlot/data/" + json_file + "').getVLuminosityBlockRange()\n")
#       outputfile.write("else:\n")
#       outputfile.write("   import PhysicsTools.PythonAnalysis.LumiList as LumiList\n")
#       if (isCRAB):
#           outputfile.write("   myLumis = LumiList.LumiList(filename = 'CalibCode/FillEpsilonPlot/data/" + json_file + "').getCMSSWString().split(',')\n")
#       else:
#           outputfile.write("   myLumis = LumiList.LumiList(filename = '" + pwd + "/../../CalibCode/FillEpsilonPlot/data/" + json_file + "').getCMSSWString().split(',')\n")
#       outputfile.write("   process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()\n")
#       outputfile.write("   process.source.lumisToProcess.extend(myLumis)\n")
    outputfile.write("\n")
    outputfile.write("process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')\n")
    outputfile.write("process.analyzerFillEpsilon.OutputDir = cms.untracked.string('" +  outputDir + "')\n")
    outputfile.write("process.analyzerFillEpsilon.OutputFile = cms.untracked.string('" + NameTag +  outputFile + "_" + str(ijob) + ".root')\n")
    outputfile.write("process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + ExternalGeometry + "')\n")
    if (isCRAB):
        outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + NameTag + calibMapName + "')\n")
        outputfile.write("process.analyzerFillEpsilon.isCRAB  = cms.untracked.bool(True)\n")
    else:
        if (SubmitFurtherIterationsFromExisting and iteration == 0):
            outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + startingCalibMap + "')\n")
        else:
            outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")
    outputfile.write("process.analyzerFillEpsilon.useEBContainmentCorrections = cms.untracked.bool(" + useEBContainmentCorrections + ")\n")
    outputfile.write("process.analyzerFillEpsilon.useEEContainmentCorrections = cms.untracked.bool(" + useEEContainmentCorrections + ")\n")
    outputfile.write("process.analyzerFillEpsilon.EBContainmentCorrections = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + EBContainmentCorrections + "')\n")
    if (new_pi0ContainmentCorrections):
	    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + new_MVAEBContainmentCorrections_01 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + new_MVAEBContainmentCorrections_02 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEEContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + new_MVAEEContainmentCorrections_01 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEEContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + new_MVAEEContainmentCorrections_02 + "')\n")
    else:
	    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEBContainmentCorrections_01 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEBContainmentCorrections_02 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEEContainmentCorrections_01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEEContainmentCorrections_01 + "')\n")
	    outputfile.write("process.analyzerFillEpsilon.MVAEEContainmentCorrections_02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEEContainmentCorrections_02 + "')\n")

    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_eta01  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEBContainmentCorrections_eta01 + "')\n")
    outputfile.write("process.analyzerFillEpsilon.MVAEBContainmentCorrections_eta02  = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + MVAEBContainmentCorrections_eta02 + "')\n")
    outputfile.write("process.analyzerFillEpsilon.Endc_x_y                        = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + Endc_x_y + "')\n")
    outputfile.write("process.analyzerFillEpsilon.EBPHIContainmentCorrections = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + EBPHIContainmentCorrections + "')\n")
    outputfile.write("process.analyzerFillEpsilon.EEContainmentCorrections    = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + EEContainmentCorrections + "')\n")
    outputfile.write("process.analyzerFillEpsilon.ContCorr_EB                 = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + EBContCorr + "')\n")
    #outputfile.write("process.analyzerFillEpsilon.json_file                   = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + json_file + "')\n")
    outputfile.write("process.analyzerFillEpsilon.HLTResults                  = cms.untracked.bool(" + HLTResults + ")\n")
    if(HLTResultsNameEB!=""):
        outputfile.write("process.analyzerFillEpsilon.HLTResultsNameEB            = cms.untracked.string('" + HLTResultsNameEB + "')\n")
    if(HLTResultsNameEE!=""):
        outputfile.write("process.analyzerFillEpsilon.HLTResultsNameEE            = cms.untracked.string('" + HLTResultsNameEE + "')\n")
    outputfile.write("process.analyzerFillEpsilon.RemoveDead_Flag             = cms.untracked.bool(" + RemoveDead_Flag + ")\n")
    outputfile.write("process.analyzerFillEpsilon.RemoveDead_Map              = cms.untracked.string('" + RemoveDead_Map + "')\n")
    if(EtaRingCalibEB):
      outputfile.write("process.analyzerFillEpsilon.EtaRingCalibEB    = cms.untracked.bool(True)\n")
    if(EtaRingCalibEE):
      outputfile.write("process.analyzerFillEpsilon.EtaRingCalibEE    = cms.untracked.bool(True)\n")
    if(SMCalibEB):
      outputfile.write("process.analyzerFillEpsilon.SMCalibEB    = cms.untracked.bool(True)\n")
    if(SMCalibEE):
      outputfile.write("process.analyzerFillEpsilon.SMCalibEE    = cms.untracked.bool(True)\n")
    if(EtaRingCalibEB or SMCalibEB or EtaRingCalibEE or SMCalibEE):
      outputfile.write("process.analyzerFillEpsilon.CalibMapEtaRing = cms.untracked.string('" + CalibMapEtaRing + "')\n")
    if(MC_Asssoc):
        outputfile.write("process.analyzerFillEpsilon.GenPartCollectionTag = cms.untracked." + genPartInputTag + "\n")
        outputfile.write("process.analyzerFillEpsilon.MC_Asssoc            = cms.untracked.bool(True)\n")
        outputfile.write("process.analyzerFillEpsilon.MC_Asssoc_DeltaR     = cms.untracked.double(0.3)\n")
    if(Are_pi0):
        outputfile.write("process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(False)\n")
    
    if(useMVAContainmentCorrections):
        outputfile.write("process.analyzerFillEpsilon.useMVAContainmentCorrections                 = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.analyzerFillEpsilon.useMVAContainmentCorrections                 = cms.untracked.bool(False)\n")

    if(new_pi0ContainmentCorrections):
         outputfile.write("process.analyzerFillEpsilon.new_pi0ContainmentCorrections                 = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.analyzerFillEpsilon.new_pi0ContainmentCorrections                 = cms.untracked.bool(False)\n")

    
    outputfile.write("process.analyzerFillEpsilon.useOnlyEEClusterMatchedWithES = cms.untracked.bool(" + useOnlyEEClusterMatchedWithES + ")\n\n")

    outputfile.write("### choosing proper input tag (recalibration module changes the collection names)\n")
    outputfile.write("if correctHits:\n")
    outputfile.write("    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')\n")
    outputfile.write("    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')\n")
    outputfile.write("else:\n")
    outputfile.write("    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked." + ebInputTag + "\n")
    outputfile.write("    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked." + eeInputTag + "\n")
    outputfile.write("process.analyzerFillEpsilon.ESRecHitCollectionTag = cms.untracked." + esInputTag + "\n")
    #outputfile.write("process.analyzerFillEpsilon.l1InputTag = cms.untracked." + l1InputTag + "\n")

    outputfile.write("process.analyzerFillEpsilon.L1TriggerTag = cms.untracked." + hltGtDigis + "\n")
    outputfile.write("process.analyzerFillEpsilon.triggerTag   = cms.untracked." + triggerTag + "\n")
    outputfile.write("process.analyzerFillEpsilon.hltL1GtObjectMap   = cms.untracked." + hltL1GtObjectMap + "\n")
    outputfile.write("process.analyzerFillEpsilon.CalibType    = cms.untracked.string('" + CalibType + "')\n")
    outputfile.write("process.analyzerFillEpsilon.CurrentIteration = cms.untracked.int32(" + str(iteration) + ")\n")
    if( EB_Seed_E!='' ):
        outputfile.write("process.analyzerFillEpsilon.EB_Seed_E = cms.untracked.double(" + EB_Seed_E + ")\n")
    if( useEE_EtSeed!='' ):
        outputfile.write("process.analyzerFillEpsilon.useEE_EtSeed = cms.untracked.bool(" + useEE_EtSeed + ")\n")
    if( EE_Seed_E!='' ):
        outputfile.write("process.analyzerFillEpsilon.EE_Seed_E = cms.untracked.double(" + EE_Seed_E + ")\n")
    if( EE_Seed_Et!='' ):
        outputfile.write("process.analyzerFillEpsilon.EE_Seed_Et = cms.untracked.double(" + EE_Seed_Et + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0PtCutEB_low = cms.untracked.double(" + Pi0PtCutEB_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0PtCutEB_high = cms.untracked.double(" + Pi0PtCutEB_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0PtCutEE_low = cms.untracked.double(" + Pi0PtCutEE_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0PtCutEE_high = cms.untracked.double(" + Pi0PtCutEE_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.gPtCutEB_low = cms.untracked.double(" + gPtCutEB_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.gPtCutEB_high = cms.untracked.double(" + gPtCutEB_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.gPtCutEE_low = cms.untracked.double(" + gPtCutEE_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.gPtCutEE_high = cms.untracked.double(" + gPtCutEE_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0IsoCutEB_low = cms.untracked.double(" + Pi0IsoCutEB_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0IsoCutEB_high = cms.untracked.double(" + Pi0IsoCutEB_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0IsoCutEE_low = cms.untracked.double(" + Pi0IsoCutEE_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0IsoCutEE_high = cms.untracked.double(" + Pi0IsoCutEE_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.CutOnHLTIso = cms.untracked.bool(" + CutOnHLTIso + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0HLTIsoCutEB_low = cms.untracked.double(" + Pi0HLTIsoCutEB_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0HLTIsoCutEB_high = cms.untracked.double(" + Pi0HLTIsoCutEB_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0HLTIsoCutEE_low = cms.untracked.double(" + Pi0HLTIsoCutEE_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Pi0HLTIsoCutEE_high = cms.untracked.double(" + Pi0HLTIsoCutEE_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_1_EB_low = cms.untracked.double(" +  nXtal_1_EB_low+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_1_EB_high = cms.untracked.double(" +  nXtal_1_EB_high+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_2_EB_low = cms.untracked.double(" +  nXtal_2_EB_low+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_2_EB_high = cms.untracked.double(" +  nXtal_2_EB_high+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_1_EE_low = cms.untracked.double(" +  nXtal_1_EE_low+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_1_EE_high = cms.untracked.double(" +  nXtal_1_EE_high+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_2_EE_low = cms.untracked.double(" +  nXtal_2_EE_low+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.nXtal_2_EE_high = cms.untracked.double(" +  nXtal_2_EE_high+ ")\n")
    outputfile.write("process.analyzerFillEpsilon.S4S9_EB_low = cms.untracked.double(" + S4S9_EB_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.S4S9_EB_high = cms.untracked.double(" + S4S9_EB_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.S4S9_EE_low = cms.untracked.double(" + S4S9_EE_low + ")\n")
    outputfile.write("process.analyzerFillEpsilon.S4S9_EE_high = cms.untracked.double(" + S4S9_EE_high + ")\n")
    outputfile.write("process.analyzerFillEpsilon.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
    if(len(json_file)>0):
       outputfile.write("process.analyzerFillEpsilon.JSONfile = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + json_file + "')\n")
    if GeometryFromFile:
       outputfile.write("process.analyzerFillEpsilon.GeometryFromFile = cms.untracked.bool(True)\n")
    if isMC:
       outputfile.write("process.analyzerFillEpsilon.isMC = cms.untracked.bool(True)\n")
    if MakeNtuple4optimization:
       outputfile.write("process.analyzerFillEpsilon.MakeNtuple4optimization = cms.untracked.bool(True)\n")
    if( L1TriggerInfo ):
        outputfile.write("process.analyzerFillEpsilon.L1TriggerInfo = cms.untracked.bool(True)\n")
    if not( L1Seed=='' ):
        outputfile.write("process.analyzerFillEpsilon.L1_Bit_Sele = cms.untracked.string('" + L1Seed + "')\n")
    outputfile.write("process.p = cms.Path()\n")
    outputfile.write("if useHLTFilter:\n")
    outputfile.write("    process.p *= process.AlcaP0Filter\n")
    outputfile.write("if correctHits:\n")
    outputfile.write("    print 'ADDING RECALIB RECHIT MODULE WITH PARAMETERS'\n")
    outputfile.write("    print 'ENERGY SCALE '+str(process.ecalPi0ReCorrected.doEnergyScale)\n")
    outputfile.write("    print 'INTERCALIBRATION '+str(process.ecalPi0ReCorrected.doIntercalib)\n")
    outputfile.write("    print 'LASER '+str(process.ecalPi0ReCorrected.doLaserCorrections)\n")
    outputfile.write("    process.p *= process.ecalPi0ReCorrected\n")
    if (FROMDIGI):
        outputfile.write("process.p *= process.dummyHits\n")
        if(FixGhostDigis):
            outputfile.write("process.p *= process.cleanedEcalDigis\n")
        if(MULTIFIT):
           outputfile.write("process.p *= process.ecalMultiFitUncalibRecHit\n")
        if (WEIGHTS):
           outputfile.write("process.p *= process.ecalweight\n")
        outputfile.write("process.p *= process.ecalLocalRecoSequence\n")
    outputfile.write("process.p *= process.analyzerFillEpsilon\n")

def printFitCfg( outputfile, iteration, outputDir, nIn, nFin, EBorEE, nFit ):
    outputfile.write("import FWCore.ParameterSet.Config as cms\n")
    outputfile.write("process = cms.Process('FitEpsilonPlot')\n")
    outputfile.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
    outputfile.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )\n")
    outputfile.write("process.source =   cms.Source('EmptySource')\n")
    outputfile.write("process.fitEpsilon = cms.EDAnalyzer('FitEpsilonPlot')\n")
    outputfile.write("process.fitEpsilon.OutputFile = cms.untracked.string('" + NameTag + EBorEE + "_" + str(nFit) + "_" + calibMapName + "')\n")
    outputfile.write("process.fitEpsilon.CalibType = cms.untracked.string('" + CalibType + "')\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        outputfile.write("process.fitEpsilon.OutputDir = cms.untracked.string('$TMPDIR')\n")
    else:
        outputfile.write("process.fitEpsilon.OutputDir = cms.untracked.string('" +  outputDir + "')\n")
    outputfile.write("process.fitEpsilon.CurrentIteration = cms.untracked.int32(" + str(iteration) + ")\n")
    outputfile.write("process.fitEpsilon.NInFit = cms.untracked.int32(" + str(nIn) + ")\n")
    outputfile.write("process.fitEpsilon.NFinFit = cms.untracked.int32(" + str(nFin) + ")\n")
    outputfile.write("process.fitEpsilon.EEorEB = cms.untracked.string('" + EBorEE + "')\n")
    outputfile.write("process.fitEpsilon.isNot_2010 = cms.untracked.bool(" + isNot_2010 + ")\n")
    if(Are_pi0):
        outputfile.write("process.fitEpsilon.Are_pi0 = cms.untracked.bool( True )\n")
    else:
        outputfile.write("process.fitEpsilon.Are_pi0 = cms.untracked.bool( False )\n")

    if(useMVAContainmentCorrections):
        outputfile.write("process.fitEpsilon.useMVAContainmentCorrections = cms.untracked.bool( True )\n")
    else:
        outputfile.write("process.fitEpsilon.useMVAContainmentCorrections = cms.untracked.bool( False )\n")
    
    if(new_pi0ContainmentCorrections):
        outputfile.write("process.fitEpsilon.new_pi0ContainmentCorrections = cms.untracked.bool( True )\n")
    else:
        outputfile.write("process.fitEpsilon.new_pi0ContainmentCorrections = cms.untracked.bool( False )\n")

    outputfile.write("process.fitEpsilon.StoreForTest = cms.untracked.bool( True )\n")
    outputfile.write("process.fitEpsilon.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
    if not(isCRAB): #If CRAB you have to put the correct path, and you do it on calibJobHandler.py, not on ./submitCalibration.py
        outputfile.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration) + "/" + NameTag + "epsilonPlots.root')\n")
        if (SubmitFurtherIterationsFromExisting and iteration == 0):
            outputfile.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + startingCalibMap + "')\n")
        else:
            outputfile.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")
    outputfile.write("process.p = cms.Path(process.fitEpsilon)\n")


def printSubmitFitSrc(outputfile, cfgName, source, destination, pwd, logpath):
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
        outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
        # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
        outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/{quote}' > " + logpath  + "\n")
    outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/' >> " + logpath  + "\n")
    outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
    outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
    if not(isCRAB): #If CRAB you have to put the correct path, anbd you do it on calibJobHandler.py, not on ./submitCalibration.py
       destrooplot = destination.replace("calibMap","fitRes")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
           outputfile.write("echo 'eos cp " + source + " " + destination + "' >> " + logpath  + "\n")
           outputfile.write("echo 'eos cp /tmp/Fit_Stored.root " + destrooplot + "' >> " + logpath  + "\n")
       else:
           outputfile.write("echo 'cmsStage -f " + source + " " + destination + "' >> " + logpath  + "\n")           
           outputfile.write("echo 'cmsStage -f /tmp/Fit_Stored.root " + destrooplot + "' >> " + logpath  + "\n")
       outputfile.write(myeosstage + source + " " + destination + " >> " + logpath + " 2>&1 \n")
       outputfile.write(myeosstage + " /tmp/Fit_Stored.root " + destrooplot + " >> " + logpath + " 2>&1 \n")
       outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
       outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")
       outputfile.write("echo 'rm -f /tmp/Fit_Stored.root' >> " + logpath + " \n")
       outputfile.write("rm -f /tmp/Fit_Stored.root >> " + logpath + " 2>&1 \n")

def printSubmitSrc(outputfile, cfgName, source, destination, pwd, logpath):
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n") this line produces problem when running in CMSSW_8_0_3, anyway we don't use crab
    if ( not isOtherT2 and isCRAB ):
        outputfile.write("setenv X509_USER_PROXY " + CRAB_CopyCert + "\n")
    if not(Silent):
        outputfile.write("echo 'cmsRun " + cfgName + "'\n")
        outputfile.write("cmsRun " + cfgName + "\n")
        if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
            outputfile.write("echo 'eos cp " + source + " " + destination + "'\n")
        else:
            outputfile.write("echo 'cmsStage -f " + source + " " + destination + "'\n")
        outputfile.write(myeosstage + source + " " + destination + "\n")
        outputfile.write("echo 'rm -f " + source + "'\n")
        outputfile.write("rm -f " + source + "\n")
    else:
        outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FILL_COUT:/{quote}' > " + logpath  + "\n")
        outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FILL_COUT:/' >> " + logpath  + "\n")
        outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
        outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
        if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
            outputfile.write("echo 'eos cp " + source + " " + destination + "' >> " + logpath  + "\n")
        else:
            outputfile.write("echo 'cmsStage -f " + source + " " + destination + "' >> " + logpath  + "\n")
        outputfile.write(myeosstage + source + " " + destination + " >> " + logpath + " 2>&1 \n")
        outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
        outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")

def printCrab(outputfile, iter):
    #outputfile.write("[CMSSW]\n")
    #outputfile.write("pset=fillEpsilonPlot_iter_" + str(iter) + ".py\n")
    #outputfile.write("events_per_job=" + events_per_job + "\n")
    #outputfile.write("total_number_of_events="+ total_number_of_events +"\n")
    #outputfile.write("datasetpath=" + CRAB_Data_Path + "\n")
    #outputfile.write("output_file=" + NameTag + outputFile + "_0.root\n")
    #outputfile.write("\n")
    #outputfile.write("[USER]\n")
    #outputfile.write("ui_working_dir=" + dirname + "_iter_" + str(iter) + "_CRAB\n")
    #outputfile.write("return_data=0\n")
    #outputfile.write("copy_data=1\n")
    #outputfile.write("storage_element = srm-eoscms.cern.ch\n")
    #outputfile.write("storage_path=/srm/v2/server?SFN=/eos/cms/store\n")
    #outputfile.write("user_remote_dir=" + CRAB_Storage + dirname + "/iter_" + str(iter) + "\n")
    #outputfile.write("check_user_remote_dir=0\n")
    #outputfile.write("\n")
    #outputfile.write("[CRAB]\n")
    #outputfile.write("\n")
    #outputfile.write("scheduler=remoteGlidein\n")
    #outputfile.write("jobtype=cmssw\n")
    ##outputfile.write("from CRABClient.client_utilities import getBasicConfig\n")
    ##outputfile.write("config = getBasicConfig()\n")
    outputfile.write("from CRABClient.UserUtilities import config\n")
    outputfile.write("config = config()\n")
    outputfile.write("config.General.requestName = 'CRAB_Folder'\n")
    outputfile.write("config.General.workArea = 'crab_projects'\n")
    outputfile.write("config.General.transferLogs = True\n")
    outputfile.write("config.JobType.pluginName = 'Analysis'\n")
    outputfile.write("config.JobType.psetName = 'fillEpsilonPlot_iter_" + str(iter) + ".py'\n")
    outputfile.write("config.Data.inputDataset = '" + CRAB_Data_Path + "'\n")
    outputfile.write("config.Data.splitting = 'FileBased'\n")
    outputfile.write("config.Data.unitsPerJob = " + str(unitsPerJob) + "\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" ):
       outputfile.write("config.Data.outLFNDirBase = '" + outLFN + "/iter_" + str(iter) + "'\n")
       outputfile.write("config.User.voGroup = '" + voGroup + "'\n") #Only needed from lxplus, not from m-machines
    else:
       outputfile.write("config.Data.inputDBS = 'global'\n")
       outputfile.write("config.Data.outLFNDirBase = '" + eosPath + "/" + dirname + "/" + "iter_" + str(iter) + "/'\n")
    outputfile.write("config.Site.storageSite = '" + storageSite + "'\n")
    outputfile.write("config.JobType.outputFiles = ['" + NameTag + outputFile + "_0.root']\n")
    outputfile.write("config.Data.publication = False\n")

def printCrabHadd(outputfile, iter, pwd):
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("#qsub -q localgrid@cream02 " + pwd + "/" + dirname + "/CRAB_files/HaddSendafterCrab_" + iter + ".sh\n")
       outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
       outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
       # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
       outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    else:
       outputfile.write("#bsub -q " + queueForDaemon + " 'bash " + pwd + "/" + dirname + "/CRAB_files/HaddSendafterCrab_" + iter + ".sh'\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("AddPath='putPATHhere' #Use path1~path2 if you have more folder from CRAB. The ICs will go to the 1st path\n")
    outputfile.write("AddPathOLDIter='putPATHhere' #Empty if iter 0, Additional path of previous iter, otherwise\n")
    outputfile.write("echo 'python calibJobHandler.py CRAB " + iter + " " + queue + "' $AddPath $AddPathOLDIter\n")
    outputfile.write("if [ '$AddPath' == 'putPATHhere' -o '$AddPathOLDIter' == 'putPATHhere' ]; then\n")
    outputfile.write("   echo 'Wrong Use of HaddSendafterCrab_X.sh, add the additional path of the CRAB output'\n")
    outputfile.write("else\n")
    outputfile.write("   python calibJobHandler.py CRAB " + iter + " " + queue + " $AddPath $AddPathOLDIter;\n")
    outputfile.write("fi\n")

def printParallelHadd(outputfile, outFile, list, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
       outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
       # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
       outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("echo 'hadd -f -k $TMPDIR/" + outFile + " @" + list + "'\n")
       outputfile.write("hadd -f -k $TMPDIR/" + outFile + " @" + list  + "\n")
       outputfile.write("echo 'srmcp file:///$TMPDIR/" + outFile + " " + destination + "/" + outFile + "'\n")
       outputfile.write("srmcp file:///$TMPDIR/" + outFile + " " + destination + "/" + outFile + "\n")
       outputfile.write("rm -f $TMPDIR/" + outFile + "\n")
    else:
       outputfile.write("echo 'hadd -f -k /tmp/" + outFile + " @" + list + "'\n")
       outputfile.write("hadd -f -k /tmp/" + outFile + " @" + list  + "\n")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:       
           outputfile.write("echo 'eos cp /tmp/" + outFile + " " + destination + "'\n")
       else:
           outputfile.write("echo 'cmsStage -f /tmp/" + outFile + " " + destination + "'\n")
       outputfile.write(myeosstage + "/tmp/" + outFile + " " + destination + "\n")
       outputfile.write("rm -f /tmp/" + outFile + "\n")

def printFinalHadd(outputfile, list, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("export SCRAM_ARCH=slc6_amd64_gcc491\n")
       outputfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
       # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n")
       outputfile.write("export X509_USER_PROXY=/localgrid/lpernie/x509up_u20580\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       outputfile.write("echo 'hadd -f -k $TMPDIR/" + NameTag + "epsilonPlots.root @" + list + "'\n")
       outputfile.write("hadd -f -k $TMPDIR/" + NameTag + "epsilonPlots.root @" + list  + "\n")
       outputfile.write("echo 'srmcp file:///$TMPDIR/" + NameTag + "epsilonPlots.root " + destination + "/epsilonPlots.root" + "'\n")
       outputfile.write("srmcp file:///$TMPDIR/" + NameTag + "epsilonPlots.root " + destination + "/epsilonPlots.root" + "\n")
       outputfile.write("rm -f $TMPDIR/" + NameTag + "epsilonPlots.root\n")
    else:
       outputfile.write("echo 'hadd -f -k /tmp/" + NameTag + "epsilonPlots.root @" + list + "'\n")
       outputfile.write("hadd -f -k /tmp/" + NameTag + "epsilonPlots.root @" + list  + "\n")
       if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:       
           outputfile.write("echo 'eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "'\n")
       else:
           outputfile.write("echo 'cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "'\n")
       outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
       outputfile.write("rm -f /tmp/" + NameTag + "epsilonPlots.root\n")

# following function is copied and modified just after. Keep it commented for reference for now

# def printParallelHaddFAST(outputfile, outFile, listReduced, destination, pwd, numList):
#     import os, sys, imp, re
#     CMSSW_VERSION=os.getenv("CMSSW_VERSION")
#     outputfile.write("#!/bin/bash\n")
#     if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
#          print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
#          print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
#          outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
#     else:
#          outputfile.write("cd " + pwd + "\n")
#     outputfile.write("eval `scramv1 runtime -sh`\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")
# #if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
#     if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
#         outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms/eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"" + myeosstage + "root://eoscms/eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
#     elif "cmsStage -f" in myeosstage:
#         outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
#     else:
#         outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
#     outputfile.write("files=`cat " + listReduced + "`\n")
#     outputfile.write("filesHadd=''\n")
#     outputfile.write("for file in $files;\n")
#     outputfile.write("do\n")
#     if( isCRAB ):
#         outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
#     else:
#         outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $10 }'`\n")  # since I added a directory on eos, must print 10th position, not 9th
#     outputfile.write("   echo \"-> outEncode=$( { fastHadd encode -o /tmp/${SUBSTRING}.pb /tmp/${SUBSTRING}; } 2>&1 )\"\n")
#     outputfile.write("   outEncode=$( { fastHadd encode -o /tmp/${SUBSTRING}.pb /tmp/${SUBSTRING}; } 2>&1 )\n")
#     outputfile.write('   echo "outEncode is: $outEncode"\n')
#     outputfile.write('   if [[ "$outEncode" =~ "DEBUG: Encoding" ]]\n')
#     outputfile.write('   then\n')
#     outputfile.write('      if [[ ! "$outEncode" =~ "Error" ]]\n')
#     outputfile.write('      then\n')
#     outputfile.write('        filesHadd="$filesHadd /tmp/${SUBSTRING}.pb"\n')
#     outputfile.write("      fi\n")
#     outputfile.write("   fi\n")
#     outputfile.write("done\n")
#     outputfile.write("echo \"fastHadd add -o /tmp/" + NameTag + "FinalFile.pb $filesHadd\"\n")
#     outputfile.write("fastHadd add -o /tmp/" + NameTag + "FinalFile.pb $filesHadd\n")
#     outputfile.write("fastHadd convert -o /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root /tmp/" + NameTag + "FinalFile.pb\n")
#     if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
#         outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
#     else:
#         outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
#     outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")

#############################################################

# copying and modifying previous method

def printParallelHaddFAST(outputfile, outFile, listReduced, destination, pwd, numList):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
#if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms//eos/cms\\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "root://eoscms//eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    elif "cmsStage -f" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\"  /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
    else:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\"  /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("filesHadd=''\n")
    outputfile.write("for file in $files;\n")
    outputfile.write("do\n")
    if( isCRAB ):
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
    else:
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $10 }'`\n")  # since I added a directory on eos, must print 10th position, not 9th
    outputfile.write('   filesHadd="$filesHadd /tmp/$SUBSTRING"\n')
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -k /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root $filesHadd\"\n")
    outputfile.write("hadd -k /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root $filesHadd\n")
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots_" + str(numList) + ".root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + outputFile + "_*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")


####################################

# commented as previous function. Will use modifyed version below

# def printFinalHaddFAST(outputfile, listReduced, destination, pwd):
#     import os, sys, imp, re
#     CMSSW_VERSION=os.getenv("CMSSW_VERSION")
#     outputfile.write("#!/bin/bash\n")
#     if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
#          print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
#          print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
#          outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
#     else:
#          outputfile.write("cd " + pwd + "\n")
#     outputfile.write("eval `scramv1 runtime -sh`\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")
# #if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
#     if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
#         outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms//eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"" + myeosstage + "root://eoscms//eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
#     elif "cmsStage -f" in myeosstage:
#         outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
#     else:
#         outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\"\n")
#         outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
#     outputfile.write("files=`cat " + listReduced + "`\n")
#     outputfile.write("filesHadd=''\n")
#     outputfile.write("for file in $files;\n")
#     outputfile.write("do\n")
#     if( isCRAB ):
#         outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
#     else:
#         outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $10 }'`\n")   # since I added a directory on eos, must print 10th position, not 9th
#     outputfile.write("   echo \"-> outEncode=$( { fastHadd encode -o /tmp/${SUBSTRING}.pb /tmp/${SUBSTRING}; } 2>&1 )\"\n")
#     outputfile.write("   outEncode=$( { fastHadd encode -o /tmp/${SUBSTRING}.pb /tmp/${SUBSTRING}; } 2>&1 )\n")
#     outputfile.write('   echo "outEncode is: $outEncode"\n')
#     outputfile.write('   if [[ "$outEncode" =~ "DEBUG: Encoding" ]]\n')
#     outputfile.write('   then\n')
#     outputfile.write('      if [[ ! "$outEncode" =~ "Error" ]]\n')
#     outputfile.write('      then\n')
#     outputfile.write('        filesHadd="$filesHadd /tmp/${SUBSTRING}.pb"\n')
#     outputfile.write("      fi\n")
#     outputfile.write("   fi\n")
#     outputfile.write("done\n")
#     outputfile.write("echo \"fastHadd add -o /tmp/" + NameTag + "FinalFile.pb $filesHadd\"\n")
#     outputfile.write("fastHadd add -o /tmp/" + NameTag + "FinalFile.pb $filesHadd\n")
#     outputfile.write("fastHadd convert -o /tmp/" + NameTag + "epsilonPlots.root /tmp/" + NameTag + "FinalFile.pb\n")
#     if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
#         outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
#     else:
#         outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
#     outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
#     outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")


##############################################################

def printFinalHaddFAST(outputfile, listReduced, destination, pwd):
    import os, sys, imp, re
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    if(re.match("CMSSW_5_.*_.*",CMSSW_VERSION)):
         print "WARNING!!!! ----> I'm ging to use a harcoded path: /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src"
         print "This because you are in a release CMSSW_5_*_*, that do not allow a hadd with a @file.list."
         outputfile.write("cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_4_2_4/src\n")
    else:
         outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    #outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")  # removing this line, we don't have a file named like this (we had it in another function, not anymore)
#if we leave "cmsStage -f" to cpy file from eos to /tmp, then ok, otherwise, with "eos cp" files on eos must be preceeded by "root://eoscms//eos/cms". In the lines below $0 is a file read from listreduced, which will be of the form /store/blabla/file.root 
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"eos cp root://eoscms//eos/cms\\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "root://eoscms//eos/cms\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    elif "cmsStage -f" in myeosstage:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"" + myeosstage + "\"$0 \" /tmp/\"}' " + listReduced + " | bash\n")        
    else:
        outputfile.write("echo \"Copying files locally: awk '{print \"cmsStage -f \\\"$\\0 \\\" /tmp/\"}' " + listReduced + " | bash\"\n")
        outputfile.write("awk '{print \"cmsStage -f \"$0 \" /tmp/\"}' " + listReduced + " | bash\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("filesHadd=''\n")
    outputfile.write("for file in $files;\n")
    outputfile.write("do\n")
    if( isCRAB ):
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $14 }'`\n")
    else:
        outputfile.write("   SUBSTRING=`echo ${file} | awk -F / '{ print $10 }'`\n")   # since I added a directory on eos, must print 10th position, not 9th
    outputfile.write('   filesHadd="$filesHadd /tmp/$SUBSTRING"\n')
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -k /tmp/" + NameTag + "epsilonPlots.root $filesHadd\"\n")
    outputfile.write("hadd -k /tmp/" + NameTag + "epsilonPlots.root $filesHadd\n")
    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")



def printFinalHaddRegroup(outputfile, listReduced, destination, pwd, grouping=10):
    import os, sys, imp, re, ntpath
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    #outputfile.write("rm -rf /tmp/" + NameTag + "FinalFile*\n")  # removing this line, we don't have a file named like this (we had it in another function, not anymore)
    
    fileWithList = open(listReduced,"r")
    files = fileWithList.readlines()
    idx=0
    grouped_files = []
    while len(files)>0:
        filesToMerge = files[:grouping]
        mergedfile_n = "/tmp/hadded_epsilon_"+str(idx)+".root"
        outputfile.write("echo Copying files locally\n")
        strippedFiles = []
        for f in filesToMerge:
            f = f.strip()
            if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
                outputfile.write(myeosstage + "root://eoscms//eos/cms" + f + " /tmp \n")
            elif "cmsStage -f" in myeosstage:
                outputfile.write(myeosstage + f + " /tmp \n")        
            else:
                outputfile.write("cmsStage -f " + f + " /tmp \n")
            strippedFiles.append(ntpath.basename(f))
        outputfile.write("filesHadd=\"/tmp/" + " /tmp/".join(strippedFiles) + "\"\n")
        outputfile.write("echo \"hadd -k " + mergedfile_n + " $filesHadd\"\n")
        outputfile.write("hadd -k " + mergedfile_n + " $filesHadd\n")
        outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n\n")

        grouped_files.append(mergedfile_n)
        idx += 1
        files = files[grouping:]

    outputfile.write("echo now hadding the intermediate hadded files: " + " ".join(grouped_files) + "\n")
    outputfile.write("hadd -k /tmp/" +  NameTag + "epsilonPlots.root " + " ".join(grouped_files) + "\n")

    if "/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select" in myeosstage:
        outputfile.write("echo \"eos cp /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    else:
        outputfile.write("echo \"cmsStage -f /tmp/" + NameTag + "epsilonPlots.root " + destination + "\"\n")
    outputfile.write(myeosstage + "/tmp/" + NameTag + "epsilonPlots.root " + destination + "\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "epsilonPlots*\n")
    outputfile.write("rm -rf /tmp/" + NameTag + "hadded_epsilon*\n")

