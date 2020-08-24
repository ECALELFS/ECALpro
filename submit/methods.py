from parameters import *
import os

def printFillCfg1( outputfile ):
    outputfile.write("useHLTFilter = " + useHLTFilter + "\n")
    outputfile.write("correctHits = " + correctHits + "\n\n")
    outputfile.write('import FWCore.ParameterSet.Config as cms\n')
    outputfile.write('import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi\n')
    outputfile.write("import os, sys, imp, re\n")
    outputfile.write('CMSSW_VERSION=os.getenv("CMSSW_VERSION")\n')
    if runCalibrationFromRecHits:
        outputfile.write('process = cms.Process("analyzerFillEpsilonFromRecHits")\n')
    else:
        outputfile.write('process = cms.Process("analyzerFillEpsilon")\n')
    outputfile.write('process.load("FWCore.MessageService.MessageLogger_cfi")\n\n')
    outputfile.write('process.load("Configuration.Geometry.GeometryIdeal_cff")\n')

    # if (nThread > 1 ):
    #     outputfile.write("\n")
    #     outputfile.write("process.options = cms.untracked.PSet(\n")
    #     outputfile.write("    numberOfThreads = cms.untracked.uint32( 4 ),\n")
    #     outputfile.write("    numberOfStreams = cms.untracked.uint32( 0 ),\n")
    #     outputfile.write("    sizeOfStackForThreadsInKB = cms.untracked.uint32( 10*1024 )\n")
    #     outputfile.write(")\n\n")


    outputfile.write('process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")\n')
    #outputfile.write('process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")\n')
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
        if not(PFRechitTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + PFRechitTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + PFRechitTag + "'),\n")
            outputfile.write("              connect = cms.string('" + PFRechitDB + "')\n")
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
        if not(ESIntercalibTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + ESIntercalibTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + ESIntercalibTag + "'),\n")
            outputfile.write("              connect = cms.string('" + ESIntercalibDB + "')\n")
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
        if not(EcalChannelStatusTag==''):
            outputfile.write("     cms.PSet(record = cms.string('" + EcalChannelStatusTagRecord + "'),\n")
            outputfile.write("              tag = cms.string('" + EcalChannelStatusTag + "'),\n")
            outputfile.write("              connect = cms.string('" + EcalChannelStatusDB + "')\n")
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
    outputfile.write("    ),\n")
    outputfile.write("    skipBadFiles = cms.untracked.bool(True)\n")
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
    if(not useJsonFilterInCpp and len(json_file)>0):
        outputfile.write("import FWCore.PythonUtilities.LumiList as LumiList\n")
        outputfile.write("json_file = '" + json_file + "'\n")
        if json_file.startswith('/afs/cern.ch/'):
            outputfile.write("process.source.lumisToProcess = LumiList.LumiList(filename = json_file).getVLuminosityBlockRange()\n")
        else:
            CMSSW_BASE=os.getenv("CMSSW_BASE")
            outputfile.write("process.source.lumisToProcess = LumiList.LumiList(filename = '" + CMSSW_BASE + "/src/CalibCode/FillEpsilonPlot/data/" + json_file + "').getVLuminosityBlockRange()\n")


    outputfile.write("\n")
    if justCreateRecHits:
        # here we create rechits through output module, without running fillEpsilon analyzer
        outputfile.write("## create rechits and save them\n")
        outputfile.write("outRecHitsFileName = '%s/AlCaP0_RecHitsFromDigis_%s/out_%s.root'\n" % (eosOutputPathForRecHits, dirname, str(ijob)))
        outputfile.write("process.outputALCAP0 = cms.OutputModule( 'PoolOutputModule',\n")
        outputfile.write("    fileName = cms.untracked.string( str(outRecHitsFileName) ),\n")
        outputfile.write("    fastCloning = cms.untracked.bool( False ),\n")
        outputfile.write("    dataset = cms.untracked.PSet(\n")
        outputfile.write("        filterName = cms.untracked.string( '' ),\n")
        outputfile.write("        dataTier = cms.untracked.string( 'RAW' )\n")
        outputfile.write("    ),\n")
        if filterEventsByAlCaTrigger:
            outputfile.write("    SelectEvents = cms.untracked.PSet(\n")
            outputfile.write("        SelectEvents = cms.vstring('p')\n")
            outputfile.write("    ),\n")
        outputfile.write("    outputCommands = cms.untracked.vstring( \n")
        outputfile.write("        'drop *',\n")
        if Barrel_or_Endcap == 'ONLY_BARREL':
            outputfile.write("        'keep *_ecalRecHit_EcalRecHitsEB*_*',\n")
        elif Barrel_or_Endcap == 'ONLY_ENDCAP':
            outputfile.write("        'keep *_ecalRecHit_EcalRecHitsEE*_*',\n")
        else:
            outputfile.write("        'keep *_ecalRecHit_EcalRecHitsE*_*',\n")
                         
        if not Barrel_or_Endcap == 'ONLY_BARREL':
            if not filterEventsByAlCaTrigger or (filterEventsByAlCaTrigger and 'AlCa_EcalEta' in HLTPaths):
                outputfile.write("        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*',\n")
                outputfile.write("        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*',\n")
            if not filterEventsByAlCaTrigger or (filterEventsByAlCaTrigger and 'AlCa_EcalPi0' in HLTPaths):
                outputfile.write("        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*',\n")
                outputfile.write("        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*',\n")
        if L1TriggerInfo:
            outputfile.write("        'keep *_hltGtStage2Digis_*_*',\n")
        outputfile.write("        'keep edmTriggerResults_*_*_HLT' ) # collection is recreated for this process, keep only original from HLT\n")
        outputfile.write(")\n")
    else:
        outputfile.write("process.analyzerFillEpsilon = cms.EDAnalyzer('FillEpsilonPlot')\n")
        outputfile.write("process.analyzerFillEpsilon.OutputDir = cms.untracked.string('" +  outputDir + "')\n")
        outputfile.write("process.analyzerFillEpsilon.OutputFile = cms.untracked.string('" + NameTag +  outputFile + "_" + str(ijob) + ".root')\n")
        outputfile.write("process.analyzerFillEpsilon.ExternalGeometry = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + ExternalGeometry + "')\n")
        if (isCRAB):
            outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + NameTag + calibMapName + "')\n")
            outputfile.write("process.analyzerFillEpsilon.isCRAB  = cms.untracked.bool(True)\n")
        else:
            if (SubmitFurtherIterationsFromExisting and iteration == 0):
                outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('" + startingCalibMap + "')\n")
            else:
                outputfile.write("process.analyzerFillEpsilon.calibMapPath = cms.untracked.string('" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")
            if SubmitFurtherIterationsFromExisting:
                if SystOrNot != 0:
                    outputfile.write("process.analyzerFillEpsilon.SystOrNot = cms.untracked.int32(" + str(SystOrNot) + ")\n")

        outputfile.write("process.analyzerFillEpsilon.Endc_x_y                        = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + Endc_x_y + "')\n")
        outputfile.write("process.analyzerFillEpsilon.HLTResults                  = cms.untracked.bool(" + HLTResults + ")\n")
        if(HLTResultsNameEB!=""):
            outputfile.write("process.analyzerFillEpsilon.HLTResultsNameEB            = cms.untracked.string('" + HLTResultsNameEB + "')\n")
        if(HLTResultsNameEE!=""):
            outputfile.write("process.analyzerFillEpsilon.HLTResultsNameEE            = cms.untracked.string('" + HLTResultsNameEE + "')\n")
        if RemoveSeedsCloseToDeadXtal:
            outputfile.write("process.analyzerFillEpsilon.RemoveSeedsCloseToDeadXtal             = cms.untracked.bool(True)\n")        
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
        if(Are_pi0):
            outputfile.write("process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(True)\n")
        else:
            outputfile.write("process.analyzerFillEpsilon.Are_pi0                 = cms.untracked.bool(False)\n")


        if useContainmentCorrectionsFromEoverEtrue:
            outputfile.write("process.analyzerFillEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( True )\n")
            outputfile.write("process.analyzerFillEpsilon.scalingEoverEtrueCC_g1 = cms.untracked.double("+ scalingEoverEtrueCC_g1 +")\n")
            outputfile.write("process.analyzerFillEpsilon.scalingEoverEtrueCC_g2 = cms.untracked.double("+ scalingEoverEtrueCC_g2 +")\n")
            if copyCCfileToTMP:
                copiedCCfile = str(fileEoverEtrueContainmentCorrections.split('/')[-1])
                copiedCCfile = "/tmp/" + copiedCCfile.replace(".root","_iter_{ni}_job_{nj}.root".format(ni=iteration, nj=ijob))
                # do not copy here, but inside job .sh
                #os.system("xrdcp {rf} /tmp/{rfcopy}".format(rf=fileEoverEtrueContainmentCorrections, rfcopy=copiedCCfile)
                outputfile.write("process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string(\"" + copiedCCfile + "\")\n")
            else:
                outputfile.write("process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string(\"" + fileEoverEtrueContainmentCorrections + "\")\n")
        else:
            outputfile.write("process.analyzerFillEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( False )\n")
            outputfile.write("process.analyzerFillEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string(\"\")\n")


        outputfile.write("process.analyzerFillEpsilon.useOnlyEEClusterMatchedWithES = cms.untracked.bool(" + useOnlyEEClusterMatchedWithES + ")\n\n")

        outputfile.write("### choosing proper input tag (recalibration module changes the collection names)\n")
        outputfile.write("if correctHits:\n")
        outputfile.write("    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')\n")
        outputfile.write("    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')\n")
        outputfile.write("else:\n")
        outputfile.write("    process.analyzerFillEpsilon.EBRecHitCollectionTag = cms.untracked." + ebInputTag + "\n")
        outputfile.write("    process.analyzerFillEpsilon.EERecHitCollectionTag = cms.untracked." + eeInputTag + "\n")
        outputfile.write("process.analyzerFillEpsilon.ESRecHitCollectionTag = cms.untracked." + esInputTag + "\n")

        outputfile.write("process.analyzerFillEpsilon.triggerTag   = cms.untracked." + triggerTag + "\n")
        outputfile.write("process.analyzerFillEpsilon.L1GTobjmapTag   = cms.untracked." + L1GTobjmapTag + "\n")
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
        outputfile.write("process.analyzerFillEpsilon.nXtal_1_EB_low = cms.untracked.int32(" +  nXtal_1_EB_low+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_1_EB_high = cms.untracked.int32(" +  nXtal_1_EB_high+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_2_EB_low = cms.untracked.int32(" +  nXtal_2_EB_low+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_2_EB_high = cms.untracked.int32(" +  nXtal_2_EB_high+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_1_EE_low = cms.untracked.int32(" +  nXtal_1_EE_low+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_1_EE_high = cms.untracked.int32(" +  nXtal_1_EE_high+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_2_EE_low = cms.untracked.int32(" +  nXtal_2_EE_low+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.nXtal_2_EE_high = cms.untracked.int32(" +  nXtal_2_EE_high+ ")\n")
        outputfile.write("process.analyzerFillEpsilon.S4S9_EB_low = cms.untracked.double(" + S4S9_EB_low + ")\n")
        outputfile.write("process.analyzerFillEpsilon.S4S9_EB_high = cms.untracked.double(" + S4S9_EB_high + ")\n")
        outputfile.write("process.analyzerFillEpsilon.S4S9_EE_low = cms.untracked.double(" + S4S9_EE_low + ")\n")
        outputfile.write("process.analyzerFillEpsilon.S4S9_EE_high = cms.untracked.double(" + S4S9_EE_high + ")\n")
        outputfile.write("process.analyzerFillEpsilon.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
        if(useJsonFilterInCpp and len(json_file)>0):
           if json_file.startswith('/afs/cern.ch/'): 
               outputfile.write("process.analyzerFillEpsilon.JSONfile = cms.untracked.string('" + json_file + "')\n")
           else:
               outputfile.write("process.analyzerFillEpsilon.JSONfile = cms.untracked.string('CalibCode/FillEpsilonPlot/data/" + json_file + "')\n")
        if GeometryFromFile:
           outputfile.write("process.analyzerFillEpsilon.GeometryFromFile = cms.untracked.bool(True)\n")
        if useMassInsteadOfEpsilon:
            outputfile.write("process.analyzerFillEpsilon.useMassInsteadOfEpsilon = cms.untracked.bool(True)\n")
        else:
            outputfile.write("process.analyzerFillEpsilon.useMassInsteadOfEpsilon = cms.untracked.bool(False)\n")
        if isDebug:
            outputfile.write("process.analyzerFillEpsilon.isDebug = cms.untracked.bool(True)\n")
        if isMC:
           outputfile.write("process.analyzerFillEpsilon.isMC = cms.untracked.bool(True)\n")
           outputfile.write("process.analyzerFillEpsilon.pileupSummaryTag = cms.untracked." + pileupInputTag + "\n")
           if(MC_Assoc):
               outputfile.write("process.analyzerFillEpsilon.GenPartCollectionTag = cms.untracked." + genPartInputTag + "\n")
               outputfile.write("process.analyzerFillEpsilon.MC_Assoc             = cms.untracked.bool(True)\n")
               outputfile.write("process.analyzerFillEpsilon.MC_Assoc_DeltaR      = cms.untracked.double(" + MC_Assoc_DeltaR + ")\n")        
           if isEoverEtrue:
               outputfile.write("process.analyzerFillEpsilon.isEoverEtrue = cms.untracked.bool(True)\n")
        if fillKinematicVariables:
            outputfile.write("process.analyzerFillEpsilon.fillKinematicVariables = cms.untracked.bool(True)\n")
        else:
            outputfile.write("process.analyzerFillEpsilon.fillKinematicVariables = cms.untracked.bool(False)\n")
        if MakeNtuple4optimization:
           outputfile.write("process.analyzerFillEpsilon.MakeNtuple4optimization = cms.untracked.bool(True)\n")
        if( L1TriggerInfo ):
            outputfile.write("process.analyzerFillEpsilon.L1TriggerInfo = cms.untracked.bool(True)\n")
            outputfile.write("process.analyzerFillEpsilon.L1SeedsPi0Stream = cms.untracked.string(\"" + L1SeedExpression + "\")\n")
            nSeeds = L1SeedExpression.count(" OR ") + 1 
            outputfile.write("process.analyzerFillEpsilon.nL1SeedsPi0Stream = cms.untracked.int32(" + str(nSeeds) + ")\n")
        else:
            # if L1TriggerInfo is false, pass following two variables anyway, because FIllEpsilonPlot.cc expects to get them, even though they won't be used
            outputfile.write("process.analyzerFillEpsilon.L1TriggerInfo = cms.untracked.bool(False)\n")
            outputfile.write("process.analyzerFillEpsilon.L1SeedsPi0Stream = cms.untracked.string(\"\")\n")
            outputfile.write("process.analyzerFillEpsilon.nL1SeedsPi0Stream = cms.untracked.int32(0)\n")

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
    if not justCreateRecHits:
        outputfile.write("process.p *= process.analyzerFillEpsilon\n")
    if justCreateRecHits:
        outputfile.write("process.endp = cms.EndPath(process.outputALCAP0)\n")
    else:
        outputfile.write("process.endp = cms.EndPath()\n")

def printFitCfg( outputfile, iteration, outputDir, nIn, nFin, EBorEE, nFit, justDoHistogramFolding=False ):
    if isEoverEtrue and localFolderToWriteFits:
        outputDir = outputDir.replace("/tmp",localFolderToWriteFits)
    outputfile.write("import FWCore.ParameterSet.Config as cms\n")
    outputfile.write("process = cms.Process('FitEpsilonPlot')\n")
    outputfile.write("process.load('FWCore.MessageService.MessageLogger_cfi')\n")
    outputfile.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )\n")
    outputfile.write("process.source =   cms.Source('EmptySource')\n")
    outputfile.write("process.fitEpsilon = cms.EDAnalyzer('FitEpsilonPlot')\n")
    outputfile.write("process.fitEpsilon.OutputFile = cms.untracked.string('" + NameTag + EBorEE + "_" + str(nFit) + "_" + calibMapName + "')\n")
    outputfile.write("process.fitEpsilon.CalibType = cms.untracked.string('" + CalibType + "')\n")
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
    

    if useContainmentCorrectionsFromEoverEtrue:
        outputfile.write("process.fitEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( True )\n")
        outputfile.write("process.fitEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string(\"" + fileEoverEtrueContainmentCorrections + "\")\n")
    else:
        outputfile.write("process.fitEpsilon.useContainmentCorrectionsFromEoverEtrue = cms.untracked.bool( False )\n")
        outputfile.write("process.fitEpsilon.fileEoverEtrueContainmentCorrections = cms.untracked.string(\"\")\n")

    if useMassInsteadOfEpsilon:
        outputfile.write("process.fitEpsilon.useMassInsteadOfEpsilon = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.fitEpsilon.useMassInsteadOfEpsilon = cms.untracked.bool(False)\n")
    if isEoverEtrue:
        outputfile.write("process.fitEpsilon.isEoverEtrue = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.fitEpsilon.isEoverEtrue = cms.untracked.bool(False)\n")
    outputfile.write("process.fitEpsilon.StoreForTest = cms.untracked.bool( True )\n")
    if foldInSuperModule:
        outputfile.write("process.fitEpsilon.foldInSuperModule = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.fitEpsilon.foldInSuperModule = cms.untracked.bool(False)\n")
    if useFit_RooMinuit:
        outputfile.write("process.fitEpsilon.useFit_RooMinuit = cms.untracked.bool( True )\n")        
    outputfile.write("process.fitEpsilon.Barrel_orEndcap = cms.untracked.string('" + Barrel_or_Endcap + "')\n")
    if not(isCRAB): #If CRAB you have to put the correct path, and you do it on calibJobHandler.py, not on ./submitCalibration.py
        outputfile.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('" + eosPath + "/" + dirname + "/iter_" + str(iteration) + "/" + NameTag + "epsilonPlots.root')\n")
        if (SubmitFurtherIterationsFromExisting and iteration == 0):
            outputfile.write("process.fitEpsilon.calibMapPath = cms.untracked.string('" + startingCalibMap + "')\n")
        else:
            outputfile.write("process.fitEpsilon.calibMapPath = cms.untracked.string('" + eosPath + "/" + dirname + "/iter_" + str(iteration-1) + "/" + NameTag + calibMapName + "')\n")
    if justDoHistogramFolding:
        outputfile.write("process.fitEpsilon.makeFoldedHistograms = cms.untracked.bool(True)\n")
    else:
        outputfile.write("process.fitEpsilon.makeFoldedHistograms = cms.untracked.bool(False)\n")

    outputfile.write("process.p = cms.EndPath(process.fitEpsilon)\n")


def printSubmitFitSrc(outputfile, cfgName, source, destination, pwd, logpath, justDoHistogramFolding=False):
    if isEoverEtrue and localFolderToWriteFits:        
        source = source.replace("/tmp",localFolderToWriteFits)
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    # if using RooMinuit fo the fit (obsolete according too RooFit guide), then save some pieces from stdout to check status of fit
    # this is not needed if one uses RooMinimizer (suggested option) because the printing is different
    # anyway, these prints doesn't affect the code behaviour, they just fall in the fit log file
    # we keep only the output containing 'FIT_EPSILON:'
    if useFit_RooMinuit:
        outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/{quote}' > " + logpath  + "\n")
        outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FIT_EPSILON:/ || /WITHOUT CONVERGENCE/ || /HAS CONVERGED/' >> " + logpath  + "\n")
    else:
        outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FIT_EPSILON:/{quote}' > " + logpath  + "\n")
        outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FIT_EPSILON:/' >> " + logpath  + "\n")
    # if only folding there won't be any actual output in /tmp
    if not justDoHistogramFolding:
        sourcerooplot = source.replace("calibMap","fitRes")
        destrooplot = destination.replace("calibMap","fitRes")
        outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
        outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
        outputfile.write("echo 'ls " + sourcerooplot + " >> " + logpath + " 2>&1' \n" )
        outputfile.write("ls " + sourcerooplot + " >> " + logpath + " 2>&1 \n" )
        outputfile.write("echo 'cp " + source + " " + destination + "' >> " + logpath  + "\n")           
        outputfile.write("echo 'cp " + sourcerooplot + " " + destrooplot + "' >> " + logpath  + "\n")
        outputfile.write("cp " + source + " " + destination + " >> " + logpath + " 2>&1 \n")
        outputfile.write("cp " + sourcerooplot + " " + destrooplot + " >> " + logpath + " 2>&1 \n")
        outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
        outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")
        outputfile.write("echo 'rm -f " + sourcerooplot + "' >> " + logpath + " \n")
        outputfile.write("rm -f " + sourcerooplot + " >> " + logpath + " 2>&1 \n")

def printSubmitSrc(outputfile, cfgName, source, destination, pwd, logpath):
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    # outputfile.write("source /cvmfs/cms.cern.ch/crab3/crab.sh\n") this line produces problem when running in CMSSW_8_0_3, anyway we don't use crab
    copiedCCfile = ""
    if useContainmentCorrectionsFromEoverEtrue and copyCCfileToTMP:
        # get iter and job numbers
        #cfgName is like /bla/bla/fillEps_iter_N_job_M.py, need to get iter_N_job_M
        iterJob = str(cfgName.split('/')[-1])
        iterJob = iterJob.replace("fillEps_","").replace(".py","")
        copiedCCfile = str(fileEoverEtrueContainmentCorrections.split('/')[-1])
        copiedCCfile = "/tmp/" + copiedCCfile.replace(".root","_{itj}.root".format(itj=iterJob))
        # copy file to tmp (recreate if already exists using -f)
        cpcmd = "xrdcp -f {rf} {rfcopy}".format(rf=fileEoverEtrueContainmentCorrections, rfcopy=copiedCCfile)
        #cpcmd = "xrdcp {rf} {rfcopy}".format(rf=fileEoverEtrueContainmentCorrections,rfcopy=copiedCCfile) 
        outputfile.write("echo '" + cpcmd + "'\n")
        outputfile.write(cpcmd + "\n")
    if not Silent:
        outputfile.write("echo 'cmsRun " + cfgName + "'\n")
        outputfile.write("cmsRun " + cfgName + "\n")

        outputfile.write("if test -f " + source + "; then\n")
        outputfile.write("    echo 'file exists in %s and is good, now copying to eos'\n" % source)
        outputfile.write("    echo 'cp " + source + " " + destination + "'\n")
        outputfile.write("    cp " + source + " " + destination + "\n")
        outputfile.write("    echo 'rm -f " + source + "'\n")
        outputfile.write("    rm -f " + source + "\n")
        outputfile.write("    echo 'now checking goodness of file on eos'\n")
        outputfile.write("    cmdpy='python " + pwd + "/Utilities/checkGoodnessFileEOS.py -d " + destination + "'\n")
        outputfile.write("    echo \"${cmdpy}\"\n") # note: need " and not '  here !!
        outputfile.write("    echo \"${cmdpy}\" | bash\n") # here as well
        outputfile.write("    echo ''\n")
        outputfile.write("else\n")
        outputfile.write("    echo 'file did not exist in /tmp/ (probably it was bad and deleted already'\n")
        outputfile.write("fi\n")
    else:
        outputfile.write("echo 'cmsRun " + cfgName + " 2>&1 | awk {quote}/FILL_COUT:/{quote}' > " + logpath  + "\n")
        outputfile.write("cmsRun " + cfgName + " 2>&1 | awk '/FILL_COUT:/' >> " + logpath  + "\n")
        outputfile.write("echo 'ls " + source + " >> " + logpath + " 2>&1' \n" )
        outputfile.write("ls " + source + " >> " + logpath + " 2>&1 \n" )
        outputfile.write("echo 'cp " + source + " " + destination + "' >> " + logpath  + "\n")
        outputfile.write("cp " + source + " " + destination + " >> " + logpath + " 2>&1 \n")
        outputfile.write("echo 'rm -f " + source + "' >> " + logpath + " \n")
        outputfile.write("rm -f " + source + " >> " + logpath + " 2>&1 \n")
    if len(copiedCCfile):
        outputfile.write("echo 'rm -f " + copiedCCfile + "'\n")
        outputfile.write("rm -f " + copiedCCfile + "\n")        

def printParallelHaddFAST(outputfile, outFile, listReduced, destination, pwd, numList):
    import os, sys, imp, re
    destinationWithFinalSlash = destination 
    if not destinationWithFinalSlash.endswith("/"):
        destinationWithFinalSlash = destination + "/"
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("haddstr=''\n")
    outputfile.write("for file in ${files};\n")
    outputfile.write("do\n")
    outputfile.write("   haddstr=\"${haddstr} ${file}\"\n")
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -f -k " + destinationWithFinalSlash + NameTag + "epsilonPlots_" + str(numList) + ".root ${haddstr}\"\n")
    outputfile.write("hadd -f -k " + destinationWithFinalSlash + NameTag + "epsilonPlots_" + str(numList) + ".root ${haddstr}\n")

def printFinalHaddFAST(outputfile, listReduced, destination, pwd):
    import os, sys, imp, re
    destinationWithFinalSlash = destination 
    if not destinationWithFinalSlash.endswith("/"):
        destinationWithFinalSlash = destination + "/"
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    outputfile.write("files=`cat " + listReduced + "`\n")
    outputfile.write("haddstr=''\n")
    outputfile.write("for file in ${files};\n")
    outputfile.write("do\n")
    outputfile.write('   haddstr="${haddstr} ${file}"\n')
    outputfile.write("done\n")
    outputfile.write("echo \"hadd -f -k " + destinationWithFinalSlash + NameTag + "epsilonPlots.root ${haddstr}\"\n")
    outputfile.write("hadd -f -k " + destinationWithFinalSlash + NameTag + "epsilonPlots.root ${haddstr}\n")

def printFinalHaddRegroup(outputfile, listReduced, destination, pwd, grouping=10):
    import os, sys, imp, re, ntpath
    destinationWithFinalSlash = destination 
    if not destinationWithFinalSlash.endswith("/"):
        destinationWithFinalSlash = destination + "/"
    CMSSW_VERSION=os.getenv("CMSSW_VERSION")
    outputfile.write("#!/bin/bash\n")
    outputfile.write("cd " + pwd + "\n")
    outputfile.write("eval `scramv1 runtime -sh`\n")
    fileWithList = open(listReduced,"r")
    files = fileWithList.readlines()
    idx=0
    grouped_files = []
    while len(files)>0:
        filesToMerge = files[:grouping]
        mergedfile_n = ("%s" + "hadded_epsilon_"+str(idx)+".root") % destinationWithFinalSlash
        strippedFiles = []
        for f in filesToMerge:
            f = f.strip()
            strippedFiles.append(ntpath.basename(f))
        outputfile.write(("filesHadd=\"{eos}" + " {eos}".join(strippedFiles) + "\"\n").format(eos=destinationWithFinalSlash))
        outputfile.write("echo \"hadd -f -k " + mergedfile_n + " ${filesHadd}\"\n")
        outputfile.write("hadd -f -k " + mergedfile_n + " ${filesHadd}\n")

        grouped_files.append(mergedfile_n)
        idx += 1
        files = files[grouping:]

    outputfile.write("echo \"now hadding the intermediate hadded files: " + " ".join(grouped_files) + "\"\n")
    outputfile.write("hadd -f -k " + destinationWithFinalSlash +  NameTag + "epsilonPlots.root " + " ".join(grouped_files) + "\n")
    outputfile.write("rm " + destinationWithFinalSlash + "hadded_epsilon*\n")

