#Do not modify these
nEventsPerJob      = '-1'
outputFile         = 'EcalNtp'           # without .root suffix
calibMapName       = 'calibMap.root'
GeometryFromFile   = False               # Keep that False, you want the cmssw geometry. Anyway the geometry file is needed
ExternalGeometry   = 'caloGeometry.root' 
CalibType          = 'xtal'              # Calibrating single xtals. I never try but you could calibrate EtaRing ot Trigger Towers

#Are Pi0
Are_pi0            = True               # True = using Pi0, False = using Eta
#Fold per Eta Ring
EtaRingCalibEB     = False
SMCalibEB          = False
EtaRingCalibEE     = False
SMCalibEE          = False
CalibMapEtaRing    = "CalibCode/FillEpsilonPlot/data/calibMap.root"
#PATH
eosPath = '/store/caf/user/zhicaiz'
#eosPath = '/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian'
#
#adding following variables to use commands like "eos ls" and "eos ls -l" commands instead of cmsLs.
#See also here for more details --> https://twiki.cern.ch/twiki/bin/view/CMSPublic/CERNStorageTools 
#   
myeoscmd = '/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select '  #this call directly the eos command (note that eos is an alias, see link above)
myeosls = myeoscmd + 'ls '  #to avoid use of cmsLs that is deprecated since January 2016   
myeoslsl = myeosls + '-l '
myeosmkdir = myeoscmd + 'mkdir '
#myeosstage = myeoscmd + 'cp '  
myeosstage = 'cmsStage -f '
# I called it myeosstage instead of myeoscp to remember that it substitutes cmsStage command
# as a convention, when adding commands like: command = myeoscmd + "some_option ", just leave a space AFTER the some_option, not before
# note that code used cmsStage -f, but eos cp doesn't support -f option
# also, code will copy *.root files from /tmp/ (where they are initially created) to eosPath, but eosPath must be preceeded by "root://eoscms/eos/cms" to have eos cp
# work as expected. So the destination will be root://eoscms/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/... . For this reason, we define here
#myPrefixToEosPath = 'root://eoscms//eos/cms'
myPrefixToEosPath = ''
# will modify calibJobHandler.py with this prefix to destination
#
# end of my additions
#  
#CRAB
isCRAB           = False               # If not is batch
CRAB_Data_Path   = '/SinglePion_FlatPt-1To15_AsymptNoPU/emanuele-SinglePion_FlatPt-1To15_AsymptNoPU-9709e5e865f17288f5a53621cf8e9935/USER'
CRAB_CopyCert    = '/afs/cern.ch/user/l/lpernie/private/x509up_u12147'
storageSite      = "T2_CH_CERN"
unitsPerJob = 10   #DBS File per Job
isOtherT2        = False
if(isCRAB):
   eosPath = '/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/' #For reason of space is better the group area
   if(isOtherT2):
       eosPath = '/pnfs/roma1.infn.it/data/cms/store/user/mciprian/piZero2016/'
       voGroup     = "itcms"
       storageSite = "T2_IT_Rome"
       outLFN      = "/store/user/mciprian/piZero2016/"
#MC and Selection Optimization
isMC = False
MakeNtuple4optimization = False
#InputList and Folder name
inputlist_n      = 'InputList/run272818.list' # list of input files
dirname          = 'run272818_YongCC'
Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output
#TAG, QUEUE and ITERS
NameTag          = 'testScript'                   # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'          # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'
nIterations      = 1
#N files
ijobmax          = 5                     # 5 number of files per job
nHadd            = 35                    # 35 number of files per hadd
fastHadd         = True                 # From 7_4_X we can use this faster mathod. But files have to be copied on /tmp/ to be converted in .db
if( isCRAB and isOtherT2 ):
   fastHadd      = False                 # No fastHadd on a different T2
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'
#Remove Xtral Dead
RemoveDead_Flag = "True"
RemoveDead_Map  = ""
#RemoveDead_Map  = "/afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_6_2_5/src/CalibCode/submit/AfterCalibTools/DeadXtals/plots/h_DeadXtal.root"

#L1 Bit Collection
L1TriggerInfo = False;                              # If we want to Fill the L1 Trigger Bit Histo (and if we perform the cut based on a L1Bit of L1Seed != "")
hltGtDigis = 'InputTag("simGtDigis")'               # Not used anymore in the Fill.cc -> To take the info to Fill the L1 Bit histo
triggerTag = 'InputTag("TriggerResults")'           # To run the FillEB only if the HLTName for EB is present
hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap")'   # To fill the L1 Trigger fired
L1Seed = ""                                         # You can ask that one Bit is FIRED: Ex: "L1_SingleJet16" or more complicated stuff "L1_SingleJet16 OR L1_SingleJet36"
#Association with GenPart
MC_Asssoc = False

#Seeds (Comment if you want the standard cuts ones)
EB_Seed_E    = '0.5'
useEE_EtSeed = 'False'
EE_Seed_Et   = '0.0'
EE_Seed_E    = '1.0' #1.5 for 40PU25
#Selection
CutOnHLTIso = "False"
if(Are_pi0):
   #inner barrel
   Pi0PtCutEB_low = '2.0'
   gPtCutEB_low = '1.0'
   Pi0IsoCutEB_low = '0.0'
   Pi0HLTIsoCutEB_low = "999"
   nXtal_1_EB_low = '5'
   nXtal_2_EB_low = '5'
   S4S9_EB_low = '0.9'
   #outer barrel # P. Jarry cuts
   Pi0PtCutEB_high = '3.0'
   gPtCutEB_high = '2.0'
   Pi0IsoCutEB_high = '0.0'
   Pi0HLTIsoCutEB_high = '999'
   nXtal_1_EB_high = '5'
   nXtal_2_EB_high = '5'
   S4S9_EB_high = '0.9'
   #low eta EE
   Pi0PtCutEE_low = '3.0'
   gPtCutEE_low = '0.95'
   Pi0IsoCutEE_low = '.0'
   Pi0HLTIsoCutEE_low = '999'
   nXtal_1_EE_low = '5'
   nXtal_2_EE_low = '5'
   S4S9_EE_low = '0.95'
   #high eta EE
   Pi0PtCutEE_high = '1.5'
   gPtCutEE_high = '0.65'
   Pi0IsoCutEE_high = '0.0'
   Pi0HLTIsoCutEE_high = '999'
   nXtal_1_EE_high = '5'
   nXtal_2_EE_high = '5'
   S4S9_EE_high = '0.95'
   if MakeNtuple4optimization:
      #inner barrel
      Pi0PtCutEB_low = '1'
      gPtCutEB_low = '0.4'
      Pi0IsoCutEB_low = '0.0'
      Pi0HLTIsoCutEB_low = "999"
      nXtal_1_EB_low = '0'
      nXtal_2_EB_low = '0'
      S4S9_EB_low = '0.6'
      #outer barrel
      Pi0PtCutEB_high = '1.0'
      gPtCutEB_high = '0.4'
      Pi0IsoCutEB_high = '0.0'
      Pi0HLTIsoCutEB_high = '999'
      nXtal_1_EB_high = '0'
      nXtal_2_EB_high = '0'
      S4S9_EB_high = '0.6'
      #low eta EE
      Pi0PtCutEE_low = '1.0'
      gPtCutEE_low = '0.4'
      Pi0IsoCutEE_low = '.0'
      Pi0HLTIsoCutEE_low = '999'
      nXtal_1_EE_low = '0'
      nXtal_2_EE_low = '0'
      S4S9_EE_low = '0.6'
      #high eta EE
      Pi0PtCutEE_high = '1.0'
      gPtCutEE_high = '0.4'
      Pi0IsoCutEE_high = '0.0'
      Pi0HLTIsoCutEE_high = '999'
      nXtal_1_EE_high = '0'
      nXtal_2_EE_high = '0'
      S4S9_EE_high = '0.6'
#ETA
else:
   #inner barrel
   Pi0PtCutEB_low = '1'
   gPtCutEB_low = '.4'
   Pi0IsoCutEB_low = '0.0'
   Pi0HLTIsoCutEB_low = "999"
   nXtal_1_EB_low = '0'
   nXtal_2_EB_low = '0'
   S4S9_EB_low = '0.6'
   #outer barrel
   Pi0PtCutEB_high = '1.0'
   gPtCutEB_high = '.4'
   Pi0IsoCutEB_high = '0.0'
   Pi0HLTIsoCutEB_high = '999'
   nXtal_1_EB_high = '0'
   nXtal_2_EB_high = '0'
   S4S9_EB_high = '0.6'
   #low eta EE
   Pi0PtCutEE_low = '1.0'
   gPtCutEE_low = '.4'
   Pi0IsoCutEE_low = '.0'
   Pi0HLTIsoCutEE_low = '999'
   nXtal_1_EE_low = '0'
   nXtal_2_EE_low = '0'
   S4S9_EE_low = '0.6'
   #high eta EE
   Pi0PtCutEE_high = '1.0'
   gPtCutEE_high = '0.4'
   Pi0IsoCutEE_high = '0.0'
   Pi0HLTIsoCutEE_high = '999'
   nXtal_1_EE_high = '0'
   nXtal_2_EE_high = '0'
   S4S9_EE_high = '0.6'
   if MakeNtuple4optimization:
      #inner barrel
      Pi0PtCutEB_low = '1'
      gPtCutEB_low = '.4'
      Pi0IsoCutEB_low = '0.0'
      Pi0HLTIsoCutEB_low = "999"
      nXtal_1_EB_low = '0'
      nXtal_2_EB_low = '0'
      S4S9_EB_low = '0.6'
      #outer barrel
      Pi0PtCutEB_high = '1.0'
      gPtCutEB_high = '.4'
      Pi0IsoCutEB_high = '0.0'
      Pi0HLTIsoCutEB_high = '999'
      nXtal_1_EB_high = '0'
      nXtal_2_EB_high = '0'
      S4S9_EB_high = '0.6'
      #low eta EE
      Pi0PtCutEE_low = '1.0'
      gPtCutEE_low = '.4'
      Pi0IsoCutEE_low = '.0'
      Pi0HLTIsoCutEE_low = '999'
      nXtal_1_EE_low = '0'
      nXtal_2_EE_low = '0'
      S4S9_EE_low = '0.6'
      #high eta EE
      Pi0PtCutEE_high = '1.0'
      gPtCutEE_high = '0.4'
      Pi0IsoCutEE_high = '0.0'
      Pi0HLTIsoCutEE_high = '999'
      nXtal_1_EE_high = '0'
      nXtal_2_EE_high = '0'
      S4S9_EE_high = '0.6'
#containment corrections
useEBContainmentCorrections = 'True'
useEEContainmentCorrections = 'True'
EBContainmentCorrections = 'totNewPi0TupleMB_fillingTot.fittedcorrectionsEB.root'
MVAEBContainmentCorrections_01 = 'JOSH_MVA_pi01_Mediumtrain.root'
MVAEBContainmentCorrections_02 = 'JOSH_MVA_pi02_Mediumtrain.root'
MVAEEContainmentCorrections_01 = 'JOSH_MVA_pi01_Mediumtrain_EE.root'
MVAEEContainmentCorrections_02 = 'JOSH_MVA_pi02_Mediumtrain_EE.root'

useMVAContainmentCorrections = False
new_pi0ContainmentCorrections = False
new_MVAEBContainmentCorrections_01 = 'new_JOSH_MVA_pi01_Mediumtrain.root'
new_MVAEBContainmentCorrections_02 = 'new_JOSH_MVA_pi02_Mediumtrain.root'
new_MVAEEContainmentCorrections_01 = 'new_JOSH_MVA_pi01_Mediumtrain_EE.root'
new_MVAEEContainmentCorrections_02 = 'new_JOSH_MVA_pi02_Mediumtrain_EE.root'

MVAEBContainmentCorrections_eta01 = 'JOSH_MVA_eta1_Mediumtrain.root'
MVAEBContainmentCorrections_eta02 = 'JOSH_MVA_eta2_Mediumtrain.root'
Endc_x_y = 'Endc_x_y_ring.txt'
EBPHIContainmentCorrections = 'correctionsEB_PHI.root'
EEContainmentCorrections = 'totNewPi0TupleMB_fillingTot.fittedcorrectionsEE.root'
EBContCorr = 'correctionsEB.root'

# preshower
useOnlyEEClusterMatchedWithES = 'True'

#-----------------------------------------------------------------------------------
laserTagRecord='';laserTag='';laserDB=''
alphaTagRecord2='';alphaTag2='';alphaDB2=''
GeVTagRecord='';GeVTag='';GeVDB=''

######################################################################
# Now decomment the part that correspond to data you want to run on. #
######################################################################

##2015C AlCaP0 RAW
isMC               = False
isNot_2010         = 'True'                                    # Fit Parameter Range
HLTResults         = 'True'                                    # Fill the EB(EE) histos only is Eb()ee is fired: it uses GetHLTResults(iEvent, HLTResultsNameEB.Data() );
json_file          = 'json_DCSONLY.txt' if isMC==False else ''            #/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/json_ecalonly/
overWriteGlobalTag = False                                     # Allow to overwrite AlphaTag, Laser correction etc
doEnenerScale      = 'False'
doIC               = 'False'                                   # Member of Recalibration Module
doLaserCorr        = "False"
hltGtDigis         = "InputTag('simGtDigis')"        # Not used in the Fill.cc   
triggerTag         = 'InputTag("TriggerResults")'    # Run Fill EB only if the HLTPaths for EB(ee) exist. In this sample also extist InputTag('simGtDigis','','HLT')
hltL1GtObjectMap   = 'InputTag("hltL1GtObjectMap")'
useHLTFilter       = "True" if isMC==False else "False"                                  # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
correctHits        = 'False'
globaltag          = '80X_dataRun2_Prompt_v8' if isMC==False else '80X_mcRun2_asymptotic_v5' #old is GR_P_V56
globaltag_New      = True
FROMDIGI           = True
DigiCustomization  = False   # keep this False since CMSSW_7_4_15, there is a module in CMSSW providing the bunchSpacing.  ===> NEW - 03/05/2016 - : can set it True because to run (at least) on data, that introduces --> outputfile.write("process.ecalMultiFitUncalibRecHit.algoPSet.useLumiInfoRunHeader = False\n") <-- in fillEpsilonPlot*.py file, which is needed to run without errors, but it also add another line to activate process.ecalMultiFitUncalibRecHit.algoPSet.activeBXs, so keep False for now
MULTIFIT           = True;   # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
is50ns             = False      # If DigiCustomization and MULTIFIT is True
WEIGHTS            = False;   # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
if(Are_pi0):                                           # Member of Recalibration Module
   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES')"
   HLTPaths='AlCa_EcalPi0E*'                        # HLT Name to ask before running the event. It can contain a *.
   HLTResultsNameEB   = 'AlCa_EcalPi0EB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
   HLTResultsNameEE   = 'AlCa_EcalPi0EE'
else:
   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonlyRegional','etaEcalRecHitsES')"
   HLTPaths='AlCa_EcalEtaE*' #AlCa_EcalEtaEBonly_LowPU_v1 AlCa_EcalEtaEEonly_LowPU_v1
   HLTResultsNameEB   = 'AlCa_EcalEtaEB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
   HLTResultsNameEE   = 'AlCa_EcalEtaEE'
if(FROMDIGI):
   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
   if(Are_pi0): 
      EBdigi = 'InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis")'
      EEdigi = 'InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis")'
   else:
      EBdigi = 'InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis")'
      EEdigi = 'InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis")'
else:
   if isMC:
      ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","RECO")'
      eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","RECO")'
      esInputTag = 'InputTag("ecalPreshowerRecHit","EcalRecHitsES","RECO")'
   else:
      if(Are_pi0):
         ebInputTag = 'InputTag("hltAlCaPi0EBUncalibrator","pi0EcalRecHitsEB")'
         eeInputTag = 'InputTag("hltAlCaPi0EEUncalibrator","pi0EcalRecHitsEE")'
      else:
         ebInputTag = 'InputTag("hltAlCaEtaEBUncalibrator","etaEcalRecHitsEB")'
         eeInputTag = 'InputTag("hltAlCaEtaEEUncalibrator","etaEcalRecHitsEE")'
if isMC:
   MC_Asssoc = True
   genPartInputTag = 'InputTag("genParticles","")'
      
#<<<<<<< HEAD

##2015B AlCaP0 RAW
#isMC               = False
#isNot_2010         = 'True'                                    # Fit Parameter Range
#HLTResults         = 'True'                                    # Fill the EB(EE) histos only is Eb()ee is fired: it uses GetHLTResults(iEvent, HLTResultsNameEB.Data() );
#json_file          = 'goodrunlist_json2015Bred.txt'            #/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/json_ecalonly/
#overWriteGlobalTag = False                                     # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale      = 'False'
#doIC               = 'False'                                   # Member of Recalibration Module
#doLaserCorr        = "False"
#hltGtDigis         = "InputTag('simGtDigis','','HLT')"        # Not used in the Fill.cc   
#triggerTag         = 'InputTag("TriggerResults","","HLT")'    # Run Fill EB only if the HLTPaths for EB(ee) exist. In this sample also extist InputTag('simGtDigis','','HLT')
#hltL1GtObjectMap   = 'InputTag("hltL1GtObjectMap","","HLT")'
#useHLTFilter       = "True"                                   # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits        = 'False'
#globaltag          = '74X_dataRun2_Prompt_v0' #old is GR_P_V56
#globaltag_New      = True
#FROMDIGI           = True
#DigiCustomization  = False
#MULTIFIT           = True;   # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
#is50ns             = True      # If DigiCustomization and MULTIFIT is True
#WEIGHTS            = False;   # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
#if(Are_pi0):                                           # Member of Recalibration Module
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES','HLT')"
#   HLTPaths='AlCa_EcalPi0E*'                        # HLT Name to ask before running the event. It can contain a *.
#   HLTResultsNameEB   = 'AlCa_EcalPi0EB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalPi0EE'
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonlyRegional','etaEcalRecHitsES','HLT')"
#   HLTPaths='AlCa_EcalEtaE*' #AlCa_EcalEtaEBonly_LowPU_v1 AlCa_EcalEtaEEonly_LowPU_v1
#   HLTResultsNameEB   = 'AlCa_EcalEtaEB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalEtaEE'
#if(FROMDIGI):
#   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
#   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
#   if(Are_pi0): 
#      EBdigi = 'InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","HLT")'
#      EEdigi = 'InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","HLT")'
#   else:
#      EBdigi = 'InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis","HLT")'
#      EEdigi = 'InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis","HLT")'
#else:
#   if(Are_pi0):
#      ebInputTag = 'InputTag("hltAlCaPi0EBUncalibrator","pi0EcalRecHitsEB","HLT")'
#      eeInputTag = 'InputTag("hltAlCaPi0EEUncalibrator","pi0EcalRecHitsEE","HLT")'
#   else:
#      ebInputTag = 'InputTag("hltAlCaEtaEBUncalibrator","etaEcalRecHitsEB","HLT")'
#      eeInputTag = 'InputTag("hltAlCaEtaEEUncalibrator","etaEcalRecHitsEE","HLT")'

###2015A AlCaP0 RAW
#isMC               = False
#isNot_2010         = 'True'                                    # Fit Parameter Range
#HLTResults         = 'True'                                    # Fill the EB(EE) histos only is Eb()ee is fired: it uses GetHLTResults(iEvent, HLTResultsNameEB.Data() );
#json_file          = ''                                        #/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/json_ecalonly/
#overWriteGlobalTag = False                                     # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale      = 'False'
#doIC               = 'False'                                   # Member of Recalibration Module
#doLaserCorr        = "False"
#hltGtDigis         = "InputTag('simGtDigis','','HLT')"         # Not used in the Fill.cc   
#triggerTag         = 'InputTag("TriggerResults","","HLT")'     # Run Fill EB only if the HLTPaths for EB(ee) exist. In this sample also extist InputTag('simGtDigis','','HLT')
#hltL1GtObjectMap   = 'InputTag("hltL1GtObjectMap","","HLT")'
#useHLTFilter       = "True"                                    # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits        = 'False'
#globaltag          = 'GR_P_V55::All'
#globaltag_New      = False
#FROMDIGI           = True
#DigiCustomization  = False
#MULTIFIT           = True;     # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
#is50ns             = True      # If DigiCustomization and MULTIFIT is True
#WEIGHTS            = False;    # Choose WEIGHTS or MULTIFIT (MULTIFIT is standard)
#if(Are_pi0):                                           # Member of Recalibration Module
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES','HLT')"
#   HLTPaths='AlCa_EcalPi0E*'                        # HLT Name to ask before running the event. It can contain a *.
#   HLTResultsNameEB   = 'AlCa_EcalPi0EB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalPi0EE'
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonlyRegional','etaEcalRecHitsES','HLT')"
#   HLTPaths='AlCa_EcalEtaE*' #AlCa_EcalEtaEBonly_LowPU_v1 AlCa_EcalEtaEEonly_LowPU_v1
#   HLTResultsNameEB   = 'AlCa_EcalEtaEB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalEtaEE'
#if(FROMDIGI):
#   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
#   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
#   if(Are_pi0): 
#      EBdigi = 'InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","HLT")'
#      EEdigi = 'InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","HLT")'
#   else:
#      EBdigi = 'InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis","HLT")'
#      EEdigi = 'InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis","HLT")'
#else:
#   if(Are_pi0):
#      ebInputTag = 'InputTag("hltAlCaPi0EBUncalibrator","pi0EcalRecHitsEB","HLT")'
#      eeInputTag = 'InputTag("hltAlCaPi0EEUncalibrator","pi0EcalRecHitsEE","HLT")'
#   else:
#      ebInputTag = 'InputTag("hltAlCaEtaEBUncalibrator","etaEcalRecHitsEB","HLT")'
#      eeInputTag = 'InputTag("hltAlCaEtaEEUncalibrator","etaEcalRecHitsEE","HLT")'

###2015A Commissioning MinBias AOD
#isMC               = False
#isNot_2010         = 'True'                                    # Fit Parameter Range
#HLTResults         = 'False'                                    # Fill the EB(EE) histos only is Eb()ee is fired: it uses GetHLTResults(iEvent, HLTResultsNameEB.Data() );
#json_file          = ''
#overWriteGlobalTag = False                                     # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale      = 'False'
#doIC               = 'False'                                   # Member of Recalibration Module
#doLaserCorr        = "False"
#hltGtDigis         = "InputTag('simGtDigis','','TEST')"        # Not used in the Fill.cc   
#triggerTag         = 'InputTag("TriggerResults","","RECO")'    # Run Fill EB only if the HLTPaths for EB(ee) exist. In this sample also extist InputTag('simGtDigis','','HLT')
#hltL1GtObjectMap   = 'InputTag("hltL1GtObjectMap","","TEST")'
#useHLTFilter       = "False"                                    # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits        = 'False'
#globaltag_New      = True
#globaltag          = 'GR_H_V58C' #GR_P_V55:All GR_H_V58C
#FROMDIGI           = False
#DigiCustomization  = False
#is50ns = True #If DigiCustomization is True
#if(Are_pi0):                                           # Member of Recalibration Module
#   esInputTag = "InputTag('reducedEcalRecHitsES','','RECO')"
#   HLTPaths='AlCa_EcalPi0E*'                        # HLT Name to ask before running the event. It can contain a *.
#   HLTResultsNameEB   = 'AlCa_EcalPi0EB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalPi0EE'
#else:
#   esInputTag = "InputTag('reducedEcalRecHitsES','','RECO')"
#   HLTPaths='AlCa_EcalEtaE*' #AlCa_EcalEtaEBonly_LowPU_v1 AlCa_EcalEtaEEonly_LowPU_v1
#   HLTResultsNameEB   = 'AlCa_EcalEtaEB'            # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalEtaEE'
#if(FROMDIGI):
#   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
#   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
#   if(Are_pi0): 
#      EBdigi = 'InputTag("selectDigi","selectedEcalEBDigiCollection","RECO")'
#      EEdigi = 'InputTag("selectDigi","selectedEcalEBDigiCollection","RECO")'
#   else:
#      EBdigi = 'InputTag("selectDigi","selectedEcalEBDigiCollection","RECO")'
#      EEdigi = 'InputTag("selectDigi","selectedEcalEBDigiCollection","RECO")'
#else:
#   ebInputTag = 'InputTag("reducedEcalRecHitsEB","","RECO")'
#   eeInputTag = 'InputTag("reducedEcalRecHitsEE","","RECO")'

###2015 Commissioning MinBias Josh
#isMC               = False
#isNot_2010         = 'True'                                   # Fit Parameter Range
#HLTResults         = 'True'                                   # Fill the EB(EE) histos only is Eb()ee is fired: it uses GetHLTResults(iEvent, "AlCa_EcalPi0E*");
#json_file          = ''
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale      = 'False'
#doIC               = 'False'                                  # Member of Recalibration Module
#doLaserCorr        = "False"
#hltGtDigis         = "InputTag('simGtDigis','','TEST')"
#triggerTag         = 'InputTag("TriggerResults","","TEST")'
#hltL1GtObjectMap   = 'InputTag("hltL1GtObjectMap","","TEST")'
#useHLTFilter       = "True"                                   # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits        = 'False'
#globaltag          = 'GR_P_V54::All' #GR_P_V55 soon!
#FROMDIGI           = False
#DigiCustomization  = False
#is50ns = True #If DigiCustomization is True
#if(Are_pi0):                                           # Member of Recalibration Module
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES','TEST')"
#   HLTPaths='AlCa_EcalPi0E*' #AlCa_EcalPi0EBonly_LowPU_v1 AlCa_EcalPi0EEonly_LowPU_v1
#   HLTResultsNameEB   = 'AlCa_EcalPi0EB'                         # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalPi0EE'                                  
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonlyRegional','etaEcalRecHitsES','TEST')"
#   HLTPaths='AlCa_EcalEtaE*' #AlCa_EcalEtaEBonly_LowPU_v1 AlCa_EcalEtaEEonly_LowPU_v1
#   HLTResultsNameEB   = 'AlCa_EcalEtaEB'                         # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#   HLTResultsNameEE   = 'AlCa_EcalEtaEE'                                  
#if(FROMDIGI):
#   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
#   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
#   if(Are_pi0): 
#      EBdigi = 'InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","TEST")'
#      EEdigi = 'InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","TEST")'
#   else:
#      EBdigi = 'InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis","TEST")'
#      EEdigi = 'InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis","TEST")'
#else:
#   if(Are_pi0):
#      ebInputTag = 'InputTag("hltAlCaPi0EBUncalibrator","pi0EcalRecHitsEB","TEST")'
#      eeInputTag = 'InputTag("hltAlCaPi0EEUncalibrator","pi0EcalRecHitsEE","TEST")'
#   else:
#      ebInputTag = 'InputTag("hltAlCaEtaEBUncalibrator","etaEcalRecHitsEB","TEST")'
#      eeInputTag = 'InputTag("hltAlCaEtaEEUncalibrator","etaEcalRecHitsEE","TEST")'

##2012D
####/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012D.txt'
#isNewTag=True
#HLTResults = 'True'
#isNot_2010 = 'True' #Just for the fit, put true 
#useHLTFilter="False" #Should be True
#correctHits='False' #Should be True
#overWriteGlobalTag = False #Should be True
#if not(isNewTag):
#   globaltag='GR_P_V40::All'
#else:
#   globaltag='FT_R_53_V21::All'
#if(Are_pi0): 
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalPi0*' 
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalEta*'
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#if not(isNewTag):
#   GeVTagRecord='EcalADCToGeVConstantRcd'
#   GeVTag='EcalADCToGeVConstant_Bon_V20111129'
#   GeVDB='frontier://FrontierProd/CMS_COND_31X_ECAL'
#   laserTagRecord='EcalIntercalibConstantsRcd'
#   laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#   laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#   alphaTagRecord2='EcalLaserAlphasRcd'
#   alphaTag2='EcalLaserAlphas_EB_sic1_btcp152_EE_sic1_btcp116'
#   alphaDB2='frontier://FrontierInt/CMS_COND_ECAL'
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20130124_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#else:
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20130130_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'


##############
#2012C
##/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
##(_2): https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
#isNewTag=True
#json_file ='goodrunlist_json2012C.txt'
#HLTResults = 'True'
#isNot_2010 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='False' #SHOULD be true
#overWriteGlobalTag = True
#if not(isNewTag):
#   globaltag='GR_P_V42::All'
#else:
#   globaltag='GR_R_70_V2::All'
#if(Are_pi0): 
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalPi0*' 
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalEta*'
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#if not(isNewTag):
#   laserTagRecord='EcalIntercalibConstantsRcd'
#   laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#   laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20121020_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#else:
#   alphaTagRecord=''
#   alphaTag=''
#   alphaDB=''
#else:
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20130130_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
##############
##2012B
##/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012C.txt'
#isNewTag=True
#HLTResults = 'True'
#isNot_2010 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#if not(isNewTag):
#   globaltag='FT_R_53_V6::All'
#else:
#   globaltag='FT_R_53_V21::All'
#if(Are_pi0): 
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalPi0*' 
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalEta*'
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#if not(isNewTag):
#   laserTagRecord='EcalIntercalibConstantsRcd'
#   laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#   laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20121020_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#else:
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20130130_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'

##############
#2012A
##/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012C.txt'
#isNewTag=True
#HLTResults = 'True'
#isNot_2010 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#if not(isNewTag):
#   globaltag='FT_R_53_V6::All'
#else:
#   globaltag='FT_R_53_V21::All'
#if(Are_pi0): 
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalPi0*' 
#else:
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES', 'HLT')"
#   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','HLT')"
#   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','HLT')"
#   HLTPaths='AlCa_EcalEta*'
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#if not(isNewTag):
#   laserTagRecord='EcalIntercalibConstantsRcd'
#   laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#   laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20121020_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#else:
#   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#   alphaTag='EcalLaserAPDPNRatios_20130130_447_p1_v2'
#   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'

##############
##2011 AlcaRAW
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt
#json_file ='goodrunlist_json2011.txt'
#HLTResults = 'False'
#isNot_2010 = 'True'
#overWriteGlobalTag = False
#doEnenerScale='False'
#doIC='False'
#doLaserCorr='True'
#ebInputTag = "InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEB', 'HLT')"
#eeInputTag = "InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEE', 'HLT')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES', 'HLT')"
#useHLTFilter = 'True'
#correctHits = 'True'
#globaltag='GR_R_42_V24::All'
#laserTagRecord='EcalLaserAPDPNRatiosRcd'
#alphaTagRecord='EcalLaserAlphasRcd'
#laserTag            = 'EcalLaserAPDPNRatios_data_20120814_2011-2012_v3'
#laserDB             = 'frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#alphaTag            = 'EcalLaserAlphas_EB_sic_btcp152_EE_sic1_btcp116'
#alphaDB             = 'frontier://FrontierPrep/CMS_COND_ECAL'
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter

##############

##2010 AlcaRECO
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/goodrunlist_json.txt
#json_file = 'goodrunlist_json2010.txt'
#HLTResults = 'False'
#isNot_2010 = 'False'
#overWriteGlobalTag = False
#doEnenerScale='False'
#doIC='False'
#doLaserCorr="True"
#ebInputTag = "InputTag('ecalPi0Corrected','pi0EcalRecHitsEB')"
#eeInputTag = "InputTag('ecalPi0Corrected','pi0EcalRecHitsEE')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES')"
#useHLTFilter = "False"
#correctHits = 'False'
#globaltag='GR_R_42_V21B::All'
#laserTagRecord='EcalLaserAPDPNRatiosRcd'
#laserTag = 'EcalLaserAPDPNRatios_data_20111122_158851_180363'
#laserDB  = 'frontier://FrontierPrep/CMS_COND_ECAL'
#alphaTagRecord='EcalLaserAlphasRcd'
#alphaTag = 'EcalLaserAlphas_lto420-620_progr_data_20111122'
#alphaDB  = 'frontier://FrontierPrep/CMS_COND_ECAL'
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter

##2010 AlcaRECO
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/goodrunlist_json.txt
#HLTResults = 'False'
#json_file = ''
#isNot_2010 = 'False'
#overWriteGlobalTag = False
#doEnenerScale='False'
#doIC='False'
#doLaserCorr="True"
#ebInputTag = "InputTag('ecalRecHit','EcalRecHitsEB','RECO')"
#eeInputTag = "InputTag('ecalRecHit','EcalRecHitsEE','RECO')"
#esInputTag = "InputTag('ecalPreshowerRecHit','EcalRecHitsES')"
#useHLTFilter = "False"
#correctHits = 'False'
#globaltag='GR_R_42_V21B::All'
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter
#isMC = True

##MC 2010 AlcaRECO
#HLTResults = 'False'
#json_file = ''
#isNot_2010 = 'False'
#overWriteGlobalTag = False
#doEnenerScale='False'
#doIC='False'
#doLaserCorr="True"
#ebInputTag = "InputTag('ecalRecHit','EcalRecHitsEB','RECO')"
#eeInputTag = "InputTag('ecalRecHit','EcalRecHitsEE','RECO')"
#esInputTag = "InputTag('ecalPreshowerRecHit','EcalRecHitsES')"
#useHLTFilter = "False"
#correctHits = 'False'
#globaltag='GR_R_42_V21B::All'
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter
#isMC = True

##MC MINBIAS_PIZERO_ALCARAW_NOL1_v2
#HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
#json_file = ''
#isNot_2010 = 'False'                                             # Fit Parameter Range
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale='False'
#doIC='False'                                                  # Member of Recalibration Module
#doLaserCorr="True"                                            # Member of Recalibration Module
#ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','TEST')"
#eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEB','TEST')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES','TEST')"
#useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits = 'False'
#globaltag='MCRUN2_74_V6A::All'
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter

## MC 40bx25 HLT ALCARAW
#HLTResults = 'False'
#json_file = ''
#isNot_2010 = 'False'
#overWriteGlobalTag = False
#doEnenerScale='False'
#doIC='False'
#doLaserCorr="False"
#ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','TEST')"
#eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEB','TEST')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES','TEST')"
#useHLTFilter = "False"
#correctHits = 'False'
#globaltag='POSTLS162_V2::All'
#if(Are_pi0):
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter

##MC MINBIAS_PIZERO_ALCARAW_NOL1_v2
#HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
#json_file = ''
#isNot_2010 = 'False'                                             # Fit Parameter Range
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale='False'
#doIC='False'                                                  # Member of Recalibration Module
#doLaserCorr="False"
#if(Are_pi0):                                           # Member of Recalibration Module
#   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB')"
#   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE')"
#   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonlyRegional','pi0EcalRecHitsES','TEST')"
#else:
#   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','TEST')"
#   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','TEST')"
#   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonlyRegional','etaEcalRecHitsES','TEST')"
#hltGtDigis = "InputTag('simGtDigis','','TEST')"
#triggerTag = 'InputTag("TriggerResults","","TEST")'
#hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap","","TEST")'
#useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits = 'False'
#globaltag='MCRUN2_74_V6A::All'
#if(Are_pi0):
#   HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#else:
#   HLTPaths='AlCa_EcalEta_*'                                     # Name of the HLT path selected with useHLTFilter
#isMC = True
#FROMDIGI = True
#is50ns = False
#if(FROMDIGI):
#   ebInputTag = 'InputTag("ecalRecHit","EcalRecHitsEB","analyzerFillEpsilon")'
#   eeInputTag = 'InputTag("ecalRecHit","EcalRecHitsEE","analyzerFillEpsilon")'
#   if(Are_pi0): 
#      EBdigi = 'InputTag("hltAlCaPi0EBRechitsToDigis","pi0EBDigis","TEST")'
#      EEdigi = 'InputTag("hltAlCaPi0EERechitsToDigis","pi0EEDigis","TEST")'
#   else:
#      EBdigi = 'InputTag("hltAlCaEtaEBRechitsToDigis","etaEBDigis","TEST")'
#      EEdigi = 'InputTag("hltAlCaEtaEERechitsToDigis","etaEEDigis","TEST")'

##MC CRAB Neutrino GUN
#HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, HLTResultsNameEB or EE);
#HLTResultsNameEB   = 'AlCa_EcalPi0EB'                         # HLT Name to ask for into the GetHLTResults (do not use name_EB* please)
#HLTResultsNameEE   = 'AlCa_EcalPi0EE'
#json_file = ''
#isNot_2010 = 'False'                                             # Fit Parameter Range
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale='False'
#doIC='False'                                                  # Member of Recalibration Module
#doLaserCorr="False"
#FROMDIGI=False
#ebInputTag = "InputTag('reducedEcalRecHitsEB','')"
#eeInputTag = "InputTag('reducedEcalRecHitsEE','')"
#esInputTag = "InputTag('reducedEcalRecHitsES','')"
#hltGtDigis = "InputTag('simGtDigis','','TEST')"
#triggerTag = 'InputTag("TriggerResults","","TEST")'
#hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap","","TEST")'
#useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits = 'False'
#globaltag='MCRUN2_74_V6A::All'
#HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#isMC = True

#Pi0Gun
#HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
#json_file = ''
#isNot_2010 = 'False'                                             # Fit Parameter Range
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale='False'
#doIC='False'                                                  # Member of Recalibration Module
#doLaserCorr="False"
#ebInputTag = "InputTag('ecalRecHit','EcalRecHitsEB')"
#eeInputTag = "InputTag('ecalRecHit','EcalRecHitsEE')"
#esInputTag = "InputTag('ecalPreshowerRecHit','EcalRecHitsES')"
#hltGtDigis = "InputTag('simGtDigis','','TEST')"
#triggerTag = 'InputTag("TriggerResults","","TEST")'
#hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap","","TEST")'
#useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits = 'False'
#globaltag='MCRUN2_74_V6A::All'
#HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
#isMC = True
#MC_Asssoc = True
#genPartInputTag = "InputTag('genParticles','')"


