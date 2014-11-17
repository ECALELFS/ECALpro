#Do not modify these
nEventsPerJob = '-1'
outputFile    = 'EcalNtp' # without .root suffix
calibMapName = 'calibMap.root'
GeometryFromFile = False
ExternalGeometry = 'caloGeometry.root'
CalibType  = 'xtal'

#Are Pi0
Are_pi0  = False # True = using Pi0, False = using Eta
#IS CRAB
isCRAB = False
CRAB_Data_Path = '/Neutrino_Pt-2to20_gun/Fall13dr-tsg_PU40bx25_POSTLS162_V2-v1/AODSIM'
CRAB_CopyCert  = '/afs/cern.ch/user/l/lpernie/private/x509up_u12147'
CRAB_Storage   = 'group/alca_ecalcalib/lpernie/' #the absence of the beginning slash is mandatory
events_per_job = '5000'
total_number_of_events = '-1'
#MC and TTree
isMC = True
MakeNtuple4optimization = True
#PATH
eosPath = '/store/caf/user/lpernie'
#eosPath = '/store/group/alca_ecalcalib/lpernie/'
if(isCRAB):
   eosPath = '/store/group/alca_ecalcalib/lpernie/'
inputlist_n      = 'ALL_MINBIAS_PIZERO_ALCARAW_NO_UNCAL.list' # list of the input files
dirname          = 'ALL_MINBIAS_ETA_ALCARAW_NO_UNCAL_01_TryL1Geom'

Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output
#TAG, QUEUE and ITERS
NameTag          = 'MC_'                # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'          # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'
nIterations = 1
#N files
ijobmax          = 5                     # 5 number of files per job
nHadd            = 35                    # 50 number of files per hadd
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'
#Remove Xtral Dead
RemoveDead_Flag = "False"
RemoveDead_Map = "/afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_6_2_5/src/CalibCode/submit/AfterCalibTools/DeadXtals/plots/h_DeadXtal.root"
#L1 Bit Collection
L1TriggerInfo = False;                          # If we want to Fill the L1 Trigger Bit Histo (and if we perform the cut based on a L1Bit of L1Seed != "")
hltGtDigis = 'InputTag("simGtDigis")'     # To take the info to Fill the L1 Bit histo
triggerTag = 'InputTag("TriggerResults")'
hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap")'
L1Seed = ""                                     # You can ask taht one Bit is FIRED: Ex: "L1_SingleJet16" or more complicated stuff "L1_SingleJet16 OR L1_SingleJet36"

#Seeds (Comment if you want the standard cuts ones)
EB_Seed_E    = '0.5'
useEE_EtSeed = 'False'
EE_Seed_Et   = '0.5'
EE_Seed_E    = '1.5'
#Selection
CutOnHLTIso = "False"
if(Are_pi0):
   #inner barrel
   Pi0PtCutEB_low = '2.1'
   gPtCutEB_low = '0.8'
   Pi0IsoCutEB_low = '0.'
   nXtal_1_EB_low = '6'
   nXtal_2_EB_low = '4'
   S4S9_EB_low = '0.7'
   #outer barrel
   Pi0PtCutEB_high = '2.1'
   gPtCutEB_high = '0.8'
   Pi0IsoCutEB_high = '0.8'
   nXtal_1_EB_high = '6'
   nXtal_2_EB_high = '4'
   S4S9_EB_high = '0.7'
   #low eta EE
   Pi0PtCutEE_low = '2.1'
   gPtCutEE_low = '0.8'
   Pi0IsoCutEE_low = '0.25'
   nXtal_1_EE_low = '6'
   nXtal_2_EE_low = '4'
   S4S9_EE_low = '0.85'   
   #high eta EE
   Pi0PtCutEE_high = '2.1'
   gPtCutEE_high = '0.8'
   Pi0IsoCutEE_high = '0.25'
   nXtal_1_EE_high = '6'
   nXtal_2_EE_high = '4'
   S4S9_EE_high = '0.85'
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
#ETA
else:
   #inner barrel
   Pi0PtCutEB_low = '3.2'
   gPtCutEB_low = '1.4'
   Pi0IsoCutEB_low = '0.'
   nXtal_1_EB_low = '6'
   nXtal_2_EB_low = '4'
   S4S9_EB_low = '0.8'
   #outer barrel
   Pi0PtCutEB_high = '3.2'
   gPtCutEB_high = '1.4'
   Pi0IsoCutEB_high = '0.8'
   nXtal_1_EB_high = '6'
   nXtal_2_EB_high = '4'
   S4S9_EB_high = '0.8'
   #low eta EE
   Pi0PtCutEE_low = '3.2'
   gPtCutEE_low = '1.4'
   Pi0IsoCutEE_low = '0.25'
   nXtal_1_EE_low = '6'
   nXtal_2_EE_low = '4'
   S4S9_EE_low = '0.85'   
   #high eta EE
   Pi0PtCutEE_high = '3.2'
   gPtCutEE_high = '1.4'
   Pi0IsoCutEE_high = '0.25'
   nXtal_1_EE_high = '6'
   nXtal_2_EE_high = '4'
   S4S9_EE_high = '0.85'
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
useEEContainmentCorrections = 'False'
EBContainmentCorrections = 'totNewPi0TupleMB_fillingTot.fittedcorrectionsEB.root'
MVAEBContainmentCorrections_01 = 'JOSH_MVA_pi01_Mediumtrain.root'
MVAEBContainmentCorrections_02 = 'JOSH_MVA_pi02_Mediumtrain.root'
MVAEEContainmentCorrections_01 = 'JOSH_MVA_pi01_Mediumtrain_EE.root'
MVAEEContainmentCorrections_02 = 'JOSH_MVA_pi02_Mediumtrain_EE.root'
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
##2012D
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012D.txt'
#isNewTag=True
#HLTResults = 'True'
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
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
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
###(_2): https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
#isNewTag=True
#json_file ='goodrunlist_json2012C.txt'
#HLTResults = 'True'
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#if not(isNewTag):
#   globaltag='GR_P_V42::All'
#else:
####   globaltag='FT_R_53_V21::All'
#   globaltag='FT_53_V21_AN6::All'
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
###else:
###   alphaTagRecord='EcalLaserAPDPNRatiosRcd'
###   alphaTag='EcalLaserAPDPNRatios_20130130_447_p1_v2'
###   alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
##############
##2012B
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012C.txt'
#isNewTag=True
#HLTResults = 'True'
#is_2011 = 'True' #Just for the fit, put true 
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
#is_2011 = 'True' #Just for the fit, put true 
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
#is_2011 = 'True'
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
#HLTPaths='AlCa_EcalPi0_*'

##############

##2010 AlcaRECO
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/goodrunlist_json.txt
#json_file = 'goodrunlist_json2010.txt'
#HLTResults = 'False'
#is_2011 = 'False'
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
#HLTPaths='AlCa_EcalPi0_*' 

##2010 AlcaRECO
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/goodrunlist_json.txt
#HLTResults = 'False'
#json_file = ''
#is_2011 = 'False'
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
#HLTPaths='AlCa_EcalPi0_*'
#isMC = True

##MC 2010 AlcaRECO
#HLTResults = 'False'
#json_file = ''
#is_2011 = 'False'
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
#HLTPaths='AlCa_EcalPi0_*'
#isMC = True

##MC MINBIAS_PIZERO_ALCARAW_NOL1_v2
#HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
#json_file = ''
#is_2011 = 'False'                                             # Fit Parameter Range
#overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
#doEnenerScale='False'
#doIC='False'                                                  # Member of Recalibration Module
#doLaserCorr="True"                                            # Member of Recalibration Module
#ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','TEST')"
#eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEB','TEST')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES','TEST')"
#useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
#correctHits = 'False'
#globaltag='MCRUN2_72_V1A::All'
#HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter

## MC 40bx25 HLT ALCARAW
#HLTResults = 'False'
#json_file = ''
#is_2011 = 'False'
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
#HLTPaths='AlCa_EcalPi0_*'

##MC MINBIAS_PIZERO_ALCARAW_NOL1_v2
HLTResults = 'False'                                          # Use the function GetHLTResults(iEvent, "AlCa_EcalPi0EBonly.*");
json_file = ''
is_2011 = 'False'                                             # Fit Parameter Range
overWriteGlobalTag = False                                    # Allow to overwrite AlphaTag, Laser correction etc
doEnenerScale='False'
doIC='False'                                                  # Member of Recalibration Module
doLaserCorr="False"
if(Are_pi0):                                           # Member of Recalibration Module
   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','TEST')"
   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','TEST')"
   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES','TEST')"
else:
   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','TEST')"
   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','TEST')"
   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES','TEST')"
hltGtDigis = "InputTag('simGtDigis','','TEST')"
triggerTag = 'InputTag("TriggerResults","","TEST")'
hltL1GtObjectMap = 'InputTag("hltL1GtObjectMap","","TEST")'
useHLTFilter = "False"                                        # Add to the path the request of a HLT path:  process.AlcaP0Filter.HLTPaths = 
correctHits = 'False'
globaltag='MCRUN2_72_V1A::All'
HLTPaths='AlCa_EcalPi0_*'                                     # Name of the HLT path selected with useHLTFilter
isMC = True

##2012 Selection
      #2012
      #Pi0PtCutEB_low = '1.'
      #Pi0PtCutEB_high = '1.'
      #Pi0PtCutEE_low = '1.'
      #Pi0PtCutEE_high = '1.'
      #gPtCutEB_low = '0.4'
      #gPtCutEB_high = '0.4'
      #gPtCutEE_low = '0.4'
      #gPtCutEE_high = '0.4'
      #Pi0IsoCutEB_low = '0.'
      #Pi0IsoCutEB_high = '0.'
      #Pi0IsoCutEE_low = '0.'
      #Pi0IsoCutEE_high = '0.'
      #nXtal_1_EB_low = '0.'
      #nXtal_1_EB_high = '0.'
      #nXtal_2_EB_low = '0.'
      #nXtal_2_EB_high = '0.'
      #nXtal_1_EE_low = '0.'
      #nXtal_1_EE_high = '0.'
      #nXtal_2_EE_low = '0.'
      #nXtal_2_EE_high = '0.'
      #S4S9_EB_low = '0.6'
      #S4S9_EB_high = '0.6'
      #S4S9_EE_low = '0.6'
      #S4S9_EE_high = '0.6'
