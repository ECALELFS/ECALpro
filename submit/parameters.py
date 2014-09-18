#Do not modify these
nEventsPerJob = '-1'
outputFile    = 'EcalNtp' # without .root suffix
calibMapName = 'calibMap.root'
ExternalGeometry = 'caloGeometry.root'
CalibType  = 'xtal'

#Are Pi0
Are_pi0  = True # True = using Pi0, False = using Eta
#IS CRAB
isCRAB = True
CRAB_Data_Path = '/Neutrino_Pt-2to20_gun/Fall13dr-tsg_PU40bx25_POSTLS162_V2-v1/AODSIM'
events_per_job = '5000'
total_number_of_events = '-1'
#MC and TTree
isMC = True
MakeNtuple4optimization = False
#PATH
eosPath = '/store/caf/user/lpernie'
#eosPath = '/store/group/alca_ecalcalib/lpernie/'
if(isCRAB):
   eosPath = '/store/group/alca_ecalcalib/lpernie/'
inputlist_n      = 'ALL_NeuPt2_20_PU40x25_reduced.list' # list of the input files
dirname          = 'ALL_NeuPt2_20_PU40x25_02'
Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output
#TAG, QUEUE and ITERS
NameTag          = 'MC13_'              # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'                 # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'
nIterations = 13
#N files
ijobmax          = 5                     # 5 number of files per job
nHadd            = 35                    # 50 number of files per hadd
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'
#Remove Xtral Dead
RemoveDead_Flag = "False"
RemoveDead_Map = "/afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_6_2_5/src/CalibCode/submit/AfterCalibTools/DeadXtals/plots/h_DeadXtal.root"
#L1Seeds
L1Seed = "" #You can ask "L1_SingleJet16" or more complicated stuff "L1_SingleJet16 OR L1_SingleJet36"

#Selection
if(Are_pi0):
   Pi0PtCutEB_low = '2.1'
   Pi0PtCutEB_high = '2.1'
   Pi0PtCutEE_low = '2.1'
   Pi0PtCutEE_high = '2.1'
   gPtCutEB_low = '0.8'
   gPtCutEB_high = '0.8'
   gPtCutEE_low = '0.8'
   gPtCutEE_high = '0.8'
   Pi0IsoCutEB_low = '0.'
   Pi0IsoCutEB_high = '0.'
   Pi0IsoCutEE_low = '0.25'
   Pi0IsoCutEE_high = '0.25'
   nXtal_1_EB_low = '6'
   nXtal_1_EB_high = '6'
   nXtal_2_EB_low = '4'
   nXtal_2_EB_high = '4'
   nXtal_1_EE_low = '6'
   nXtal_1_EE_high = '6'
   nXtal_2_EE_low = '4'
   nXtal_2_EE_high = '4'
   S4S9_EB_low = '0.7'
   S4S9_EB_high = '0.7'
   S4S9_EE_low = '0.85'
   S4S9_EE_high = '0.85'
   if(isMC and MakeNtuple4optimization):
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
      #2015
      Pi0PtCutEB_low = '1.8'
      Pi0PtCutEB_high = '1.8'
      Pi0PtCutEE_low = '3.'
      Pi0PtCutEE_high = '3.'
      gPtCutEB_low = '1.5'
      gPtCutEB_high = '1.5'
      gPtCutEE_low = '1.5'
      gPtCutEE_high = '1.5'
      Pi0IsoCutEB_low = '0.'
      Pi0IsoCutEB_high = '0.'
      Pi0IsoCutEE_low = '0.'
      Pi0IsoCutEE_high = '0.'
      nXtal_1_EB_low = '7'
      nXtal_1_EB_high = '8'
      nXtal_2_EB_low = '7'
      nXtal_2_EB_high = '7'
      nXtal_1_EE_low = '6'
      nXtal_1_EE_high = '7'
      nXtal_2_EE_low = '4'
      nXtal_2_EE_high = '2'
      S4S9_EB_low = '0.9'
      S4S9_EB_high = '0.9'
      S4S9_EE_low = '0.9'
      S4S9_EE_high = '0.9'
#ETA
else:
   Pi0PtCutEB = '3.2'
   Pi0PtCutEE = '3.2'
   gPtCutEB = '1.4'
   gPtCutEE = '1.4'
   Pi0IsoCutEB = '0.1'
   Pi0IsoCutEE = '0.25'
   nXtal_1_EB = '5'
   nXtal_2_EB = '4'
   nXtal_1_EE = '7'
   nXtal_2_EE = '5'
   S4S9_EB = '0.9'
   S4S9_EE = '0.9'
   if(isMC and MakeNtuple4optimization):
      Pi0PtCutEB = '1.'
      Pi0PtCutEE = '1.'
      gPtCutEB = '0.4'
      gPtCutEE = '0.4'
      Pi0IsoCutEB = '0.'
      Pi0IsoCutEE = '0.'
      nXtal_1_EB = '0'
      nXtal_2_EB = '0'
      nXtal_1_EE = '0'
      nXtal_2_EE = '0'
      S4S9_EB = '0.6'
      S4S9_EE = '0.6'
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
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = False
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
#is2012 = False
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = False
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
#is2012 = False
#l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
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
HLTResults = 'False'
json_file = ''
is_2011 = 'False'
is2012 = False
l1InputTag =  "InputTag('hltGtDigis','', 'HLT')"
overWriteGlobalTag = False
doEnenerScale='False'
doIC='False'
doLaserCorr="True"
ebInputTag = "InputTag('reducedEcalRecHitsEB','','RECO')"
eeInputTag = "InputTag('reducedEcalRecHitsEE','','RECO')"
esInputTag = "InputTag('reducedEcalRecHitsES','','RECO')"
useHLTFilter = "False"
correctHits = 'False'
globaltag='POSTLS162_V2::All'
HLTPaths='AlCa_EcalPi0_*'
isMC = True
