nEventsPerJob = '-1' # default -1
eosPath = '/store/group/alca_ecalcalib/lpernie/'
#eosPath = '/store/caf/user/lpernie'
outputFile    = 'EcalNtp' # without .root suffix
calibMapName = 'calibMap.root'
ExternalGeometry = 'caloGeometry.root'
CalibType  = 'xtal'
Are_pi0  = False # True = using Pi0, False = using Eta
#Pi0
if(Are_pi0):
   Pi0PtCutEB = '2.1'
   Pi0PtCutEE = '2.1'
   gPtCutEB = '0.8'
   gPtCutEE = '0.8'
   Pi0IsoCutEB = '0.'
   Pi0IsoCutEE = '0.25'
   nXtal_1_EB = '6'
   nXtal_2_EB = '4'
   nXtal_1_EE = '6'
   nXtal_2_EE = '4'
   S4S9_EB = '0.7'
   S4S9_EE = '0.85'
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

inputlist_n      = 'ALL_2012D_good.list' # list of the input files
ijobmax          = 5                     # 5 number of files per job
nHadd            = 35                    # 50 number of files per hadd
nFit             = 2000                  # number of fits done in parallel
Barrel_or_Endcap = 'ALL_PLEASE'          # Option: 'ONLY_BARREL','ONLY_ENDCAP','ALL_PLEASE'
dirname          = 'ALL_2012D_Eta_03'
Silent           = False                 # True->Fill modules is silent; False->Fill modules has a standard output
NameTag          = '2012D_'              # Tag to the names to avoid overlap
queueForDaemon   = 'cmscaf1nw'           # Option suggested: 2nw/2nd, 1nw/1nd, cmscaf1nw/cmscaf1nd... even cmscaf2nw
queue            = 'cmscaf1nd'

nIterations = 10

#-----------------------------------------------------------------------------------
alphaTagRecord2='';alphaTag2='';alphaDB2=''
GeVTagRecord='';GeVTag='';GeVDB=''

######################################################################
# Now decomment the part that correspond to data you want to run on. #
######################################################################
##2012D
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt
json_file ='goodrunlist_json2012D.txt'
HLTResults = 'True'
is2012 = True
is_2011 = 'True' #Just for the fit, put true 
useHLTFilter="True"
correctHits='True'
overWriteGlobalTag = True
globaltag='GR_P_V40::All'
if(Are_pi0): 
   esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
   ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
   eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
   HLTPaths='AlCa_EcalPi0*' 
else:
   esInputTag = "InputTag('hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES', 'HLT')"
   ebInputTag = "InputTag('hltAlCaEtaEBUncalibrator','etaEcalRecHitsEB','HLT')"
   eeInputTag = "InputTag('hltAlCaEtaEEUncalibrator','etaEcalRecHitsEE','HLT')"
   HLTPaths='AlCa_EcalEta*'
doEnenerScale='True'
doIC='True'
doLaserCorr="True"
GeVTagRecord='EcalADCToGeVConstantRcd'
GeVTag='EcalADCToGeVConstant_Bon_V20111129'
GeVDB='frontier://FrontierProd/CMS_COND_31X_ECAL'
laserTagRecord='EcalIntercalibConstantsRcd'
laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
alphaTagRecord2='EcalLaserAlphasRcd'
alphaTag2='EcalLaserAlphas_EB_sic1_btcp152_EE_sic1_btcp116'
alphaDB2='frontier://FrontierInt/CMS_COND_ECAL'
alphaTagRecord='EcalLaserAPDPNRatiosRcd'
alphaTag='EcalLaserAPDPNRatios_20130124_447_p1_v2'
alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
##############
#2012C
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt
#json_file ='goodrunlist_json2012C.txt'
#HLTResults = 'True'
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#globaltag='GR_P_V42::All'
#ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#laserTagRecord='EcalIntercalibConstantsRcd'
#laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#alphaTag='EcalLaserAPDPNRatios_20121020_447_p1_v2'
#alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#HLTPaths='AlCa_EcalPi0*' 

##############
##2012B
###/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208357_8TeV_PromptReco_Collisions12_JSON.txt
#json_file ='goodrunlist_json2012C.txt'
#HLTResults = 'True'
#is2012 = True
#is_2011 = 'True' #Just for the fit, put true 
#useHLTFilter="True"
#correctHits='True'
#overWriteGlobalTag = True
#globaltag='FT_R_53_V6::All'
#ebInputTag = "InputTag('hltAlCaPi0EBUncalibrator','pi0EcalRecHitsEB','HLT')"
#eeInputTag = "InputTag('hltAlCaPi0EEUncalibrator','pi0EcalRecHitsEE','HLT')"
#esInputTag = "InputTag('hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES', 'HLT')"
#doEnenerScale='True'
#doIC='True'
#doLaserCorr="True"
#laserTagRecord='EcalIntercalibConstantsRcd'
#laserTag = 'EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies'
#laserDB  = 'frontier://FrontierInt/CMS_COND_ECAL'
#alphaTagRecord='EcalLaserAPDPNRatiosRcd'
#alphaTag='EcalLaserAPDPNRatios_20121020_447_p1_v2'
#alphaDB='frontier://FrontierProd/CMS_COND_42X_ECAL_LAS'
#HLTPaths='AlCa_EcalPi0*'


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
