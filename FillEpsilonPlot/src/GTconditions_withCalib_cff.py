import FWCore.ParameterSet.Config as cms


GTconditions = cms.VPSet(
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
