#!/bin/bash


inputFiles=('root://cms-xrd-global.cern.ch//store/data/Run2018D/AlCaP0/RAW/v1/000/321/396/00000/248ADA80-FBA1-E811-8742-FA163E7B2F96.root')
#globaltag='105X_dataRun2_v8'
globaltag='121X_mcRun3_2021_realistic_v18'
#jsonFile='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
jsonFile=''
outputDir='./'
outputFile='AlCaP0.root'
conditionsFileName='GTconditions_withCalib_cff' ###for without calibration - conditionsFileName='GTconditions_withoutCalib_cff'
#cd /afs/cern.ch/user/e/ecalgit/pi0Monitoring4Run3/CMSSW_10_6_8/src/CalibCode/submit
cd /afs/cern.ch/user/e/ecalgit/pi0Monitoring4Run3/CMSSW_12_2_0/src/CalibCode/submit/monitoring/
eval `scramv1 runtime -sh`
cmsRun reconstructPi0_template.py inputFiles=${inputFiles} globaltag=${globaltag} jsonFile=${jsonFile} outputDir=${outputDir} outputFile=${outputFile} conditionsFileName=${conditionsFileName}
