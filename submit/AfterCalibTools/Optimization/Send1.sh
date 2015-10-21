#!/bin/bash
cd /afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/WD/CMSSW_7_4_0_pre8_March15/src/CalibCode/submit/AfterCalibTools/Optimization
eval `scramv1 runtime -sh`
root -b -q .x Step1_MakeHistos.C+
