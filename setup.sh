#!/bin/bash

RELEASE=$1
case $RELEASE in
	CMSSW_7_4_15)
		;;
	*)
		echo "[ERROR] Release not supported" >> /dev/stderr
		exit 1
		;;
esac

scram project $RELEASE
cd $RELEASE/src
eval `scramv1 runtime -sh`
git cms-init
git clone git@github.com:ECALELFS/ECALpro.git CalibCode
scram b -j16
git cms-addpkg RecoLocalCalo/EcalRecProducers
git cms-addpkg CondFormats/EcalObjects
git cms-addpkg RecoLocalCalo/EcalRecProducers
