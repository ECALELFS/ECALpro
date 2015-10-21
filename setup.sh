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
cmsenv
git cms-init
git clone git@github.com:ECALELFS/ECALpro.git CalibCode
scram b -j16
