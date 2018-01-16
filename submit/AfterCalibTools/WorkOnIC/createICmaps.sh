#!/bin/bash

# compute IC and draw maps normalized to 1 in eta rings

minArgs=3
thisScriptName=$(basename $BASH_SOURCE)

if [[ ( "$#" < ${minArgs} ) ]]; then
    echo "Usage: source ${thisScriptName} dirName nIter year drawonly[yes:0/no:1] username [options]"
    return 0
fi

dirName="$1"
nIter="$2"
year="$3"
drawonly="$4"
username="$5"
otherOptions="$6"
outDir="./"

if [[ "X${username}" == "X" ]]; then
    username="mciprian"
fi
  
if [ "${drawonly}" != "0" ]; then
    echo ""
    echo "Running MultiplyIC_txt_root.py"
    echo ""
    ./MultiplyIC_txt_root.py /afs/cern.ch/user/m/mciprian/public/ECALproTools/dump_EcalIntercalibConstants/dump_EcalIntercalibConstants__since_00271951_till_18446744073709551615.dat root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero${year}/${username}/${dirName}/iter_${nIter}/${dirName}_calibMap.root ic_${year}_${dirName}_iter${nIter}/ Absolute_IC.root ${otherOptions}
fi

echo ""

if [[ "${username}" == "mciprian" ]]; then
    outDir="${HOME}/www/pi0calib/ICplot/${dirName}/iter_${nIter}/2DMaps/ICmaps/"
    echo "I will save plots in folder: ${outDir}"
else
    echo "WARNING: I will save maps in current folder! If you want a specific path you should modify the script ${thisScriptName}"
fi 

echo "Drawing IC not normalized to 1 in etaRing"
echo ""
python DrawIC.py ic_${year}_${dirName}_iter${nIter}/IC_fromECALpro.txt -o ${outDir}notNormalized/ --max-EB 0.1 --max-EE 0.25 ${otherOptions}
echo ""
echo "Drawing IC normalized to 1 in etaRing"
echo ""
python DrawIC.py ic_${year}_${dirName}_iter${nIter}/IC_fromECALpro.txt -o ${outDir}norm1etaring/ --max-EB 0.1 --max-EE 0.25 --normalize-etaring ${otherOptions}
echo ""
cp /afs/cern.ch/user/m/mciprian/public/index.php ${outDir}notNormalized/
cp /afs/cern.ch/user/m/mciprian/public/index.php ${outDir}norm1etaring/


