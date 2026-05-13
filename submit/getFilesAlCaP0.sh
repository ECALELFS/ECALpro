#! /bin/bash

thisYear="26"   # use 16, 17, 18
dayMonthYear=`date +%d_%m_%Y`
dataset="AlCaP0"
runYear="20${thisYear}"
# eras=(B C D E F G)
eras=(D)
#JsonFilter="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions${thisYear}/13TeV/DCSOnly/json_DCSONLY.txt"   # caution to the year
JsonFilter="/afs/cern.ch/user/d/delvecch/CMSSW_15_0_18/src/CalibCode/submit/InputList/Collisions26_13p6TeV_Latest.json"
firstRunBrilcalc="" # can be empty, otherwise it is used to select a run range for the brilcalc command
lastRunBrilcalc="" # can be empty, otherwise it is used to select a run range for the brilcalc command

ecalproFolder="${CMSSW_BASE}/src/CalibCode/submit/"
outputdir="${ecalproFolder}InputList/"
# outputfile="${dataset}_Run${runYear}${dataEra}_${dayMonthYear}.list"
# outputfile="${dataset}_Run${runYear}${eras[0]}${eras[-1]}_${dayMonthYear}.list"
outputfile="${dataset}_Run${runYear}${eras[0]}_${dayMonthYear}.list"
fileList="${outputdir}${outputfile}"
brilcalcOutputFile="${ecalproFolder}brilcalc_${outputfile/list/txt}"

cd ${ecalproFolder}
eval `scramv1 runtime -sh`

host=`echo "$HOSTNAME"`
if [[ ${host} != *"lxplus"* ]]; then
    echo "Error! You must be on lxplus to get list of files. Do ssh -XY lxplus and work from a release."
    return 0
fi

echo "Creating list of file running dasgoclient"
echo "" > ${fileList}
for dataEra in ${eras[@]} 
do 
    echo "dasgoclient --query=\"file dataset=/${dataset}/Run${runYear}${dataEra}-v1/RAW\""
    dasgoclient -query="file dataset=/${dataset}/Run${runYear}${dataEra}-v1/RAW" >> ${fileList}
done
echo "List created in ${fileList}"

echo ""
echo "Filtering list of files with json ${JsonFilter}"
purifyCmd="python3 Utilities/Purify_List.py ${fileList} ${JsonFilter}"
echo "${purifyCmd}"
echo "${purifyCmd}" | bash

#echo ""
#echo "Checking integrated luminosity with brilcalc"
#brilcalcCmd="brilcalc lumi -u /fb -i ${JsonFilter} --without-checkjson --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json" # normtag_PHYSICS.json 
#if [[ "X${firstRunBrilcalc}" != "X" ]]; then
#    brilcalcCmd="${brilcalcCmd} --begin ${firstRunBrilcalc}"
#fi
#if [[ "X${lastRunBrilcalc}" != "X" ]]; then
#    brilcalcCmd="${brilcalcCmd} --end ${lastRunBrilcalc}"
#fi
#echo "${brilcalcCmd}"
#echo "Will save output in file ${brilcalcOutputFile}"
#echo "${brilcalcCmd} > ${brilcalcOutputFile}" | bash

echo ""
echo "DONE!" 
echo ""
