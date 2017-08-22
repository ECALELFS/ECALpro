#!/bin/bash                                                                        
                     
if [[ ( "$#" < 3 ) ]]; then
    thisScriptName=$(basename $BASH_SOURCE)
    echo "Usage: source ${thisScriptName} eosPath dirName iter_ini [iter_fin]"
    echo "An example of eosPath is: /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
    echo "iter_fin is optional"
    return 0
fi

eosPath="$1"
echo "Selected path: ${eosPath}"
dirName="$2"
echo "Directory: ${dirName}"                                                             
iter_ini="$3"
echo "Initial iteration: ${iter_ini}"
if [[ "X$4" != "X" ]]; then
    iter_fin="$4"
else
    iter_fin="${iter_ini}"
fi
echo "Final iteration: ${iter_fin}"

if [ ! -d "${eosPath}${dirName}/" ]; then
    echo "Error: folder ${eosPath}${dirName}/ does not exist! Abort."
    return 0
fi

tagName="${dirName}_"

host=`echo "$HOSTNAME"`

if [[ ${host} != *"lxplus"* ]]; then
  echo "Error! You must be on lxplus because you need eos to be mounted. Do ssh -XY lxplus and work from a release."
  return 0
fi

release=`echo ${CMSSW_BASE}`
if [[ ${release} != *"CMSSW"* ]]; then
  echo "Error! You must be in a release"
  return 0
fi

eval `scramv1 runtime -sh`


for i in `seq $iter_ini $iter_fin`
do

    if [ ! -d "${eosPath}${dirName}/iter_${i}/" ]; then
	echo "Error: folder ${eosPath}${dirName}/iter_${i}/ does not exist! Abort."
	return 0
    fi


    barrel_list=""
    for j in {0..30}
    do
	file="${eosPath}${dirName}/iter_${i}/${tagName}Barrel_${j}_fitRes.root"
	barrel_list="${barrel_list} ${file}" 
    done
    
    endcap_list=""
    for j in {0..7}
    do  
        file="${eosPath}${dirName}/iter_${i}/${tagName}Endcap_${j}_fitRes.root"
        endcap_list="${endcap_list} ${file}"
    done

    echo "iteration ${i}: going to merge files on EOS ..."	
    echo "... for barrel ..."
    mergedFile="${eosPath}${dirName}/iter_${i}/${tagName}AllBarrel_fitRes.root"
    hadd -f -k $mergedFile $barrel_list
    echo "... and now for endcap ..."
    mergedFile="${eosPath}${dirName}/iter_${i}/${tagName}AllEndcap_fitRes.root"
    hadd -f -k $mergedFile $endcap_list
    echo  ""

done

