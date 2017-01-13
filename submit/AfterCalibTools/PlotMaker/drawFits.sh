#!/bin/bash                                                                        
                                                                                  
iter_ini=7
iter_fin=7  # it is included in sequence below                                   

eosPrefix="root://eoscms//eos/cms"                                        
wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"                             
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
dirName="AlcaP0_2016_json3p99fb_weight_8iter_noCC"
tagName="AlcaP0_2016_json3p99fb_weight_8iter_noCC_"

skipHaddFile="no"   # yes or no, decide if you need to merge files on eos, if not, skipping saves a lot of time

for i in `seq $iter_ini $iter_fin`
do

    if [[ "${skipHaddFile}" == "no" ]]; then
	barrel_list=""
	for j in {0..30}
	do
	    file="${eosPrefix}${eosPath}${dirName}/iter_${i}/${tagName}Barrel_${j}_fitRes.root"
	    barrel_list="${barrel_list} ${file}" 
	done
	
	endcap_list=""
	for j in {0..7}
	do  
            file="${eosPrefix}${eosPath}${dirName}/iter_${i}/${tagName}Endcap_${j}_fitRes.root"
            endcap_list="${endcap_list} ${file}"
	done

	echo "Going to merge files on EOS ..."	
	echo "... for barrel ..."
	mergedFile="${tagName}Barrel_fitRes.root"
	hadd -f $mergedFile $barrel_list
	eos cp $mergedFile ${eosPrefix}${eosPath}${dirName}/iter_${i}/
	echo "... and now for endcap ..."
	mergedFile="${tagName}Endcap_fitRes.root"
	hadd -f $mergedFile $endcap_list
	eos cp $mergedFile ${eosPrefix}${eosPath}${dirName}/iter_${i}/
    fi	

    iterNumber="iter_$i"
    echo  "iter_$i"
    root -l -b -q 'drawFits.C+("'$wwwPath'","'$eosPath'","'$dirName'","'$iterNumber'","'$tagName'")'
done

