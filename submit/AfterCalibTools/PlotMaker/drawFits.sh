#!/bin/bash                                                                        
                                                                                  
iter_ini=3
iter_fin=3  # it is included in sequence below                                   

eosPrefix="root://eoscms//eos/cms"                                        
wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"                             
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
dirName="AlcaP0_2016EF_mar2017newCond_reg2012_ext1"
tagName="AlcaP0_2016EF_mar2017newCond_reg2012_ext2_"

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
	rm $mergedFile
	echo "... and now for endcap ..."
	mergedFile="${tagName}Endcap_fitRes.root"
	hadd -f $mergedFile $endcap_list
	eos cp $mergedFile ${eosPrefix}${eosPath}${dirName}/iter_${i}/
	rm $mergedFile
    fi	

    iterNumber="iter_$i"
    echo  "iter_$i"
    root -l -b -q 'drawFits.C+("'$wwwPath'","'$eosPath'","'$dirName'","'$iterNumber'","'$tagName'")'
done

