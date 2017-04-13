#! /bin/bash

iter_ini=0
iter_fin=4  # it is included in sequence below                                                                           

eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
dirName="AlcaP0_2016CD_mar2017newCond_reg2012_fromIter1Run2016CD"

# you can use "epsilonPlots_" as pattern to delete all directory with the mass distributions. The ending underscore prevents the merged "*epsilonPlots.root" file
# from being deleted as well (you might want to keep it)

pattern="epsilonPlots_"
#pattern="EcalNtp"  # use it with grep to select which file to remove
# use following string to test if eos directory exists: we use a regular expression to test whether this string is in the output of "eos ls ..." 
noDirFound="No such file or directory" 

for i in `seq $iter_ini $iter_fin`
do
    eos_ls_output=`eos ls ${eosPath}${dirName}/iter_${i}`
    echo "Testing existence of ${eosPath}${dirName}/iter_${i}"

    if [[ ${eos_ls_output} =~ ${noDirFound} ]]; then
	echo "Directory ${eosPath}${dirName}/iter_${i} not found!"
    else

	echo "Ok, directory exists :)"
	filesToRemove=`eos ls ${eosPath}${dirName}/iter_${i} | grep ${pattern}`
	if [ "${filesToRemove}" == "" ]; then
	    echo "No files in ${eosPath}${dirName}/iter_${i} matching '${pattern}'"
	else 
	    echo "Removing files matching '${pattern}' in ${eosPath}${dirName}/iter_${i}"
	    for thisfile in $filesToRemove
	    do
		eos rm ${eosPath}${dirName}/iter_${i}/${thisfile}
	    done	
	fi
	
    fi

done

echo "THE END!"