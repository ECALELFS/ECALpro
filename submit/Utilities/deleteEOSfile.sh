#! /bin/bash

# now EOS is mounted, but you must be on lxplus

host=`echo "$HOSTNAME"`
if [[ ${host} != *"lxplus"* ]]; then
    echo "Error! You must be on lxplus to use this script. Do ssh -XY lxplus and work from a release."
    return 0
fi

iter_ini=4
iter_fin=4  # it is included in sequence below                                                                          

eosPath="/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/"
dirName="AlCaEta_2016_ULrereco_ext1_fromIter4"

# you can use "epsilonPlots_" as pattern to delete all directory with the mass distributions. The ending underscore prevents the merged "*epsilonPlots.root" file
# from being deleted as well (you might want to keep it)

patterns=( "EcalNtp" "epsilonPlots_" )

#pattern="epsilonPlots_"
#pattern="EcalNtp"  # use it with grep to select which file to remove
# use following string to test if eos directory exists: we use a regular expression to test whether this string is in the output of "eos ls ..." 
noDirFound="No such file or directory" 

echo ""

for pattern in "${patterns[@]}"
do

    echo "==================================="
    echo "Pattern --> ${pattern}"
    echo "==================================="

    for i in `seq $iter_ini $iter_fin`
    do

    	thisFolder="${eosPath}${dirName}/iter_${i}"
    	echo "Testing existence of ${eosPath}${dirName}/iter_${i}"

    	if [ ! -d "${thisFolder}" ]; then
    	    echo "WARNING: no folder named ${thisFolder}"
    	else

    	    echo "Ok, directory exists :)"
    	    filesToRemove=`ls ${thisFolder} | grep ${pattern}`
    	    if [ "${filesToRemove}" == "" ]; then
    		echo "No files in ${thisFolder} matching '${pattern}'"
    	    else 
		nFilesToRemove=`ls ${thisFolder} | grep ${pattern} | wc -l`
    		echo "Removing ${nFilesToRemove} files matching '${pattern}' in ${thisFolder}"
    		for thisfile in $filesToRemove
    		do
    		    rm ${thisFolder}/${thisfile}
    		done	
    	    fi
	    
    	fi

    	echo ""

    done

done

echo "THE END!"
