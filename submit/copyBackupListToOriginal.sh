#! /bin/bash

# this script is useful when for some reason the Hadd josb started before all the Fill ones were over
# in this case, the list of EcalNtp files passed to the hadd jobs are screwed up
# therefore, we restore the original situation

# specify the folder and the iteration number to be affected
# note that some files are removed, so be careful!

folder="AlCaP0_2016_ULrereco_from0"
iter="4"

path="${folder}/src/hadd/"
fullpath="$PWD/${folder}/src/hadd/"
iter_match="iter_${iter}_"

#echo "${fullpath}"
backupFiles=`ls ${fullpath} | grep ${iter_match} | grep _backup.list`
prunedFiles=`ls ${fullpath} | grep ${iter_match} | grep _pruned.list`

for f in ${backupFiles};
do
    echo "Copying ${path}${f} to ${path}${f/_backup/}"
    cp ${path}${f} ${path}${f/_backup/}  # remove _backup, so copy to the original file
    echo "Removing ${path}${f}"
    rm ${path}${f}                      # remove backup
done


for f in ${prunedFiles};
do 
    echo "Removing file ${path}${f}"
    rm ${path}${f}                      # remove pruned file (no longer needed)
done
