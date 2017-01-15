#!/bin/bash

eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/" 
dirName="AlcaP0_2016H_23SeptReReco_EBonly_reg2012"                            # dirname (see CalibCode/submit/parameters.py)  
iter_number="7"                                                          # number n of iterations (iter_0 to iter_{n-1})
tagName="AlcaP0_2016H_23SeptReReco_EBonly_reg2012_"                           # TagName (see CalibCode/submit/parameters.py)  

# will copy output here, if directory exists
wwwTargetDir="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/${dirName}/TestConvergence/"               

nJump=1
# leave extension as "noExtension" in you don't need to add additional steps that start from the one above
# format is newDirName,newIterNumber,newTagName
extension="noExtension"
#extension="AlcaP0_2016_json3p99fb_weight_extV2,4,AlcaP0_2016_json3p99fb_weight_extV2_:AlcaP0_2016_json3p99fb_weight_extV2_4more,4,AlcaP0_2016_json3p99fb_weight_extV2_4more_"
detectorToSkip="EE"   # detectorToSkip = "" to skip nothing, "EB" to skip EB, "EE" to skip EE

if [ "${extension}" != "noExtension" ]
then
    echo "Extension added: will store plots in --> ${wwwTargetDir}extension/ (if it exists)"
    wwwTargetDir=${wwwTargetDir}extension/
fi

for option in "$@";
do
    if [ "$option" = "-noEB" ]; then
        detectorToSkip="EB"
    elif [ "$option" = "-noEE" ]; then
        detectorToSkip="EE"
    fi
done


#compile Convergence.C
g++ -Wall -pedantic -lm -o Convergence Convergence.C `$ROOTSYS/bin/root-config --cflags --libs` 

if [ $? -ne 0 ]
then
    echo "Compiling Convergence.C : failed!"
    return 0
fi
echo "Compiling Convergence.C : success :)"

echo "Now runnning Convergence.C"
./Convergence $eosPath $dirName $iter_number $tagName $nJump $extension $detectorToSkip

if [ $? -ne 0 ]
then
    echo "An error occurred! Exit"
    return 0
fi
# copy output to wwwTargetDir if it exists and remove local directory, otherwise just keep in local
test -d ${wwwTargetDir}/ && cp ./plot_${dirName}/* ${wwwTargetDir}/ && rm -r ./plot_${dirName}/

echo "THE END!"
