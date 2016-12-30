#!/bin/bash

eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/" 
dirName="AlcaP0_2016_json3p99fb_weight_8iter_noCC"                            # dirname (see CalibCode/submit/parameters.py)  
iter_number="8"                                                          # number n of iterations (iter_0 to iter_{n-1})
tagName="AlcaP0_2016_json3p99fb_weight_8iter_noCC_"                           # TagName (see CalibCode/submit/parameters.py)  

# will copy output here, if directory exists
wwwTargetDir="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/${dirName}/TestConvergence/"               

nJump=1
# leave extension as "noExtension" in you don't need to add additional steps that start from the one above
# format is newDirName,newIterNumber,newTagName
extension="noExtension"
#extension="AlcaP0_2016_json3p99fb_weight_ext,4,AlcaP0_2016_json3p99fb_weight_ext_"

if [ "${extension}" != "noExtension" ]
then
    echo "Extension added: will store plots in --> ${wwwTargetDir}extension/ (if it exists)"
    wwwTargetDir=${wwwTargetDir}extension/
fi

#compile Convergence.C
g++ -Wall -pedantic -lm -o Convergence Convergence.C `$ROOTSYS/bin/root-config --cflags --libs` 

if [ $? -ne 0 ]
then
    echo "Compiling Convergence.C : failed!"
    return 0
fi
echo "Compiling Convergence.C : success :)"

echo "Now runnning Convergence.C"
./Convergence $eosPath $dirName $iter_number $tagName $nJump $extension

if [ $? -ne 0 ]
then
    echo "An error occurred! Exit"
    return 0
fi
# copy output to wwwTargetDir if it exists, otherwise just keep in local
test -d ${wwwTargetDir}/ && cp ./plot_${dirName}/* ${wwwTargetDir}/

echo "THE END!"
