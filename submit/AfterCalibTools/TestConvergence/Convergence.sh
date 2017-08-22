#!/bin/bash

eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/" 
#eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/" 
#eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/" 
dirName="AlCaP0_IC2017_upTo31July2017_noCC"                            # dirname (see CalibCode/submit/parameters.py)  
iter_number="5"                                                          # number n of iterations (iter_0 to iter_{n-1})
tagName="${dirName}_"                           # TagName (see CalibCode/submit/parameters.py)  

# will copy output here, if directory exists
wwwTargetDir="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/${dirName}/TestConvergence/"               

nJump=1
# leave extension as "noExtension" in you don't need to add additional steps that start from the one above
# format is newDirName_ext1,newIterNumber_ext1,newTagName_ext1:newDirName_ext2,newIterNumber_ext2,newTagName_ext2 and so on (different extensions separated by : )
extension="noExtension"
#extension="AlcaP0_2016CD_mar2017newCond_reg2012_fromIter1Run2016CD,5,AlcaP0_2016CD_mar2017newCond_reg2012_fromIter1Run2016CD_"
detectorToSkip="no"   # detectorToSkip = "no" to skip nothing, "EB" to skip EB, "EE" to skip EE

for option in "$@";
do
    if [ "$option" = "-noEB" ]; then
        detectorToSkip="EB"
    elif [ "$option" = "-noEE" ]; then
        detectorToSkip="EE"
    elif [ "$option" = "-noext" ]; then
	extension="noExtension"
    fi
done


if [ "${extension}" != "noExtension" ]
then
    echo "Extension added: will store plots in --> ${wwwTargetDir}extension/ (if it exists)"
    wwwTargetDir=${wwwTargetDir}extension
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
./Convergence $eosPath $dirName $iter_number $tagName $nJump $extension $detectorToSkip

if [ $? -ne 0 ]
then
    echo "An error occurred! Exit"
    return 0
fi
# copy output to wwwTargetDir if it exists and remove local directory, otherwise just keep in local
test -d ${wwwTargetDir}/ && cp ./plot_${dirName}/* ${wwwTargetDir}/ && rm -r ./plot_${dirName}/ 


echo "THE END!"
