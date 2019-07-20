#!/bin/bash   

currentPath="$PWD"

wwwPath="/afs/cern.ch/user/m/mciprian/www/"  # your afs path to directory associated to your website
basePath="${wwwPath}pi0calib/ICplot_Legacy/"        # here you will create the dirname folder (you can choose the name you want, add more folders or simply use one)
baseDir="AlCaP0_2018_ULrereco_1every2_ext1_fromIter6"       # dirname (could use any name, but better to stick with dirname in parameters.py
iter_ini=1 # first iteration to use (in general it would be 0)                                               
iter_fin=1  # last iteration to use: it is included in sequence below (if you did n iterations, this should be n-1)                  

for i in `seq $iter_ini $iter_fin`
do
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEp
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEm
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Endcap
    cp ${wwwPath}index.php ${basePath}${baseDir}/iter_${i}/2DMaps/Barrel/
    cp ${wwwPath}index.php ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEp/
    cp ${wwwPath}index.php ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEm/
    cp ${wwwPath}index.php ${basePath}${baseDir}/iter_${i}/fitResPlots/Barrel/
    cp ${wwwPath}index.php ${basePath}${baseDir}/iter_${i}/fitResPlots/Endcap/
done

mkdir -p ${basePath}${baseDir}/TestConvergence/extension   # use extension directory to store plots when you add extensions to a set of ICs
cp ${wwwPath}index.php ${basePath}${baseDir}/TestConvergence/extension/
cp ${wwwPath}index.php ${basePath}${baseDir}/TestConvergence/

# cd $wwwPath
# ./copyphp.sh   # you should already have this script to be able to see files in your website
# cd $currentPath
