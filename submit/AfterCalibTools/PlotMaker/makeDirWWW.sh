#!/bin/bash   

currentPath="$PWD"

wwwPath="/afs/cern.ch/user/m/mciprian/www/"  # your afs path to directory associated to your website
basePath="${wwwPath}pi0calib/ICplot/"        # here you will create the dirname folder (you can choose the name you want, add more folders or simply use one)
baseDir="AlcaP0_2016_json3p99fb_weight_8iter_noCC"       # dirname (could use any name, but better to stick with dirname in parameters.py)
iter_ini=0  # first iteration to use (in general it would be 0)                                                                                   
iter_fin=7  # last iteration to use: it is included in sequence below (if you did n iterations, this should be n-1)                                                  

for i in `seq $iter_ini $iter_fin`
do
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEp
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEm
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Endcap
done

mkdir -p ${basePath}${baseDir}/TestConvergence/extension   # use extension directory to store plots when you add extensions to a set of ICs

cd $wwwPath
./copyphp.sh   # you should already have this script to be able to see files in your website
cd $currentPath