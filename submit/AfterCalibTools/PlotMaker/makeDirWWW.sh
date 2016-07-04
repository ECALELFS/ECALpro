#!/bin/bash   

currentPath="$PWD"

wwwPath="/afs/cern.ch/user/m/mciprian/www/"  # your afs path to directory associated to your website
basePath="${wwwPath}pi0calib/ICplot/"        # here you will create the dirname folder (you can choose the name you want, add more folders or simply use one)
baseDir="AlcaP0_2016_json3p99fb_V2ext"       # dirname (could use any name, but better to stick with dirname in parameters.py)
nIter=4 # directory named from 0 to nIter

for i in `seq 0 $nIter`
do
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEp
    mkdir -p ${basePath}${baseDir}/iter_${i}/2DMaps/Endcap/EEm
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Barrel
    mkdir -p ${basePath}${baseDir}/iter_${i}/fitResPlots/Endcap
done

mkdir -p ${basePath}${baseDir}/TestConvergence

cd $wwwPath
./copyphp.sh   # you should already have this script to be able to see files in your website
cd $currentPath