#!/bin/bash                                                                        
                                                                                  
iter_ini=0
iter_fin=0  # it is included in sequence below                                   

eosPrefix="root://eoscms//eos/cms"                                        
wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"                             
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
dirName="AlCaEta_Run2017_23June2017_7xtal_v3"
tagName="${dirName}_"

useMergedFitFile=false # when true, no need to specify a file index
BarrelOrEndcap="Barrel"  # Barrel, Endcap
fileIndex=15 # index for EB goes from 0 to 30 and for EE it goes from 0 to 7
#fileIndex=3
#BarrelOrEndcap="Endcap"  # Barrel, Endcap

nFitToPlot=150  # there are at most 2000 plots in each file

for i in `seq $iter_ini $iter_fin`
do

    outputDIR=${wwwPath}${dirName}"/iter_${i}/fitResPlots/${BarrelOrEndcap}/"
    echo "iter_${i}"
    if [ "$useMergedFitFile" = true ]; then
	file="${eosPrefix}${eosPath}${dirName}/iter_${i}/${tagName}All${BarrelOrEndcap}_fitRes.root" 
    else
	file="${eosPrefix}${eosPath}${dirName}/iter_${i}/${tagName}${BarrelOrEndcap}_${fileIndex}_fitRes.root" 
    fi
    echo "file --> ${file}"
    echo "output directory --> ${outputDIR}"
    root -l -b -q 'drawFitsSingleFile.C+("'${file}'","'${BarrelOrEndcap}'","'${outputDIR}'",'${nFitToPlot}')'

done

