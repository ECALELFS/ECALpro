#!/bin/bash                                                                        
                                                                                  
iter_ini=4
iter_fin=4  # it is included in sequence below                                   

eosPrefix="root://eoscms//eos/cms"                                        
wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"                             
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
#dirName="AlcaP0_Run2016G_sel16_reg12"
dirName="AlcaP0_Run2016G_sel17optim_reg12"
tagName="${dirName}_"

useMergedFitFile=false # when true, no need to specify a file index
BarrelOrEndcap="Barrel"  # Barrel, Endcap
fileIndex=0 # index for EB goes from 0 to 30 and for EE it goes from 0 to 7
#fileIndex=3
#BarrelOrEndcap="Endcap"  # Barrel, Endcap

nFitsToPlot=10  # there are at most 2000 plots in each file

# if fitIndexToPlot>=0, look for this specific index in the rooplots in the file and just plot that one
# if negative just plot nFistToPlot plots
#fitIndexToPlot="54049"  
fitIndexToPlot="-1"  

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
    root -l -b -q 'drawFitsSingleFile.C+("'${file}'","'${BarrelOrEndcap}'","'${outputDIR}'",'${nFitsToPlot}','${fitIndexToPlot}')'

done

