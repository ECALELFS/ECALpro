#!/bin/bash                                                                        
                                                                                  
iter_ini=4
iter_fin=4  # it is included in sequence below                                   

eosPrefix="root://eoscms//eos/cms"                                        
wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"                             
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
dirName="AlcaP0_2016CD_mar2017newCond_reg2012_fromIter1Run2016CD"
tagName="AlcaP0_2016CD_mar2017newCond_reg2012_fromIter1Run2016CD_"

BarrelOrEndcap="Endcap"  # Barrel, Endcap
fileIndex=2 # index for EB goes from 0 to 30 and for EE it goes from 0 to 7
file="${eosPrefix}${eosPath}${dirName}/iter_${i}/${tagName}${BarrelOrEndcap}_${fileIndex}_fitRes.root" 

nFitToPlot=10  # there are at most 2000 plots in each file

for i in `seq $iter_ini $iter_fin`
do

    outputDIR=${wwwPath}${dirName}"/iter_${i}/fitResPlots/${BarrelOrEndcap}/"
    echo "iter_${i}"
    echo "file --> ${file}"
    echo "output directory --> ${outputDIR}"
    root -l -b -q 'drawFitsSingleFile.C+("'${file}'","'${BarrelOrEndcap}'","'${outputDIR}'",'${nFitToPlot}')'

done

