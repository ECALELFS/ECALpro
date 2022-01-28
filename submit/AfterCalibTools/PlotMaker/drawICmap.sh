#!/bin/bash                                                                        

# usage:  source drawICmap.sh [options]
# available options
# -noEB  -->  skip EB
# -noEE  -->  skip EE
                                                                                  
# iter_ini=0
# iter_fin=7  # it is included in sequence below                                   
                                                                     
# path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
# dirName="AlcaP0_2016_json2p07fb"
# tagName="AlcaP0_2016_json2p07fb_"

iter_ini=0
iter_fin=6  # it is included in sequence below                     

wwwPath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"
eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
#eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/"
#eosPath="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/"
dirName="AlCaP0_IC2017_upTo21September2017"
tagName="${dirName}_"

ECALdetToSkip=""

for option in "$@";
do
    if [ "$option" = "-noEB" ]; then
	ECALdetToSkip="EB"
    elif [ "$option" = "-noEE" ]; then
        ECALdetToSkip="EE"
    fi
done


for i in `seq $iter_ini $iter_fin`
do
    iterNumber="iter_$i"
    echo  "iter_$i"
    root -l -b -q 'drawICmap.C+("'$wwwPath'","'$eosPath'","'$dirName'","'$iterNumber'","'$tagName'","'$ECALdetToSkip'")'
done

