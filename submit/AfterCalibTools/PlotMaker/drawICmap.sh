#!/bin/bash                                                                        
                                                                                  
# iter_ini=0
# iter_fin=7  # it is included in sequence below                                   
                                                                     
# path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
# dirName="AlcaP0_2016_json2p07fb"
# tagName="AlcaP0_2016_json2p07fb_"

#Emanuele
iter_ini=0
iter_fin=13  # it is included in sequence below                     

path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/"
dirName="pi0data_2015"
tagName="pi0data_2015"


for i in `seq $iter_ini $iter_fin`
do
    iterNumber="iter_$i"
    echo  "iter_$i"
    root -l -b -q 'drawICmap.C+("'$path'","'$dirName'","'$iterNumber'","'$tagName'")'
done

