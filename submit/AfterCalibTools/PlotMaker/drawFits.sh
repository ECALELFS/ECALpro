#!/bin/bash                                                                        
                                                                                  
iter_ini=3
iter_fin=3  # it is included in sequence below                                   
                                                                     
path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"
dirName="AlcaP0_2016_json2p07fb"
tagName="AlcaP0_2016_json2p07fb_"

for i in `seq $iter_ini $iter_fin`
do
    iterNumber="iter_$i"
    echo  "iter_$i"
    root -l -b -q 'drawFits.C+("'$path'","'$dirName'","'$iterNumber'","'$tagName'")'
done

