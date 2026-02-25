#! /bin/bash

# index=`root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,1)' | grep int | awk '{ print $2 }'`
# echo "${index}"

# at the borders of SM and modules
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,1)'
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,65)'
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,85)'

# at the borders of SM but not modules
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,15)'
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,55)'
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(1,75)'

# last crsytal in eta, but far from SM border
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,85)' 

# far from SM border and at ieta=84, where there seems to be a boundary for photon2
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(4,84)' 

# far from gaps
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,15)' 
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,35)' 
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,55)' 
root -l -b -q 'getFitIndex_from_iphiix_ietaiy.C+(10,75)' 


