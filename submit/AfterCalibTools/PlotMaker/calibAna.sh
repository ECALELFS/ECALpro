#!/bin/bash   

# launch as --> source calibAna.sh
# if you just want to compile and end script execution, use option -c

echo
echo "Compiling ..."

# to comile using root libraries I created an alias
#alias rootlib='$ROOTSYS/bin/root-config --cflags --libs'
# to be used like --> g++ -o file.exe file.C `rootlib`
# otherwise just use the following

g++ -Wall -pedantic -lm -o calibAnaEB calibAnaEB.C `$ROOTSYS/bin/root-config --cflags --libs`
if [ $? -ne 0 ]
then
    echo "Compilation of CalibAnaEB.C failed!"
    return 0
fi
echo "Compilation of CalibAnaEB.C succeeded :)"


g++ -Wall -pedantic -lm -o calibAnaEE calibAnaEE.C `$ROOTSYS/bin/root-config --cflags --libs`
if [ $? -ne 0 ]
then
    echo "Compilation of CalibAnaEE.C failed!"
    return 0
fi
echo "Compilation of CalibAnaEE.C succeeded :)"


for option in "$@";
do
    if [ $option = "-c" ]
    then
    # if option -c passed, we just want to complile and exit script                                                                                                     
        return 0
    fi
done

echo "Now launching executable ..."


iter_ini=4  # first iteration to use
iter_fin=4  # last iteration to use: it is included in sequence below
path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/"  # path to directory on eos
dirName="AlcaP0_2016_json3p99fb_V2ext"                            # dirname (see CalibCode/submit/parameters.py)
tagName="AlcaP0_2016_json3p99fb_V2ext_"                           # TagName (see CalibCode/submit/parameters.py) 

# iter_ini=0                                                          
# iter_fin=7  # it is included in sequence below                                                                                                                      
# path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"                                                                                                     
# dirName="AlcaP0_2016_json2p07minus0p8fb"                                               
# tagName="AlcaP0_2016_json2p07minus0p8fb_"

for i in `seq $iter_ini $iter_fin` 
do
    echo  "iter_$i --> EB"
    ./calibAnaEB $path $dirName $i $tagName
    echo  "iter_$i --> EE+"
    ./calibAnaEE EEp $path $dirName $i $tagName
    echo  "iter_$i --> EE-"
    ./calibAnaEE EEm $path $dirName $i $tagName
done
