#!/bin/bash   

# launch as --> source calibAna.sh
# if you just want to compile and end script execution, use option -c

##############################################
# THINGS TO CHANGE BEFORE LAUNCHING SCRIPT
##############################################


iter_ini=6                                                        # first iteration to use
iter_fin=6                                                        # last iteration to use: it is included in sequence below
path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/"  # path to directory on eos
dirName="AlcaP0_2016_json3p99fb_V2ext"                            # dirname (see CalibCode/submit/parameters.py)
tagName="AlcaP0_2016_json3p99fb_V2ext_"                           # TagName (see CalibCode/submit/parameters.py) 

# iter_ini=0                                                          
# iter_fin=7  # it is included in sequence below                                                                                                                      
# path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"                                                                                                     
# dirName="AlcaP0_2016_json2p07minus0p8fb"                                               
# tagName="AlcaP0_2016_json2p07minus0p8fb_"

wwwBasePath="/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/"   # base directory where plots are stored (other directories are created inside this path) 

# list of source files to compile (except main) without extension (which is supposed to be .C in the following)
sourceFileList=(calibAnaEcal_base calibAnaEcal calibAnaEcalEB calibAnaEcalEE)  
mainSourceFile="main"

###############################################
###############################################

###############################################
#flags used in the script to decide what to do (can be set by user through options passed)
onlyCompile=false
doNotCompile=false
skipEB=false
skipEEp=false
skipEEm=false

#########################################
#parsing options
for option in "$@";
do
    if [ $option = "-c" ]; then
        onlyCompile=true
    elif [ $option = "-nc" ]; then
        doNotCompile=true
    elif [ $option = "-noEB" ]; then
	skipEB=true
    elif [ $option = "-noEEp" ]; then
	skipEEp=true
    elif [ $option = "-noEEm" ]; then
	skipEEm=true
    fi
done

########################################

olist=""
sourceFile=""

if [ "$doNotCompile" = false ]; then

    echo
    echo "Compiling ..."

# to compile using root libraries I created an alias
#alias rootlib='$ROOTSYS/bin/root-config --cflags --libs'
# to be used like --> g++ -o file.exe file.C `rootlib`
# otherwise just use "$ROOTSYS/bin/root-config --cflags --libs" as below
    
    for sourceFile in ${sourceFileList[@]}
    do
	g++ -Wall -pedantic -c ${sourceFile}.C `$ROOTSYS/bin/root-config --cflags --libs`
	if [ $? -ne 0 ]
	then
	    echo "Compiling ${sourceFile}.C : failed!"
	    return 0
	fi
	echo "Compiling ${sourceFile}.C : success :)"
	olist="${olist} ${sourceFile}.o"
    done
    

    echo "Now linking ..."

    g++ -Wall -pedantic -lm -o ${mainSourceFile} ${mainSourceFile}.cc ${olist} `$ROOTSYS/bin/root-config --cflags --libs`
    if [ $? -ne 0 ]
    then
	echo "Compiling ${mainSourceFile}.cc : failed!"
	return 0
    fi
    echo "Compiling ${mainSourceFile}.cc : success :)"
    olist="${olist} ${mainSourceFile}.o"

else
    echo "Options -nc passed: not compiling"
fi

if [ "$onlyCompile" = true ]; then
    # if option -c passed, we just want to complile and exit script              
    return 0
fi


echo "Now launching executable ..."

#EBorEE="EB"
# BarrelOrEndcap=""
# if [ $EBorEE = "EB" ]
# then
#     BarrelOrEndcap="Barrel"
# else
#     BarrelOrEndcap="Endcap"
# fi

for i in `seq $iter_ini $iter_fin` 
do
    if [ "$skipEB" = false ]; then
	echo  "iter_$i --> EB"
	wwwPath="${wwwBasePath}${dirName}/iter_${i}/2DMaps/Barrel/"
	./$mainSourceFile $path $dirName $i $tagName EB $wwwPath
    fi
    if [ "$skipEEp" = false ]; then
	echo  "iter_$i --> EE+"
	wwwPath="${wwwBasePath}${dirName}/iter_${i}/2DMaps/Endcap/EEp/"
	./$mainSourceFile $path $dirName $i $tagName EEp $wwwPath
    fi
    if [ "$skipEEm" = false ]; then
	echo  "iter_$i --> EE-"
	wwwPath="${wwwBasePath}${dirName}/iter_${i}/2DMaps/Endcap/EEm/"
	./$mainSourceFile $path $dirName $i $tagName EEm $wwwPath
    fi
done
