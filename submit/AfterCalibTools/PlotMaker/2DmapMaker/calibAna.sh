#!/bin/bash   

# launch as --> source calibAna.sh [OPTIONS]
# if you just want to compile and end script execution, use option -c

##############################################
# THINGS TO CHANGE BEFORE LAUNCHING SCRIPT
##############################################


iter_ini=6                                                     # first iteration to use
iter_fin=6                                                      # last iteration to use: it is included in sequence below
#path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/emanuele/"  # path to directory on eos
#path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/"  # path to directory on eos
#path="/store/group/dpg_ecal/alca_ecalcalib/piZero2016/zhicaiz/"  # path to directory on eos
path="/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"  # path to directory on eos
#dirName="AlcaP0_Run2016G_sel16_reg12"                            # dirname (see CalibCode/submit/parameters.py)
dirName="AlCaP0_Run2017_C_CCiter0"
#tagName="AlCaP0_Run2017B_3July_upToRun297723_ext1_"                           # TagName (see CalibCode/submit/parameters.py) 
tagName="${dirName}_"

Pi0orEta="Pi0"  # possible options are Pi0 and Eta, axis ranges are set a little differently

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
    elif [ $option = "--nc" ]; then
        doNotCompile=true
    elif [ $option = "--noEB" ]; then
	skipEB=true
    elif [ $option = "--noEEp" ]; then
	skipEEp=true
    elif [ $option = "--noEEm" ]; then
	skipEEm=true
    elif [ $option = "--pi0" ]; then
	Pi0orEta="Pi0"
    elif [ $option = "--eta" ]; then
	Pi0orEta="Eta"
    fi
done

########################################

olist=""
sourceFile=""

if [ "$doNotCompile" = false ]; then

    echo
    echo "Compiling ..."

# to compile using root libraries I created an alias
# alias rootlib='$ROOTSYS/bin/root-config --cflags --libs'
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
	mkdir -p ${wwwPath}
	cp /afs/cern.ch/user/m/mciprian/public/index.php {wwwPath}
	./$mainSourceFile $path $dirName $i $tagName EB $wwwPath $Pi0orEta
    fi
    if [ "$skipEEp" = false ]; then
	echo  "iter_$i --> EE+"
	wwwPath="${wwwBasePath}${dirName}/iter_${i}/2DMaps/Endcap/EEp/"
	mkdir -p ${wwwPath}
	cp /afs/cern.ch/user/m/mciprian/public/index.php {wwwPath}
	./$mainSourceFile $path $dirName $i $tagName EEp $wwwPath $Pi0orEta
    fi
    if [ "$skipEEm" = false ]; then
	echo  "iter_$i --> EE-"
	wwwPath="${wwwBasePath}${dirName}/iter_${i}/2DMaps/Endcap/EEm/"
	mkdir -p ${wwwPath}
	cp /afs/cern.ch/user/m/mciprian/public/index.php {wwwPath}
	./$mainSourceFile $path $dirName $i $tagName EEm $wwwPath $Pi0orEta
    fi
done
