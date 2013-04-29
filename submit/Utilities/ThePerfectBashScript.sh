#!bin/bash
#Script by Luca Pernie'
#bash ThePerfectBashScript.sh help -> for help

# Function
usage (){
echo "--->The correct usage is: $>bash $0 param1. param2. ecc..."
echo "---->To have a lessons try:    man_arg: to know how to use the argument in bash."
echo "---->                          man_oper: to know how are the operator in bash."
echo "---->                          man_perl: to know how to use perl language and awk"
echo "---->To make something useful: iamroot: to know if you are a root user"
echo "---->                          changeFiles-directory-worlds1-worlds2: to chanhe all the 'worlds1' with 'worlds2' in all the file 'in directory'"
echo "---->                          removeEOS-group_Eos-dirName-TagFiles-iterMin-iterMax: to remove all Useless-files in group_Eos/dirName'"
echo "---->                          CopyEosFolder-group-DirFromToCopy-NewDirectory: make a cosy of a EOS folder with a different name."
echo "---->                          RemoveHaddfailed-group-Directory-TagFiles: Remove from Hadd list the corrupt files."
echo "---->                          ResendFillfailed-group-Directory-TagFiles-queue: Resubmit Fill jobs failed."
echo "---->                          EventInFile-filePath: Total number of events in list of file."
echo "---->If you are sad use:       IamSad"
echo "--------------"
}
Smile (){
echo "          --------             "
echo "         -        -            "
echo "        -   0  0   -           "
echo "       -            -    ______________________   "
echo "      -      ||      -  |  _      _  _  _   _  |  "
echo "      -              -  | |  |_| |_ |_ |_| |_  |  "
echo "       -    ____    -  <  |_ | | |_ |_ | \ __| |  "
echo "        -  (____)  -    |______________________|  "
echo "         -        -   "
echo "          --------    "
}
argum (){
echo 'How To Use Argument:'
echo '----> $* = list of all arg.'
echo '----> $0 = script name.'
echo '----> $1 = first arg.'
echo '----> $# = numb of arg.'
echo '----> ${!#} = last arg.'
echo "--------------"
}
operator (){
echo '----> -a = AND; -o = OR;'
echo '----> -gt = >; -ge = >=; -lt = <; -le = <=; -ne = !=; eq = ==;'
echo '----> -e = exist; -O = if you are the own; -f = if is a file; -d = if is a directory; -h = if is a link;'
echo '----> -s = if the dim. is>0; -z = empty string; -n = not empty string;'
echo '----> -r,w,x = if you have the permissions of reading, writing, execute;'
echo '---->  file1 -nt(-ot) file2 = if file1 is more(less) recent;'
echo "--------------"
}
man_perl (){
echo '----> PERL:'
echo '       Sobstitute: name =~ s/london/London/g; Opt: g: globally, i: case-insensitive'
echo '       Match:      name =~ m/sought_text/;'
echo '       NB: sed -e "s/london/London/g" in bash script. If /is part of the name, use \/ or #'
echo '       Special caracter:  *, ., ^ (beginning line), $ (end line), +, [aA]'
echo '       Special caracter: * (match 0 or more), + (match 1 or more_, ? (match 0 or 1), {n} (match n times), {n,} (match more than n), {n,m} (match n up to m times)'
echo '       Special caracter2: \w (alphanumeric plus "_"), \W (the rest), \s (whitespace), \S (the rest), \d (digit), \D, \t (tab), \n, \r (return), \e (escape)'
echo ''
echo '----> awk:'
echo '       es: ls cartella | awk '{print $5 " other things to add to the line"}' >> File.txt'
echo '       es: ls files_list | awk '{print "mv "$1" "$1".new"}' | sh'
echo "--------------"
}

iamroot (){
ROOT_UID=0
E_NONROOT=67 #codice di uscuta non root
if [ "$UID" -ne "$ROOT_UID" ]; then echo "You are a root user (UID = $UID).";
else echo "You are NOT a root user (UID = $UID)."; exit $E_NONROOT;
fi
echo "--------------"
}
SobstituteWords (){
E_SOB=69
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
dir=${array[1]}
s1=${array[2]}
s2=${array[3]}

if [ ${#array[@]} -ne "4" ]; then echo "Wrong Use of 'change'"; usage; exit $E_SOB;
fi
echo "I will sobstitute \"$s1\" with \"$s2\" in all \"$dir\" files... continue? [y-n]"
read var1
if [ "$var1" == "y" ]; then 
   dir=${dir}/*
   for file in $dir;
   do
      echo "Changing ${file}"
      sed -i "s/$s1/$s2/g" ${file}
   done
fi
echo "--------------"
}
removeEOS (){
E_EOS=70
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
group=${array[1]}
dir=${array[2]}
folder=${group}${dir}
tag=${array[3]}
iter_min=${array[4]}
iter_max=${array[5]}

if [ ${#array[@]} -ne "6" ]; then echo "Wrong Use of 'RemoveEOS'"; usage; exit $E_SOB;
fi
echo "I will REMOVE files in \"$folder\" from iter \"$iter_min\" to iter \"$iter_max\"... continue? [y-n]"
read var1
if [ "$var1" == "y" ]; then
   for Iter in $(eval echo "{$iter_min..$iter_max}")
   do
     for Barrel in {0..30}
     do
         echo "${folder}/iter_${Iter}/${tag}Barrel_${Barrel}_calibMap.root"
         cmsRm ${folder}/iter_${Iter}/${tag}Barrel_${Barrel}_calibMap.root
     done
     for Endcap in {0..7}
     do
         echo "${folder}/iter_${Iter}/${tag}Endcap_${Endcap}_calibMap.root"
         cmsRm ${folder}/iter_${Iter}/${tag}Endcap_${Endcap}_calibMap.root
     done
   
     for EcalNtp in {0..765}
     do
         echo "${folder}/iter_${Iter}/${tag}EcalNtp_${EcalNtp}.root"
         cmsRm ${folder}/iter_${Iter}/${tag}EcalNtp_${EcalNtp}.root
     done
     for Eps in {0..22}
     do
         echo "${folder}/iter_${Iter}/${tag}epsilonPlots_${Eps}.root"
         cmsRm ${folder}/iter_${Iter}/${tag}epsilonPlots_${Eps}.root
     done  
   done
fi
echo "--------------"
}
CopyEosFolder (){
E_CopyEosFold=71
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
group=${array[1]}
dir=${array[2]}
folder=${group}${dir}
newdir=${array[3]}
Newfolder=${group}${newdir}

if [ ${#array[@]} -ne "4" ]; then echo "Wrong Use of 'CopyEosFolder'"; usage; exit $E_CopyEosFold;
fi
eos mkdir $Newfolder
files=`cmsLs ${folder} | awk '{print $5}'`
echo $files

for file in ${files}
do
    echo "cmsStage -f $file /tmp"
    cmsStage -f $file /tmp
    echo "cmsStage -f /tmp/$file $Newfolder"
    cmsStage -f /tmp/$file $Newfolder
done
echo "--------------"
}

RemoveHaddfailed (){
echo "Removing corrupted files from hadd lists."
E_RemoveHadd=72
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
group=${array[1]}
dir=${array[2]}
folder=${group}${dir}
Tag=${array[3]}
mimSize=90000.
Here=`pwd`
if [ ${#array[@]} -ne "4" ]; then echo "Wrong Use of 'RemoveHaddfailed'"; usage; exit $E_RemoveHadd;
fi
List=`cmsLs -l ${folder} | awk '{print $2 "-" $5}' | grep "EcalNtp"`
read -a arrayList <<<$List

maxValue=$(echo "scale=0; ${#arrayList[@]}-1" | bc) #scale=numero cifre decimali
for file in $(eval echo "{0..${maxValue}}")
do
    Myfile=$(echo ${arrayList[$file]} | tr "-" " ")
    read -a arrayFinal <<<$Myfile
    if [ 1 -eq `echo "${arrayFinal[0]} < ${mimSize}" | bc` ]
    then
       BadList=$(echo ${arrayFinal[1]} | tr "/" " ")
       read -a arrayBad <<<$BadList
       #echo "grep ${arrayBad[0]}/${arrayBad[1]}/${arrayBad[2]}/${arrayBad[3]}//${arrayBad[4]}/${arrayBad[5]}/${arrayBad[6]} ${Here}/${arrayBad[4]}/src/hadd/*iter_0_*.list"
       Location=`grep ${arrayBad[0]}/${arrayBad[1]}/${arrayBad[2]}/${arrayBad[3]}//${arrayBad[4]}/${arrayBad[5]}/${arrayBad[6]} ${Here}/${arrayBad[4]}/src/hadd/*iter_0_*.list`
       if [ ! "${Location}" != "${p/[Rr]oot/}" ]; then continue
       fi
       Location=$(echo ${Location} | tr ":" " ")
       read -a arrayLocation <<<$Location
       #echo "cat ${arrayLocation[0]} | grep ${arrayBad[0]}/${arrayBad[1]}/${arrayBad[2]}/${arrayBad[3]}//${arrayBad[4]}/${arrayBad[5]}/${arrayBad[6]}"
       lineToRemove=`cat ${arrayLocation[0]} | grep ${arrayBad[0]}/${arrayBad[1]}/${arrayBad[2]}/${arrayBad[3]}//${arrayBad[4]}/${arrayBad[5]}/${arrayBad[6]}`
       echo "I will remove $lineToRemove in  ${arrayLocation[0]}"
       `cp ${arrayLocation[0]} ${arrayLocation[0]}.bkp`
       lineToRemove=$(echo ${lineToRemove} | tr "/" " ")
       read -a arrayRem <<<$lineToRemove
       sed "/${arrayRem[0]}\/\/${arrayRem[1]}\/\/${arrayRem[2]}\/${arrayRem[3]}\/${arrayRem[4]}\/${arrayRem[5]}\/${arrayRem[6]}\/${arrayRem[7]}\/\/${arrayRem[8]}\/${arrayRem[9]}\/${arrayRem[10]}/d" ${arrayLocation[0]} > ${arrayLocation[0]}.tmp
       mv ${arrayLocation[0]}.tmp ${arrayLocation[0]}
    fi
done
echo "--------------"
}

ResendFillfailed (){
E_Resend=73
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
group=${array[1]}
dir=${array[2]}
folder=${group}${dir}
Tag=${array[3]}
queue=${array[4]}
mimSize=90000.
Here=`pwd`
if [ ${#array[@]} -ne "5" ]; then echo "Wrong Use of 'ResendFillfailed'"; usage; exit $E_Resend;
fi
List=`cmsLs -l ${folder} | awk '{print $2 "-" $5}' | grep "EcalNtp"`
read -a arrayList <<<$List

maxValue=$(echo "scale=0; ${#arrayList[@]}-1" | bc) #scale=numero cifre decimali
for file in $(eval echo "{0..${maxValue}}")
do
    Myfile=$(echo ${arrayList[$file]} | tr "-" " ")
    read -a arrayFinal <<<$Myfile
    if [ 1 -eq `echo "${arrayFinal[0]} < ${mimSize}" | bc` ]
    then
       BadList=$(echo ${arrayFinal[1]} | tr "/" " ")
       read -a arrayBad <<<$BadList
       Location=`grep ${arrayBad[0]}/${arrayBad[1]}/${arrayBad[2]}/${arrayBad[3]}//${arrayBad[4]}/${arrayBad[5]}/${arrayBad[6]} ${Here}/${arrayBad[4]}/src/Fill/*${arrayBad[5]}*.sh`
       Location=$(echo ${Location} | tr ":" " ")
       read -a arrayLocation <<<$Location
       echo "Is in: ${arrayLocation[0]}"
       echo "bsub -q ${queue} -o ${Here}/${arrayBad[4]}/log/fillEpsilonPlot_${arrayBad[5]}_job_Resub.log ${arrayLocation[0]}"
       `bsub -q ${queue} -o ${Here}/${arrayBad[4]}/log/fillEpsilonPlot_${arrayBad[5]}_job_Resub.log ${arrayLocation[0]}`
    fi
done
echo "--------------"
}

EventInFile (){
E_EventInFile=74
param=$(echo $1 | tr "-" " ")
read -a array <<<$param
list=${array[1]}
if [ ${#array[@]} -ne "2" ]; then echo "Wrong Use of 'EventInFile'"; usage; exit $E_EventInFile;
fi
List_array=`cat $list`
Tot_Ev=0
for Line in ${List_array}
do
   Events=`edmFileUtil $Line | grep events | awk '{print $6}'`
   echo $Events
   Tot_Ev=$(echo "scale=0; $Tot_Ev+$Events" | bc)
done
echo "The total numbers of events is: $Tot_Ev"
echo "--------------"
}

#-----------------------------------------------------------------------------------
echo "Hello, I'm \"$0\": The first bash-script who teach you how to make a script in bash!"
echo "->Use 'bash ThePerfectBashScript.sh help' to see what I can do for you"
echo "*******----START----*******"

# Exit codes, ID codes
nArg_min=1
E_nArg=68
if [ "$#" -lt "$nArg_min" ]; then exit $E_nArg;
fi

# Info.
curr_dir=`pwd`
nomeutente=`whoami`
today=`date`

# Reading Parameters
for p in $*;
  do
  if [ "$p" == "help" ]; then usage;
  elif [ "$p" != "${p/[Mm]an/}" -a "$p" != "${p/[Aa]rg/}" ]; then argum;
  elif [ "$p" != "${p/[Mm]an/}" -a "$p" != "${p/[Oo]per/}" ]; then operator;
  elif [ "$p" != "${p/[Aa]m/}" -a "$p" != "${p/[Rr]oot/}" ]; then iamroot;
  elif [ "$p" != "${p/[Mm]an/}" -a "$p" != "${p/[Pp]erl/}" ]; then man_perl;
  elif [ "$p" != "${p/[Cc]hangeFiles/}" ]; then SobstituteWords $p;
  elif [ "$p" != "${p/[Rr]emove/}" -a "$p" != "${p/[Ee][Oo][Ss]/}" ]; then removeEOS $p;
  elif [ "$p" == "IamSad" ]; then Smile;
  elif [ "$p" != "${p/[Cc]opy/}" -a "$p" != "${p/[Ee][Oo][Ss]/}" -a "$p" != "${p/[Ff]older/}" ]; then CopyEosFolder $p;
  elif [ "$p" != "${p/[Rr]/emove}" -a "$p" != "${p/[Hh]add/}" -a "$p" != "${p/[Ff]ailed/}" ]; then RemoveHaddfailed $p;
  elif [ "$p" != "${p/[Rr]/esend}" -a "$p" != "${p/[Ff]ill/}" -a "$p" != "${p/[Ff]ailed/}" ]; then ResendFillfailed $p;
  elif [ "$p" != "${p/[Ee]/vent}" -a "$p" != "${p/[Ii][Nn]/}" -a "$p" != "${p/[Ff]ile/}" ]; then EventInFile $p;
  else echo "Sorry but $p is not a valid option... use help to see the option I have"
  fi
done

echo "*******-----END-----*******"
exit
