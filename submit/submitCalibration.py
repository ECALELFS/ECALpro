#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ): # Beacause in IIHE the pwd give a link to the area, and you don't want that
    pwd         = os.getenv('PWD')
else:
    pwd         = os.getcwd()

#-------- check if you have right access to queues --------#
checkAccessToQueues = subprocess.Popen(['bjobs'], stderr=subprocess.PIPE, shell=True);
output = checkAccessToQueues.communicate()[1]
if(output.find('command not found')==-1):
    print "[calib] Correct setup for batch submission"
else:
    print "[calib] Missing access to queues"
    if not( isCRAB and storageSite=="T2_BE_IIHE" ):
       sys.exit(1)

#-------- create folders --------#

workdir = pwd+'/'+dirname
cfgFillPath = workdir + '/cfgFile/Fill'
cfgFitPath  = workdir + '/cfgFile/Fit'
cfgHaddPath  = workdir + '/src/hadd'
srcPath  = workdir + '/src'

print "[calib] Creating local folders (" + dirname + ")"
folderCreation = subprocess.Popen(['mkdir -p ' + workdir], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/cfgFile/'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/CRAB_files/'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + cfgFillPath], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + cfgFitPath], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/log/'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath + '/Fill'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath + '/Fit'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath + '/hadd'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()

print "[calib] Storing parameter.py for future reference"
CopyParam = subprocess.Popen(['cp  parameters.py ' + workdir], stdout=subprocess.PIPE, shell=True);
CopyParam.communicate()

if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
   print "[calib] Creating folders on PNFS"
   folderCreation = subprocess.Popen(['srmmkdir srm://maite.iihe.ac.be:8443' + eosPath + '/' + dirname ], stdout=subprocess.PIPE, shell=True);
   folderCreation.communicate()
else:
   print "[calib] Creating folders on EOS"
   folderCreation = subprocess.Popen(['cmsMkdir ' + eosPath + '/' + dirname ], stdout=subprocess.PIPE, shell=True);
   folderCreation.communicate()

for iter in range(nIterations):
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       print "[calib]  ---  srmmkdir " + eosPath + '/' + dirname + '/iter_' + str(iter)
       folderCreation = subprocess.Popen(['srmmkdir srm://maite.iihe.ac.be:8443' + eosPath + '/' + dirname + '/iter_' + str(iter)], stdout=subprocess.PIPE, shell=True);
       folderCreation.communicate()
    else:
       print "[calib]  ---  cmsMkdir " + eosPath + '/' + dirname + '/iter_' + str(iter)
       folderCreation = subprocess.Popen(['cmsMkdir ' + eosPath + '/' + dirname + '/iter_' + str(iter)], stdout=subprocess.PIPE, shell=True);
       folderCreation.communicate()

#-------- fill cfg files --------#
if( isCRAB ):
    print "--------------This inter-calibration will use CRAB: Good Luck!------------------------"

# open list of input files
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = inputlist_f.readlines()

print "[calib] Total number of files to be processed: " , len(inputlistbase_v)
print "[calib] Creating cfg Files"

for iter in range(nIterations):
    print "[calib]  '-- Fill::Iteration " + str(iter)
    # copy by value and not by reference
    inputlist_v = inputlistbase_v[:]
    ijob=0

    # Creating different list for hadd
    NrelJob = float(len(inputlist_v)) / float(ijobmax)
    if( float(int(NrelJob) - NrelJob) < 0. ):
        NrelJob = int(NrelJob) + 1
    Nlist_flo = float(NrelJob/nHadd) + 1.
    Nlist = int(Nlist_flo)

    haddSrc_n_s = list()
    haddSrc_f_s = list()

    print "[calib]  '-- Hadd::Number of hadd tasks: " + str(Nlist) + "  (" + str(nHadd) + " files per task)"

    haddSrc_final_n_s = srcPath + "/hadd/hadd_iter_" + str(iter) + "_final.list"
    haddSrc_final_f_s = open(  haddSrc_final_n_s, 'w')
    for num_list in range(Nlist):
        haddSrc_n_s.append( srcPath + "/hadd/hadd_iter_" + str(iter) + "_step_" + str(num_list)+ ".list")
        haddSrc_f_s.append( open(  haddSrc_n_s[num_list], 'w') )
        if(fastHadd):
            fileToAdd_final_n_s = eosPath + '/' + dirname + '/iter_' + str(iter) + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
        else:
            fileToAdd_final_n_s = 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iter) + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
        for nj in range(nHadd):
            nEff = num_list*nHadd+nj
            if(fastHadd):
                fileToAdd_n_s = eosPath + '/' + dirname + '/iter_' + str(iter) + '/' + NameTag + outputFile + '_' + str(nEff) + '.root\n'
            else:
                fileToAdd_n_s = 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iter) + '/' + NameTag + outputFile + '_' + str(nEff) + '.root\n'
            if(nEff < NrelJob) :
                haddSrc_f_s[num_list].write(fileToAdd_n_s)
        haddSrc_final_f_s.write(fileToAdd_final_n_s)
        haddSrc_f_s[num_list].close()
    haddSrc_final_f_s.close()

    # create Hadd cfg file
    dest = eosPath + '/' + dirname + '/iter_' + str(iter) + '/'
    for num_list in range(Nlist):
        hadd_cfg_n = cfgHaddPath + "/HaddCfg_iter_" + str(iter) + "_job_" + str(num_list) + ".sh"
        hadd_cfg_f = open( hadd_cfg_n, 'w' )
        HaddOutput = NameTag + "epsilonPlots_" + str(num_list) + ".root"
        if(fastHadd):
            printParallelHaddFAST(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd, num_list )
        else:
            printParallelHadd(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd )
        hadd_cfg_f.close()
        changePermission = subprocess.Popen(['chmod 777 ' + hadd_cfg_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()
    # print Final hadd
    Fhadd_cfg_n = cfgHaddPath + "/Final_HaddCfg_iter_" + str(iter) + ".sh"
    Fhadd_cfg_f = open( Fhadd_cfg_n, 'w' )
    if(fastHadd):
        printFinalHaddFAST(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd )
    else:
        printFinalHadd(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd )
    Fhadd_cfg_f.close()
    # loop over the whole list
    while (len(inputlist_v) > 0):

        # create cfg file
        fill_cfg_n = cfgFillPath + "/fillEpsilonPlot_iter_" + str(iter) + "_job_" + str(ijob) + ".py"
        #print "writing " + fill_cfg_n + " ..."
        fill_cfg_f = open( fill_cfg_n, 'w' )

        # print first part of the cfg file
        printFillCfg1( fill_cfg_f )
        # loop over the names of the input files to be put in a single cfg
        lastline = min(ijobmax,len(inputlist_v)) - 1
        for line in range(min(ijobmax,len(inputlist_v))):
            ntpfile = inputlist_v.pop(0)
            ntpfile = ntpfile.rstrip()
            if ntpfile != '':
                if(line != lastline):
                    fill_cfg_f.write("        '" + ntpfile + "',\n")
                else:
                    fill_cfg_f.write("        '" + ntpfile + "'\n")

        # print the last part of the cfg file
        if( isCRAB ):
            printFillCfg2( fill_cfg_f, pwd, iter , "", ijob )
        else: 
            printFillCfg2( fill_cfg_f, pwd, iter , "/tmp/", ijob )
        fill_cfg_f.close()

        # print source file for batch submission of FillEpsilonPlot task
        fillSrc_n = srcPath + "/Fill/submit_iter_" + str(iter) + "_job_" + str(ijob) + ".sh"
        fillSrc_f = open( fillSrc_n, 'w')
        source_s = NameTag +outputFile + "_" + str(ijob) + ".root"
        destination_s = eosPath + '/' + dirname + '/iter_' + str(iter) + "/" + source_s
        logpathFill = pwd + "/" + dirname + "/log/" + "fillEpsilonPlot_iter_" + str(iter) + "_job_" + str(ijob) + ".log"
        printSubmitSrc(fillSrc_f, fill_cfg_n, "/tmp/" + source_s, destination_s , pwd, logpathFill)
        fillSrc_f.close()

        # make the source file executable
        changePermission = subprocess.Popen(['chmod 777 ' + fillSrc_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()

        ijob = ijob+1

njobs = ijob

#-------- fit cfg files --------#
    # Fit parallelized
nEB = 61199/nFit
if (61199%nFit != 0) :
    nEB = int(nEB) +1
nEE = 14647/nFit
if (14647%nFit != 0) :
    nEE = int(nEE) +1

print '[calib] Splitting Fit Task: ' + str(nEB) + ' jobs on EB, ' + str(nEE) + ' jobs on EE'
#print 'I will submit ' + str(nEB) + ' jobs to fit the Barrel'
#print 'I will submit ' + str(nEE) + ' jobs to fit the Endcap'
print '[calib] Creating Fit cfg files'
inListB = list()
finListB = list()
inListE = list()
finListE = list()
for tmp in range(nEB):
    inListB.append( 2000*tmp )
    finListB.append( 2000*tmp+1999 )
for tmp in range(nEE):
    inListE.append( 2000*tmp )
    finListE.append( 2000*tmp+1999 )
    # cfg
for iter in range(nIterations):
    print "[calib]  '-- Fit::Iteration " + str(iter)
    for nFit in range(nEB):
        # create cfg file
        fit_cfg_n = cfgFitPath + "/fitEpsilonPlot_EB_" + str(nFit) + "_iter_" + str(iter) + ".py"
        fit_cfg_f = open( fit_cfg_n, 'w' )

        # print the cfg file
        printFitCfg( fit_cfg_f , iter, "/tmp",inListB[nFit],finListB[nFit],"Barrel",nFit)
        fit_cfg_f.close()

        # print source file for batch submission of FitEpsilonPlot task
        fitSrc_n = srcPath + "/Fit/submit_EB_" + str(nFit) + "_iter_" + str(iter) + ".sh"
        fitSrc_f = open( fitSrc_n, 'w')
        destination_s = eosPath + '/' + dirname + '/iter_' + str(iter) + "/" + NameTag + "Barrel_" + str(nFit)+ "_" + calibMapName
        logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EB_" + str(nFit) + "_iter_" + str(iter) + ".log"
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            printSubmitFitSrc(fitSrc_f, fit_cfg_n, "$TMPDIR/" + NameTag + "Barrel_" + str(nFit) + "_" + calibMapName, destination_s, pwd, logpath)
        else:
            printSubmitFitSrc(fitSrc_f, fit_cfg_n, "/tmp/" + NameTag + "Barrel_" + str(nFit) + "_" + calibMapName, destination_s, pwd, logpath)
        fitSrc_f.close()

        # make the source file executable
        changePermission = subprocess.Popen(['chmod 777 ' + fitSrc_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()

    for nFit in range(nEE):
        # create cfg file
        fit_cfg_n = cfgFitPath + "/fitEpsilonPlot_EE_" + str(nFit) + "_iter_" + str(iter) + ".py"
        fit_cfg_f = open( fit_cfg_n, 'w' )

        # print the cfg file
        printFitCfg( fit_cfg_f , iter, "/tmp",inListE[nFit],finListE[nFit],"Endcap",nFit)

        fit_cfg_f.close()

        # print source file for batch submission of FitEpsilonPlot task
        fitSrc_n = srcPath + "/Fit/submit_EE_" + str(nFit) + "_iter_" + str(iter) + ".sh"
        fitSrc_f = open( fitSrc_n, 'w')
        destination_s = eosPath + '/' + dirname + '/iter_' + str(iter) + "/" + NameTag + "Endcap_" + str(nFit) + "_" + calibMapName
        logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EE_" + str(nFit) + "_iter_" + str(iter) + ".log"
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            printSubmitFitSrc(fitSrc_f, fit_cfg_n, "$TMPDIR/" + NameTag + "Endcap_" + str(nFit)+ "_" + calibMapName, destination_s, pwd, logpath)
        else:
            printSubmitFitSrc(fitSrc_f, fit_cfg_n, "/tmp/" + NameTag + "Endcap_" + str(nFit)+ "_" + calibMapName, destination_s, pwd, logpath)
        fitSrc_f.close()

        # make the source file executable
        changePermission = subprocess.Popen(['chmod 777 ' + fitSrc_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()

### setting environment
env_script_n = workdir + "/submit.sh"
env_script_f = open(env_script_n, 'w')
env_script_f.write("#!/bin/bash\n")
env_script_f.write("cd " + pwd + "\n")
env_script_f.write("ulimit -c 0\n")
#if(is2012):
#   env_script_f.write("export SCRAM_ARCH=slc5_amd64_gcc462\n")
#else:
#   env_script_f.write("export SCRAM_ARCH=slc5_amd64_gcc434\n")

env_script_f.write("eval `scramv1 runtime -sh`\n")
env_script_f.write( "python " + pwd + "/calibJobHandler.py " + str(njobs) + " " + queue + "\n")
env_script_f.write( "rm -rf " + pwd + "/core.*")
env_script_f.close()

# make the source file executable
changePermission = subprocess.Popen(['chmod 777 ' + env_script_n], stdout=subprocess.PIPE, shell=True);
debugout = changePermission.communicate()

if( isCRAB ):
    for iter in range(nIterations):
        CRAB1_n =  workdir + "/CRAB_files/crab_eos_" + str(iter) + ".py"
        CRAB1_f = open( CRAB1_n, 'w' )
        printCrab( CRAB1_f, iter)
        CopyFill = subprocess.Popen(['cp ' +cfgFillPath + '/fillEpsilonPlot_iter_' + str(iter) + '_job_0.py ' + workdir + '/CRAB_files/fillEpsilonPlot_iter_' + str(iter) + '.py' ], stdout=subprocess.PIPE, shell=True);
        CopyFill.communicate()
        CrabSendHadd_n =  workdir + "/CRAB_files/HaddSendafterCrab_" + str(iter) + ".sh"
        CrabSendHadd_f = open( CrabSendHadd_n, 'w' )
        printCrabHadd( CrabSendHadd_f, str(iter), pwd)
    #Instructions
    print "---------------------------------"
    print "Here it is how it works with CRAB:"
    print "--> 1) You will run the crab_eos_0.cfg I wrote for you in: " + workdir + "/CRAB_files/crab_eos.cfg: \n  --->crab submit -c crab_eos_X.py"
    print "--> 2) When all the outputs are on EOS you will launch the second part of the script to do the HADD and the FIT with the command:\n  --->bsub -q " + queueForDaemon + " 'bash " + workdir + "/CRAB_files/HaddSendafterCrab_XXX.sh'"
    print "--> 3) Once it has finished you will re-run CRAB importing the constant you produced" #!!! this part is not clear.
    print "--> 4) Then you repeat all these steps for all the iterations you need. Good luck."
    # in the futur launch a script that send crab automatically

else:
    # configuring calibration handler
    print "[calib] Number of jobs created = " + str(njobs)
    print "[calib] Submitting calibration handler"
    submit_s = 'bsub -q ' + queueForDaemon + ' -o ' + workdir + '/calibration.log "source ' + env_script_n + '"'
    print "[calib]  '-- " + submit_s
    
    # submitting calibration handler
    submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
    output = (submitJobs.communicate()[0]).splitlines()
    print "[calib]  '-- " + output[0]
    
    #    print "usage thisPyton.py pwd njobs queue"
