#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

from optparse import OptionParser                                                                                                                   
                                                                                          
parser = OptionParser(usage="%prog [options]")    
parser.add_option("-c", "--create",           dest="create", action="store_true", default=False, help="Do not submit the jobs, only create the subfolders")
parser.add_option("-l", "--daemon-local",     dest="daemonLocal", action="store_true", default=False, help="Do not submit a job to manage the daemon, do it locally")
(options, args) = parser.parse_args()
pwd = os.getcwd()

if not options.create:
    #-------- check if you have right access to queues --------#
    checkAccessToQueues = subprocess.Popen(['bjobs'], stderr=subprocess.PIPE, shell=True);
    output = checkAccessToQueues.communicate()[1]
    if(output.find('command not found')==-1):
        print "[calib] Correct setup for batch submission"
    else:
        print "[calib] Missing access to queues"
        if not( isCRAB and storageSite=="T2_BE_IIHE" ):
            sys.exit(1)
#-------- copy regression files -------#
if ContainmentCorrection == '2017reg':
	os.system("cp /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/zhicaiz/GBRForest_2017/* ../FillEpsilonPlot/data/")

#-------- create folders --------#

workdir = pwd+'/'+dirname
condordir = workdir + "/condor_files/"
cfgFillPath = workdir + '/cfgFile/Fill'
cfgFitPath  = workdir + '/cfgFile/Fit'
cfgHaddPath  = workdir + '/src/hadd'
srcPath  = workdir + '/src'

print "[calib] Creating local folders (" + dirname + ")"
folderCreation = subprocess.Popen(['mkdir -p ' + workdir], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + condordir], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
for it in range(nIterations):
    folderCreation = subprocess.Popen(['mkdir -p ' + condordir + "/iter_" + str(it)], stdout=subprocess.PIPE, shell=True);
    folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/cfgFile/'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/CRAB_files/'], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
for it in range(nIterations):
    folderCreation = subprocess.Popen(['mkdir -p ' + cfgFillPath + '/iter_' + str(it)], stdout=subprocess.PIPE, shell=True);
    folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + cfgFitPath], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
# created in calibJobHandlerCondor.py to follow the flow from local folders
#for it in range(nIterations):    
#    for dirtype in ["Fill", "Hadd", "Final_Hadd", "Fit"]:
#        folderCreation = subprocess.Popen(['mkdir -p ' + workdir + '/log/' + dirtype + '/iter_' + str(it)], stdout=subprocess.PIPE, shell=True);
#        folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
for it in range(nIterations):
    folderCreation = subprocess.Popen(['mkdir -p ' + srcPath + '/Fill/iter_' + str(it)], stdout=subprocess.PIPE, shell=True);
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
   folderCreation = subprocess.Popen(['mkdir -p ' + eosPath + '/' + dirname ], stdout=subprocess.PIPE, shell=True);
   folderCreation.communicate()

for it in range(nIterations):
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
       print "[calib]  ---  srmmkdir " + eosPath + '/' + dirname + '/iter_' + str(it)
       folderCreation = subprocess.Popen(['srmmkdir srm://maite.iihe.ac.be:8443' + eosPath + '/' + dirname + '/iter_' + str(it)], stdout=subprocess.PIPE, shell=True);
       folderCreation.communicate()
    else:
       print "[calib]  ---  mkdir " + eosPath + '/' + dirname + '/iter_' + str(it)
       folderCreation = subprocess.Popen(['mkdir ' + eosPath + '/' + dirname + '/iter_' + str(it)], stdout=subprocess.PIPE, shell=True);
       folderCreation.communicate()

#-------- fill cfg files --------#
if( isCRAB ):
    print "--------------This inter-calibration will use CRAB: Good Luck!------------------------"

# open list of input files
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = [x for x in inputlist_f.readlines() if not x.lstrip().startswith('#')]  # do not consider commented line

print "[calib] Total number of files to be processed: " , len(inputlistbase_v)
print "[calib] Creating cfg Files"

ijob = 0
for it in range(nIterations):
    print "[calib]  '-- Fill::Iteration " + str(it)
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

    haddSrc_final_n_s = srcPath + "/hadd/hadd_iter_" + str(it) + "_final.list"
    haddSrc_final_f_s = open(  haddSrc_final_n_s, 'w')
    for num_list in range(Nlist):
        haddSrc_n_s.append( srcPath + "/hadd/hadd_iter_" + str(it) + "_step_" + str(num_list)+ ".list")
        haddSrc_f_s.append( open(  haddSrc_n_s[num_list], 'w') )
        fileToAdd_final_n_s = eosPath + '/' + dirname + '/iter_' + str(it) + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
        for nj in range(nHadd):
            nEff = num_list*nHadd+nj
            fileToAdd_n_s = eosPath + '/' + dirname + '/iter_' + str(it) + '/' + NameTag + outputFile + '_' + str(nEff) + '.root\n'
            if(nEff < NrelJob) :
                haddSrc_f_s[num_list].write(fileToAdd_n_s)
        haddSrc_final_f_s.write(fileToAdd_final_n_s)
        haddSrc_f_s[num_list].close()
    haddSrc_final_f_s.close()

    # create Hadd cfg file
    dest = eosPath + '/' + dirname + '/iter_' + str(it) + '/'
    for num_list in range(Nlist):
        hadd_cfg_n = cfgHaddPath + "/HaddCfg_iter_" + str(it) + "_job_" + str(num_list) + ".sh"
        hadd_cfg_f = open( hadd_cfg_n, 'w' )
        HaddOutput = NameTag + "epsilonPlots_" + str(num_list) + ".root"
        printParallelHaddFAST(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd, num_list )
        hadd_cfg_f.close()
        changePermission = subprocess.Popen(['chmod 777 ' + hadd_cfg_n], stdout=subprocess.PIPE, shell=True);
        debugout = changePermission.communicate()
    # print Final hadd
    Fhadd_cfg_n = cfgHaddPath + "/Final_HaddCfg_iter_" + str(it) + ".sh"
    Fhadd_cfg_f = open( Fhadd_cfg_n, 'w' )
    printFinalHaddRegroup(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd )
    Fhadd_cfg_f.close()
    # loop over the whole list
    while (len(inputlist_v) > 0):

        # create cfg file
        fill_cfg_n = cfgFillPath + "/iter_" + str(it) + "/fillEps_iter_" + str(it) + "_job_" + str(ijob) + ".py"
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
                prefixSourceFileToUse = ""
                if prefixSourceFile not in ntpfile:
                    prefixSourceFileToUse = prefixSourceFile
                if(line != lastline):
                    fill_cfg_f.write("        '" + prefixSourceFileToUse + ntpfile + "',\n")
                else:
                    fill_cfg_f.write("        '" + prefixSourceFileToUse + ntpfile + "'\n")

        # print the last part of the cfg file
        if( isCRAB ):
            printFillCfg2( fill_cfg_f, pwd, it , "", ijob )
        else: 
            printFillCfg2( fill_cfg_f, pwd, it , "/tmp/", ijob )
        fill_cfg_f.close()

        # print source file for batch submission of FillEpsilonPlot task
        fillSrc_n = srcPath + "/Fill/iter_" + str(it) + "/submit_iter_" + str(it) + "_job_" + str(ijob) + ".sh"
        fillSrc_f = open( fillSrc_n, 'w')
        source_s = NameTag +outputFile + "_" + str(ijob) + ".root"
        destination_s = eosPath + '/' + dirname + '/iter_' + str(it) + "/" + source_s
        logpathFill = pwd + "/" + dirname + "/log/iter_" + str(it) + "/fillEps_iter_" + str(it) + "_job_" + str(ijob) + ".log"
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
for it in range(nIterations):
    print "[calib]  '-- Fit::Iteration " + str(it)
    for nFit in range(nEB):
        # create cfg file
        fit_cfg_n = cfgFitPath + "/fitEpsilonPlot_EB_" + str(nFit) + "_iter_" + str(it) + ".py"
        fit_cfg_f = open( fit_cfg_n, 'w' )

        # print the cfg file
        printFitCfg( fit_cfg_f , it, "/tmp",inListB[nFit],finListB[nFit],"Barrel",nFit)
        fit_cfg_f.close()

        # print source file for batch submission of FitEpsilonPlot task
        fitSrc_n = srcPath + "/Fit/submit_EB_" + str(nFit) + "_iter_" + str(it) + ".sh"
        fitSrc_f = open( fitSrc_n, 'w')
        destination_s = eosPath + '/' + dirname + '/iter_' + str(it) + "/" + NameTag + "Barrel_" + str(nFit)+ "_" + calibMapName
        logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EB_" + str(nFit) + "_iter_" + str(it) + ".log"
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
        fit_cfg_n = cfgFitPath + "/fitEpsilonPlot_EE_" + str(nFit) + "_iter_" + str(it) + ".py"
        fit_cfg_f = open( fit_cfg_n, 'w' )

        # print the cfg file
        printFitCfg( fit_cfg_f , it, "/tmp",inListE[nFit],finListE[nFit],"Endcap",nFit)

        fit_cfg_f.close()

        # print source file for batch submission of FitEpsilonPlot task
        fitSrc_n = srcPath + "/Fit/submit_EE_" + str(nFit) + "_iter_" + str(it) + ".sh"
        fitSrc_f = open( fitSrc_n, 'w')
        destination_s = eosPath + '/' + dirname + '/iter_' + str(it) + "/" + NameTag + "Endcap_" + str(nFit) + "_" + calibMapName
        logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EE_" + str(nFit) + "_iter_" + str(it) + ".log"
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
env_script_f.write("eval `scramv1 runtime -sh`\n")
env_script_f.write( "python " + pwd + "/calibJobHandlerCondor.py " + str(njobs) + " " + queue + "\n")
env_script_f.write( "rm -rf " + pwd + "/core.*\n")
env_script_f.close()

# make the source file executable
changePermission = subprocess.Popen(['chmod 777 ' + env_script_n], stdout=subprocess.PIPE, shell=True);
debugout = changePermission.communicate()

dummy_exec = open(condordir+'/dummy_exec_daemon.sh','w')
dummy_exec.write('#!/bin/bash\n')
dummy_exec.write('bash $*\n')
dummy_exec.close()

condor_file_name = condordir+'/condor_submit_daemon.condor'
condor_file = open(condor_file_name,'w')
# line 'next_job_start_delay = 1' not needed here
condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {ld}/$(ProcId).log
Output     = {ld}/$(ProcId).out
Error      = {ld}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
request_memory = 4000
+MaxRuntime = 604800
+JobBatchName = "ecalpro_daemon"\n
'''.format(de=os.path.abspath(dummy_exec.name), ld=os.path.abspath(condordir), here=os.environ['PWD'] ) )

condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(env_script_n)))
condor_file.close()


# configuring calibration handler

submit_s = 'condor_submit {cdf} '.format(cdf=condor_file_name)
if not options.create:
    print "[calib] Number of jobs created = " + str(njobs)
    print "[calib] Submitting calibration handler"
    # submitting calibration handler
    if options.daemonLocal:
        print "[calib]  '-- source " + os.path.abspath(env_script_n)
        os.system("source " + os.path.abspath(env_script_n))
    else:        
        print "[calib]  '-- " + submit_s
        submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
        output = (submitJobs.communicate()[0]).splitlines()
        print "[calib]  '-- " + output[0]

    #    print "usage thisPyton.py pwd njobs queue"
else:
    print "options -c was given: jobs are not submitted, but all folders and files were created normally. You can still do local tests."
    print "To run the whole code use the following command."
    print submit_s


