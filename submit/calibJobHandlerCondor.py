#!/usr/bin/env python

import subprocess, time, sys, os
from ROOT import *
from methods import *
from datetime import datetime

from optparse import OptionParser

# def checkNjobsCondor(noDaemon=True):
#     checkJobs = subprocess.Popen(['condor_q'], stdout=subprocess.PIPE, shell=True);
#     datalines = (checkJobs.communicate()[0]).splitlines()
#     nRetjobs = 0
#     for l in datalines:
#         if all(x in l for x in ["jobs", "completed", "removed"]):
#             #print l
#             nRetjobs = int(l.split("jobs")[0].strip()) # note, one is the daemon, so subtract 1
#             if noDaemon: nRetjobs = nRetjobs - 1
#     return nRetjobs

def renewTokenAFS(daemonLocal=False, infile=""):
    if daemonLocal and len(infile):
        print "Renewing AFS token"
        os.system("cat {infile} | kinit".format(infile=infile))
        os.system("kinit -R")

def checkNjobsCondor(grepArg="ecalpro"):

    # condor_q might fail from time to time, so check whether it is working by asking condor_q first
    # this should return some lines like:
    # -- Schedd: pccmsrm31.dyndns.cern.ch : <128.141.87.213:9618?... @ 02/05/19 15:06:41
    # OWNER BATCH_NAME      SUBMITTED   DONE   RUN    IDLE   HOLD  TOTAL JOB_IDS
    #
    # 0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended

    checkJobs = subprocess.Popen(['condor_q'], stdout=subprocess.PIPE, shell=True);
    tmp = str(checkJobs.communicate()[0])
    if all(x in tmp for x in ["OWNER", "BATCH_NAME", "SUBMITTED", "jobs", "completed", "removed"]):   # a very dumb check, I know
        checkJobs = subprocess.Popen(['condor_q -af JobBatchName| grep {gr} | wc -l'.format(gr=grepArg)], stdout=subprocess.PIPE, shell=True);
        nRetjobs = checkJobs.communicate()[0]
        # check nRetjobs is not a null string, but a number (either 0 or something else)
        if len(nRetjobs) and nRetjobs != "\n" and nRetjobs.replace('\n','').isdigit():
            return int(nRetjobs)
        else:
            print "Error 1: output = " + str(nRetjobs)
            return 9999
    else:
        # something wrong with condor, return a dummy positive number to tell the code to keep checking
        print "Error 2: output = " + tmp
        return 9999 


# helper function to save some lines, the file is not opened not closed here, this must be handled outside
def writeCondorSubmitBase(condor_file="", dummy_exec_name="", logdir="", jobBatchName="undefined", memory=4000, maxtime=86400):    
    condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {ld}/$(ProcId).log
Output     = {ld}/$(ProcId).out
Error      = {ld}/$(ProcId).error
getenv      = True
environment = "LS_SUBCWD={here}"
next_job_start_delay = 1
request_memory = {mem}
#requirements = (OpSysAndVer =?= "SLCern6")
+MaxRuntime = {time}
+JobBatchName = "{jbn}"
'''.format(de=os.path.abspath(dummy_exec_name), ld=os.path.abspath(logdir), here=os.environ['PWD'], jbn=jobBatchName, mem=memory, time=maxtime ) )
    if os.environ['USER'] in ['mciprian']:
        # mydate = datetime.today()
        # month = int(mydate.month)
        # year  = int(mydate.year)
        # if month == 10 and year == 2019:
        #     pass
        #     condor_file.write('+AccountingGroup = "group_u_CMS.u_zh.priority"\n\n')
        # else:
        #     condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
        condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
    else:
        condor_file.write('\n')

parser = OptionParser(usage="%prog [options]")
parser.add_option(  "--recover-fill",     dest="recoverFill", action="store_true", default=False, help="When resubmitting calibration from hadd, first try to recover failed fills")
parser.add_option("-l", "--daemon-local",     dest="daemonLocal", action="store_true", default=False, help="Tells this script if the daemon is running locally (needed to renew the AFS token)")
parser.add_option("-t", "--token-file", dest="tokenFile",  type="string", default="", help="File needed to renew token (when daemon running locally)")
parser.add_option("--min-efficiency-recover-fill",   dest="minEfficiencyToRecoverFill",   type="float", default=0.97, help="Tolerance of EcalNtp loss. Require fraction of good EcalNtp abive this number to skip recover");
(options, args) = parser.parse_args()

mode = str(args[0])
pwd         = os.getcwd()
num         = 2

print "sys.argv = " + "  ".join(str(x) for x in sys.argv)
print "args = " + "  ".join(str(x) for x in args)
#print "mode = ", mode
print "Nargv = %d   Nargs = %d" % (len(sys.argv), len(args))

if ( mode.find('BATCH_RESU')==-1 ): # Batch system
   if len(args) != 2:
       print "usage thisPyton.py nITER queue"
       sys.exit(1)
elif ( mode.find('BATCH_RESU') != -1 ):                                   # Batch Resubmission
    if len(args) != 4:
       print "usage thisPyton.py BATCH_RESU nITER queue nJobs"
       sys.exit(1)
#Selec what mode you are running
RunCRAB = False; RunBatch = True; RunResub = True;
if ( mode.find('BATCH_RESU') != -1 ):
     RunBatch = False; RunResub = True;
else:
     RunBatch = True; RunResub = False;
ONLYHADD = False; ONLYFINHADD = False; ONLYFIT = False; ONLYMERGEFIT = False
if ( mode.find('ONLYHADD') != -1 ):
     ONLYHADD = True;
if ( mode.find('ONLYFINALHADD') != -1 ):
     ONLYFINHADD = True;
if ( mode.find('ONLYFIT') != -1 ):
     ONLYFIT = True;
if ( mode.find('ONLYMERGEFIT') != -1 ):
     ONLYMERGEFIT = True;

Add_path = ''
Add_pathOLDIter = ''
ListPaths = []
if ( RunResub ):
    njobs = int(sys.argv[4])
    queue = sys.argv[3]
    nIterations = nIterations - int(sys.argv[2])
else:
    njobs = int(sys.argv[1])
    queue = sys.argv[2]

outputdir = pwd+'/'+dirname
condorPath = outputdir + '/condor_files/'
logPath = outputdir + '/log'
srcPath  = outputdir + '/src'
cfgHaddPath  = outputdir + '/src/hadd'

# To compute the num of hadd
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = inputlist_f.readlines()

for iters in range(nIterations):

    if ( RunResub ):
        iters = iters + int(sys.argv[2])

    condordir = condorPath + '/iter_' + str(iters) 
    if not os.path.exists(condordir): os.makedirs(condordir)
    if not os.path.exists(logPath): os.makedirs(logPath)

    if ( not ONLYHADD and not ONLYFIT and not ONLYFINHADD and not ONLYMERGEFIT):

        renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
        # PREPARE CONDOR FILES FOR FILL AT ITER iters
        dummy_exec = open(condordir+'/dummy_exec_fill.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('bash $*\n')
        dummy_exec.close()

        logdir    = logPath    + '/Fill/iter_' + str(iters) 
        if not os.path.exists(logdir): os.makedirs(logdir)
        condor_file_name = condordir+'/condor_submit_fill.condor'
        condor_file = open(condor_file_name,'w')
        mymaxtimeFill = int(1. + ijobmax / 8.) * 12 * 3600
        if mymaxtimeFill > 48 * 3600:
            mymaxtimeFill = 48 * 3600
        writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Fill", 
                              memory=2500, maxtime=mymaxtimeFill) # this does not close the file

        print "\n*******  ITERATION " + str(iters) + "/" + str(nIterations-1) + "  *******"
        print "Submitting " + str(njobs) + " jobs"
        for ijob in range(njobs):
            #In case you want the stat. syst
            if ( mode.find('BATCH_RESU_SYST_1') != -1 ):
                 env_script_n = open(outputdir + "/cfgFile/Fill/iter_" + str(iters) + "/fillEps_iter_" + str(iters) + "_job_" + str(ijob) + ".py", 'a')
                 SystParamLine = 'process.analyzerFillEpsilon.SystOrNot = cms.untracked.double(1)\n'
                 env_script_n.write(SystParamLine)
                 env_script_n.close()
            if ( mode.find('BATCH_RESU_SYST_2') != -1 ):
                 env_script_n = open(outputdir + "/cfgFile/Fill/iter_" + str(iters) + "/fillEps_iter_" + str(iters) + "_job_" + str(ijob) + ".py", 'a')
                 SystParamLine = 'process.analyzerFillEpsilon.SystOrNot = cms.untracked.double(2)\n'
                 env_script_n.write(SystParamLine)
                 env_script_n.close()
            # preparing submission of filling tasks
            fill_src_n = srcPath + "/Fill/iter_" + str(iters) + "/submit_iter_"     + str(iters) + "_job_" + str(ijob) + ".sh"
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(fill_src_n)))
            #print '\n[job #' + str(ijob) + '] :: ' + submit_s

        condor_file.close()

        submit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        print ">>> Running --> " + submit_s
        # actually submitting filling tasks
        submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
        output = submitJobs.communicate()
        print "Out: " + str(output)

        time.sleep(15)    
        nFilljobs = checkNjobsCondor("ecalpro_Fill")
        print "There are {n} jobs for Fill part".format(n=nFilljobs)
        
        print 'Waiting for filling jobs to be finished...'
        # Daemon cheking running jobs
        nCheck = 0
        while nFilljobs > 0 :
            sleeptime = 900
            time.sleep(sleeptime)
            nFilljobs = checkNjobsCondor("ecalpro_Fill")
            print "I still see {n} jobs for Fill part".format(n=nFilljobs)
            checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
            datalines2 = (checkJobs2.communicate()[0]).splitlines()
            nCheck += 1
            if nCheck * sleeptime > 43200: 
                renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
                nCheck = 0

        print 'Done with the Fill part'

        # print "="*30
        # print "This is a test: execution will end here"
        # print "="*30
        # quit()

        ##########
        # only for ntuples, resubmit failed *EcalNtp*.root jobs (max number of resubmission is hardcoded, currently it is only 2 in order not to waste too much time)
        ##########
        # if MakeNtuple4optimization:

        #     NtpRecoveryAttempt = 0
        #     goodNtp = 0
        #     while goodNtp < njobs and NtpRecoveryAttempt < 2:
                
        #         logdir    = logPath    + '/Fill/iter_{it}_recovery_{nr}'.format(it=str(iters), nr=str(NtpRecoveryAttempt))
        #         if not os.path.exists(logdir): os.makedirs(logdir)
        #         condor_file_name = condordir+'/condor_submit_fill_recovery_{nr}.condor'.format(nr=str(NtpRecoveryAttempt))
        #         condor_file = open(condor_file_name,'w')
        #         writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Fill_recovery", memory=2500, maxtime=43200)  
        #         goodNtp = 0
        #         for ih in range(njobs):
        #             eosFile = eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + NameTag + "EcalNtp_" + str(ih) + ".root"
        #             testNtpFile_s = 'ls -l ' + eosFile   # eos is now mounted on lxplus
        #             print "checking the presence and the sanity of EcalNtp file: " + eosFile
        #             testNtpFile = subprocess.Popen([testNtpFile_s], stdout=subprocess.PIPE, shell=True);
        #             output = testNtpFile.communicate()[0]
        #             fsize = 0
        #             if len(output)>0:
        #                 print "output = ",output
        #                 fsize = int(output.split()[4])
        #             # I expect about some MB, so ask at least 100kB
        #             if len(output)==0 or fsize<100000:
        #                 print "The file " + eosFile + " is not present, or empty. Resubmitting ..."
        #                 Ntp_src_n = srcPath + "/Fill/iter_" + str(iters) + "/submit_iter_" + str(iters) + "_job_" + str(ijob) + ".sh"
        #                 condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Ntp_src_n)))
        #             else: goodNtp += 1
                    
        #         condor_file.close()
        #         Ntpsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        #         # actually submitting recovery tasks
        #         subJobs = subprocess.Popen([Ntpsubmit_s], stdout=subprocess.PIPE, shell=True);
        #         outJobs = subJobs.communicate()
        #         print outJobs

        #         time.sleep(15)
        #         nFilljobs = checkNjobsCondor("ecalpro_Fill_recovery")
        #         print "There are {n} jobs for Fill_recovery part".format(n=nFilljobs)

        #         print 'Waiting for filling jobs to be finished...'
        #         # Daemon cheking running jobs
        #         print "Checking recovery of Ntp ..."
        #         while nFilljobs > 0 :
        #             time.sleep(900)
        #             nFilljobs = checkNjobsCondor("ecalpro_Fill_recovery")
        #             print "I still see {n} jobs for Fill_recovery part".format(n=nFilljobs)
        #             checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
        #             datalines2 = (checkJobs2.communicate()[0]).splitlines()

        #         NtpRecoveryAttempt += 1
        #         print 'Done with Ntp recovery n.' + str(NtpRecoveryAttempt)

#     if MakeNtuple4optimization:
#         print """MakeNtuple4optimization is set to True in parameters.py
# Code will stop know before adding the *EcalNtp*.root files.
# It is better that you run on all the output files using a TChain. Indeed, these are big files, and the hadd part is slow and the jobs can fail in producing the output. 
# """
#         print "Done with iteration " + str(iters)
#         quit()

    #HADD for batch and CRAB, if you do not want just the finalHADD or the FIT
    if ( not ONLYFIT and not ONLYFINHADD and not ONLYMERGEFIT):
        
        #if MakeNtuple4optimization:
        if options.recoverFill:

            renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
            NtpRecoveryAttempt = 0
            goodNtp = 0
            while goodNtp < njobs and NtpRecoveryAttempt < 1:
                
                logdir    = logPath    + '/Fill/iter_{it}_recovery_{nr}'.format(it=str(iters), nr=str(NtpRecoveryAttempt))
                if not os.path.exists(logdir): os.makedirs(logdir)
                condor_file_name = condordir+'/condor_submit_fill_recovery_{nr}.condor'.format(nr=str(NtpRecoveryAttempt))
                condor_file = open(condor_file_name,'w')
                mymaxtimeFill = int(1. + ijobmax / 8.) * 12 * 3600
                if mymaxtimeFill > 48 * 3600:
                    mymaxtimeFill = 48 * 3600
                writeCondorSubmitBase(condor_file, condordir+'/dummy_exec_fill.sh', logdir, "ecalpro_Fill_recovery", 
                                      memory=2500, maxtime=mymaxtimeFill)  
                goodNtp = 0
                for ih in range(njobs):
                    eosFile = eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + NameTag + "EcalNtp_" + str(ih) + ".root"
                    Ntp_src_n = srcPath + "/Fill/iter_" + str(iters) + "/submit_iter_" + str(iters) + "_job_" + str(ih) + ".sh"
                    print "checking the presence and the sanity of EcalNtp file: " + eosFile
                    filesize=0
                    if os.path.exists(eosFile): filesize = os.path.getsize(eosFile)
                    if filesize>100000:
                        tf = TFile.Open("root://eoscms/"+eosFile)
                        if not tf or tf.IsZombie(): 
                            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Ntp_src_n)))
                        elif tf.TestBit(TFile.kRecovered):                    
                            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Ntp_src_n)))
                            tf.Close()
                        else:
                            goodNtp += 1
                            tf.Close()
                    else:
                        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Ntp_src_n)))
                                            
                condor_file.close()
                print "Found {n}/{ntot} good EcalNtp files.".format(n=goodNtp,ntot=njobs)
                nGoodOverTot = float(goodNtp)/float(njobs)
                if nGoodOverTot < options.minEfficiencyToRecoverFill:

                    print "Resubmitting failed fill jobs."
                    Ntpsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
                    # actually submitting recovery tasks
                    subJobs = subprocess.Popen([Ntpsubmit_s], stdout=subprocess.PIPE, shell=True);
                    outJobs = subJobs.communicate()
                    print outJobs

                    time.sleep(30)
                    nFilljobs = checkNjobsCondor("ecalpro_Fill_recovery")
                    print "There are {n} jobs for Fill_recovery part".format(n=nFilljobs)

                    print 'Waiting for filling jobs to be finished...'
                    # Daemon cheking running jobs
                    print "Checking recovery of Ntp ..."
                    nCheck = 0
                    while nFilljobs > 0 :
                        sleeptime = 900
                        time.sleep(sleeptime)
                        nFilljobs = checkNjobsCondor("ecalpro_Fill_recovery")
                        print "I still see {n} jobs for Fill_recovery part ({nr})".format(n=nFilljobs,nr=NtpRecoveryAttempt)
                        checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
                        datalines2 = (checkJobs2.communicate()[0]).splitlines()
                        nCheck += 1
                        if nCheck * sleeptime > 43200: 
                            renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
                            nCheck = 0


                    print 'Done with Ntp recovery n.' + str(NtpRecoveryAttempt)
                else:
                    print "Fraction of EcalNtp root file was {perc:.1%} (tolerance was set to {tol:.1%}".format(perc=nGoodOverTot, tol=options.minEfficiencyToRecoverFill)
                    print "Fill recovery was not attempted."

                NtpRecoveryAttempt += 1
        #  END OF FILL RECOVERY
        #########################

        renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
        # PREPARE CONDOR FILES FOR HADD AT ITER iters
        dummy_exec = open(condordir+'/dummy_exec_hadd.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('bash $*\n')
        dummy_exec.close()

        logdir    = logPath    + '/Hadd/iter_' + str(iters) 
        if not os.path.exists(logdir): os.makedirs(logdir)
        condor_file_name = condordir+'/condor_submit_hadd.condor'
        condor_file = open(condor_file_name,'w')
        writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Hadd",memory=2000, maxtime=5000) # this does not close the file

        print 'Now adding files...'
        Nlist = 0
        if not( RunCRAB ):
           inputlist_v = inputlistbase_v[:]
           NrelJob = float(len(inputlist_v)) / float(ijobmax)
           if( float(int(NrelJob) - NrelJob) < 0. ):
               NrelJob = int(NrelJob) + 1
           Nlist_flo = float(NrelJob/nHadd) + 1.
           Nlist = int(Nlist_flo)
        print "Number of Hadd in parallel: " + str(Nlist)
        for nHadds in range(Nlist):
            Hadd_src_n = srcPath + "/hadd/HaddCfg_iter_" + str(iters) + "_job_" + str(nHadds) + ".sh"
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Hadd_src_n)))
            #Before each HADD we need to check if the all the files in the list are present
            #BUT we do that only if you are working on batch
            if not( RunCRAB ):
               Grepcommand = "grep -i list " + Hadd_src_n + " | grep -v echo | grep -v bash | awk '{print $2}'"
               myGrep = subprocess.Popen([Grepcommand], stdout=subprocess.PIPE, shell=True )
               FoutGrep = myGrep.communicate()
               # FoutGrep is something like the following
               # ('/afs_path_to_dirName/src/hadd/hadd_iter_XXX_step_YYY.list`\n', None)
               # we want to keep /afs_path_to_dirName/src/hadd/hadd_iter_XXX_step_YYY.list
               # removing (' and `\n', None)
               FoutGrep_2 = str(FoutGrep)[2:]
               FoutGrep_2 = str(FoutGrep_2)[:-11]
               #print 'Checking ' + str(FoutGrep_2)
               #Chech The size for each line
               f = open( str(FoutGrep_2) )
               lines = f.readlines()
               f.close()
               # create backup of original list of files
               fbckp = open( FoutGrep_2.replace(".list","_backup.list"), "w")
               for l in lines:
                   fbckp.write(l)
               fbckp.close()

               # creating a list of EcalNtp files that are actually present on eos. The new list will be overwritten on the original one (which was backuped)
               newlines = []
               for filetoCheck in lines:
                   if not os.path.exists(filetoCheck.strip()):
                       #print 'HADD::MISSING: {f}'.format(f=os.path.basename(filetoCheck))
                       continue
                   else:
                       filesize = os.path.getsize(filetoCheck.strip())
                       #If is corrupted (size too small), remove it from the list
                       if( filesize<100000 ):
                           #print 'HADD::Bad size {size} for: {f}'.format(size=filesize, f=os.path.basename(filetoCheck))
                           continue
                       else:
                           # at this point the file should be good, but let's check if there are no recovered keys
                           #open and check there are no recovered keys: in this case remove these files from the list, otherwise hadd might fail
                           tf = TFile.Open("root://eoscms/"+filetoCheck.strip())
                           if not tf or tf.IsZombie(): continue
                           if tf.TestBit(TFile.kRecovered):
                               #print "HADD::Attemp to recover file {f}".format(f=filetoCheck.strip())
                               tf.Close()
                               continue
                           tf.Close()
                   newlines.append(filetoCheck)

               #moving the .list to the correct one
               if( len(newlines) ):
                   prunedfile = FoutGrep_2.replace(".list","_pruned.list")
                   fprun = open(prunedfile,"w")
                   for l in newlines:
                       fprun.write(l)
                   fprun.close()
                   MoveComm = "cp " + prunedfile + " " + str(FoutGrep_2)
                   MoveC = subprocess.Popen([MoveComm], stdout=subprocess.PIPE, shell=True);
                   mvOut = MoveC.communicate()
                   #print "Some files were removed in " + str(FoutGrep_2)
                   #print "Copied " + prunedfile + " into " + str(FoutGrep_2)
            #End of the check, sending the job
            print "Preparing job to hadd files in list number " + str(nHadds) + "/" + str(Nlist - 1)  #nHadds goes from 0 to Nlist -1

        condor_file.close()
        Hsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        print ">>> Running --> " + Hsubmit_s
        subJobs = subprocess.Popen([Hsubmit_s], stdout=subprocess.PIPE, shell=True);
        outJobs = subJobs.communicate()

        print str(outJobs)
        time.sleep(15)
        nHaddjobs = checkNjobsCondor("ecalpro_Hadd")
        print "There are {n} jobs for Hadd part".format(n=nHaddjobs)

        print 'Waiting for all the hadd...'
        # Daemon cheking running jobs
        while nHaddjobs > 0 :
            time.sleep(300)
            nHaddjobs = checkNjobsCondor("ecalpro_Hadd")
            print "I still see {n} jobs for Hadd part".format(n=nHaddjobs)
            #checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
            #datalines2 = (checkJobs2.communicate()[0]).splitlines()
        print 'Done with various hadd'

        # Check if all the hadds are there and files are not empty
        goodHadds = 0
        for ih in range(Nlist):
            eosFile = eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + NameTag + "epsilonPlots_" + str(ih) + ".root"
            filesize=0
            if os.path.exists(eosFile): filesize = os.path.getsize(eosFile)
            if filesize>100000:
                tf = TFile.Open("root://eoscms/"+eosFile)
                if not tf or tf.IsZombie(): continue
                if not tf.TestBit(TFile.kRecovered):                    
                    goodHadds += 1
                tf.Close()

        if goodHadds == Nlist:
            print "All files are present: going to hadd them"
        else:
            print "Some files are missing or empty or some keys were recovered: going to recover them"

        HaddRecoveryAttempt = 0
        # this recovery is often unsuccessful: if it failed the first time, there is a high chance that the problem is persistent
        # better to try only a second time, not more
        while goodHadds < Nlist and HaddRecoveryAttempt < 1:

            logdir    = logPath    + '/Hadd/iter_{it}_recovery_{nr}'.format(it=str(iters), nr=str(HaddRecoveryAttempt))
            if not os.path.exists(logdir): os.makedirs(logdir)
            condor_file_name = condordir+'/condor_submit_hadd_recovery_{nr}.condor'.format(nr=str(HaddRecoveryAttempt))
            condor_file = open(condor_file_name,'w')
            writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Hadd_recovery",memory=2000, maxtime=5000)  

            print "Trying to recover failed hadd. Attempt n." + str(HaddRecoveryAttempt)
            goodHadds = 0
            for ih in range(Nlist):
                eosFile = eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + NameTag + "epsilonPlots_" + str(ih) + ".root"
                filesize=0
                if os.path.exists(eosFile): 
                    filesize = os.path.getsize(eosFile)
                # should be of the order of some MB, so ask at least 100 kB (if empty, it is about 1.1 kB)
                Hadd_src_n = srcPath + "/hadd/HaddCfg_iter_" + str(iters) + "_job_" + str(ih) + ".sh"
                if filesize<100000:
                    #print "The file " + eosFile + " is not present, or empty. Redoing hadd..."
                    condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Hadd_src_n)))
                    #print Hadd_src_n                
                else: 
                    tf = TFile.Open("root://eoscms/"+eosFile)
                    if not tf or tf.IsZombie():
                        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Hadd_src_n)))
                        continue # must not close file
                    if tf.TestBit(TFile.kRecovered):
                        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(Hadd_src_n)))
                    else:
                        goodHadds += 1
                    tf.Close()
                        

            print "Good files: {num}/{den}".format(num=goodHadds,den=Nlist)
            if goodHadds == Nlist: 
                print "All files recovered successfully!"
                break   
            # inside this loop, it can be that goodHadds == Nlist, because when we enter we resubmit some jobs for sure, and at the next attempt 
            # we still don't know if they were all successful: if so, we exit the loop

            condor_file.close()
            Hsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
            print ">>> Running --> " + Hsubmit_s
            # actually submitting recovery tasks
            subJobs = subprocess.Popen([Hsubmit_s], stdout=subprocess.PIPE, shell=True);
            outJobs = subJobs.communicate()
            print outJobs

            time.sleep(15)
            nHaddjobs = checkNjobsCondor("ecalpro_Hadd_recovery")
            print "There are {n} jobs for Hadd_recovery part".format(n=nHaddjobs)

            print 'Waiting for Hadd recovery jobs to be finished...'
            # Daemon cheking running jobs
            print "Checking recovery of Hadd ..."
            nCheck = 0
            while nHaddjobs > 0 :
                sleeptime = 300
                time.sleep(sleeptime)
                nHaddjobs = checkNjobsCondor("ecalpro_Hadd_recovery")
                print "I still see {n} jobs for Hadd_recovery part".format(n=nHaddjobs)
                #checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
                #datalines2 = (checkJobs2.communicate()[0]).splitlines()
                nCheck += 1
                if nCheck * sleeptime > 43200: 
                    renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
                    nCheck = 0

            HaddRecoveryAttempt += 1
     
            print 'Done with hadd recovery'


    if ( not ONLYFIT and not ONLYMERGEFIT):
        renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
        print 'Now The Final One...'
        
        # PREPARE CONDOR FILES FOR FINAL HADD AT ITER iters
        dummy_exec = open(condordir+'/dummy_exec_finalHadd.sh','w')
        dummy_exec.write('#!/bin/bash\n')
        dummy_exec.write('bash $*\n')
        dummy_exec.close()

        logdir    = logPath    + '/Final_Hadd/iter_' + str(iters) 
        if not os.path.exists(logdir): os.makedirs(logdir)
        condor_file_name = condordir+'/condor_submit_finalHadd.condor'
        condor_file = open(condor_file_name,'w')
        writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_FinalHadd", memory=2000, maxtime=5000) # this does not close the file
        
        FHadd_src_n = srcPath + "/hadd/Final_HaddCfg_iter_" + str(iters) + ".sh"
        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(FHadd_src_n)))        
        condor_file.close()
        
        FHsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        print ">>> Running --> " + FHsubmit_s
        FsubJobs = subprocess.Popen([FHsubmit_s], stdout=subprocess.PIPE, shell=True);
        FoutJobs = FsubJobs.communicate()
        print FoutJobs
        time.sleep(5)

        nFinHaddjobs = checkNjobsCondor("ecalpro_FinalHadd")     
        print "There are {n} jobs for FinalHadd part".format(n=nFinHaddjobs)
        print 'Waiting for the Final hadd...'
        # Daemon cheking running jobs
        while nFinHaddjobs > 0 :
            time.sleep(120)
            nFinHaddjobs = checkNjobsCondor("ecalpro_FinalHadd")
            print "I still see {n} jobs for FinalHadd part".format(n=nFinHaddjobs)
        print 'Done with final hadd'


    if MakeNtuple4optimization:
        # it actually stopped already before hadding files
        print """MakeNtuple4optimization is set to True in parameters.py
From the current behaviour of FillEpsilonPlot.cc code (version 11/06/2017), it means the histogram used to do the fit for 
each crystal are not saved and therefore the Fit part will crash because these histograms will not be found in '*epsilonPlots.root' file.
Code will stop know, since it is assumed that if you are optimizing selection then the Fit part is not needed (and you don't need further iterations)
If this is not the case, modify FillEpsilonPlot.cc
"""
        print "Done with iteration " + str(iters)
        quit()

    # if running with flag to fold SM, first do that part. It is not needed to run a job, can be done locally in few minutes
    # in case we only need to merge fit, this is likely not needed, so skip it

    if foldInSuperModule and not ONLYMERGEFIT:
        # check if the file is already present, in which case this step can be skipped
        histograms_foldedInSM_exists = False
        hFoldFile = eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'histograms_foldedInSM.root'
        if not os.path.isfile(hFoldFile):
            srcFold_n = srcPath + "/Fit/submit_justFoldSM_iter_" + str(iters) + ".sh"
            Fitsubmit_s = "{src}".format(src=srcFold_n)
            print ">>> Running --> " + Fitsubmit_s      
            FsubJobs = subprocess.Popen([Fitsubmit_s], stdout=subprocess.PIPE, shell=True);
            FoutJobs = FsubJobs.communicate()  # do I really need to communicate to stdout?
            print FoutJobs        
            #os.system(Fitsubmit_s)    

    # N of Fit to send
    nEBindependentXtals = 1699 if foldInSuperModule else 61199
    nEB = nEBindependentXtals/nFit
    if (nEBindependentXtals%nFit != 0):
        nEB = int(nEB) +1
    nEE = 14647/nFit
    if (14647%nFit != 0) :
        nEE = int(nEE) +1

    if Barrel_or_Endcap == "ONLY_ENDCAP": nEB = 0
    if Barrel_or_Endcap == "ONLY_BARREL": nEE = 0
        
    # For final hadd
    ListFinalHaddEB = list()
    ListFinalHaddEE = list()

    # PREPARE CONDOR FILES FOR FIT AT ITER iters
    dummy_exec = open(condordir+'/dummy_exec_fit.sh','w')
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('bash $*\n')
    dummy_exec.close()
    
    logdir    = logPath    + '/Fit/iter_' + str(iters) 
    if not os.path.exists(logdir): os.makedirs(logdir)
    condor_file_name = condordir+'/condor_submit_fit.condor'
    condor_file = open(condor_file_name,'w')
    writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Fit", memory=2000, maxtime=86400) # this does not close the file

    # preparing submission of fit tasks (EB)
    if (not ONLYMERGEFIT): print 'Submitting ' + str(nEB) + ' jobs to fit the Barrel'
    for inteb in range(nEB):
        fit_src_n = srcPath + "/Fit/submit_EB_" + str(inteb) + "_iter_"     + str(iters) + ".sh"
        ListFinalHaddEB.append(eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName )
        if (not ONLYMERGEFIT):
            print 'About to EB fit:'
            print eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName            
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(fit_src_n)))

    # preparing submission of fit tasks (EE)
    if (not ONLYMERGEFIT): print 'Submitting ' + str(nEE) + ' jobs to fit the Endcap'
    for inte in range(nEE):        
        fit_src_n = srcPath + "/Fit/submit_EE_" + str(inte) + "_iter_"     + str(iters) + ".sh"
        ListFinalHaddEE.append(eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName)
        if (not ONLYMERGEFIT):
            print 'About to EE fit:'
            print eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(fit_src_n)))

    condor_file.close()
            
    if (not ONLYMERGEFIT ):

        renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
        Fitsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        print ">>> Running --> " + Fitsubmit_s
        FsubJobs = subprocess.Popen([Fitsubmit_s], stdout=subprocess.PIPE, shell=True);
        FoutJobs = FsubJobs.communicate()
        print FoutJobs
        time.sleep(5)

        nFitjobs = checkNjobsCondor("ecalpro_Fit")     
        print "There are {n} jobs for Fit part".format(n=nFitjobs)
        print 'Waiting for fit jobs to be finished...'
        # Daemon cheking running jobs
        while nFitjobs > 0 :
            time.sleep(300)
            nFitjobs = checkNjobsCondor("ecalpro_Fit")
            print "I still see {n} jobs for Fit part".format(n=nFitjobs)
        print "Done with fitting! Now we have to merge all fits in one Calibmap.root"

    # check all fits are there (this check is made only once at the moment, if the file is still missing the code will exit later)
    # this check should be made until all fits are present, otherwise the code crashes, butif it keeps happening there might be something serious to solve
    # so better to have the code crash and investigate locally


    renewTokenAFS(daemonLocal=options.daemonLocal, infile=options.tokenFile) 
    allFitsGood = True
    fit_src_toResub = []
    for inteb in range(nEB):
        fit_src_n = srcPath + "/Fit/submit_EB_" + str(inteb) + "_iter_"     + str(iters) + ".sh"
        thisfile = eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName
        #thisfile_f = TFile.Open(thisfile)
        if not os.path.isfile(thisfile):
            print "Will resubmit missing file {f}".format(f=thisfile)
            allFitsGood = False
            fit_src_toResub.append(fit_src_n)
    for inte in range(nEE):        
        fit_src_n = srcPath + "/Fit/submit_EE_" + str(inte) + "_iter_"     + str(iters) + ".sh"
        thisfile = eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName
        #thisfile_f = TFile.Open(thisfile)
        #if not thisfile_f:
        if not os.path.isfile(thisfile):
            print "Will resubmit missing file {f}".format(f=thisfile)
            allFitsGood = False
            fit_src_toResub.append(fit_src_n)

    if not allFitsGood:
        logdir    = logPath    + '/Fit_recovery/iter_' + str(iters) 
        if not os.path.exists(logdir): os.makedirs(logdir)
        condor_file_name = condordir+'/condor_submit_fit_recovery.condor'
        condor_file = open(condor_file_name,'w')
        writeCondorSubmitBase(condor_file, dummy_exec.name, logdir, "ecalpro_Fit_recovery", memory=2000, maxtime=86400) # this does not close the file
        for fit in fit_src_toResub:
            condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(fit)))               
        condor_file.close()

        Fitsubmit_s = "condor_submit {cfn}".format(cfn=condor_file_name)
        print ">>> Running --> " + Fitsubmit_s
        FsubJobs = subprocess.Popen([Fitsubmit_s], stdout=subprocess.PIPE, shell=True);
        FoutJobs = FsubJobs.communicate()
        print FoutJobs
        time.sleep(5)

        nFitjobs = checkNjobsCondor("ecalpro_Fit_recovery")     
        print "There are {n} jobs for Fit recovery part".format(n=nFitjobs)
        print 'Waiting for fit jobs to be finished...'
        # Daemon cheking running jobs
        while nFitjobs > 0 :
            time.sleep(60)
            nFitjobs = checkNjobsCondor("ecalpro_Fit_recovery")
            print "I still see {n} jobs for Fit recovery part".format(n=nFitjobs)
        print "Done with fitting recovery! Now we have to merge all fits in one Calibmap.root"

    # Merge Final CalibMap

    #Run only on EB or EE if needed
    ListFinalHadd = list()
    if Barrel_or_Endcap=='ONLY_BARREL':
       ListFinalHadd = ListFinalHaddEB
    if Barrel_or_Endcap=='ONLY_ENDCAP':
       ListFinalHadd = ListFinalHaddEE
    if (Barrel_or_Endcap=='ALL_PLEASE'):
       ListFinalHadd = ListFinalHaddEB
       ListFinalHadd = ListFinalHadd + ListFinalHaddEE

    from PhysicsTools.PythonAnalysis import *
    gSystem.Load("libFWCoreFWLite.so")
    #AutoLibraryLoader.enable()
    FWLiteEnabler.enable()
    finalCalibMapFileName = eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + Add_path + "/" + NameTag + calibMapName
    f = TFile.Open(finalCalibMapFileName, 'recreate')
    if not f:
        print "WARNING in calibJobHandlerCondor.py: file '" + finalCalibMapFileName +  "' not opened correctly. Quitting ..."
        quit()
    else:
        f.cd()

    for n_repeat in range(2):

        # when isEoverEtrue is True we have two maps for the two photons (for each detector), otherwise we only have one 
        # the things of the second photon are named adding a suffix "_g2"

        if n_repeat == 1 and not isEoverEtrue:
            continue

        # suffix not used at the moment, I just change the name below
        # map_suffix = ""
        # if isEoverEtrue:
        #     if n_repeat == 1:
        #         map_suffix = "_g2"
        # else: 
        #     if n_repeat == 1:
        #     continue

        # Create 2 struct objects for EB and 2 for EE, but only once in case of MC and isEoverEtrue == True
        if n_repeat == 0:
            if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
               gROOT.ProcessLine(\
                 "struct EBStruct{\
                   Int_t rawId_;\
                   Int_t hashedIndex_;\
                   Int_t ieta_;\
                   Int_t iphi_;\
                   Int_t iSM_;\
                   Int_t iMod_;\
                   Int_t iTT_;\
                   Int_t iTTeta_;\
                   Int_t iTTphi_;\
                   Int_t iter_;\
                   Double_t coeff_;\
                   Double_t Signal_;\
                   Double_t Backgr_;\
                   Double_t Chisqu_;\
                   Double_t Ndof_;\
                   Double_t fit_mean_;\
                   Double_t fit_mean_err_;\
                   Double_t fit_sigma_;\
                   Double_t fit_Snorm_;\
                   Double_t fit_b0_;\
                   Double_t fit_b1_;\
                   Double_t fit_b2_;\
                   Double_t fit_b3_;\
                   Double_t fit_Bnorm_;\
                 };")
               gROOT.ProcessLine(\
                 "struct EB1Struct{\
                   Int_t rawId;\
                   Int_t hashedIndex;\
                   Int_t ieta;\
                   Int_t iphi;\
                   Int_t iSM;\
                   Int_t iMod;\
                   Int_t iTT;\
                   Int_t iTTeta;\
                   Int_t iTTphi;\
                   Int_t iter;\
                   Double_t coeff;\
                   Double_t Signal;\
                   Double_t Backgr;\
                   Double_t Chisqu;\
                   Double_t Ndof;\
                   Double_t fit_mean;\
                   Double_t fit_mean_err;\
                   Double_t fit_sigma;\
                   Double_t fit_Snorm;\
                   Double_t fit_b0;\
                   Double_t fit_b1;\
                   Double_t fit_b2;\
                   Double_t fit_b3;\
                   Double_t fit_Bnorm;\
                 };")

            if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
               gROOT.ProcessLine(\
                 "struct EEStruct{\
                   Int_t ix_;\
                   Int_t iy_;\
                   Int_t zside_;\
                   Int_t sc_;\
                   Int_t isc_;\
                   Int_t ic_;\
                   Int_t iquadrant_;\
                   Int_t hashedIndex_;\
                   Int_t iter_;\
                   Double_t coeff_;\
                   Double_t Signal_;\
                   Double_t Backgr_;\
                   Double_t Chisqu_;\
                   Double_t Ndof_;\
                   Double_t fit_mean_;\
                   Double_t fit_mean_err_;\
                   Double_t fit_sigma_;\
                   Double_t fit_Snorm_;\
                   Double_t fit_b0_;\
                   Double_t fit_b1_;\
                   Double_t fit_b2_;\
                   Double_t fit_b3_;\
                   Double_t fit_Bnorm_;\
                 };")
               gROOT.ProcessLine(\
                 "struct EE1Struct{\
                   Int_t ix;\
                   Int_t iy;\
                   Int_t zside;\
                   Int_t sc;\
                   Int_t isc;\
                   Int_t ic;\
                   Int_t iquadrant;\
                   Int_t hashedIndex;\
                   Int_t iter;\
                   Double_t coeff;\
                   Double_t Signal;\
                   Double_t Backgr;\
                   Double_t Chisqu;\
                   Double_t Ndof;\
                   Double_t fit_mean;\
                   Double_t fit_mean_err;\
                   Double_t fit_sigma;\
                   Double_t fit_Snorm;\
                   Double_t fit_b0;\
                   Double_t fit_b1;\
                   Double_t fit_b2;\
                   Double_t fit_b3;\
                   Double_t fit_Bnorm;\
                 };")
               
        if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
            s = EBStruct()
            s1 = EB1Struct()
        if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
            t = EEStruct()
            t1 = EE1Struct()

        if f.IsOpen():
            f.cd()
        else:
            print "ERROR: it seems the output file '" + finalCalibMapFileName + "' is no longer opened! n_repeat = %d" % n_repeat
            quit()

        if isEoverEtrue and n_repeat == 1:
            calibMap_EB = TH2F("calibMap_EB_g2", "EB calib coefficients: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
            calibMap_EEm = TH2F("calibMap_EEm_g2", "EE- calib coefficients", 100,0.5,100.5,100,0.5,100.5)
            calibMap_EEp = TH2F("calibMap_EEp_g2", "EE+ calib coefficients", 100,0.5,100.5,100,0.5,100.5)
            TreeEB = TTree("calibEB_g2", "Tree of EB Inter-calibration constants")
            TreeEE = TTree("calibEE_g2", "Tree of EE Inter-calibration constants")
        else:
            calibMap_EB = TH2F("calibMap_EB", "EB calib coefficients: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
            calibMap_EEm = TH2F("calibMap_EEm", "EE- calib coefficients", 100,0.5,100.5,100,0.5,100.5)
            calibMap_EEp = TH2F("calibMap_EEp", "EE+ calib coefficients", 100,0.5,100.5,100,0.5,100.5)
            TreeEB = TTree("calibEB", "Tree of EB Inter-calibration constants")
            TreeEE = TTree("calibEE", "Tree of EE Inter-calibration constants")

        if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
           TreeEB.Branch('rawId_'      , AddressOf(s,'rawId_'),'rawId_/I')
           TreeEB.Branch('hashedIndex_', AddressOf(s,'hashedIndex_'),'hashedIndex_/I')
           TreeEB.Branch('ieta_'       , AddressOf(s,'ieta_'),'ieta_/I')
           TreeEB.Branch('iphi_'       , AddressOf(s,'iphi_'),'iphi_/I')
           TreeEB.Branch('iSM_'        , AddressOf(s,'iSM_'),'iSM_/I')
           TreeEB.Branch('iMod_'       , AddressOf(s,'iMod_'),'iMod_/I')
           TreeEB.Branch('iTT_'        , AddressOf(s,'iTT_'),'iTT_/I')
           TreeEB.Branch('iTTeta_'     , AddressOf(s,'iTTeta_'),'iTTeta_/I')
           TreeEB.Branch('iTTphi_'     , AddressOf(s,'iTTphi_'),'iTTphi_/I')
           TreeEB.Branch('iter_'       , AddressOf(s,'iter_'),'iter_/I')
           TreeEB.Branch('coeff_'      , AddressOf(s,'coeff_'),'coeff_/F')
           TreeEB.Branch('Chisqu_'     , AddressOf(s,'Chisqu_'),'Chisqu_/F')
           TreeEB.Branch('Ndof_'       , AddressOf(s,'Ndof_'),'Ndof_/F')
           TreeEB.Branch('fit_mean_'   , AddressOf(s,'fit_mean_'),'fit_mean_/F')
           TreeEB.Branch('fit_mean_err_'   , AddressOf(s,'fit_mean_err_'),'fit_mean_err_/F')
           TreeEB.Branch('fit_sigma_'  , AddressOf(s,'fit_sigma_'),'fit_sigma_/F')
           if not isEoverEtrue:
               TreeEB.Branch('Signal_'     , AddressOf(s,'Signal_'),'Signal_/F')
               TreeEB.Branch('Backgr_'     , AddressOf(s,'Backgr_'),'Backgr_/F')
               TreeEB.Branch('fit_Snorm_'  , AddressOf(s,'fit_Snorm_'),'fit_Snorm_/F')
               TreeEB.Branch('fit_b0_'     , AddressOf(s,'fit_b0_'),'fit_b0_/F')
               TreeEB.Branch('fit_b1_'     , AddressOf(s,'fit_b1_'),'fit_b1_/F')
               TreeEB.Branch('fit_b2_'     , AddressOf(s,'fit_b2_'),'fit_b2_/F')
               TreeEB.Branch('fit_b3_'     , AddressOf(s,'fit_b3_'),'fit_b3_/F')
               TreeEB.Branch('fit_Bnorm_'  , AddressOf(s,'fit_Bnorm_'),'fit_Bnorm_/F')

    
        if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
           TreeEE.Branch('ix_'         , AddressOf(t,'ix_'),'ix_/I')
           TreeEE.Branch('iy_'         , AddressOf(t,'iy_'),'iy_/I')
           TreeEE.Branch('zside_'      , AddressOf(t,'zside_'),'zside_/I')
           TreeEE.Branch('sc_'         , AddressOf(t,'sc_'),'sc_/I')
           TreeEE.Branch('isc_'        , AddressOf(t,'isc_'),'isc_/I')
           TreeEE.Branch('ic_'         , AddressOf(t,'ic_'),'ic_/I')
           TreeEE.Branch('iquadrant_'  , AddressOf(t,'iquadrant_'),'iquadrant_/I')
           TreeEE.Branch('hashedIndex_', AddressOf(t,'hashedIndex_'),'hashedIndex_/I')
           TreeEE.Branch('iter_'       , AddressOf(t,'iter_'),'iter_/I')
           TreeEE.Branch('coeff_'      , AddressOf(t,'coeff_'),'coeff_/F')
           TreeEE.Branch('Chisqu_'     , AddressOf(t,'Chisqu_'),'Chisqu_/F')
           TreeEE.Branch('Ndof_'       , AddressOf(t,'Ndof_'),'Ndof_/F')
           TreeEE.Branch('fit_mean_'   , AddressOf(t,'fit_mean_'),'fit_mean_/F')
           TreeEE.Branch('fit_mean_err_'   , AddressOf(t,'fit_mean_err_'),'fit_mean_err_/F')
           TreeEE.Branch('fit_sigma_'  , AddressOf(t,'fit_sigma_'),'fit_sigma_/F')
           if not isEoverEtrue:
               TreeEE.Branch('Signal_'     , AddressOf(t,'Signal_'),'Signal_/F')
               TreeEE.Branch('Backgr_'     , AddressOf(t,'Backgr_'),'Backgr_/F')
               TreeEE.Branch('fit_Snorm_'  , AddressOf(t,'fit_Snorm_'),'fit_Snorm_/F')
               TreeEE.Branch('fit_b0_'     , AddressOf(t,'fit_b0_'),'fit_b0_/F')
               TreeEE.Branch('fit_b1_'     , AddressOf(t,'fit_b1_'),'fit_b1_/F')
               TreeEE.Branch('fit_b2_'     , AddressOf(t,'fit_b2_'),'fit_b2_/F')
               TreeEE.Branch('fit_b3_'     , AddressOf(t,'fit_b3_'),'fit_b3_/F')
               TreeEE.Branch('fit_Bnorm_'  , AddressOf(t,'fit_Bnorm_'),'fit_Bnorm_/F')

        # print "Printing list of files on eos ..."
        # print "############################"
        # cmdEosLs = 'ls ' + eosPath + '/' + dirname + '/iter_' + str(iters) + "/"
        # eosFileList = subprocess.Popen([cmdEosLs], stdout=subprocess.PIPE, shell=True);
        # print eosFileList.communicate()
        # print "############################"

        print "Going to merge the following files ..."
        for thisfile_s in ListFinalHadd:
            thisfile_s = thisfile_s.rstrip()
            print "file --> " + str(thisfile_s)
            thisfile_f = TFile.Open(thisfile_s)
            if not thisfile_f: 
                print "Error in calibJobHandlerCondor.py --> file not found" 
                quit()
            #Taking Interval and EB or EE
            h_Int = thisfile_f.Get("hint")
            #
            #sys.exit()
            init = h_Int.GetBinContent(1)
            finit = h_Int.GetBinContent(2)
            EEoEB = h_Int.GetBinContent(3)

            #TTree
            if EEoEB == 0:
               if isEoverEtrue and n_repeat == 1:
                   thisTree = thisfile_f.Get("calibEB_g2")
               else:
                   thisTree = thisfile_f.Get("calibEB")
               thisTree.SetBranchAddress( 'rawId',AddressOf(s1,'rawId'));
               thisTree.SetBranchAddress( 'hashedIndex',AddressOf(s1,'hashedIndex'));
               thisTree.SetBranchAddress( 'ieta',AddressOf(s1,'ieta'));
               thisTree.SetBranchAddress( 'iphi',AddressOf(s1,'iphi'));
               thisTree.SetBranchAddress( 'iSM',AddressOf(s1,'iSM'));
               thisTree.SetBranchAddress( 'iMod',AddressOf(s1,'iMod'));
               thisTree.SetBranchAddress( 'iTT',AddressOf(s1,'iTT'));
               thisTree.SetBranchAddress( 'iTTeta',AddressOf(s1,'iTTeta'));
               thisTree.SetBranchAddress( 'iTTphi',AddressOf(s1,'iTTphi'));
               thisTree.SetBranchAddress( 'iter',AddressOf(s1,'iter'));
               thisTree.SetBranchAddress( 'coeff',AddressOf(s1,'coeff'));
               thisTree.SetBranchAddress( 'Chisqu',AddressOf(s1,'Chisqu'));
               thisTree.SetBranchAddress( 'Ndof',AddressOf(s1,'Ndof'));
               thisTree.SetBranchAddress( 'fit_mean',AddressOf(s1,'fit_mean'));
               thisTree.SetBranchAddress( 'fit_mean_err',AddressOf(s1,'fit_mean_err'));
               thisTree.SetBranchAddress( 'fit_sigma',AddressOf(s1,'fit_sigma'));
               if not isEoverEtrue:
                   thisTree.SetBranchAddress( 'Signal',AddressOf(s1,'Signal'));
                   thisTree.SetBranchAddress( 'Backgr',AddressOf(s1,'Backgr'));
                   thisTree.SetBranchAddress( 'fit_Snorm',AddressOf(s1,'fit_Snorm'));
                   thisTree.SetBranchAddress( 'fit_b0',AddressOf(s1,'fit_b0'));
                   thisTree.SetBranchAddress( 'fit_b1',AddressOf(s1,'fit_b1'));
                   thisTree.SetBranchAddress( 'fit_b2',AddressOf(s1,'fit_b2'));
                   thisTree.SetBranchAddress( 'fit_b3',AddressOf(s1,'fit_b3'));
                   thisTree.SetBranchAddress( 'fit_Bnorm',AddressOf(s1,'fit_Bnorm'));
               for ntre in range(thisTree.GetEntries()):
                   thisTree.GetEntry(ntre);
                   if (ntre>=init and ntre<=finit):
                       s.rawId_ = s1.rawId
                       s.hashedIndex_ = s1.hashedIndex
                       s.ieta_ = s1.ieta
                       s.iphi_ = s1.iphi
                       s.iSM_ = s1.iSM
                       s.iMod_ = s1.iMod
                       s.iTT_ = s1.iTT
                       s.iTTeta_ = s1.iTTeta
                       s.iTTphi_ = s1.iTTphi
                       s.iter_ = s1.iter
                       s.coeff_ = s1.coeff
                       s.Chisqu_ = s1.Chisqu
                       s.Ndof_ = s1.Ndof
                       s.fit_mean_ = s1.fit_mean
                       s.fit_mean_err_ = s1.fit_mean_err
                       s.fit_sigma_ = s1.fit_sigma
                       if not isEoverEtrue:
                           s.Signal_ = s1.Signal
                           s.Backgr_ = s1.Backgr
                           s.fit_Snorm_ = s1.fit_Snorm
                           s.fit_b0_ = s1.fit_b0
                           s.fit_b1_ = s1.fit_b1
                           s.fit_b2_ = s1.fit_b2
                           s.fit_b3_ = s1.fit_b3
                           s.fit_Bnorm_ = s1.fit_Bnorm
                       TreeEB.Fill()
            else:
               if isEoverEtrue and n_repeat == 1:
                   thisTree = thisfile_f.Get("calibEE_g2")
               else:
                   thisTree = thisfile_f.Get("calibEE")
               thisTree.SetBranchAddress( 'ix',AddressOf(t1,'ix'));
               thisTree.SetBranchAddress( 'iy',AddressOf(t1,'iy'));
               thisTree.SetBranchAddress( 'zside',AddressOf(t1,'zside'));
               thisTree.SetBranchAddress( 'sc',AddressOf(t1,'sc'));
               thisTree.SetBranchAddress( 'isc',AddressOf(t1,'isc'));
               thisTree.SetBranchAddress( 'ic',AddressOf(t1,'ic'));
               thisTree.SetBranchAddress( 'iquadrant',AddressOf(t1,'iquadrant'));
               thisTree.SetBranchAddress( 'hashedIndex',AddressOf(t1,'hashedIndex'));
               thisTree.SetBranchAddress( 'iter',AddressOf(t1,'iter'));
               thisTree.SetBranchAddress( 'coeff',AddressOf(t1,'coeff'));
               thisTree.SetBranchAddress( 'Chisqu',AddressOf(t1,'Chisqu'));
               thisTree.SetBranchAddress( 'Ndof',AddressOf(t1,'Ndof'));
               thisTree.SetBranchAddress( 'fit_mean',AddressOf(t1,'fit_mean'));
               thisTree.SetBranchAddress( 'fit_mean_err',AddressOf(t1,'fit_mean_err'));
               thisTree.SetBranchAddress( 'fit_sigma',AddressOf(t1,'fit_sigma'));
               if not isEoverEtrue:
                   thisTree.SetBranchAddress( 'Signal',AddressOf(t1,'Signal'));
                   thisTree.SetBranchAddress( 'Backgr',AddressOf(t1,'Backgr'));
                   thisTree.SetBranchAddress( 'fit_Snorm',AddressOf(t1,'fit_Snorm'));
                   thisTree.SetBranchAddress( 'fit_b0',AddressOf(t1,'fit_b0'));
                   thisTree.SetBranchAddress( 'fit_b1',AddressOf(t1,'fit_b1'));
                   thisTree.SetBranchAddress( 'fit_b2',AddressOf(t1,'fit_b2'));
                   thisTree.SetBranchAddress( 'fit_b3',AddressOf(t1,'fit_b3'));
                   thisTree.SetBranchAddress( 'fit_Bnorm',AddressOf(t1,'fit_Bnorm'));
               for ntre in range(thisTree.GetEntries()):
                   thisTree.GetEntry(ntre);
                   if (ntre>=init and ntre<=finit):
                       t.ix_ = t1.ix
                       t.iy_ = t1.iy
                       t.zside_ = t1.zside
                       t.sc_ = t1.sc
                       t.isc_ = t1.isc
                       t.ic = t1.ic
                       t.iquadrant_ = t1.iquadrant
                       t.hashedIndex_ = t1.hashedIndex
                       t.iter_ = t1.iter
                       t.coeff_ = t1.coeff
                       t.Chisqu_ = t1.Chisqu
                       t.Ndof_ = t1.Ndof
                       t.fit_mean_ = t1.fit_mean
                       t.fit_mean_err_ = t1.fit_mean_err
                       t.fit_sigma_ = t1.fit_sigma
                       if not isEoverEtrue:
                           t.Signal_ = t1.Signal
                           t.Backgr_ = t1.Backgr
                           t.fit_Snorm_ = t1.fit_Snorm
                           t.fit_b0_ = t1.fit_b0
                           t.fit_b1_ = t1.fit_b1
                           t.fit_b2_ = t1.fit_b2
                           t.fit_b3_ = t1.fit_b3
                           t.fit_Bnorm_ = t1.fit_Bnorm
                       TreeEE.Fill()
            #TH2
            if isEoverEtrue and n_repeat == 1:
                thisHistoEB = thisfile_f.Get("calibMap_EB_g2")
                thisHistoEEm = thisfile_f.Get("calibMap_EEm_g2")
                thisHistoEEp = thisfile_f.Get("calibMap_EEp_g2")
            else:
                thisHistoEB = thisfile_f.Get("calibMap_EB")
                thisHistoEEm = thisfile_f.Get("calibMap_EEm")
                thisHistoEEp = thisfile_f.Get("calibMap_EEp")
            if EEoEB == 0:
               MaxEta = 85
               Init = int(init)
               Fin = int(finit+1)
               for nFitB in range(Init,Fin):
                   if foldInSuperModule:  
                       # in this case Init,Fin span a region within 0-1700, 
                       # and the map is basically in a single SM (but repreated in all the others as well)
                       if nFitB < 1700:
                           for ism in range(1,37):
                               myRechit = EBDetId(ism, nFitB+1, 1)
                               bin_x = myRechit.ieta()+MaxEta+1
                               bin_y = myRechit.iphi()
                               value = thisHistoEB.GetBinContent(bin_x,bin_y)
                               calibMap_EB.SetBinContent(bin_x,bin_y,value)
                   else:
                       if nFitB < 61200:
                           myRechit = EBDetId( EBDetId.detIdFromDenseIndex(nFitB) )
                           bin_x = myRechit.ieta()+MaxEta+1
                           bin_y = myRechit.iphi()
                           value = thisHistoEB.GetBinContent(bin_x,bin_y)
                           calibMap_EB.SetBinContent(bin_x,bin_y,value)
            else :
               Init1 = int(init)
               Fin1 = int(finit+1)
               for nFitE in range(Init1,Fin1):
                   if nFitE < 14648:
                      myRechitE = EEDetId( EEDetId.detIdFromDenseIndex(nFitE) )
                      if myRechitE.zside() < 0 :
                         value = thisHistoEEm.GetBinContent(myRechitE.ix(),myRechitE.iy())
                         calibMap_EEm.SetBinContent(myRechitE.ix(),myRechitE.iy(),value)
                      if myRechitE.zside() > 0 :
                         value = thisHistoEEp.GetBinContent(myRechitE.ix(),myRechitE.iy())
                         calibMap_EEp.SetBinContent(myRechitE.ix(),myRechitE.iy(),value)

            thisfile_f.Close()

        # write objects to file before going to next objects
        f.cd()
        f.Write()


##########################################
##########################################


    # f.cd()
    # f.Write()
    f.Close()

    print "Done with iteration " + str(iters)
    if( ONLYHADD or ONLYFINHADD or ONLYFIT or ONLYMERGEFIT):
       mode = "BATCH_RESU"
       ONLYHADD = False; ONLYFINHADD = False; ONLYFIT=False; ONLYMERGEFIT=False

print "---THE END---"
