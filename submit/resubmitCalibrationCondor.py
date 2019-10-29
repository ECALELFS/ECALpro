#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *
from datetime import datetime

from optparse import OptionParser                                                                      
                                       
                                                                                          
parser = OptionParser(usage="%prog [options]")    
parser.add_option("-l", "--daemon-local",     dest="daemonLocal", action="store_true", default=True, help="Do not submit a job to manage the daemon, do it locally")
parser.add_option(      "--recover-fill",     dest="recoverFill", action="store_true", default=False, help="When resubmitting calibration from hadd, first try to recover failed fills")
parser.add_option("-t", "--token-file", dest="tokenFile",  type="string", default="", help="File needed to renew token (when daemon running locally)")
parser.add_option("--min-efficiency-recover-fill",   dest="minEfficiencyToRecoverFill",   type="float", default=0.97, help="Tolerance of EcalNtp loss. Require fraction of good EcalNtp abive this number to skip recover");
(options, args) = parser.parse_args()

if len(args) != 7:
    print str(len(args)) + "is a wrong number of arguments (" + str(len(args)) +" given, while it should be 7)."
    print "./resubmitCalibrationCondor.py iteration_to_resume isSystematicError(0,1,2) JustHADD(True,False) JustFINALHADD(True,False) JustFIT(True,False) JustMergeFIT(True,False) nJobs+1(goes from j=0 to n<YOUR_Number)"
    print "   where 0=no syst, just normal calib, 1 only even events, 2 odd events(just for last iter)"
    sys.exit(1)

iteration_to_resume = int(args[0])
SystParam           = int(args[1])
onlyHadd            = str(args[2])
onlyFinalHadd       = str(args[3]) 
OnlyFIT             = str(args[4])
OnlyMergeFIT        = str(args[5])
nJobs               = str(args[6])
pwd                 = os.getcwd()

#print str(args[6])
#quit()

workdir = pwd+'/'+dirname
condordir = workdir + "/condor_files/"

Mode = "BATCH_RESU"
if SystParam != 0 :
    Mode = Mode + "_SYST_" + str(SystParam)
if onlyHadd=="True" :
    Mode = Mode + "_ONLYHADD"
if onlyFinalHadd=="True" :
    Mode = Mode + "_ONLYFINALHADD"
if OnlyFIT=="True" :
    Mode = Mode + "_ONLYFIT"
if OnlyMergeFIT=="True" :
    Mode = Mode + "_ONLYMERGEFIT"

### setting environment
env_script_n = workdir + "/resubmit.sh"
env_script_f = open(env_script_n, 'w')
env_script_f.write("#!/bin/bash\n")
env_script_f.write("cd " + pwd + "\n")
env_script_f.write("ulimit -c 0\n")
env_script_f.write("eval `scramv1 runtime -sh`\n")
pycmd =  "python " + pwd + "/calibJobHandlerCondor.py " + Mode + " " + str(iteration_to_resume) + " " + queue + " " + str(nJobs)
if options.recoverFill: pycmd += " --recover-fill "
if options.daemonLocal: pycmd += " --daemon-local "
if options.tokenFile: pycmd += " --token-file {tf}".format(tf=options.tokenFile)
if options.minEfficiencyToRecoverFill >= 0.0:
        pycmd += (" --min-efficiency-recover-fill " + str(options.minEfficiencyToRecoverFill))

print pycmd
env_script_f.write(pycmd + "\n")
env_script_f.write("rm -rf " + pwd + "/core.*\n")
env_script_f.close()

# dummy_exec = open(condordir+'/dummy_exec_daemon.sh','w')
# dummy_exec.write('#!/bin/bash\n')
# dummy_exec.write('bash $*\n')
# dummy_exec.close()
dummy_exec_name = condordir+'/dummy_exec_daemon.sh'

condor_file_name = condordir+'/condor_resubmit_daemon.condor'
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
+JobBatchName = "ecalpro_daemon"
'''.format(de=os.path.abspath(dummy_exec_name), ld=os.path.abspath(condordir), here=os.environ['PWD'] ) )
if os.environ['USER'] in ['mciprian']:
    # mydate = datetime.today()
    # month = int(mydate.month)
    # year  = int(mydate.year)
    # if month == 10 and year == 2019:
    #     condor_file.write('+AccountingGroup = "group_u_CMS.u_zh.priority"\n\n')
    # else:
    #     condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
    condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
else:
    condor_file.write('\n')


condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(env_script_n)))
condor_file.close()


# configuring calibration handler
submit_s = 'condor_submit {cdf} '.format(cdf=condor_file_name)

if  options.daemonLocal:
    print "[resubmit]  '-- source " + os.path.abspath(env_script_n)
    os.system("source " + os.path.abspath(env_script_n))
else:        
    print "[resubmit] Resubmitting calibration handler"
    print "[resubmit]  '-- " + submit_s
    submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
    output = (submitJobs.communicate()[0]).splitlines()
    print "[resubmit]  '-- " + output[0]


# for interactive debugging
#os.system("bash "+env_script_n)
