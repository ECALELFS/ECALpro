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
parser.add_option("-i", "--iteration", dest="iteration",  type="int",     default=0,   help="Iteration to start from, usually 0 unless resubmitting jobs")
parser.add_option("-n", "--njobs", dest="njobs",  type="int",     default=0,   help="Number of jobs")
#parser.add_option("-q", "--queue", dest="queue",  type="string",     default="",   help="Queue for running jobs")
parser.add_option("-r", "--run", dest="run",  type="string",     default="",   help="Specify where to start from when resubmitting jobs [hadd,finalhadd,fit,mergefit]")
parser.add_option("-s", "--syst", dest="syst",  type="int",     default=0,   help="To run on all events (0, default), only even (1) or odd (2) events ")
(options, args) = parser.parse_args()

if not options.run:
    print("Must specify option -r [hadd,finalhadd,fit,mergefit] when using --resubmit. Abort")
    sys.exit(1)
if options.run not in ["hadd","finalhadd","fit","mergefit"]:
    print("Option -r requires one of [hadd,finalhadd,fit,mergefit], while '%s' was given. Abort" % options.run)

pwd                 = os.getcwd()
workdir = pwd+'/'+dirname
condordir = workdir + "/condor_files/"

### setting environment
env_script_n = workdir + "/resubmit.sh"
env_script_f = open(env_script_n, 'w')
env_script_f.write("#!/bin/bash\n")
env_script_f.write("cd " + pwd + "\n")
env_script_f.write("ulimit -c 0\n")
env_script_f.write("eval `scramv1 runtime -sh`\n")
pycmd =  "python3 " + pwd + "/calibJobHandlerCondor.py --resubmit -i {i} -r {r} -n {n} -s {s} ".format(i=options.iteration,n=options.njobs,r=options.run,s=options.syst)
if options.recoverFill: pycmd += " --recover-fill "
if options.daemonLocal: pycmd += " --daemon-local "
if options.tokenFile: pycmd += " --token-file {tf}".format(tf=options.tokenFile)
if options.minEfficiencyToRecoverFill >= 0.0:
        pycmd += (" --min-efficiency-recover-fill " + str(options.minEfficiencyToRecoverFill))

print(pycmd)
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
condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = True
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
    condor_file.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')
else:
    condor_file.write('\n')


condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(env_script_n)))
condor_file.close()


# configuring calibration handler
submit_s = 'condor_submit {cdf} '.format(cdf=condor_file_name)

if  options.daemonLocal:
    print("[resubmit]  '-- source " + os.path.abspath(env_script_n))
    os.system("source " + os.path.abspath(env_script_n))
else:        
    print("[resubmit] Resubmitting calibration handler")
    print("[resubmit]  '-- " + submit_s)
    submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
    output = (submitJobs.communicate()[0]).splitlines()
    print("[resubmit]  '-- " + output[0])


# for interactive debugging
#os.system("bash "+env_script_n)
