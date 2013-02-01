#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

if len(sys.argv) != 5:
    print "./resubmitCalibration.py iteration_to_resume isSystematicError(0,1,2) onlyFIT(True,False) onlyFinalHadd(True,False)"
    print "   where 0=no syst, just normal calib, 1 only even events, 2 odd events(just for last iter)"
    sys.exit(1)

iteration_to_resume = int(sys.argv[1])
SystParam           = int(sys.argv[2])
OnlyFIT             = str(sys.argv[3])
onlyFinalHadd       = str(sys.argv[4])
pwd                 = os.getcwd()


#-------- check if you have right access to queues --------#
checkAccessToQueues = subprocess.Popen(['bjobs'], stderr=subprocess.PIPE, shell=True);
output = checkAccessToQueues.communicate()[1]
if(output.find('command not found')==-1):
    print "[resubmit] Correct setup for batch submission"
else:
    print "[resubmit] Missing access to queues"
    sys.exit(1)

workdir = pwd+'/'+dirname

### setting environment
env_script_n = workdir + "/resubmit.sh"
env_script_f = open(env_script_n, 'w')
env_script_f.write("#!/bin/bash\n")
env_script_f.write("cd " + pwd + "\n")
if(is2012):
   env_script_f.write("export SCRAM_ARCH=slc5_amd64_gcc462\n")
else:
   env_script_f.write("export SCRAM_ARCH=slc5_amd64_gcc434\n")
env_script_f.write("eval `scramv1 runtime -sh`\n")
#f=os.popen("wc " + str(inputlist_n) + " | awk '{print $1}'" )
env_script_f.write(pwd + "/calibJobHandler-resubmission.py " + pwd + " " + queue + " " + str(iteration_to_resume) + " " + str(SystParam) + " " + str(OnlyFIT) + " " + str(onlyFinalHadd) + "\n")
#env_script_f.write(pwd + "/calibJobHandler-resubmission.py " + pwd + " " + queue + " " + str(iteration_to_resume) + " " + str(f.readline())  + "\n")
env_script_f.close()

# configuring calibration handler
print "[resubmit] Iteration to resume = " + str(iteration_to_resume)
print "[resubmit] Submitting calibration handler"
submit_s = "bsub -q " + queueForDaemon + " -o " + workdir + "/resume-calibration.log source " + env_script_n
print "[resubmit]  '-- " + submit_s

# submitting calibration handler
submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
output = (submitJobs.communicate()[0]).splitlines()
print "[resubmit]  '-- " + output[0]
