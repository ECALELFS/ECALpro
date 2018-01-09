#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

if len(sys.argv) != 7:
    print str(sys.argv) + "is a wrong number of srguments (" + str(len(sys.argv)) +" given, while it sould be 7)."
    print "./resubmitCalibration.py iteration_to_resume isSystematicError(0,1,2) JustHADD(True,False) JustFINALHADD(True,False) JustFIT(True,False) nJobs+1(goes from j=0 to n<YOUR_Number)"
    print "   where 0=no syst, just normal calib, 1 only even events, 2 odd events(just for last iter)"
    sys.exit(1)

iteration_to_resume = int(sys.argv[1])
SystParam           = int(sys.argv[2])
onlyHadd            = str(sys.argv[3])
onlyFinalHadd       = str(sys.argv[4]) #To be implemented
OnlyFIT             = str(sys.argv[5])
nJobs               = str(sys.argv[6])
pwd                 = os.getcwd()

workdir = pwd+'/'+dirname

Mode = "BATCH_RESU"
if SystParam != 0 :
    Mode = Mode + "_SYST_" + str(SystParam)
if onlyHadd=="True" :
    Mode = Mode + "_ONLYHADD"
if onlyFinalHadd=="True" :
    Mode = Mode + "_ONLYFINALHADD"
if OnlyFIT=="True" :
    Mode = Mode + "_ONLYFIT"

### setting environment
env_script_n = workdir + "/resubmit.sh"
env_script_f = open(env_script_n, 'w')
env_script_f.write("#!/bin/bash\n")
env_script_f.write("cd " + pwd + "\n")
env_script_f.write("ulimit -c 0\n")
env_script_f.write("eval `scramv1 runtime -sh`\n")
print "python " + pwd + "/calibJobHandler.py " + Mode + " " + str(iteration_to_resume) + " " + queue + " " + str(nJobs)
env_script_f.write("python " + pwd + "/calibJobHandler.py " + Mode + " " + str(iteration_to_resume) + " " + queue + " " + str(nJobs) + "\n")
env_script_f.write("rm -rf " + pwd + "/core.*\n")
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

# for interactive debugging
#os.system("bash "+env_script_n)
