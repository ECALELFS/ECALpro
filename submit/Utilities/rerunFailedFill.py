#!/usr/bin/env python 

import subprocess, time, sys, os, string
from ROOT import *

# example:
# python Utilities/rerunFailedFill.py -e /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/ -d AlCaP0_AllRun2017_condor_fixEBm16 -i 5 --useLSF -q cmscaf1nd --remove-zombie -p
# add -p to check what will be done, without running jobs

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")    
# dummy option
parser.add_option("-e", "--eosdir", dest="eosdir",  type="string", default="", help="Eos path with all calibration folders")
parser.add_option("-d", "--dirname", dest="dirname",  type="string", default="", help="Name of folder on eos (will have to be the same as the local one to retrieve the script to run)")
parser.add_option("-i", "--iter", dest="iter",  type="int", default=-1, help="Number of iteration to match the iter_<n> folder (n >= 0)")
parser.add_option("--useLSF", dest="useLSF",  action="store_true", default=False, help="Run jobs using LSF instead of condor")
parser.add_option("-q", "--queue", dest="queue",  type="string", default="cmscaf1nd", help="Name of queue for job submission (for LSF)")
parser.add_option("-p", "--pretend", dest="pretend",  action="store_true", default=False, help="Just print commands, do not submit jobs")
# not needed, can take it from folder
#parser.add_option("-n", "--n-tot", dest="nTot",  type="int", default=-1, help="Total number of expected fill files")
parser.add_option(       "--remove-zombie", dest="removeZombie",  action="store_true", default=False, help="Remove zombie file before submitting new jobs")
parser.add_option(       "--check-zombie", dest="checkZombie",  action="store_true", default=False, help="Only check for zombies and exit")
(options, args) = parser.parse_args()

if len(options.eosdir) == 0:
    print "Need the path to eos with --eosdir <path>"
    quit()
    
if len(options.dirname) == 0:
    print "Need the name of the folder with --dirname <name>"
    quit()

if options.iter < 0:
    print "Need the iteration number with --iter <number>"
    quit()

# if options.nTot < 0:
#     print "Need the total number of fill files with --n-tot <number>"
#     quit()

if not options.useLSF:
    print "Warning: at the moment only usage with LSF job submission is supported. Use --useLSF"
    quit()
    

hostname = os.environ['HOSTNAME']
if "lxplus" not in hostname:
    print "Warning: you need to be on lxplus to use this script"
    quit()

eosdir = options.eosdir
dirname = options.dirname
niter = int(options.iter)
#ntot = int(options.nTot)
ntot = 0 # count the value below

mypath = os.environ['CMSSW_BASE'] + "/src/CalibCode/submit/"
jobdir = "{p}{dn}/src/Fill/iter_{i}/".format(p=mypath,dn=dirname,i=niter)
logdir = "{p}{dn}/log/rerunFailedFill/iter_{i}/".format(p=mypath,dn=dirname,i=niter)
if not os.path.exists(logdir):
    os.makedirs(logdir)

if not eosdir.endswith("/"): eosdir += "/"
fulleosfolder = "{eos}{dn}/iter_{i}/".format(eos=eosdir,dn=dirname,i=niter)

files = [os.path.join(fulleosfolder, f) for f in os.listdir(fulleosfolder) if (os.path.isfile(os.path.join(fulleosfolder, f)) and "EcalNtp" in f)]

for f in os.listdir(jobdir):
    if not f.endswith('.sh'): continue
    else: ntot += 1
print "ntot = %d" % ntot

isGoodFile = {}
for i in range(ntot):
    isGoodFile[int(i)] = False

count = 0
goodfiles = []
zombiefiles = []
#nZombie = 0
for f in sorted(files):
    #base = os.path.basename(f)
    sys.stdout.write('File {num}/{tot}   \r'.format(num=count,tot=ntot-1))
    sys.stdout.flush()
    count += 1

    if os.path.getsize > 30000000:  # expect about 75 MB, so ask at least 20
        # at this point the file should be good, but let's check if there are no recovered keys                                                           
        #open and check there are no recovered keys: in this case remove these files from the list, otherwise hadd might fail                             
        tf = TFile.Open("root://eoscms/"+f)        
        #if not tf:             
        #    continue
        if not tf or tf.IsZombie(): 
            #nZombie += 1
            zombiefiles.append(f)
            continue
        if not tf.TestBit(TFile.kRecovered):
            goodfiles.append(f)
            base = os.path.basename(f)
            nfile = (base.split(".root")[0]).split('_')[-1]  # name is like AlCaP0_AllRun2017_condor_fixEBm16_EcalNtp_97.root, need to get 97
            isGoodFile[int(nfile)] = True
        tf.Close()

print "I see {n} good EcalNtp files".format(n=len(goodfiles))
#print "There were {n} zombie EcalNtp files {text}".format(n=nZombie,text= "(removed)" if options.removeZombie else "(to be removed)")
print "There were {n} zombie EcalNtp files".format(n=len(zombiefiles))
if options.checkZombie:
    print "Printing list of zombies"
    for f in sorted(zombiefiles):
        print f
    quit()

if options.removeZombie:
    print "### Removing zombies"
    for f in sorted(zombiefiles):
        cmd = "rm " + f
        if options.pretend:            
            print cmd
        else:
            os.system(cmd)

nJobToRun = 0
for key in isGoodFile:
    if not isGoodFile[key]: nJobToRun += 1

for f in os.listdir(jobdir):
    if not f.endswith('.sh'): continue
    jobN = (f.split(".sh")[0]).split("_")[-1] # name is like submit_iter_5_job_14.sh, need to take 14
    if not isGoodFile[int(jobN)]:
        cmd = "bsub -q {q} -oo {ld}/{jn}.log {jd}{job}".format(q=options.queue, ld=logdir, jn=jobN, jd=jobdir, job=f)
        if options.pretend:
            #pass
            print cmd
        else:
            os.system(cmd)

print "Submitted {n} jobs".format(n=nJobToRun)
print ""
print ""
