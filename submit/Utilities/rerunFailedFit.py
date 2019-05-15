#!/usr/bin/env python 

import subprocess, time, sys, os, string
from ROOT import *

# example:
# python Utilities/rerunFailedFill.py -e /eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/ -d AlCaP0_AllRun2017_condor_fixEBm16 -i 5 --useLSF -q cmscaf1nd --remove-zombie -p --detector EB
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
parser.add_option(       "--remove-zombie", dest="removeZombie",  action="store_true", default=False, help="Remove zombie (or recovered) files before submitting new jobs")
parser.add_option(       "--check-zombie", dest="checkZombie",  action="store_true", default=False, help="Only check for zombies (or recovered keys) and exit")
parser.add_option(       "--remove-allBad-exit", dest="removeAllBadExit",  action="store_true", default=False, help="Check for zombies or bad files, remove them and exit")
parser.add_option(       "--detector", dest="detector",  type="string", default="all", help="Specify whether to check only EB or only EE or all (default)")
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
ntotEB = 0 # count the value below
ntotEE = 0

mypath = os.environ['CMSSW_BASE'] + "/src/CalibCode/submit/"
jobdir = "{p}{dn}/src/Fit/".format(p=mypath,dn=dirname,i=niter)
logdir = "{p}{dn}/log/rerunFailedFit/iter_{i}/".format(p=mypath,dn=dirname,i=niter)
if not os.path.exists(logdir):
    os.makedirs(logdir)

if not eosdir.endswith("/"): eosdir += "/"
fulleosfolder = "{eos}{dn}/iter_{i}/".format(eos=eosdir,dn=dirname,i=niter)

detToCheck = ["Barrel", "Endcap"]
if options.detector == "EB": detToCheck = ["Barrel"]
elif options.detector == "EE": detToCheck = ["Endcap"]
files = [os.path.join(fulleosfolder, f) for f in os.listdir(fulleosfolder) if ((os.path.isfile(os.path.join(fulleosfolder, f)) and "calibMap" in f and any (det in f for det in detToCheck)))]

detIdMatch = ""
if any(options.detector == x for x in ["EB","EE"]):
    detIdMatch = options.detector

for f in os.listdir(jobdir):
    if not f.endswith('.sh'): continue
    if not f.startswith('submit_'): continue
    #if detIdMatch != "" and not detIdMatch in f: continue
    #else:
    if "EB" in f:
        if options.detector == "EE": continue
        else: ntotEB += 1
    else:
        if options.detector == "EB": continue
        else: ntotEE += 1
print "ntotEB = %d %s" % (ntotEB, "(skipping EB)" if options.detector == "EE" else "")
print "ntotEE = %d %s" % (ntotEE, "(skipping EE)" if options.detector == "EB" else "") 

ntot = ntotEB + ntotEE

isGoodFile = {}
for i in range(ntotEB):
    isGoodFile["EB_"+str(i)] = False
for i in range(ntotEE):
    isGoodFile["EE_"+str(i)] = False

count = 0
goodfiles = []
zombiefiles = []
recoveredfiles = []
#nZombie = 0
for f in sorted(files):
    #base = os.path.basename(f)
    sys.stdout.write('File {num}/{tot}   \r'.format(num=count,tot=ntot-1))
    sys.stdout.flush()
    count += 1

    if os.path.getsize > 100000:  # expect about 500 kB, so ask at least 100k
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
            # name is like pi0CC_2018_EoverEtrue_foldSM_nFit50_onlyEB_Barrel_10_calibMap.root, need to get Barrel_10 (or Endcap)
            nfile = (base.split(dirname+"_")[1]).split('_calibMap.root')[0]  
            isGoodFile[nfile.replace("Barrel","EB").replace("Endcap","EE")] = True
        else:
            recoveredfiles.append(f)
        tf.Close()

print "I see {n} good calibMap files".format(n=len(goodfiles))
#print "There were {n} zombie calibMap files {text}".format(n=nZombie,text= "(removed)" if options.removeZombie else "(to be removed)")
print "There were {n} zombie calibMap files".format(n=len(zombiefiles))
print "There were {n} recovered calibMap files".format(n=len(recoveredfiles))
if options.checkZombie:
    print "Printing list of zombies"
    for f in sorted(zombiefiles):
        print f
    for f in sorted(recoveredfiles):
        print f
    quit()

if options.removeAllBadExit:
    print "### Removing all bad files (zombies and recovered ones)"
    for f in sorted(zombiefiles):
        cmd = "rm " + f
        if options.pretend:            
            print cmd
        else:
            os.system(cmd)
    for f in sorted(recoveredfiles):
        cmd = "rm " + f
        if options.pretend:            
            print cmd
        else:
            os.system(cmd)
    quit()


if options.removeZombie:
    print "### Removing zombies"
    for f in sorted(zombiefiles):
        cmd = "rm " + f
        if options.pretend:            
            print cmd
        else:
            os.system(cmd)
    print "### Removing recovered keys"
    for f in sorted(recoveredfiles):
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
    if not f.startswith('submit_'): continue
    # name is like submit_EB_10_iter_0.sh, need to take EB_10
    jobN = (f.split("_iter")[0]).split("submit_")[-1] 
    if "EB" in jobN and options.detector == "EE": continue
    if "EE" in jobN and options.detector == "EB": continue
    if "justFoldSM" in jobN: continue
    if not isGoodFile[jobN]:
        cmd = "bsub -q {q} -oo {ld}/{jn}.log {jd}{job}".format(q=options.queue, ld=logdir, jn=jobN, jd=jobdir, job=f)
        if options.pretend:
            #pass
            print cmd
        else:
            os.system(cmd)

print "Submitted {n} jobs".format(n=nJobToRun)
print ""
print ""
