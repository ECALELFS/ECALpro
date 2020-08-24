#!/usr/bin/env python

## check if a file is good on eos, and delete it if it is not good
# the file name should be of the form /eos/cms/... or root://eoscms//eos/cms
## a file is bad if:
## 1) the size is below a minimum threshold (default is 1 MB)
## 2) the kRecovered bit is true
## 3) it it is a zombie or empty (in which case it would probably fail the size check already)

import subprocess, time, sys, os, string
import ROOT

from optparse import OptionParser                                                    

parser = OptionParser(usage="%prog [options] filename")    
parser.add_option("-d", "--delete", dest="delete",action="store_true", default=False, help="Delete bad files")
parser.add_option("-s", "--size", dest="sizeThreshold", type="int", default=1024, help="Threshold for good size of the file in kB, default is 1024 = 1 MB")
(options, args) = parser.parse_args()

if len(args)<1:
    parser.print_usage()
    quit()

print ""
eosFile = args[0]
eosFileNameToOpen = eosFile
if eosFileNameToOpen.startswith("/eos/cms"):
    eosFileNameToOpen = "root://eoscms/" + eosFile
elif not eosFileNameToOpen.startswith("root://eoscms//eos/cms"):
    print "Error: file name in input not valid, must start with /eos/cms or root://eoscms//eos/cms"
    quit()

isGood = True
sizeThreshold = 1024 * options.sizeThreshold

filesize=0
if os.path.exists(eosFile): 
    filesize = os.path.getsize(eosFile)
if filesize < sizeThreshold:
    print "%%% size is too small"
    isGood = False
else:
    tf = ROOT.TFile.Open(eosFileNameToOpen)
    if not tf or tf.IsZombie(): 
        print "%%% file is zombie"
        isGood = False
    elif tf.TestBit(ROOT.TFile.kRecovered):                    
        print "%%% file was recovered"
        isGood = False
        tf.Close()
                    
if isGood:
    print ">>> Check successful: file is good on EOS"
else:
    print "#### File is bad or non existing."
    if options.delete: 
        print "### Will be deleted if existing"
        if not ROOT.gSystem.AccessPathName(eosFile):
            # file exists, let's delete it
            # gSystem.Exec("rm " + eosFile)  
            ROOT.gSystem.Unlink(eosFile) # this works also for non-Unix systems, just in case

print ""

