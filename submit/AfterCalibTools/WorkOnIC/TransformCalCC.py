#!/usr/bin/env python

import ROOT
#from sys import argv
from math import fabs
from string import split
import operator
import xml.dom.minidom
import subprocess, time, sys, os
from ROOT import *
from PhysicsTools.PythonAnalysis import *

gSystem.Load("libFWCoreFWLite.so")
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(11111)

def usage():
    print "Usage: -----> python TransformCalCC.py fileCal.txt isEB output.txt (isEB=True/False)"
#python TransformCalCC.py /afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/FT_R_53_V21/Caltech/pi_EB.txt 
#True Compare_2012C_EcalProWithCalCC_Caltech/Cal/outputEB.txt
###START
print "-----------------START-------------------"
if (len(sys.argv) != 4):
    usage()
    sys.exit(1)
try:
    #Input
    file1 = sys.argv[1]
    isEEB = sys.argv[2]
    output = sys.argv[3]
    print "Comparing: " +str(file1)
    print "isEB: " +str(isEEB)
    print "Saving into: " +str(output)
    file1_r = open(file1,'r')
    file1_v = file1_r.readlines()
    file1_v1 = file1_v[:]
    outputfile = open( output, 'w' )
except:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

#Read IC
IC_EB_1=list()
IC_EE_1=list()
print 'Starting reading files!'
for nXtal in range(len(file1_v)):
    List1=file1_v1[nXtal].split(' ')
    if(str(isEEB)=='True'): #iEta iPhi iC
       if(float(List1[2])!=-1):
          if(int(List1[0])<85):
             outputfile.write( str(int(List1[0])-85 ) + " " + str(int(List1[1])+1) + " 0 " + str(float(List1[2]) ) + " 0\n")
          else:
             outputfile.write( str(int(List1[0])-84 ) + " " + str(int(List1[1])+1) + " 0 " + str(float(List1[2]) ) + " 0\n")
       else:
          if(int(List1[0])<85):
             outputfile.write( str(int(List1[0])-85 ) + " " + str(int(List1[1])+1) + " 0 " + str(float(List1[2]) ) + " 999.\n")
          else:
             outputfile.write( str(int(List1[0])-84 ) + " " + str(int(List1[1])+1) + " 0 " + str(float(List1[2]) ) + " 999.\n")
    else:      #iZ iX iY iC
       if(int(List1[0])==0):
          if(float(List1[3])!=-1):
             outputfile.write( str(int(List1[1]) ) + " " + str(int(List1[2])) + " -1 " + str(float(List1[3]) ) + " 0\n")
          else:
             outputfile.write( str(int(List1[1]) ) + " " + str(int(List1[2])) + " -1 " + str(float(List1[3]) ) + " 999.\n")
       else:
          if(float(List1[3])!=-1):
             outputfile.write( str(int(List1[1]) ) + " " + str(int(List1[2])) + " 1 " + str(float(List1[3]) ) + " 0\n")
          else:
             outputfile.write( str(int(List1[1]) ) + " " + str(int(List1[2])) + " 1 " + str(float(List1[3]) ) + " 999.\n")
#OutPut
outputfile.close()
