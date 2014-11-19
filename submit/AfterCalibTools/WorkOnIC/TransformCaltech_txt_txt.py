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
    print "Usage: -----> python TransformCaltech_txt_txt.py CaltechOld.txt CaltechNew.txt"
#python TransformCaltech_txt_txt.py True Calthech_2012C_NEWTAG/pi_EB.txt Calthech_2012C_NEWTAG/pi_EB_new.txt
#python TransformCaltech_txt_txt.py False Calthech_2012C_NEWTAG/pi_EE.txt Calthech_2012C_NEWTAG/pi_EE_new.txt

###START
print "-----------------START-------------------"
if (len(sys.argv) != 4):
    usage()
    sys.exit(1)
#Input
isEB  = str(sys.argv[1])
file1 = sys.argv[2]
file2 = sys.argv[3]
print "Transforming: " +str(file1)
print "in: " +str(file2)
if(isEB=='True'):
   print "is Barrel!"
else:
   print "is Endcap!"
file1_r = open(file1,'r')
file2_r = open( file2, 'w' )
file1_v = file1_r.readlines()
file1_v1 = file1_v[:]
#Reading Files
print 'Starting reading files!'
for nXtal in range(len(file1_v)):
    List1=file1_v1[nXtal].split(' ')
    if(isEB=='True'):
       iEta=int(List1[0])
       iPhi=int(List1[1])
       iC = float(List1[2]) #eta Phi IC
       print str(iEta) + " " + str(iPhi) + " " + str(iC)
       iPhi=iPhi+1
       if(iEta<=84):
          iEta=iEta-85
       else:
          iEta=iEta-84
       print str(iEta) + " " + str(iPhi) + " " + str(iC)
       if(iC!=-1.):
          file2_r.write( str(iEta) + " " + str(iPhi) + " 0 " + str(iC) + " 0\n" )
       else:
          file2_r.write( str(iEta) + " " + str(iPhi) + " 0 -1. 999.\n" )
    else:
       iZ = int(List1[0]);
       iX = int(List1[1]);
       iY = int(List1[2]);
       iC = float(List1[3]);
       print str(iX) + " " + str(iY) + " " + str(iZ) + " " + str(iC)
       if(iZ==0):
          iZ=-1
       else:
          iZ=1
       print str(iX) + " " + str(iY) + " " + str(iZ) + " " + str(iC)
       if(iC!=-1.):
          file2_r.write( str(iX) + " " + str(iY) + " " + str(iZ) + " " + str(iC) + " 0\n" )
       else:
          file2_r.write( str(iX) + " " + str(iY) + " " + str(iZ) + " -1. 999.\n" )
file2_r.close()

print "End pf the Job"
