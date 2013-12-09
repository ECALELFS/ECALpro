#!/usr/bin/env python

#from sys import argv
from math import fabs

import operator
import xml.dom.minidom

import subprocess, time, sys, os

from ROOT import gROOT, gStyle, gSystem, TCanvas, TH1F, TH2F, TFile
from PhysicsTools.PythonAnalysis import *
gSystem.Load("libFWCoreFWLite.so")
#AutoLibraryLoader.enable()

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(11111)

def usage():
    print "Usage: -----> python MoltiplicateIC_xml_txt.py pathXml pathTXT Output"
# python MoltiplicateIC_xml_txt.py InputFile/EcalIntercalibConstants_GR10_H_V6.xml 2010/IC.txt 2010/ABSIC_ResidsualFree.root
def readCellId(node):
    ieta =-1
    iphi =-1
    ix   =-1
    iy   =-1
    ixSC =-1
    iySC =-1
    zside=-1

    ieta_str =str(node.getAttribute("iEta"))
    iphi_str =str(node.getAttribute("iPhi"))
    ix_str   =str(node.getAttribute("ix"))
    iy_str   =str(node.getAttribute("iy"))
    zside_str=str(node.getAttribute("zside"))

    if ( not str(ieta_str)== "" and  not str(iphi_str)== ""):
        return ("EB",str(ieta_str),str(iphi_str))
    elif ( not str(ix_str)== "" and  not str(iy_str)== ""):
        return ("EE",str(ix_str),str(iy_str),str(zside_str))
#Division
def calibMapFromXML(input):
    doc=xml.dom.minidom.parse(input)
    for xtal in doc.getElementsByTagName("cell"):
        id=readCellId(xtal)
        for value in xtal.childNodes:  
            if (value.nodeType == value.ELEMENT_NODE and value.nodeName == "Value" ):
                if (id[0]=="EB"):
                    Founded=False
                    for find_EB1 in range(len(IC_EB_1)):
                        if(int(IC_EB_1[find_EB1][0])==int(id[1]) and int(IC_EB_1[find_EB1][1])==int(id[2]) ):
                           Founded=True
                           NewIcFile.write( str(IC_EB_1[find_EB1][0]) + ' ' + str(IC_EB_1[find_EB1][1]) + ' 0 ' + str(float(IC_EB_1[find_EB1][2])*float(value.childNodes[0].nodeValue)) + ' 0\n'  )
                           print "Found EB: " + str(IC_EB_1[find_EB1][0]) + ' ' + str(IC_EB_1[find_EB1][1]) + ' : ' + str(float(IC_EB_1[find_EB1][2])*float(value.childNodes[0].nodeValue)) + ' ('+str(IC_EB_1[find_EB1][2])+'/'+str(value.childNodes[0].nodeValue)
                           break
                    if not(Founded): print "WARNING!! IC not founded, EB " + str(id[1]) + " " + str(id[2]) 
                if (id[0]=='EE'):
                    Founded=False
                    for find_EE1 in range(len(IC_EE_1)):
                        if(int(IC_EE_1[find_EE1][0])==int(id[1]) and int(IC_EE_1[find_EE1][1])==int(id[2]) and str(IC_EE_1[find_EE1][2])==str(id[ int(len(id))-1 ])):
                           Founded=True
                           NewIcFile.write( str(IC_EE_1[find_EE1][0]) + ' ' + str(IC_EE_1[find_EE1][1]) + ' ' + str(IC_EE_1[find_EE1][2]) + ' ' + str(float(IC_EE_1[find_EE1][3])*float(value.childNodes[0].nodeValue)) + ' 0\n' )
                           print "Found EE: " + str(IC_EE_1[find_EE1][0]) + ' ' + str(IC_EE_1[find_EE1][1]) + ' ' + str(IC_EE_1[find_EE1][2]) + ' : ' + str(float(IC_EE_1[find_EE1][3])*float(value.childNodes[0].nodeValue)) + ' ('+str(IC_EE_1[find_EE1][3])+'/'+str(value.childNodes[0].nodeValue)
                           break
                    if not(Founded): print "WARNING!! IC not founded, EE " + str(id[1]) + " " + str(id[2]) + " " + str(id[ int(len(id))-1 ])

###START
print "START"
if (len(sys.argv) != 4):
    usage()
    sys.exit(1)
try:
    #Input
    pathXml = sys.argv[1]
    pathTXT = sys.argv[2]
    Output  = sys.argv[3]
    print "Dividing " +str(pathTXT)
    print "For " +str(pathXml)
    fileXml = open(pathXml,'r')
    file1_r = open(pathTXT,'r')
    file1_v = file1_r.readlines()
    file1_v1 = file1_v[:]
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
    if(nXtal<61200):
       IC_EB_1.append([List1[0],List1[1],List1[3]]) #eta Phi IC
    else:
       IC_EE_1.append([List1[0],List1[1],List1[2],List1[3]]) #iX iY iZ IC
#Output
NewIcFile = open(str(Output), "w")
#Divide IC
calibMapFromXML(fileXml)
print 'Finish...'
#Final Writing
NewIcFile.close()
