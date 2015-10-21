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
    print "Usage: -----> python MoltiplicateIC.py pathXml, pathTH2, Output"
# python MoltiplicateIC.py InputFile/EcalIntercalibConstants_GR10_H_V6.xml 2010/IC_ResidsualFree.root 2010/ABSIC_ResidsualFree.root
# python MoltiplicateIC.py InputFile/EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies_2012.xml
#                          root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_Neutrino_Pt2to20_AVE40BX25_FoldEtaRing_03/iter_7/calibMap.root
#                          ALL_Neutrino_Pt2to20_AVE40BX25_FoldEtaRing_03/Absolute_IC_AlphaStudies_2012.root
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

def calibMapFromXML(input, EBIC, NewcalibMap_EB, EEpIC, NewcalibMap_EEp, EEmIC, NewcalibMap_EEm):
    doc=xml.dom.minidom.parse(input)
    for xtal in doc.getElementsByTagName("cell"):
        id=readCellId(xtal)
        for value in xtal.childNodes:
            
            if (value.nodeType == value.ELEMENT_NODE and value.nodeName == "Value" ):
                if (id[0]=="EB"):
                    OldEBIC = EBIC.GetBinContent( int(id[1])+86 , int(id[2]) )
                    OldCoef_EB.Fill(float(value.childNodes[0].nodeValue))
                    #print 'IC EB: ' + str(id[1]) + ' ' + str(id[2]) + ' -> ' + str(OldEBIC) + ' * ' + str(float(value.childNodes[0].nodeValue)) + ' = ' + str(OldEBIC * float(value.childNodes[0].nodeValue) )
                    if(OldEBIC!=1): NewcalibMap_EB.SetBinContent(int(id[1])+86 , int(id[2]), float(OldEBIC * float(value.childNodes[0].nodeValue)) )
                    else:           NewcalibMap_EB.SetBinContent(int(id[1])+86 , int(id[2]), float(OldEBIC))
                    if(float(value.childNodes[0].nodeValue)>2 or float(value.childNodes[0].nodeValue)<0.5 ): print 'WARNING EB ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(value.childNodes[0].nodeValue)
                    if( float(OldEBIC)>2 or float(OldEBIC)<0.5 ): print 'WARNING EB OldEBIC ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(OldEBIC)
                if (id[0]=='EE'):
                    if (str(id[ int(len(id))-1 ])=='-1' ):
                        OldEEmIC = EEmIC.GetBinContent( int(id[1]) , int(id[2]) )
                        OldCoef_EEm.Fill(float(value.childNodes[0].nodeValue))
                        if(OldEEmIC!=1): NewcalibMap_EEm.SetBinContent(int(id[1]) , int(id[2]), float(OldEEmIC * float(value.childNodes[0].nodeValue)) )
                        else:            NewcalibMap_EEm.SetBinContent(int(id[1]) , int(id[2]), float(OldEEmIC))
                        if(float(value.childNodes[0].nodeValue)>2 or float(value.childNodes[0].nodeValue)<0.5 ): print 'WARNING EEm ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(value.childNodes[0].nodeValue)
                        if( float(OldEEmIC)>2 or float(OldEEmIC)<0.5 ): print 'WARNING EEm OldEEmIC ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(OldEEmIC)
                    if (str(id[ int(len(id))-1 ])=='1' ):
                        OldEEpIC = EEpIC.GetBinContent( int(id[1]) , int(id[2]) )
                        OldCoef_EEp.Fill(float(value.childNodes[0].nodeValue))
                        if(OldEEpIC!=1): NewcalibMap_EEp.SetBinContent(int(id[1]) , int(id[2]), float(OldEEpIC * float(value.childNodes[0].nodeValue)) )
                        else:            NewcalibMap_EEp.SetBinContent(int(id[1]) , int(id[2]), float(OldEEpIC))
                        if( float(value.childNodes[0].nodeValue)>2 or float(value.childNodes[0].nodeValue)<0.5 ): print 'WARNING EEp ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(value.childNodes[0].nodeValue)
                        if( float(OldEEpIC)>2 or float(OldEEpIC)<0.5 ): print 'WARNING EEp OldEEpIC ' + str(id[1]) + ' ' + str(id[2]) + ' ' + str(OldEEpIC)

###START
print "START"
if (len(sys.argv) != 4):
    usage()
    sys.exit(1)
try:
    #Input
    pathXml = sys.argv[1]
    pathTH2 = sys.argv[2]
    Output  = sys.argv[3]
    print "Reading file " +str(pathXml)
    print "Reading TH2F in " +str(pathTH2)
    fileXml = open(pathXml,'r')
    fileTH2 = ROOT.TFile.Open(pathTH2)
    EBIC    = fileTH2.Get('calibMap_EB')
    EEmIC   = fileTH2.Get('calibMap_EEm')
    EEpIC   = fileTH2.Get('calibMap_EEp')
except:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)
#outPut
f = TFile(Output, 'recreate')
NewcalibMap_EB = TH2F("calibMap_EB", "EB calib coefficients: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
OldCoef_EB = TH1F("OldCoef_EB", "OldCoeff EB", 100, 0.2, 1.8)
NewcalibMap_EEm = TH2F("calibMap_EEm", "EE- calib coefficients", 100,0.5,100.5,100,0.5,100.5)
OldCoef_EEm= TH1F("OldCoef_EEm","OldCoeff EEm",100, 0.2, 1.8)
NewcalibMap_EEp = TH2F("calibMap_EEp", "EE+ calib coefficients", 100,0.5,100.5,100,0.5,100.5)
OldCoef_EEp= TH1F("OldCoef_EEp","OldCoeff EEp",100, 0.2, 1.8)
#Extraction xml
calibMapFromXML(fileXml, EBIC, NewcalibMap_EB, EEpIC, NewcalibMap_EEp, EEmIC, NewcalibMap_EEm)

print 'Finish...'

#Final Writing
f.cd()
f.Write()
f.Close()
