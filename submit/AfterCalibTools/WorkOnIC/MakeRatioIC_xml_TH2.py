#!/usr/bin/env python

import ROOT
#from sys import argv
from math import fabs

import operator
import xml.dom.minidom

import subprocess, time, sys, os

from ROOT import *
from PhysicsTools.PythonAnalysis import *
gSystem.Load("libFWCoreFWLite.so")
#AutoLibraryLoader.enable()

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(11111)

def usage():
    print "Usage: -----> python MoltiplicateIC.py pathXml (vecchie IC), pathTH2 (mediate a 1 glob.), pathTH2_w (mediate a 1 per etaRing), OutputName"
# python MakeRatioIC_xml_TH2.py InputFile/EcalIntercalibConstants_V20120620_piZPhiSEtaScale2012_IOV2_AlphaStudies.xml 2012B/IC_SetAverageTo1_mol_global.root 2012B/IC_SetAverageTo1_mol.root
#2012B/Ratio_MineLast_Older.root
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

def calibMapFromXML(fileXml, EBIC, EEpIC, EEmIC, myList, myList_tot, EBIC_w):
   doc=xml.dom.minidom.parse(fileXml)
   for xtal in doc.getElementsByTagName("cell"):
       id=readCellId(xtal)
       for value in xtal.childNodes:
           
           if (value.nodeType == value.ELEMENT_NODE and value.nodeName == "Value" ):
               #print str(id[0]) + " and " + str(id[1]) + " " + str(id[2]) + "( " + str(len(id)) + ")"
               if (id[0]=="EB"):
                   MyEBIC = EBIC.GetBinContent( int(id[1])+86 , int(id[2]) )
                   OLDEBIC = float(value.childNodes[0].nodeValue)
                   OldCoef_EB.Fill(OLDEBIC)
                   #print str(MyEBIC/OLDEBIC)
                   if(MyEBIC!=1):
                       myList_tot[abs(int(id[1])+86)]=myList_tot[abs(int(id[1])+86)]+1.
                       myList[abs(int(id[1])+86)]=myList[abs(int(id[1])+86)]+OLDEBIC
                       Ratio_EB.Fill(OLDEBIC/MyEBIC)
                       Ratio_EBMap.SetBinContent( int(id[2]), int(id[1])+86, OLDEBIC/MyEBIC)
                       t.phi=int(id[2])
                       t.eta=int(id[1]);
                       t.iC_my=float(MyEBIC)
                       t.iC_OLD=(OLDEBIC)
                       mytree.Fill()
               if (id[0]=='EE'):
                   if (str(id[ int(len(id))-1 ])=='-1' ):
                       MyEEmIC = EEmIC.GetBinContent( int(id[1]) , int(id[2]) )
                       OLDEEmIC = float(value.childNodes[0].nodeValue)
                       if(MyEEmIC!=1):
                         Ratio_EEm.Fill(OLDEEmIC/MyEEmIC)
                         Ratio_EEmMap.SetBinContent( int(id[2]), int(id[1]), OLDEEmIC/MyEEmIC)
                   if (str(id[ int(len(id))-1 ])=='1' ):
                       MyEEpIC = EEpIC.GetBinContent( int(id[1]) , int(id[2]) )
                       OLDEEpIC = float(value.childNodes[0].nodeValue)
                       if(MyEEpIC!=1):
                         Ratio_EEp.Fill(OLDEEpIC/MyEEpIC)
                         Ratio_EEpMap.SetBinContent( int(id[2]), int(id[1]), OLDEEpIC/MyEEpIC)
   for a in range(len(myList)):
     if(myList_tot[a]!=0):
       print str(a) + ") " + str(myList[a]) + " & " + str(myList_tot[a]) + " ->" + str(myList[a]/myList_tot[a])
       myList[a]=1./float(myList[a]/myList_tot[a])
   for xtal in doc.getElementsByTagName("cell"):
       id=readCellId(xtal)
       for value in xtal.childNodes:
           if (value.nodeType == value.ELEMENT_NODE and value.nodeName == "Value" ):
               if (id[0]=="EB"):
                   MyEBIC = EBIC_w.GetBinContent( int(id[1])+86 , int(id[2]) )
                   OLDEBIC = float(value.childNodes[0].nodeValue)
                   OldCoef_EB_w.Fill(float(OLDEBIC*myList[abs(int(id[1])+86)]))
                   if(MyEBIC!=1):
                       Ratio_EB_w.Fill((OLDEBIC*myList[abs(int(id[1])+86)])/MyEBIC)
                       Ratio_EBMap_w.SetBinContent( int(id[2]), int(id[1])+86, (OLDEBIC*myList[abs(int(id[1])+86)])/MyEBIC)
                   if(MyEBIC!=1 ):
                     print "Ecco: " + str(id[2]) + "  " + str(int(id[1])+86) + " : " + str(MyEBIC) + "  " + str((OLDEBIC*myList[abs(int(id[1])+86)]))
                     s.phi=int(id[2])
                     s.eta=int(id[1]);
                     s.iC_my=float(MyEBIC)
                     s.iC_OLD=(OLDEBIC*myList[abs(int(id[1])+86)])
                     mytree_w.Fill()

###START
print "START"
if (len(sys.argv) != 5):
    usage()
    sys.exit(1)
try:
    #Input
    pathXml = sys.argv[1]
    pathTH2 = sys.argv[2]
    pathTH2_w = sys.argv[3]
    Output  = sys.argv[4]
    print "Reading OLD IC in: " +str(pathXml)
    print "Reading TH2F (av. 1. globally) in: " +str(pathTH2)
    print "Reading TH2F (av. 1 in etaRing) in: " +str(pathTH2_w)
    fileXml = open(pathXml,'r')
    fileTH2 = ROOT.TFile.Open(pathTH2)
    fileTH2_w = ROOT.TFile.Open(pathTH2_w)
    EBIC    = fileTH2.Get('calibMap_EB')
    EEmIC   = fileTH2.Get('calibMap_EEm')
    EEpIC   = fileTH2.Get('calibMap_EEp')
    EBIC_w    = fileTH2_w.Get('calibMap_EB')
    EEmIC_w    = fileTH2_w.Get('calibMap_EEm')
    EEpIC_w    = fileTH2_w.Get('calibMap_EEp')
except:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)
#outPut
f = TFile(Output, 'recreate')
OldCoef_EB = TH1F("OldCoef_EB", "OldCoeff EB", 100, 0.2, 1.8)
OldCoef_EEm= TH1F("OldCoef_EEm","OldCoeff EEm",100, 0.2, 1.8)
OldCoef_EEp = TH1F("OldCoef_EEp","OldCoeff EEp",100, 0.2, 1.8)
OldCoef_EB_w = TH1F("OldCoef_EB_w", "OldCoeff EB", 100, 0.2, 1.8)

Ratio_EB    = TH1F("Ratio_EB","Ratio EB",100, 0.9, 1.1)
Ratio_EBMap = TH2F("Ratio_EBMap","Ratio EB Map #phi on x, #eta on y", 360,0.5,360.5,171,-85.5,85.5)
Ratio_EEm    = TH1F("Ratio_EEm","Ratio EEm",100, 0.9, 1.1)
Ratio_EEmMap = TH2F("Ratio_EEmMap","Ratio EEm Map #phi on x, #eta on y", 100,0.5,100.5, 100,0.5,100.5)
Ratio_EEp    = TH1F("Ratio_EEp","Ratio EEp",100, 0.9, 1.1)
Ratio_EEpMap = TH2F("Ratio_EEpMap","Ratio EEp Map #phi on x, #eta on y", 100,0.5,100.5, 100,0.5,100.5)
Ratio_EB_w    = TH1F("Ratio_EB_w","Ratio EB",100, 0.9, 1.1)
Ratio_EBMap_w = TH2F("Ratio_EBMap_w","Ratio EB Map #phi on x, #eta on y", 360,0.5,360.5,171,-85.5,85.5)
#Tree
gROOT.ProcessLine(\
      "struct EBStruct{\
        Int_t phi;\
        Int_t eta;\
        Double_t iC_my;\
        Double_t iC_OLD;\
};")
gROOT.ProcessLine(\
      "struct EBStruct_w{\
        Int_t phi;\
        Int_t eta;\
        Double_t iC_my;\
        Double_t iC_OLD;\
};")
t=EBStruct()
s=EBStruct_w()
mytree = ROOT.TTree("IC_golabally1", "IC_golabally1")
mytree.Branch('phi', AddressOf(t,'phi'), 'phi/I')
mytree.Branch('eta', AddressOf(t,'eta'), 'eta/I')
mytree.Branch('iC_my', AddressOf(t,'iC_my'), 'iC_my/D')
mytree.Branch('iC_OLD', AddressOf(t,'iC_OLD'), 'iC_OLD/D')
mytree_w = ROOT.TTree("IC_etaRing1", "IC_etaRing1")
mytree_w.Branch('phi', AddressOf(s,'phi'), 'phi/I')
mytree_w.Branch('eta', AddressOf(s,'eta'), 'eta/I')
mytree_w.Branch('iC_my', AddressOf(s,'iC_my'), 'iC_my/D')
mytree_w.Branch('iC_OLD', AddressOf(s,'iC_OLD'), 'iC_OLD/D')

#Extraction xml
myList=[]
for a1 in xrange(172):
  myList.append([])
  myList[a1]=0.
myList_tot=[]
for a2 in xrange(172):
  myList_tot.append([])
  myList_tot[a2]=0.

calibMapFromXML(fileXml, EBIC, EEpIC, EEmIC, myList, myList_tot, EBIC_w)

print 'Finish...'

#Final Writing
f.cd()
f.Write()
f.Close()
