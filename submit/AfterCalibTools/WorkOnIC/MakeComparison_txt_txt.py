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
    print "Usage: -----> python MakeComparison_txt_txt.py fileMINE.txt fileCALT.txt output"
#python MakeComparison_txt_txt.py /afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/lpernie/IC_2012C_MVAcorrResid_13Iter_MediatedGlobally.txt
#/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/2012C_NewTags_V2GOOD/interCalibConstants.combinedPi0EtaAllPeriod.EBandEE_2012C_newLaserTag_run197770-203755_GOOD.txt
#Compare2012C 2012C_GlobalMY_Vs_Caltech.root

#python MakeComparison_txt_txt.py /afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/lpernie/IC_2012B_MVAcorrResid_11Iter_MediatedGlobally.txt
#/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/2012B_NewTags_V2GOOD/interCalibConstants.2012B.EBandEE.GOOD.txt
#Compare2012B 2012B_GlobalMY_Vs_Caltech.root

#python MakeComparison_txt_txt.py /afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/lpernie/IC_2012D_MVAcorrResid_13Iter_MediatedGlobally.txt
#/afs/cern.ch/cms/CAF/CMSALCA/ALCA_ECALCALIB/energy-calibration-repository/piZero/2012D_NewTags_V2GOOD/interCalibConstants.combined_2012D_EBandEE_GOOD.txt
#Compare2012D 2012D_GlobalMY_Vs_Caltech.root

###START
print "-----------------START-------------------"
if (len(sys.argv) != 5):
    usage()
    sys.exit(1)
try:
    #Input
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    print "Comparing: " +str(file1)
    print "with: " +str(file2)
    file1_r = open(file1,'r')
    file2_r = open(file2,'r')
    file1_v = file1_r.readlines()
    file2_v = file2_r.readlines()
    file1_v1 = file1_v[:]
    file2_v1 = file2_v[:]
    print 'Lunghezza files: ' + str(len(file1_v)) + ' vs ' + str(len(file2_v))
    if not(len(file1_v)==len(file2_v)):
       print 'WARNING, FILE WITH DIFFERENT LINES'
       sys.exit(1)
except:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)

#Read IC
IC_EB_1=list()
IC_EB_2=list()
IC_EE_1=list()
IC_EE_2=list()
print 'Starting reading files!'
for nXtal in range(len(file1_v)):
    List1=file1_v1[nXtal].split(' ')
    List2=file2_v1[nXtal].split(' ')
    if(nXtal<61200):
       IC_EB_1.append([List1[0],List1[1],List1[3]]) #eta Phi IC
       IC_EB_2.append([List2[0],List2[1],List2[3]])
       #print 'EB: ' + str(nXtal) + ': ' + str(IC_EB_1[nXtal][0]) + ' ' + str(IC_EB_1[nXtal][1]) + ' ' + str(IC_EB_1[nXtal][2])
    else:
       IC_EE_1.append([List1[0],List1[1],List1[2],List1[3]]) #iX iY iZ IC
       IC_EE_2.append([List2[0],List2[1],List2[2],List2[3]])
       #print 'EE1: ' + str(nXtal-61200) + ': ' + str(IC_EE_1[nXtal-61200][0]) + ' ' + str(IC_EE_1[nXtal-61200][1]) + ' ' + str(IC_EE_1[nXtal-61200][2]) + ' ' + str(IC_EE_1[nXtal-61200][3])
       #print 'EE2: ' + str(nXtal-61200) + ': ' + str(IC_EE_2[nXtal-61200][0]) + ' ' + str(IC_EE_2[nXtal-61200][1]) + ' ' + str(IC_EE_2[nXtal-61200][2]) + ' ' + str(IC_EE_2[nXtal-61200][3])
#Out Put
f = TFile(str(sys.argv[3]) + "/" + str(sys.argv[4]) , 'recreate')
text_file = open(str(sys.argv[3]) + "/Different_xTal.txt", "w")
Ratio1D_EB = TH1F("Ratio1D_EB", "Ratio IC EB", 100, 0.95, 1.05)
Ratio1D_EB_eta1 = TH1F("Ratio1D_EB_eta1", "Ratio IC EB", 100, 0.95, 1.05)
Ratio1D_EBm = TH1F("Ratio1D_EBm", "Ratio IC EB -", 100, 0.95, 1.05)
Ratio1D_EBp = TH1F("Ratio1D_EBp", "Ratio IC EB +", 100, 0.95, 1.05)
Ratio1D_EBm_eta1 = TH1F("Ratio1D_EBm_eta1", "Ratio IC EB_eta1 -", 100, 0.95, 1.05)
Ratio1D_EBp_eta1 = TH1F("Ratio1D_EBp_eta1", "Ratio IC EB_eta1 +", 100, 0.95, 1.05)
Ratio1D_EEm = TH1F("Ratio1D_EEm", "Ratio IC EEm", 100, 0.7, 1.3)
Ratio1D_EEp = TH1F("Ratio1D_EEp", "Ratio IC EEp", 100, 0.7, 1.3)
Ratio_EBMap = TH2F("Ratio_EBMap","Ratio EB Map #phi on x, #eta on y", 360,0.5,360.5,171,-85.5,85.5)
Ratio_EEpMap = TH2F("Ratio_EEpMap","Ratio EEp Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
Ratio_EEmMap = TH2F("Ratio_EEmMap","Ratio EEm Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
#to be mediated
cal_EBMap1 = TH2F("calibMap_EB_1","Ratio EB Map #eta on x, #phi on y",171,-85.5,85.5, 360,0.5,360.5)
cal_EEpMap1 = TH2F("calibMap_EEp_1","Ratio EEp Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
cal_EEmMap1 = TH2F("calibMap_EEm_1","Ratio EEm Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
cal_EBMap2 = TH2F("calibMap_EB_2","Ratio EB Map #eta on x, #phi on y",171,-85.5,85.5, 360,0.5,360.5)
cal_EEpMap2 = TH2F("calibMap_EEp_2","Ratio EEp Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
cal_EEmMap2 = TH2F("calibMap_EEm_2","Ratio EEm Map iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
#Tree
gROOT.ProcessLine(\
      "struct EBStruct{\
        Int_t phi;\
        Int_t eta;\
        Double_t iC_1;\
        Double_t iC_2;\
};")
t=EBStruct()
gROOT.ProcessLine(\
      "struct EEStruct{\
        Int_t iX;\
        Int_t iY;\
        Int_t iZ;\
        Int_t etaRing;\
        Double_t iC_1;\
        Double_t iC_2;\
};")
s=EEStruct()
mytree = ROOT.TTree("EB_IC_comaprison", "EB_IC_comaprison")
mytree.Branch('phi', AddressOf(t,'phi'), 'phi/I')
mytree.Branch('eta', AddressOf(t,'eta'), 'eta/I')
mytree.Branch('iC_1', AddressOf(t,'iC_1'), 'iC_1/D')
mytree.Branch('iC_2', AddressOf(t,'iC_2'), 'iC_2/D')
mytreeEE = ROOT.TTree("EE_IC_comaprison", "EE_IC_comaprison")
mytreeEE.Branch('iX', AddressOf(s,'iX'), 'iX/I')
mytreeEE.Branch('iY', AddressOf(s,'iY'), 'iY/I')
mytreeEE.Branch('iZ', AddressOf(s,'iZ'), 'iZ/I')
mytreeEE.Branch('etaRing', AddressOf(s,'etaRing'), 'etaRing/I')
mytreeEE.Branch('iC_1', AddressOf(s,'iC_1'), 'iC_1/D')
mytreeEE.Branch('iC_2', AddressOf(s,'iC_2'), 'iC_2/D')

#Compare EB
print 'Starting with EB (' + str(len(IC_EB_1)) + ' xtal)'
if not(len(IC_EB_1)==len(IC_EB_1)):
   print 'WARNING, EB1 != EB2'
   sys.exit(1)
for find_EB1 in range(len(IC_EB_1)):
    finded=False
    for find_EB2 in range(len(IC_EB_2)):
        if(IC_EB_1[find_EB1][0]==IC_EB_2[find_EB2][0] and IC_EB_1[find_EB1][1]==IC_EB_2[find_EB2][1]):
           #print 'Eta: ' + str(IC_EB_1[find_EB1][0]) + ', Phi = ' + str(IC_EB_1[find_EB1][1]) + '| Eta: ' + str(IC_EB_2[find_EB2][0]) + ', Phi = ' + str(IC_EB_2[find_EB2][1])
           t.eta=int(IC_EB_1[find_EB1][0]); t.phi=int(IC_EB_1[find_EB1][1])
           t.iC_1=float(IC_EB_1[find_EB1][2]); t.iC_2=float(IC_EB_2[find_EB2][2])
           mytree.Fill()
           if(float(IC_EB_1[find_EB1][2])!=-1. and float(IC_EB_2[find_EB2][2])!=-1. ):
              print 'Eta: ' + str(IC_EB_1[find_EB1][0]) + ', Phi = ' + str(IC_EB_1[find_EB1][1]) + ': ' + str(IC_EB_1[find_EB1][2]) + ' vs ' + str(IC_EB_2[find_EB2][2])
              Ratio1D_EB.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( abs(int(IC_EB_1[find_EB1][0]))<45. ): Ratio1D_EB_eta1.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( int(IC_EB_1[find_EB1][0])<0 ): Ratio1D_EBm.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( int(IC_EB_1[find_EB1][0])<0 and abs(int(IC_EB_1[find_EB1][0]))<45. ): Ratio1D_EBm_eta1.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( int(IC_EB_1[find_EB1][0])>0 ): Ratio1D_EBp.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( int(IC_EB_1[find_EB1][0])>0 and abs(int(IC_EB_1[find_EB1][0]))<45. ): Ratio1D_EBp_eta1.Fill(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              if( abs(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))>1.02 ): text_file.write( "iEta: " + str(int(IC_EB_1[find_EB1][0])+86) + " iPhi " + (str(IC_EB_1[find_EB1][1])) + " IC1/IC2: " + str(abs(float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))) + "\n")
              Ratio_EBMap.SetBinContent(int(IC_EB_1[find_EB1][1]), int(IC_EB_1[find_EB1][0])+86, float(IC_EB_1[find_EB1][2])/float(IC_EB_2[find_EB2][2]))
              cal_EBMap1.SetBinContent(int(IC_EB_1[find_EB1][0])+86, int(IC_EB_1[find_EB1][1]), float(IC_EB_1[find_EB1][2]))
              cal_EBMap2.SetBinContent(int(IC_EB_2[find_EB2][0])+86, int(IC_EB_2[find_EB2][1]), float(IC_EB_2[find_EB2][2]))
           else:
              cal_EBMap1.SetBinContent(int(IC_EB_1[find_EB1][0])+86, int(IC_EB_1[find_EB1][1]), 1.)
              cal_EBMap2.SetBinContent(int(IC_EB_2[find_EB2][0])+86, int(IC_EB_2[find_EB2][1]), 1.)
              print 'No: ' + str(IC_EB_1[find_EB1][2]) + ' ' + str(IC_EB_2[find_EB2][2])
           finded=True
           break
    if not(finded):
        print 'WARNING, an xtal in file1 is not in file2!!!'
#-----------------------
#Compare EE
IC1_bad=0
IC2_bad=0
print 'and then with EE (' + str(len(IC_EE_1)) + ' xtal)'
if not(len(IC_EE_1)==len(IC_EE_1)):
   print 'WARNING, EE1 != EE2'
   sys.exit(1)
for find_EE1 in range(len(IC_EE_1)):
    finded=False
    for find_EE2 in range(len(IC_EE_2)):
        if(IC_EE_1[find_EE1][0]==IC_EE_2[find_EE2][0] and IC_EE_1[find_EE1][1]==IC_EE_2[find_EE2][1] and IC_EE_1[find_EE1][2]==IC_EE_2[find_EE2][2]):
           s.iX=int(IC_EE_1[find_EE1][0]); s.iY=int(IC_EE_1[find_EE1][1]); s.iZ=int(IC_EE_1[find_EE1][2])
           s.iC_1=float(IC_EE_1[find_EE1][3]); s.iC_2=float(IC_EE_2[find_EE2][3])
           mytreeEE.Fill()
           if(float(IC_EE_1[find_EE1][3])!=-1. and float(IC_EE_2[find_EE2][3])!=-1. ):
              #print 'iX: ' + str(IC_EE_1[find_EE1][0]) + ', iY: ' + str(IC_EE_1[find_EE1][1]) + ', iZ: ' + str(IC_EE_1[find_EE1][2]) + ': ' + str(IC_EE_1[find_EE1][3]) + ' vs ' + str(IC_EE_2[find_EE2][3])
              if(int(IC_EE_1[find_EE1][2])==-1):
                 Ratio1D_EEm.Fill(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))
                 Ratio_EEmMap.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))
                 cal_EEmMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3]))
                 cal_EEmMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), float(IC_EE_2[find_EE2][3]))
                 if(abs(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))>1.7): text_file.write( "iX: " + str(IC_EE_1[find_EE1][0] + " iY " + (IC_EE_1[find_EE1][1]) + " iZ=-1 IC1/IC2: " + str(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))) + "\n")
              if(int(IC_EE_1[find_EE1][2])==1):
                 Ratio1D_EEp.Fill(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))
                 Ratio_EEpMap.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))
                 cal_EEpMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3]))
                 cal_EEpMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), float(IC_EE_2[find_EE2][3]))
                 if(abs(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))>1.7): text_file.write( "iX: " + str(IC_EE_1[find_EE1][0] + " iY " + (IC_EE_1[find_EE1][1]) + " iZ=1 IC1/IC2: " + str(float(IC_EE_1[find_EE1][3])/float(IC_EE_2[find_EE2][3]))) + "\n")
           else:
              if(int(IC_EE_1[find_EE1][2])==1): 
                 if(float(IC_EE_1[find_EE1][3])!=-1.): cal_EEpMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3]))
                 else:                                 cal_EEpMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), 1.)
                 if(float(IC_EE_2[find_EE1][3])!=-1.): cal_EEpMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), float(IC_EE_2[find_EE1][3]))
                 else:                                 cal_EEpMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), 1.)
              if(int(IC_EE_1[find_EE1][2])==-1):
                 if(float(IC_EE_1[find_EE1][3])!=-1.): cal_EEmMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), float(IC_EE_1[find_EE1][3]))
                 else:                                 cal_EEmMap1.SetBinContent(int(IC_EE_1[find_EE1][0]), int(IC_EE_1[find_EE1][1]), 1.)
                 if(float(IC_EE_2[find_EE1][3])!=-1.): cal_EEmMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), float(IC_EE_2[find_EE1][3]))
                 else                                : cal_EEmMap2.SetBinContent(int(IC_EE_2[find_EE2][0]), int(IC_EE_2[find_EE2][1]), 1.)
              #print 'No: ' + str(IC_EE_1[find_EE1][0]) + ' ' + str(IC_EE_1[find_EE1][1]) + '' + str(IC_EE_1[find_EE1][2]) + ' ' + str(IC_EE_2[find_EE2][2])
              if(float(IC_EE_1[find_EE1][3])==-1.): IC1_bad=IC1_bad+1
              if(float(IC_EE_2[find_EE1][3])==-1.): IC2_bad=IC2_bad+1
           finded=True
           break
    if not(finded):
        print 'WARNING, an xtal in file1 is not in file2!!!'
print 'Now saving files.'
print 'IC1 has no: ' + str(IC1_bad) + " xtal"
print 'IC2 has no: ' + str(IC2_bad) + " xtal"

print "------------------END--------------------"
#Final Writing
f.cd()
f.Write()
text_file.close()
f.Close()
