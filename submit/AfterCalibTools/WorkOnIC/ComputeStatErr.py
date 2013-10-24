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
    print "Usage: -----> python ComputeStatErr.py File_even.root File_odd.root Output.root"
#python ComputeStatErr.py /store/group/alca_ecalcalib/lpernie/ALL_2012D_NormSelect_NormFit_Merged_01/iter_13_even/2012Dmerg_calibMap.root
#/store/group/alca_ecalcalib/lpernie/ALL_2012D_NormSelect_NormFit_Merged_01/iter_13_odd/2012Dmerg_calibMap.root 2012D/Error_Stat_2012D.root

###START
print "START"
if (len(sys.argv) != 4):
    usage()
    sys.exit(1)
try:
    #Input
    pathTH2_eve = 'root://eoscms//eos/cms' + sys.argv[1]
    pathTH2_odd = 'root://eoscms//eos/cms' + sys.argv[2]
    Output  = sys.argv[3]
    print "Reading TH2F with Even IC in: " +str(pathTH2_eve)
    print "Reading TH2F with Odd IC in: " +str(pathTH2_odd)
    fileTH2_eve = ROOT.TFile.Open(pathTH2_eve)
    fileTH2_odd = ROOT.TFile.Open(pathTH2_odd)
    EBIC_eve    = fileTH2_eve.Get('calibMap_EB')
    EEmIC_eve   = fileTH2_eve.Get('calibMap_EEm')
    EEpIC_eve   = fileTH2_eve.Get('calibMap_EEp')
    EBIC_odd    = fileTH2_odd.Get('calibMap_EB')
    EEmIC_odd   = fileTH2_odd.Get('calibMap_EEm')
    EEpIC_odd   = fileTH2_odd.Get('calibMap_EEp')
except:
    print "ERROR: unknown options in argument %s" % sys.argv[1:]
    usage()
    sys.exit(1)
#outPut
f = TFile(Output, 'recreate')
Sig_EB = TH2F("Sig_EB","EB IC Sigma #eta on x, #phi on y",171,-85.5,85.5, 360,0.5,360.5)
Sig_EEm = TH2F("Sig_EEm","EEm IC Sigma iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
Sig_EEp = TH2F("Sig_EEp","EEp IC Sigma iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
#Tree
gROOT.ProcessLine(\
      "struct EBStruct{\
        Double_t iC_eve;\
        Double_t iC_odd;\
        Int_t phi;\
        Int_t eta;\
};")
gROOT.ProcessLine(\
      "struct EEStruct{\
        Double_t iC_eve;\
        Double_t iC_odd;\
        Int_t iX;\
        Int_t iY;\
        Int_t iZ;\
};")
t=EBStruct()
s=EEStruct()
mytree_EB = ROOT.TTree("IC_StatErr_EB", "IC_StatErr_EB")
mytree_EB.Branch('phi', AddressOf(t,'phi'), 'phi/I')
mytree_EB.Branch('eta', AddressOf(t,'eta'), 'eta/I')
mytree_EB.Branch('iC_eve', AddressOf(t,'iC_eve'), 'iC_eve/D')
mytree_EB.Branch('iC_odd', AddressOf(t,'iC_odd'), 'iC_odd/D')
mytree_EE = ROOT.TTree("IC_StatErr_EE", "IC_StatErr_EE")
mytree_EE.Branch('iX', AddressOf(s,'iX'), 'iX/I')
mytree_EE.Branch('iY', AddressOf(s,'iY'), 'iY/I')
mytree_EE.Branch('iZ', AddressOf(s,'iZ'), 'iZ/I')
mytree_EE.Branch('iC_eve', AddressOf(s,'iC_eve'), 'iC_eve/D')
mytree_EE.Branch('iC_odd', AddressOf(s,'iC_odd'), 'iC_odd/D')

#Now we really start
#EB
print 'Let\'s start with EB'
MIN_ETA=1
MAX_ETA=85
MAX_PHI=360
for iEta in range(MIN_ETA,2*MAX_ETA+2):
   for iPhi in range(1,MAX_PHI):
       ic_even = EBIC_eve.GetBinContent(iEta,iPhi)
       ic_odd  = EBIC_odd.GetBinContent(iEta,iPhi)
       if( iEta!=MAX_ETA+1 ):
           Sig_EB.SetBinContent(iEta,iPhi, ic_even-ic_odd)
           print str(ic_even) + "  " + str(ic_odd)
           t.iC_eve=float(ic_even); t.iC_odd=float(ic_odd); t.eta=int(iEta-86); t.phi=int(iPhi);
           mytree_EB.Fill()
#EEm
print 'Now EE-...'
for ix in range(1,101):
   for iy in range(1,101):
       ic_even = EEmIC_eve.GetBinContent(ix,iy)
       ic_odd  = EEmIC_odd.GetBinContent(ix,iy)
       if( ic_even!=0 and ic_odd!=0 ): 
           Sig_EEm.SetBinContent(ix,iy, ic_even-ic_odd)
           s.iC_eve=float(ic_even); s.iC_odd=float(ic_odd); s.iX=int(ix); s.iY=int(iy); s.iZ=int(-1);
           mytree_EE.Fill()
#EEp
print 'En finally EE+'
for ix in range(1,101):
   for iy in range(1,101):
       ic_even = EEpIC_eve.GetBinContent(ix,iy)
       ic_odd  = EEpIC_odd.GetBinContent(ix,iy)
       if( ic_even!=0 and ic_odd!=0 ):
           Sig_EEp.SetBinContent(ix,iy, ic_even-ic_odd)
           s.iC_eve=float(ic_even); s.iC_odd=float(ic_odd); s.iX=int(ix); s.iY=int(iy); s.iZ=int(1);
           mytree_EE.Fill()

print 'Finish!!!'

#Final Writing
f.cd()
f.Write()
f.Close()
