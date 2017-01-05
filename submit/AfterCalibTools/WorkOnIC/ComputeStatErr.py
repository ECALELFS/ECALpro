#!/usr/bin/env python

import ROOT
#from sys import argv
from math import fabs, sqrt

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
#python ComputeStatErr.py /store/group/alca_ecalcalib/lpernie/ALL_2012C_pi0_NewTag_01/iter_10_odd/2012C_calibMap.root
#/store/group/alca_ecalcalib/lpernie/ALL_2012C_pi0_NewTag_01/iter_10_even/2012C_calibMap.root Error_2012C_pi0_Stat_newTag/2012C_calibMap_pi0StatErr.root

def FindEtaBin(ix, iy, iz, GeoFile_v1):
    eta=-10
    for find_eta1 in range(len(eta_1)):
       #print 'eta1: ' + str(eta_1[find_eta1][0]) + ' '  + str(eta_1[find_eta1][1]) + ' ' + str(eta_1[find_eta1][2])
       #print 'eta2: ' + str(ix) + ' '  + str(iy) + ' ' + str(iz)
       if(int(eta_1[find_eta1][0])==ix and int(eta_1[find_eta1][1])==iy and int(eta_1[find_eta1][2])==iz):
           #print 'etaFinding!!!: ' + str(eta_1[find_eta1][0]) + ' '  + str(eta_1[find_eta1][1]) + ' ' + str(eta_1[find_eta1][2])
           eta=float(eta_1[find_eta1][3])
           break
    return float(eta)


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
Sig_EB            = TH2F("Sig_EB","EB IC Sigma #eta on x, #phi on y",171,-85.5,85.5, 360,0.5,360.5)
Sig_EEm           = TH2F("Sig_EEm","EEm IC Sigma iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
Sig_EEp           = TH2F("Sig_EEp","EEp IC Sigma iX on x, iY on y", 100,0.5,100.5, 100,0.5,100.5)
Sig_EB_Diff_iEta  = TH2F("Sig_EB_Diff_iEta","#iEta on x, #IC-IC2 on y",171,-85.5,85.5, 100,-0.02,0.02)
Sig_EEm_Diff_iEta = TH2F("Sig_EEm_Diff_iEta","#eta on x, #IC-IC2 on y",15,1.5,3, 100,-0.05,0.05)
Sig_EEp_Diff_iEta = TH2F("Sig_EEp_Diff_iEta","#eta on x, #IC-IC2 on y",15,1.5,3, 100,-0.05,0.05)
Sig_EE_Diff_iEta  = TH2F("Sig_EE_Diff_iEta","#eta on x, #IC-IC2 on y",15,1.5,3, 100,-0.05,0.05)
EB_StatErr        = TH1F("EB_StatErr","#eta on x, Stat. Err. on y",171,-85.5,85.5)
EB_StatErr.GetXaxis().SetTitle("iEta"); EB_StatErr.GetXaxis().SetTitle("#sigma");
EEm_StatErr       = TH1F("EEm_StatErr","#eta on x, Stat. Err. on y",15,1.5,3)
EEm_StatErr.GetXaxis().SetTitle("#eta"); EEm_StatErr.GetXaxis().SetTitle("#sigma");
EEp_StatErr       = TH1F("EEp_StatErr","#eta on x, Stat. Err. on y",15,1.5,3)
EEp_StatErr.GetXaxis().SetTitle("#eta"); EEp_StatErr.GetXaxis().SetTitle("#sigma");
#eta
GeoFile = "../../../../CalibCode/submit/common/geometry_ietaix_iphiiy_0iz_eta.dat"
GeoFile_r = open(GeoFile,'r')
GeoFile_v = GeoFile_r.readlines()
GeoFile_v1 = GeoFile_v[:]
eta_1=list()
for nXtal in range(len(GeoFile_v)):
    List_eta=GeoFile_v1[nXtal].split(' ')
    if(nXtal>=61200):
        eta_1.append([List_eta[0],List_eta[1],List_eta[2],List_eta[3]]) #iX iY iZ eta
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
        Float_t Eta;\
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
mytree_EE.Branch('Eta', AddressOf(s,'Eta'), 'Eta/F')
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
           Sig_EB_Diff_iEta.Fill(iEta-86,ic_even-ic_odd)
           print str(ic_even) + "  " + str(ic_odd) + " " + str(iEta) 
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
           Sig_EEm_Diff_iEta.Fill(fabs(float(FindEtaBin(int(ix),int(iy),-1,GeoFile_v1))),ic_even-ic_odd)
           Sig_EE_Diff_iEta.Fill(fabs(float(FindEtaBin(int(ix),int(iy),-1,GeoFile_v1))),ic_even-ic_odd)
           s.iC_eve=float(ic_even); s.iC_odd=float(ic_odd); s.iX=int(ix); s.iY=int(iy); s.iZ=int(-1); s.Eta=float(FindEtaBin(int(ix),int(iy),-1,GeoFile_v1));
           mytree_EE.Fill()
#EEp
print 'And finally EE+'
for ix in range(1,101):
   for iy in range(1,101):
       ic_even = EEpIC_eve.GetBinContent(ix,iy)
       ic_odd  = EEpIC_odd.GetBinContent(ix,iy)
       if( ic_even!=0 and ic_odd!=0 ):
           Sig_EEp.SetBinContent(ix,iy, ic_even-ic_odd)
           Sig_EEp_Diff_iEta.Fill(fabs(float(FindEtaBin(int(ix),int(iy),1,GeoFile_v1))),ic_even-ic_odd)
           Sig_EE_Diff_iEta.Fill(fabs(float(FindEtaBin(int(ix),int(iy),1,GeoFile_v1))),ic_even-ic_odd)
           s.iC_eve=float(ic_even); s.iC_odd=float(ic_odd); s.iX=int(ix); s.iY=int(iy); s.iZ=int(1);  s.Eta=float(FindEtaBin(int(ix),int(iy),1,GeoFile_v1));
           mytree_EE.Fill()

print 'Now Histos with the Stat. Unc. Look for: EB_StatErr, EEm_StatErr and EEp_StatErr'
for i in range(Sig_EB_Diff_iEta.GetNbinsX()):
   h1 = Sig_EB_Diff_iEta.ProjectionY("",i,i)
   EB_StatErr.SetBinContent(i,h1.GetRMS()*sqrt(2)/2)
for j in range(Sig_EEm_Diff_iEta.GetNbinsX()):
   h2 = Sig_EEm_Diff_iEta.ProjectionY("",j,j)
   EEm_StatErr.SetBinContent(j,h2.GetRMS()*sqrt(2)/2)
for k in range(Sig_EEp_Diff_iEta.GetNbinsX()):
   h3 = Sig_EEp_Diff_iEta.ProjectionY("",k,k)
   EEp_StatErr.SetBinContent(k,h3.GetRMS()*sqrt(2)/2)

print 'Finish!!!'

#Final Writing
f.cd()
f.Write()
f.Close()
