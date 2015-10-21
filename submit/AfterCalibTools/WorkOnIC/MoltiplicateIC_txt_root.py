#!/usr/bin/env python
from math import fabs
import operator
import xml.dom.minidom
import subprocess, time, sys, os, optparse
from ROOT import *
from PhysicsTools.PythonAnalysis import *
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(11111)

#python MoltiplicateIC_txt_root.py Original_IC/2015A_BOFF_dump_EcalIntercalibConstants__since_00239580_till_00251003.dat
#root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_2015A_RAW_RECHIT_SMIC_estimTime_01/iter_8/2015A_calibMap.root 2015A_Iter9 IC_Plotted.root --SystErr FIT

#python MoltiplicateIC_txt_root.py Original_IC/2015B_BON_dump_EcalIntercalibConstants__since_00251004_till_18446744073709551615.dat
#root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_2015B_Multifit_01/iter_13/2015B_calibMap.root 2015B_CMSSW746_GTGR_P_V56_Iter13_Error
#IC_Plotted.root --SystErr FIT

def usage():
    print "Usage: -----> python MoltiplicateIC_txt_root.py Ori_IC_TXT MY_IC_ROOT OutputFolder RootFile"

def MoltiplicateICFromTXT():
#Read the IC in the txt file (they should be never zero) and I multiply them for my IC. If mine are 1, I moltiply them anyway.
    for iEB in range(len(IC_EB_1)):
        OricalibMap_EB.SetBinContent( int(IC_EB_1[iEB][0]) + 86 , int(IC_EB_1[iEB][1]), float(IC_EB_1[iEB][2]) )
        OriCoef_EB.Fill( float(IC_EB_1[iEB][2]) )
        myIC = EBIC.GetBinContent( int(IC_EB_1[iEB][0]) + 86 , int(IC_EB_1[iEB][1]) )
        if(float(myIC)==0.):
            print "MoltiplicateICFromTXT::WARNING, my IC is Zero in EB"
        if(float(myIC)==1. and float(IC_EB_1[iEB][2])!=1.):
            name = str(int(IC_EB_1[iEB][0]) + 86) + "_" + str(IC_EB_1[iEB][1])
            ListBadFromMe_EB.append(name)
        newIC = float(IC_EB_1[iEB][2])*float(myIC)
        NewcalibMap_EB.SetBinContent( int(IC_EB_1[iEB][0]) + 86 , int(IC_EB_1[iEB][1]), float(newIC) )
    for iEE in range(len(IC_EE_1)):
        if( float(IC_EE_1[iEE][2]) < 0 ):
            OricalibMap_EEm.SetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]), float(IC_EE_1[iEE][3]) )
            OriCoef_EEm.Fill( float(IC_EE_1[iEE][3]) )
            myIC = EEmIC.GetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]) )
            if(float(myIC)==0.):
                print "MoltiplicateICFromTXT::WARNING, my IC is Zero in EEm"
            if(float(myIC)==1. and float(IC_EE_1[iEE][3])!=1.):
                name = str(IC_EE_1[iEE][0]) + "_" + str(IC_EE_1[iEE][1])
                ListBadFromMe_EEm.append(name)
            newIC = float(IC_EE_1[iEE][3])*float(myIC)
            NewcalibMap_EEm.SetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]), float(newIC) )
        if( float(IC_EE_1[iEE][2]) > 0 ):
            OricalibMap_EEp.SetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]), float(IC_EE_1[iEE][3]) )
            OriCoef_EEp.Fill( float(IC_EE_1[iEE][3]) )
            myIC = EEpIC.GetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]) )
            if(float(myIC)==0.):
                print "MoltiplicateICFromTXT::WARNING, my IC is Zero in EEp"
            if(float(myIC)==1. and float(IC_EE_1[iEE][3])!=1.):
                name = str(IC_EE_1[iEE][0]) + "_" + str(IC_EE_1[iEE][1])
                ListBadFromMe_EEp.append(name)
            newIC = float(IC_EE_1[iEE][3])*float(myIC)
            NewcalibMap_EEp.SetBinContent( int(IC_EE_1[iEE][0]), int(IC_EE_1[iEE][1]), float(newIC) )

def WriteTXT(hEB,hEEm,hEEp,name,errorType,whichIC):
    outputfile = open( name, 'wr+' )
    outputfile.write("#iEta(ix) iPhi(iy) 0(iZ) IC error chi2\n")
    #EB
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            if( ((ieta+1)-86) != 0 ) :
                IC = hEB.GetBinContent( ieta+1, iphi+1 )
                if( IC==1. ):
                    outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC,6)) + " 999. 0\n")
                else:
                    Name = str(ieta+1) + "_" + str(iphi+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EB) ):
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC,6)) + " 999. 0\n")
                    else:
                        index = str(((ieta+1)-86)) + "_" + str(iphi+1)
                        Error = 0
                        if( SystE=="FIT" and whichIC=="mine" ): Error = SystEBList[index]
                        if( SystE=="FIT" and whichIC=="abs" ): Error = SystEBList[index] * OricalibMap_EB.GetBinContent(ieta+1, iphi+1)
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC,6)) + " " + str(Error) + " " + str(Chi2EBList[index]) + "\n")
    #EEm
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEm.GetBinContent( ix+1, iy+1 )
            if( IC!=0 ):
                if( IC==1. ):
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC,6)) + " 999. 0\n")
                else:
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEm) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC,6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_-1"
                        Error = 0
                        if( SystE=="FIT" and whichIC=="mine" ): Error = SystEEList[index]
                        if( SystE=="FIT" and whichIC=="abs" ): Error = SystEEList[index] * OricalibMap_EEm.GetBinContent(ix+1, iy+1)
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC,6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    #EEp
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEp.GetBinContent( ix+1, iy+1 )
            if( IC!=0 ):
                if( IC==1. ):
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC,6)) + " 999. 0\n")
                else:
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEp) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC,6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_1"
                        Error = 0
                        if( SystE=="FIT" and whichIC=="mine" ): Error = SystEEList[index]
                        if( SystE=="FIT" and whichIC=="abs" ): Error = SystEEList[index] * OricalibMap_EEp.GetBinContent(ix+1, iy+1)
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC,6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    outputfile.close()
    if( int(len(open(name).readlines())-1) != int(TotalIC) ):
        print "WARNING: Final IC has a number of lines different from the number of lines of the original IC!!!"
        print str(int(len(open(name).readlines()))) + " vs " + str(TotalIC)

def AverageGlobally(hEB,hEEm,hEEp,name,errorType):
    outputfile = open( name, 'wr+' )
    outputfile.write("#iEta(ix) iPhi(iy) 0(iZ) IC error chi2\n")
    #EB
    IC_tmp=0.
    IC_tot=0.
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            if( ((ieta+1)-86) != 0 ) :
                IC = hEB.GetBinContent( ieta+1, iphi+1 )
                if( IC!=1. ):
                    IC_tmp+=IC
                    IC_tot+=1
                    if( IC==0 ):
                        print "WARNING IC==0 in EB. This should not happen!!!"
    IC_tmp/=IC_tot
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            if( ((ieta+1)-86) != 0 ) :
                IC = hEB.GetBinContent( ieta+1, iphi+1 )
                if( IC==1. ):
                    NewcalibMap_Glob1_EB.SetBinContent( ieta+1, iphi+1, 1. )
                    outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 1. 999. 0\n")
                else:
                    NewcalibMap_Glob1_EB.SetBinContent( ieta+1, iphi+1, IC/IC_tmp )
                    Name = str(ieta+1) + "_" + str(iphi+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EB) ):
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC/IC_tmp,6)) + " 999. 0\n")
                    else:
                        index = str(((ieta+1)-86)) + "_" + str(iphi+1)
                        Error = 0
                        if (SystE=="FIT"): Error = SystEBList[index] * OricalibMap_EB.GetBinContent(ieta+1, iphi+1) / IC_tmp
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC/IC_tmp,6)) + " " + str(Error) + " " + str(Chi2EBList[index]) + "\n")
    #EEm
    IC_tmp=0.
    IC_tot=0.
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEm.GetBinContent( ix+1, iy+1 )
            if( IC!=0 and IC!=1. ):
                IC_tmp+=IC
                IC_tot+=1
    IC_tmp/=IC_tot
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEm.GetBinContent( ix+1, iy+1 )
            if( IC!=0):
                if( IC==1):
                    NewcalibMap_Glob1_EEm.SetBinContent( ix+1, iy+1, 1. )
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 1. 999. 0\n")
                else:
                    NewcalibMap_Glob1_EEm.SetBinContent( ix+1, iy+1, IC/IC_tmp )
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEm) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC/IC_tmp,6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_-1"
                        Error = 0
                        if (SystE=="FIT"): Error = SystEEList[index] * OricalibMap_EEm.GetBinContent(ix+1, iy+1) / IC_tmp
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC/IC_tmp,6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    #EEp
    IC_tmp=0.
    IC_tot=0.
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEp.GetBinContent( ix+1, iy+1 )
            if( IC!=0 and IC!=1. ):
                IC_tmp+=IC
                IC_tot+=1
    IC_tmp/=IC_tot
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEp.GetBinContent( ix+1, iy+1 )
            if( IC!=0):
                if( IC==1):
                    NewcalibMap_Glob1_EEp.SetBinContent( ix+1, iy+1, 1. )
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 1. 999. 0\n")
                else:
                    NewcalibMap_Glob1_EEp.SetBinContent( ix+1, iy+1, IC/IC_tmp )
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEp) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC/IC_tmp,6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_1"
                        Error = 0
                        if (SystE=="FIT"): Error = SystEEList[index] * OricalibMap_EEp.GetBinContent(ix+1, iy+1) / IC_tmp
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC/IC_tmp,6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    outputfile.close()
    if( int(len(open(name).readlines())-1) != int(TotalIC)):
        print "WARNING: Final IC has a number of lines different from the number of lines of the original IC!!!"
        print str(int(len(open(name).readlines()))) + " vs " + str(TotalIC)

def AverageEtaRing(hEB,hEEm,hEEp,name,EtaList,errorType):
    outputfile = open( name, 'wr+' )
    outputfile.write("#iEta(ix) iPhi(iy) 0(iZ) IC error chi2\n")
    #EB
    IC_tmp = [0.] * 171
    IC_tot = [0.] * 171
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            if( ((ieta+1)-86) != 0 ) :
                IC = hEB.GetBinContent( ieta+1, iphi+1 )
                if( IC!=1. ):
                    if( IC==0 ):
                        print "WARNING IC==0 in EB. This should not happen!!!"
                    IC_tmp[ieta]+=IC
                    IC_tot[ieta]+=1
    for ieta in range(171):
        if(IC_tot[ieta]>0):
            IC_tmp[ieta]/=IC_tot[ieta]
        else:
            IC_tmp[ieta]=1.
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            if( ((ieta+1)-86) != 0 ) :
                IC = hEB.GetBinContent( ieta+1, iphi+1 )
                if( IC==1. ):
                    NewcalibMap_EtaR1_EB.SetBinContent( ieta+1, iphi+1, 1. )
                    outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 1. 999. 0\n")
                else:
                    NewcalibMap_EtaR1_EB.SetBinContent( ieta+1, iphi+1, IC/IC_tmp[ieta] )
                    Name = str(ieta+1) + "_" + str(iphi+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EB) ):
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC/IC_tmp[ieta],6)) + " 999. 0\n")
                    else:
                        index = str(((ieta+1)-86)) + "_" + str(iphi+1)
                        Error = 0
                        if (SystE=="FIT"): Error = SystEBList[index] * OricalibMap_EB.GetBinContent(ieta+1, iphi+1) / IC_tmp[ieta]
                        outputfile.write( str(((ieta+1)-86)) + " " + str(iphi+1) + " 0 " + str(round(IC/IC_tmp[ieta],6)) + " " + str(Error) + " " + str(Chi2EBList[index]) + "\n")
    #EEm
    IC_tmp = [0.] * 39
    IC_tot = [0.] * 39
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEm.GetBinContent( ix+1, iy+1 )
            if( IC!=0 and IC!=1. ):
                if not( (str(ix+1) + "_" + str(iy+1)) in EtaList):
                    "NOT PRESENT!!!"
                Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                IC_tmp[int(Ring)]+=IC
                IC_tot[int(Ring)]+=1
    for iRing in range(39):
        if(IC_tot[iRing]>0):
            IC_tmp[iRing]/=IC_tot[iRing]
        else:
            IC_tmp[iRing]=1
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEm.GetBinContent( ix+1, iy+1 )
            if( IC!=0):              
                if( IC==1):
                    NewcalibMap_EtaR1_EEm.SetBinContent( ix+1, iy+1, 1. )
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 1. 999. 0\n")
                else:
                    Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                    NewcalibMap_EtaR1_EEm.SetBinContent( ix+1, iy+1, IC/IC_tmp[Ring] )
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEm) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC/IC_tmp[Ring],6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_-1"
                        Error = 0
                        if (SystE=="FIT"): Error = SystEEList[index] * OricalibMap_EEm.GetBinContent(ix+1, iy+1) / IC_tmp[Ring]
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " -1 " + str(round(IC/IC_tmp[Ring],6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    #EEp
    IC_tmp = [0.] * 39
    IC_tot = [0.] * 39
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEp.GetBinContent( ix+1, iy+1 )
            if( IC!=0 and IC!=1. ):
                if not( (str(ix+1) + "_" + str(iy+1)) in EtaList):
                    "NOT PRESENT!!!"
                Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                IC_tmp[Ring]+=IC
                IC_tot[Ring]+=1
    for iRing in range(39):
        if(IC_tot[iRing]>0):
            IC_tmp[iRing]/=IC_tot[iRing]
        else:
            IC_tmp[iRing]=1
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            IC = hEEp.GetBinContent( ix+1, iy+1 )
            if( IC!=0):              
                if( IC==1):
                    NewcalibMap_EtaR1_EEp.SetBinContent( ix+1, iy+1, 1. )
                    outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 1. 999. 0\n")
                else:
                    Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                    NewcalibMap_EtaR1_EEp.SetBinContent( ix+1, iy+1, IC/IC_tmp[Ring] )
                    Name = str(ix+1) + "_" + str(iy+1)
                    if( errorType=="ErrorFromMyIC" and  Name in set(ListBadFromMe_EEp) ):
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC/IC_tmp[Ring],6)) + " 999. 0\n")
                    else:
                        index = str(ix+1) + "_" + str(iy+1) + "_1"
                        Error = 0
                        if (SystE=="FIT"): Error = SystEEList[index] * OricalibMap_EEp.GetBinContent(ix+1, iy+1) / IC_tmp[Ring]
                        outputfile.write( str(ix+1) + " " + str(iy+1) + " 1 " + str(round(IC/IC_tmp[Ring],6)) + " " + str(Error) + " " + str(Chi2EEList[index]) + "\n")
    outputfile.close()
    if( int(len(open(name).readlines())-1) != int(TotalIC)):
        print "WARNING: Final IC has a number of lines different from the number of lines of the original IC!!!"
        print str(int(len(open(name).readlines()))) + " vs " + str(TotalIC)

def TEST_average():
    IC_tmp=0.
    IC_tot=0.
    IC_tmpEta = [0.] * 171
    IC_totEta = [0.] * 171
    for ieta in range(171):#[0,170]
        for iphi in range(360):#[0,359]
            ICG = NewcalibMap_Glob1_EB.GetBinContent( ieta+1, iphi+1 )
            ICE = NewcalibMap_EtaR1_EB.GetBinContent( ieta+1, iphi+1 )
            if( ICG!=1. ):
                IC_tmp+=ICG
                IC_tot+=1
            if( ICE!=1. ):
                IC_tmpEta[ieta]+=ICE
                IC_totEta[ieta]+=1
    print "Test the EB IC: Globally are " + str(IC_tmp/IC_tot)
    for iEta in range(171):
        if( fabs(float(IC_tmpEta[iEta]/IC_totEta[iEta]-1.)) > 0.0001 and ((iEta+1)-86) != 0 ):
            print "WARNING::In EB and iEta " + str(iEta) + " are " + str(IC_tmpEta[iEta]) + " / " + str(IC_totEta[iEta]) + " = " + str(IC_tmpEta[iEta]/IC_totEta[iEta]) 
    IC_tmp_m=0.
    IC_tot_m=0.
    IC_tmp_p=0.
    IC_tot_p=0.
    IC_tmpEta_m = [0.] * 39 
    IC_totEta_m = [0.] * 39
    IC_tmpEta_p = [0.] * 39
    IC_totEta_p = [0.] * 39
    for ix in range(100):#[0,99]
        for iy in range(100):#[0,99]
            ICmG = NewcalibMap_Glob1_EEm.GetBinContent( ix+1, iy+1 )
            ICpG = NewcalibMap_Glob1_EEp.GetBinContent( ix+1, iy+1 )
            ICmE = NewcalibMap_EtaR1_EEm.GetBinContent( ix+1, iy+1 )
            ICpE = NewcalibMap_EtaR1_EEp.GetBinContent( ix+1, iy+1 )
            if( ICmG!=1. and ICmG!=0 ):
                IC_tmp_m+=ICmG
                IC_tot_m+=1
            if( ICpG!=1. and ICpG!=0 ):
                IC_tmp_p+=ICpG
                IC_tot_p+=1
            if( ICmE!=1. and ICmE!=0 ):
                Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                IC_tmpEta_m[Ring]+=ICmE
                IC_totEta_m[Ring]+=1
            if( ICpE!=1. and ICpE!=0 ):
                Ring = EtaList[str(ix+1) + "_" + str(iy+1)]
                IC_tmpEta_p[Ring]+=ICpE
                IC_totEta_p[Ring]+=1
    print "Test the EEm IC: Globally are " + str(IC_tmp_m) + " / " + str(IC_tot_m) + " = " + str(IC_tmp_m/IC_tot_m)
    print "Test the EEp IC: Globally are " + str(IC_tmp_p) + " / " + str(IC_tot_p) + " = " + str(IC_tmp_p/IC_tot_p)
    for iRing in range(39):
        if(IC_totEta_m[iRing]>0):
            if( fabs(float(IC_tmpEta_m[iRing]/IC_totEta_m[iRing]-1.)) > 0.0001 ):
                print "WARNING::In EEm and iRing " + str(iRing) + " are " + str(IC_tmpEta_m[iRing]) + " / " + str(IC_totEta_m[iRing]) + " = " + str(IC_tmpEta_m[iRing]/IC_totEta_m[iRing])
        else:
            print "No Ring " + str(iRing)
        if(IC_totEta_p[iRing]>0):
            if( fabs(float(IC_tmpEta_p[iRing]/IC_totEta_p[iRing]-1.)) > 0.0001 ):
                print "WARNING::In EEp and iRing " + str(iRing) + " are " + str(IC_tmpEta_p[iRing]) + " / " + str(IC_totEta_p[iRing]) + " = " + str(IC_tmpEta_p[iRing]/IC_totEta_p[iRing])
        else:
            print "No Ring " + str(iRing)

print "STARTING"
Usage = """python MoltiplicateIC_txt_root.py Original_IC/2015A_BOFF_dump_EcalIntercalibConstants__since_00239580_till_00251003.dat
root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_2015A_RAW_RECHIT_SMIC_estimTime_01/iter_8/2015A_calibMap.root 2015A Absolute_IC.root
--SystErr FIT"""

if (len(sys.argv) < 5):
    usage(); sys.exit(1)
try:
    parser = optparse.OptionParser(Usage)
    parser.add_option("-S", "--SystErr",
                       help="Systematic error to be applied to the ICs. FIT means take the error form the mean param. of the fit.",
                       type="string",
                       default="none",
                       dest="SystErr")
    options, args = parser.parse_args()
    SystE = options.SystErr
    pathTXT1 = str(sys.argv[1])
    pathTH2  = str(sys.argv[2])
    OutputF  = str(sys.argv[3])
    Output   = str(sys.argv[4])
    folderCreation = subprocess.Popen(['mkdir -p ' + OutputF], stdout=subprocess.PIPE, shell=True);
    folderCreation.communicate()
    print "My IC are taken from: " + pathTH2
    print "The original IC are from: " + pathTXT1
    if not (os.path.isfile(pathTXT1)):
        print str(pathTXT1) + " doesn't exist!"
    file1_r  = open(pathTXT1,'r')
    file1_v  = file1_r.readlines()
    file1_v1 = file1_v[:]
    fileTH2  = ROOT.TFile.Open(pathTH2)
    EBIC     = fileTH2.Get('calibMap_EB')
    EEmIC    = fileTH2.Get('calibMap_EEm')
    EEpIC    = fileTH2.Get('calibMap_EEp')
    TreeEB   = fileTH2.Get("calibEB")
    TreeEE   = fileTH2.Get("calibEE")
    #Read EtaRing
    Endc_x_y_ring="../../../FillEpsilonPlot/data/Endc_x_y_ring.txt"
    print "The File to di the iRing Map is: " + str(Endc_x_y_ring)
    if not (os.path.isfile(Endc_x_y_ring)):
        print str(Endc_x_y_ring) + " doesn't exist!"
    EtaRing_r  = open(Endc_x_y_ring,'r')
    EtaRing_v  = EtaRing_r.readlines()
    EtaRing_v1 = EtaRing_v[:]
except:
    print "ERROR: something wrong with the arguments: \n %s" % sys.argv[1:]
    usage(); sys.exit(1)
#Read Ori IC
print 'Reading the original IC'
IC_EB_1=list(); IC_EE_1=list()
TotalIC = len(file1_v)
for nXtal in range(len(file1_v)):
    List1=file1_v1[nXtal].split(' ')
    if(nXtal<61200):
       IC_EB_1.append([List1[0],List1[1],List1[3]]) #eta Phi IC
    else:
       IC_EE_1.append([List1[0],List1[1],List1[2],List1[3]]) #iX iY iZ IC
print 'Now Storing EtaRing Map'
EtaList={}; maxRing=0
for nRing in range(len(EtaRing_v)):
     ListEta=EtaRing_v1[nRing].split(' ')
     EtaList[ str(int(ListEta[0])+1) + "_" + str(int(ListEta[1])+1) ] = int(ListEta[3])
     if(float(ListEta[3]) > float(maxRing)):
         maxRing = ListEta[3]
print "The max iRing is: " + str(maxRing).rstrip()
#Tree
print "Reading syst error from TTree"
SystEBList={}; Chi2EBList={}
gROOT.ProcessLine(\
    "struct MyStructEB{\
        Float_t fit_mean_err_;\
        Float_t Chisqu_;\
        Int_t   ieta_;\
        Int_t   iphi_;\
    };")
sEB = MyStructEB()
TreeEB.SetBranchAddress('fit_mean_err_',AddressOf(sEB,'fit_mean_err_'));
TreeEB.SetBranchAddress('Chisqu_',AddressOf(sEB,'Chisqu_'));
TreeEB.SetBranchAddress('ieta_',AddressOf(sEB,'ieta_'));
TreeEB.SetBranchAddress('iphi_',AddressOf(sEB,'iphi_'));
for nT in range(TreeEB.GetEntries()):
    TreeEB.GetEntry(nT);
    name = str(sEB.ieta_) + "_" + str(sEB.iphi_)
    SystEBList[ str(name) ] = sEB.fit_mean_err_
    Chi2EBList[ str(name) ] = sEB.Chisqu_
SystEEList={}; Chi2EEList={};
gROOT.ProcessLine(\
    "struct MyStructEE{\
        Float_t fit_mean_err_;\
        Float_t Chisqu_;\
        Int_t   ix_;\
        Int_t   iy_;\
        Int_t   zside_;\
    };")
sEE = MyStructEE()
TreeEE.SetBranchAddress('fit_mean_err_',AddressOf(sEE,'fit_mean_err_'));
TreeEE.SetBranchAddress('Chisqu_',AddressOf(sEE,'Chisqu_'));
TreeEE.SetBranchAddress('ix_',AddressOf(sEE,'ix_'));
TreeEE.SetBranchAddress('iy_',AddressOf(sEE,'iy_'));
TreeEE.SetBranchAddress('zside_',AddressOf(sEE,'zside_'));
for nT in range(TreeEE.GetEntries()):
    TreeEE.GetEntry(nT);
    name = str(sEE.ix_) + "_" + str(sEE.iy_) + "_" + str(sEE.zside_)
    SystEEList[ str(name) ] = sEE.fit_mean_err_
    Chi2EEList[ str(name) ] = sEE.Chisqu_

#Create Histos
print "Creating Histos"
f = TFile(OutputF + "/" + Output, 'recreate')
MyIC_EB               = fileTH2.Get('calibMap_EB')
NewcalibMap_EB        = TH2F("Abs_CalibMap_EB", "Absolute EB IC: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
OricalibMap_EB        = TH2F("OricalibMap_EB", "EB IC from GT: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
OriCoef_EB            = TH1F("OriCoef_EB", "IC EB from GT", 100, 0.2, 1.8)
NewcalibMap_Glob1_EB  = TH2F("NewcalibMap_Glob1_EB", "Absolute EB IC Globally to 1: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
NewcalibMap_EtaR1_EB  = TH2F("NewcalibMap_EtaR1_EB", "Absolute EB IC EtaRing to 1: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
MyIC_EEm              = fileTH2.Get('calibMap_EEm')
NewcalibMap_EEm       = TH2F("Abs_calibMap_EEm", "Absolute EEm IC", 100,0.5,100.5,100,0.5,100.5)
OricalibMap_EEm       = TH2F("OricalibMap_EEm", "EEm IC from GT", 100,0.5,100.5,100,0.5,100.5)
OriCoef_EEm           = TH1F("OriCoef_EEm","EEm IC from GT",100, 0.2, 1.8)
NewcalibMap_Glob1_EEm = TH2F("NewcalibMap_Glob1_EEm", "Absolute EEm IC Globally to 1", 100,0.5,100.5,100,0.5,100.5)
NewcalibMap_EtaR1_EEm = TH2F("NewcalibMap_EtaR1_EEm", "Absolute EEm IC EtaRing to 1", 100,0.5,100.5,100,0.5,100.5)
MyIC_EEp              = fileTH2.Get('calibMap_EEp')
NewcalibMap_EEp       = TH2F("Abs_calibMap_EEp", "Absolute EEp IC", 100,0.5,100.5,100,0.5,100.5)
OricalibMap_EEp       = TH2F("OricalibMap_EEp", "EEp IC from GT", 100,0.5,100.5,100,0.5,100.5)
OriCoef_EEp           = TH1F("OriCoef_EEp","EEp IC from GT",100, 0.2, 1.8)
NewcalibMap_Glob1_EEp = TH2F("NewcalibMap_Glob1_EEp", "Absolute EEp IC Globally to 1", 100,0.5,100.5,100,0.5,100.5)
NewcalibMap_EtaR1_EEp = TH2F("NewcalibMap_EtaR1_EEp", "Absolute EEp IC EtaRing to 1", 100,0.5,100.5,100,0.5,100.5)

#Moltiplicate IC
print 'Executing MoltiplicateICFromTXT'
ListBadFromMe_EB=list(); ListBadFromMe_EEm=list(); ListBadFromMe_EEp=list();
MoltiplicateICFromTXT()
#Write txt
print 'Executing WriteTXT1 for IC_fromECALpro.txt'
name = OutputF + "/IC_fromECALpro.txt"
WriteTXT(EBIC,EEmIC,EEpIC,name,"none","mine")
print 'Executing WriteTXT1 for IC_fromECALpro_Absolute.txt'
name = OutputF + "/IC_fromECALpro_Absolute.txt"
WriteTXT(NewcalibMap_EB,NewcalibMap_EEm,NewcalibMap_EEp,name,"ErrorFromMyIC","abs") #ErrorFromMyIC does that if I have no IC, you place the Original IC with 999. error.
#Average to 1 Globally
print 'Executing AverageGlobally'
name = OutputF + "/IC_fromECALpro_Absolute_Global1.txt"
AverageGlobally(NewcalibMap_EB,NewcalibMap_EEm,NewcalibMap_EEp,name,"ErrorFromMyIC")
#Average to 1 PerEtaRing
print 'Executing AverageEtaRing'
name = OutputF + "/IC_fromECALpro_Absolute_EtaRing1.txt"
AverageEtaRing(NewcalibMap_EB,NewcalibMap_EEm,NewcalibMap_EEp,name,EtaList,"ErrorFromMyIC")
#Test the average
print 'Now the final test...'
TEST_average()
print 'THE END!'
f.cd()
f.Write()
f.Close()
