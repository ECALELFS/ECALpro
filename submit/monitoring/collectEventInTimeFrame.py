#!/usr/bin/env python3

import sys
import ROOT
from ROOT import gStyle
from ROOT import gROOT
import math
import os
import re

import uproot
import numpy as np
import pandas as pd

from optparse import OptionParser
parser = OptionParser()

ROOT.gROOT.SetBatch(True)


nbins = 120
xmin = 0.
xmax = 0.3


def collectEvents(inputFileNameList, inputTreeName, outputFileName, nPi0EB, nPi0EE):
    fout = ROOT.TFile("%s"%outputFileName,"RECREATE")

    #print(pd.__version__)
    df_Init = uproot.iterate({inputFile : inputTreeName for inputFile in inputFileNameList}, library="pd")
    df_list = [data for data in df_Init]
    dfFull = pd.concat(df_list)

    ###First sort the df by time 
    df = dfFull.sort_values(by='event_time', ascending=True)
    ###group by events collected on the same day, month, year 
    df = df.groupby([dfFull.event_year,dfFull.event_month,dfFull.event_day,dfFull.isPi0EB])
    


    nPi0InHist = []  ### number of pi0's we need in each region
    nEvents = [] ## number of rows in each group
    nChunks = [] ####number of chunks formed
    group_size = df.size() ### number of rows in each group
    ngroups = df.ngroups ### numnber of groups formed
    #for ig in range(0,ngroups):
    ig = 0
    for key, item in df:
        
        if(item.all().isPi0EB==0):
            nPi0InHist.append(nPi0EE)
        if(item.all().isPi0EB==1):
            nPi0InHist.append(nPi0EB)
        
        nEvents.append(group_size.values[ig])
        print(group_size.values[ig])
        ig = ig+1 
 
    for ig in range(0,ngroups):
        print(nPi0InHist[ig])
    
    
    list_df = []
    ig = 0
    for key, item in df:
        n = nEvents[ig]/nPi0InHist[ig]  ### number of groups needed = #events in that group/#pi0's we want
        if(n<1):
            n = 1
        nChunks.append(int(n))
        print("Number of events in this group and number of pi0s needed and hence number of group %f : %f : %f"%(nEvents[ig],nPi0InHist[ig],n))
        list_df.append(np.array_split(df.get_group(key), n))
        ig = ig + 1
    
    ####Now loop over the groups, inside the group, loop over chunks and then fill them in the same histogram
    ig = 0
    histmap = {}
    for key, item in df:
        region = "EB"
        if(item.all().isPi0EB==0):
            region = "EE"
        year = df.get_group(key).event_year.iloc[0]  ###first element of this group will have the same year,month, day as the rest
        month = df.get_group(key).event_month.iloc[0]
        day = df.get_group(key).event_day.iloc[0]
        print(year)
        n = nChunks[ig]
        for i in range(0,n): ###number of time chunks in that group
        ###define histogram for each time/nEV chunk in that group
            time_of_first_ev = list_df[ig][i].event_time.iloc[0] ##time of the first event in that chunk
            histName = "pi0_mass_%d_%d_%d_%f_%s"%(year,month,day,time_of_first_ev,region)
        
            print(histName)
            histmap[histName] = ROOT.TH1F(histName,histName,nbins,xmin,xmax)
        
            nPi0s = len(list_df[ig][i])
            print("Number of pi0s in this chunk %d"%nPi0s)
            for ipi0 in range(0,nPi0s): ###number of pi0 in that time chunk
                histmap[histName].Fill(list_df[ig][i].pi0_mass.iloc[ipi0])
        ig = ig + 1        ### group index
    fout.cd()
    fout.Write()    



def main():
    parser.add_option("-f", "--inFile", dest="inputFileNameList",
                      help="Path to input rootfiles/*root, e.g. data/*root", 
                      default="data/*root")
    
    parser.add_option("-i", "--inTree", dest="inputTreeName",
                      help="Input tree name", 
                      default="monitoring")
    

    parser.add_option("-o", "--outFile", dest="outputFileName",
                      help="Path to output file", 
                      default="hist_pi0mass.root")
    
    parser.add_option("-B", "--nPi0EB", dest="nPi0EB",
                      help="Number of pi0s in the EB", 
                      default=250000)

    parser.add_option("-E", "--nPi0EE", dest="nPi0EE",
                      help="Number of pi0s in the EE", 
                      default=100000)


    (options, args) = parser.parse_args()
    
    collectEvents(options.inputFileNameList.split(','), options.inputTreeName, options.outputFileName, options.nPi0EB, options.nPi0EE)
    

    
if __name__ == "__main__":
    main()
