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
import awkward as ak

from optparse import OptionParser
parser = OptionParser()

ROOT.gROOT.SetBatch(True)


nbins = 120
xmin = 0.
xmax = 0.3


def collectEvents(inputFileNameList, inputTreeName, outputPath, nPi0EB, nPi0EE):
    """
    The splitting logic is very simple: collect enough pi0 candidates before declaring a
    now group. This of course works well when processing data from single fills, if processing more
    fills together one migh end up merging events with potentially large time gaps in between.
    """
    ROOT.gROOT.LoadMacro(os.environ['CMSSW_BASE']+'/src/CalibCode/submit/monitoring/fitMass.C')

    #print(pd.__version__)
    data = uproot.lazy({inputFile : inputTreeName for inputFile in inputFileNameList}, 
                       filter_name=["event_time", "isPi0EB", "pi0_mass"])

    index = ak.argsort(data.event_time)
    for f in data.fields:
        data[f] = data[f][index]

    data_part = {'EB' : data[data.isPi0EB==1],
                 'EE' : data[data.isPi0EB!=1]}

    outdata = []    

    #split data in chuncks and fill histograms
    #math.ceil is needed in order to collect the last period with potentially not enough pi0
    for part, dset in data_part.items():
        for ick, chunck in enumerate(np.array_split(dset, math.ceil(len(dset)/nPi0EB))):
            print(f'{part}: filling chunk {ick} with {len(chunck)} pi0 candidates')
            # fill mass histogram
            histName = f'pi0_mass_{ick}'
            hist = ROOT.TH1F(histName,histName,nbins,xmin,xmax)
            times = []
            for t,m,_ in chunck:
                times.append(t)
                hist.Fill(m)    
            mean_time = np.mean(times)

            # fit
            print("fitting")
            res = ROOT.fitMass(hist, ick, str(os.path.dirname(outputPath))+'/', part=="EB")

            # output data: mean time, EB/EE, mass, mass_unc
            outdata.append(f'{mean_time} {part} {res[0]} {res[1]}')

    with open(outputPath, 'w') as outf:
        outf.write('# avg-unix-time partition fit-mass fit-mass-unc\n')
        outf.write('\n'.join(outdata))


def main():
    parser.add_option("-f", "--inFile", dest="inputFileNameList",
                      help="Path to input rootfiles/*root, e.g. data/*root", 
                      default="data/*root")
    
    parser.add_option("-i", "--inTree", dest="inputTreeName",
                      help="Input tree name", 
                      default="monitoring")    

    parser.add_option("-o", "--outFile", dest="outFile",
                      help="Path to output file", 
                      default="./")
    
    parser.add_option("-B", "--nPi0EB", dest="nPi0EB",
                      help="Number of pi0s in the EB", 
                      default=250000)

    parser.add_option("-E", "--nPi0EE", dest="nPi0EE",
                      help="Number of pi0s in the EE", 
                      default=100000)


    (options, args) = parser.parse_args()
    
    collectEvents(inputFileNameList=options.inputFileNameList.split(','), 
                  inputTreeName=options.inputTreeName, 
                  outputPath=os.path.abspath(options.outFile), 
                  nPi0EB=options.nPi0EB, 
                  nPi0EE=options.nPi0EE)
    
    
if __name__ == "__main__":
    main()
