#! /usr/bin/env python
import os,sys
from math import *
import numpy as np
import ROOT as rt
rt.gROOT.SetBatch(True)

class CalibSmearer:
    
    def __init__(self,inputfile,options):
        self.mean = options.bias if hasattr(options,"bias") else 0
        self.sigma = options.smearing if hasattr(options,"smearing") else 0
        self.loadCalibMaps(inputfile)

    def loadCalibMaps(self,rtfile):
        print "smearing calibrations in file ",rtfile
        tf = rt.TFile.Open(rtfile)
        self.maps = {}
        self.maps["EB"] = tf.Get("calibMap_EB").Clone()
        self.maps["EEm"] = tf.Get("calibMap_EEm").Clone()
        self.maps["EEp"] = tf.Get("calibMap_EEp").Clone()
        for k,m in self.maps.iteritems(): m.SetDirectory(None)
        tf.Close()
    
    def smearOne(self,histo):
        coeffs = rt.TH1F("coeffs_%s" % histo.GetName(),"",1000,0.80,1.20)
        rnd = rt.TRandom3()
        for xb in range(1,histo.GetNbinsX()+1):
            for yb in range(1,histo.GetNbinsY()+1):
                val = histo.GetBinContent(xb,yb)
                valnew = val * (1 + rnd.Gaus(float(self.mean),float(self.sigma))) if val!=1 else val
                histo.SetBinContent(xb,yb,valnew)
                coeffs.Fill(valnew)
        return coeffs

    def writeMaps(self,options):
        outputfile = options.output
        tf = rt.TFile.Open(outputfile,"recreate")
        for k,v in self.maps.iteritems():
            print "smearing calib for ",k
            c = self.smearOne(v)
            v.SetDirectory(tf)
            tf.WriteTObject(v)
            c.SetDirectory(tf)
            tf.WriteTObject(c)
        tf.Close()

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] calibMap.root ")
    parser.add_option("-b","--bias", dest="bias", default=0, help="bias of the calibrations")
    parser.add_option("-s","--smearing", dest="smearing", default=0, help="smearing of the calibrations")
    parser.add_option("-o","--output-file", dest="output", type="string", default="smearedCalibMap.root", help="modified calibration map")

    (options, args) = parser.parse_args()
    
    cs = CalibSmearer(args[0],options)
    cs.writeMaps(options)

