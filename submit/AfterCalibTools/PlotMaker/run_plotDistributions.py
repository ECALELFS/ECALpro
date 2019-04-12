#!/bin/env python

# need to be on lxplus, so to have EOS mounted

import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

script = "plotDistributions.py"
folder = "AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC"
itern = "0"
inputfile = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/{fld}/iter_{it}/{fld}_epsilonPlots.root".format(fld=folder, it=itern)
outdir="/afs/cern.ch/user/m/mciprian/www/pi0calib/test/plotDistributions/{fld}_iter{it}/".format(fld=folder, it=itern)
otherOptions = " -l 9.8 -e 13 "

regions = ["1EB", "2EB", "1EE", "2EE"]

scaleH2 = {"1EB" : -4,
           "2EB" : -4,
           "1EE" : -3,
           "2EE" : -2
           }

title = {"1EB" : "EB: |#eta| < 1.0",
         "2EB" : "EB: |#eta| > 1.0",
         "1EE" : "EE: |#eta| < 1.8",
         "2EE" : "EE: |#eta| > 1.8"
         }

print ""


for region in regions:
    print ""
    print "region: " + region
    print "-"*30
    cmdToRun = 'python {scr} {infile}'.format(scr=script, infile=inputfile)
    cmdToRun = cmdToRun + ' "g1pt_afterCuts_region{reg}::leading #gamma" "g2pt_afterCuts_region{reg}::trailing #gamma" '.format(reg=region)
    cmdToRun = cmdToRun + ' -o {out} -c "ptGamma_region{reg}" --xAxisTitle "photon p_{{T}} [GeV]" {opt} '.format(out=outdir, reg=region, opt=otherOptions)
    cmdToRun = cmdToRun + ' --scale-hist2 {scale} -t "{title}" '.format(scale=scaleH2[region], title=title[region])

    os.system(cmdToRun)

print ""
print "THE END"
print ""
