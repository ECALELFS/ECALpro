#!/bin/env python

# need to be on lxplus, so to have EOS mounted

import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

script = "plotDistributions.py"
isPi0 = True
folder = "AlCaP0_AllRun2018_1fileEvery5_testCC2018"
itern = "0"
inputfile = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/{fld}/iter_{it}/{fld}_epsilonPlots.root".format(fld=folder, it=itern)
outdir="/afs/cern.ch/user/m/mciprian/www/pi0calib/test/plotDistributions/{fld}_iter{it}/".format(fld=folder, it=itern)
otherOptions = " -l 41 -e 13 "


regions = ["1EB", "2EB", "1EE", "2EE"]

scaleH2 = {"1EB" : -4,   
           "2EB" : -5, 
           "1EE" : -3, 
           "2EE" : -3
           }
# eta
# scaleH2 = {"1EB" : -2, #-4,  
#            "2EB" : -2, #-5,
#            "1EE" : -2, #-3,
#            "2EE" : -3
#           }

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

    print "-"*30
    cmdToRun = 'python {scr} {infile}'.format(scr=script, infile=inputfile)
    cmdToRun = cmdToRun + ' "g1Nxtal_afterCuts_region{reg}::leading #gamma" "g2Nxtal_afterCuts_region{reg}::trailing #gamma" '.format(reg=region)
    cmdToRun = cmdToRun + ' -o {out} -c "nXtalGamma_region{reg}" --xAxisTitle "number of crystals in 3x3" {opt} '.format(out=outdir, reg=region, opt=otherOptions)
    cmdToRun = cmdToRun + ' -t "{title}" '.format(title=title[region])
    os.system(cmdToRun)

    print "-"*30
    cmdToRun = 'python {scr} {infile}'.format(scr=script, infile=inputfile)
    cmdToRun = cmdToRun + ' "pi0pt_afterCuts_region{reg}" "PTPI0" '.format(reg=region)
    cmdToRun = cmdToRun + ' -o {out} -c "ptPi0_region{reg}" --xAxisTitle "{mes} p_{{T}} [GeV]" {opt} '.format(out=outdir, reg=region, 
                                                                                                              mes="#pi^{0}" if isPi0 else "#eta^{0}",
                                                                                                              opt=otherOptions)
    cmdToRun = cmdToRun + ' -t "{title}" '.format(title=title[region])

    os.system(cmdToRun)


print ""
print "THE END"
print ""
