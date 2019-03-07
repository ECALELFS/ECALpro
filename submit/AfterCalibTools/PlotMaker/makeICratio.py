import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

outpath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/ratioIC/"
#outdir = "AlCaP0_AllRun2017_condor_iter1__Over__AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6_iter6/"
outdir = "AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6_iter6__Over__AlCaP0_AllRun2017_condor_fixEBm16_iter2/"
canvasSuffix = "ratioIC"
label1 = "Prompt 2018 (9.8 fb^{-1 })"
label2 = "Legacy 2017 (40 fb^{-1 })"

f1 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6/iter_6/2DMaps/ICmaps/IC_work/calibrationMaps.root"
f2 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/AlCaP0_AllRun2017_condor_fixEBm16/iter_2/2DMaps/ICmaps/IC_work/calibrationMaps.root"

n1 = { "EB"  : "calibMap_EB",
       #"EEp" : "calibMap_EEp",
       #"EEm" : "calibMap_EEm"
       }

n2 = { "EB"  : "calibMap_EB",
       "EEp" : "calibMap_EEp",
       "EEm" : "calibMap_EEm"
       }

rangeIC = { "EB"  : [0.95, 1.05],
            "EEp" : [0.9, 1.1],
            "EEm" : [0.9, 1.1]
            }

numDet = { "EB"  : 0,
           "EEp" : 1,
           "EEm" : 2
           }

outfull = outpath + outdir

print "="*30

for det in n1:

    cmd = "root -l -b -q 'makeICratio.C++(\"{out}\",\"{cs}\",\"{f1}\",\"{f2}\",\"{n1}\",\"{n2}\",{minIC},{maxIC},{nDet},\"{l1}\",\"{l2}\")'".format(out=outfull, 
                                                                                                                                                    cs=canvasSuffix, 
                                                                                                                                                    f1=f1, f2=f2, 
                                                                                                                                                    n1=n1[det], n2=n2[det], 
                                                                                                                                                    minIC=rangeIC[det][0], 
                                                                                                                                                    maxIC=rangeIC[det][1],
                                                                                                                                                    nDet=numDet[det],
                                                                                                                                                    l1=label1,
                                                                                                                                                    l2=label2)
    print "-"*30
    print ">>> " + det
    print cmd
    os.system(cmd)


print "="*30
print "THE END"
print "="*30


