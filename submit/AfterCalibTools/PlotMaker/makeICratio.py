import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

outpath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/ratioIC/"
#outdir = "AlCaP0_AllRun2017_condor_iter1__Over__AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6_iter6/"
outdir = "TESTmemoryOptimWithTH2_AlCaP0_2018_ULrereco_thirdThirdNotUsedFor2018IC_fromIC2018_iter0__Over_AlCaP0_2018_ULrereco_thirdThirdNotUsedFor2018IC_fromIC2018_iter0/"
canvasSuffix = "ratioIC"
#label1 = "2018 UL (half 2018, ~30 fb^{-1 })"
#label2 = "2018 test (other half, ~30 fb^{-1 })"
label1 = "2018, third ~10 fb^{-1 } not used for UL (Memory Optim)"
label2 = "2018, third ~10 fb^{-1 } not used for UL"

f1 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/TESTmemoryOptimWithTH2_AlCaP0_2018_ULrereco_thirdThirdNotUsedFor2018IC_fromIC2018/iter_0/2DMaps/ICmaps/IC_work/calibrationMaps.root"
f2 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/AlCaP0_2018_ULrereco_thirdThirdNotUsedFor2018IC_fromIC2018/iter_0/2DMaps/ICmaps/IC_work/calibrationMaps.root"

n1 = { "EB"  : "calibMap_EB",
       "EEp" : "calibMap_EEp",
       "EEm" : "calibMap_EEm"
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


