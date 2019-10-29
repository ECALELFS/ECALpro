import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

foldername = "AlCaP0_2016_ULrereco_from0_ext1_fromIter6"
niter = 0       # generally it starts from 0
eosPi0Folder = "piZero_Run2"
excludeMod2EBm16 = True
excludeMod4EBm06 = True if "2016" in foldername else False # might also have iphi in 41-60

outpath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/"
outfull = outpath + foldername + "/iter_" + str(niter) + "/2DMaps/ICmaps/IC_work/"

f1 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/{eosfd}/mciprian/{fd}/iter_{n}/{fd}_calibMap.root".format(eosfd=eosPi0Folder,
                                                                                                                           fd=foldername, 
                                                                                                                           n=niter)

print "="*30

cmd = "root -l -b -q 'makeICmapAndDump.C++(\"{out}\",\"{f1}\",\"dumpIC_norm1etaRing.dat\",\"calibMap_EB\"".format(out=outfull, 
                                                                                                                  f1=f1)
                                                                                           
cmd += ", 0.95, 1.05,true,false,0,{ex1},{ex2})'".format(ex1="true" if excludeMod2EBm16 else "false",
                                                        ex2="true" if excludeMod4EBm06 else "false")


print "-"*30
print cmd
os.system(cmd)


print "="*30
print "THE END"
print "="*30


