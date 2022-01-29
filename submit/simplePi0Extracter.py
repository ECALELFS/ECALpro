import subprocess, time, sys, os
from methods import *
from datetime import datetime

from optparse import OptionParser                              


fill_cfg_n = "try.py"
fill_cfg_f = open( fill_cfg_n, 'w' )

pwd = os.getcwd()

'''
###build Pi0
printFillCfg1(fill_cfg_f)
fill_cfg_f.write("        '" + "'\n")
printFillCfg2( fill_cfg_f, pwd, 0 , "/tmp/", 0 )
fill_cfg_f.close() 
'''


###fit pi0
fit_cfg_n = "fitEpsilonPlot_justFoldSM_try.py"
fit_cfg_f = open( fit_cfg_n, 'w' )
# print the cfg file
printFitCfg( fit_cfg_f , 0, "",0,0,"Barrel",0,justDoHistogramFolding=True)
fit_cfg_f.close()

