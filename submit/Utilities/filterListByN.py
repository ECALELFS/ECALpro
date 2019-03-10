#!/usr/bin/env python

import subprocess, time, sys, os, string
import ROOT

from optparse import OptionParser                                                    

parser = OptionParser(usage="%prog [options] filelist.txt")    
parser.add_option("-n", "--filter_every-n", dest="filterEveryN",  type="int",     default=3,   help="Will take 1 file every N specified by this option")
(options, args) = parser.parse_args()

if len(args)<1:
    parser.print_usage()
    quit()

print ""
fin = args[0]

fileListDir = os.path.dirname(fin)
if len(fileListDir):
   fileListDir += "/"
else:
   fileListDir = "./"

base = os.path.basename(fin).split('.')[0]
ext = os.path.basename(fin).split('.')[1]
fout = fileListDir + base + "_1every{n}.".format(n=options.filterEveryN) + ext

Filelist_f = open( fin )
Filelistbase_v = Filelist_f.readlines()
Filelist_f.close()

NEW_f = open( fout, 'w' )
NEW_f.write("# created filtering {f}: 1 file every {n}\n".format(f=args[0],n=options.filterEveryN))
ifile = -1 # keep -1, it becomes 0 for the first not-commented line
nCopied = 0

for f in Filelistbase_v:
    if f.startswith('#'): continue
    ifile += 1 # so the very first value is 0
    if ifile % options.filterEveryN: continue
    nCopied += 1
    NEW_f.write(f)
NEW_f.close()

print "-"*30
print "Copied {n} files into {f}".format(n=nCopied, f=fout)
print "-"*30
