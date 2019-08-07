#!/usr/bin/env python

import subprocess, time, sys, os, string
import ROOT

from optparse import OptionParser                                                    

parser = OptionParser(usage="%prog [options] filelist.txt")    
parser.add_option("-n", "--filter_every-n", dest="filterEveryN",  type="int",     default=3,   help="Will take 1 file every N specified by this option")
parser.add_option("-i", "--invert", dest="invert",  action="store_true", default=False, help="Instead of taking 1 file every N, keep the other N-1 files")
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
fout = fileListDir + base + "_{i}every{n}{pfx}.".format(i=(options.filterEveryN-1) if options.invert else 1,
                                                        n=options.filterEveryN,
                                                        pfx="_inverted" if options.invert else "")
fout +=  ext

Filelist_f = open( fin )
Filelistbase_v = Filelist_f.readlines()
Filelist_f.close()

NEW_f = open( fout, 'w' )
NEW_f.write("# created filtering {f}: {i} file every {n}\n".format(f=args[0],
                                                                   i=(options.filterEveryN-1) if options.invert else 1,
                                                                   n=options.filterEveryN))
ifile = -1 # keep -1, it becomes 0 for the first not-commented line
nCopied = 0

for f in Filelistbase_v:
    if f.startswith('#'): continue
    ifile += 1 # so the very first value is 0
    if options.invert:
        if ifile % options.filterEveryN: 
            NEW_f.write(f)
            nCopied += 1
        else:
            continue
    else:
        if ifile % options.filterEveryN:
            continue
        else:
            NEW_f.write(f)
            nCopied += 1

NEW_f.close()

print "-"*30
print "Copied {n} files into {f}".format(n=nCopied, f=fout)
print "-"*30
