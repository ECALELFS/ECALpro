#!/usr/bin/env python

import subprocess, time, sys, os, string

#######
#fileRun: list of RUN number ex: 197601/n 197602/n 197603/n...
#fileList: list of files
#fileJson: JSON file
#fileNEW: OutPut: the original list cleaned removing the file not in the JSON file
######

#file name
fileRun = 'ALL_2012D_RUN.list'
fileList = 'ALL_2012D.list'
fileJson = 'common/goodrunlist_json2012D.txt'
fileNEW = 'ALL_2012D_good.list'
#open
Runlist_f = open( fileRun )
Filelist_f = open( fileList )
Jsonlist_f = open( fileJson )
NEW_f = open( fileNEW, 'w' )
#Read
Runlistbase_v = Runlist_f.readlines()
Filelistbase_v = Filelist_f.readlines()
Jsonlistbase_v = Jsonlist_f.readlines()

for Nline in range(len(Runlistbase_v)):
  IsThere=False
  line = Runlistbase_v[Nline]
  newLine = line.strip('\n')
  newLine = '"' + str(newLine) + '":'
  print "Look For: " + str(newLine)
  for NlineJson in range(len(Jsonlistbase_v)):
      if( string.find(str(Jsonlistbase_v[NlineJson]),str(newLine))>0 ): IsThere=True

  if(IsThere):
     print "There is!"
     NEW_f.write(Filelistbase_v[Nline])
  else:
     print "There isn't."

print "---THE END---"
