#!/usr/bin/env python

import subprocess, time, sys, os, string

#######
#fileList: list of files
#fileJson: JSON file
#fileNEW: OutPut: the original list cleaned removing the file not in the JSON file
######

if len(sys.argv)<3:
   print "Usage: Purify_List.py filelist.txt json.txt"
   exit(0)

#file name
fileList = sys.argv[1]
if not( os.path.isfile(fileList) ):
   print "WARNING!!! " + str(fileList) + " not found!"
fileJson = sys.argv[2]
if not( os.path.isfile(fileJson) ):
   print "WARNING!!! " + str(fileJson) + " not found!"

fileNEW = "purified_"+os.path.basename(fileList)

if ( os.path.isfile(fileNEW) ):
   os.remove(fileNEW)
#open
Filelist_f = open( fileList )
Jsonlist_f = open( fileJson )
NEW_f = open( fileNEW, 'w' )
#Read
Filelistbase_v = Filelist_f.readlines()
Jsonlistbase_v = Jsonlist_f.readlines()

for Nline in range(len(Filelistbase_v)):
  IsThere=False
  line = Filelistbase_v[Nline]
  num = line.index('000') #assume .../v1/000/251/028/...
  newLine  = line[int(num+4):int(num+7)]
  newLine += line[int(num+8):int(num)+11]
  #print "Look For: " + str(newLine)
  for NlineJson in range(len(Jsonlistbase_v)):
      JsonLine = str(Jsonlistbase_v[NlineJson]).strip('\n')
      if( string.find(str(JsonLine),str(newLine))>0 ): IsThere=True
  if(IsThere):
     #print "There is!"
     NEW_f.write(Filelistbase_v[Nline])
  #else:
     #print "There isn't."

print "---THE END---"
