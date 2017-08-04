#!/usr/bin/env python

import subprocess, time, sys, os, string

#######
#fileList: list of files
#fileJson: JSON file
#fileNEW: OutPut: the original list cleaned removing the file not in the JSON file
######

#file name
#fileList = '../InputList/data_HLTPhysics1_Run2016H-v1_run283685_RAW.txt'
fileList = '../InputList/AlCaP0_2017_upTo31July2017.list'
if not( os.path.isfile(fileList) ):
   print "WARNING!!! " + str(fileList) + " not found!"
#fileJson = '../../FillEpsilonPlot/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
fileJson = '../../FillEpsilonPlot/data/json_DCSONLY.txt'
if not( os.path.isfile(fileJson) ):
   print "WARNING!!! " + str(fileJson) + " not found!"
#fileNEW = '../InputList/data_HLTPhysics1_Run2016H-v1_run283685_RAW_purified_with_lastJSON2016.list'
fileNEW = '../InputList/AlCaP0_2017_upTo31July2017_purified_json_DCSONLY.list'
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
