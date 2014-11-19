#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

if len(sys.argv) != 7:
    print "usage thisPyton.py pwd queue iter-to-resubmit Systparam onlyFIT onlyFinalHADD"
    sys.exit(1)

pwd = sys.argv[1]
queue = sys.argv[2]
Iter_toResub = int(sys.argv[3])
SystParam=int(sys.argv[4])
onlyFIT=str(sys.argv[5])
onlyFinalHadd=str(sys.argv[6])

outputdir = pwd+'/'+dirname
logPath = outputdir + '/log'
srcPath  = outputdir + '/src'

# To compute the num of hadd
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = inputlist_f.readlines()
#N real Job
njobs = float(len(inputlistbase_v[:])) / float(ijobmax)
if( float(int(njobs) - njobs) < 0. ):
    njobs = int(njobs) + 1

for iters in range(Iter_toResub, nIterations):
   if(onlyFIT=='False' and onlyFinalHadd=='False'):
      print "\n*******  ITERATION " + str(iters) + "/" + str(nIterations-1) + "  *******"
      for ijob in range(njobs):

          # preparing submission of filling tasks
          env_script_n = open(outputdir + "/cfgFile/Fill/fillEpsilonPlot_iter_"     + str(iters) + "_job_" + str(ijob) + ".py", 'a')
          SystParamLine = 'process.analyzerFillEpsilon.SystOrNot = cms.untracked.double(' + str(SystParam) + ')\n'
          env_script_n.write(SystParamLine)
          env_script_n.close()
          fill_log_n = logPath + "/fillEpsilonPlot_iter_" + str(iters) + "_job_" + str(ijob) + ".log"
          fill_src_n = srcPath + "/Fill/submit_iter_"     + str(iters) + "_job_" + str(ijob) + ".sh"
          if not(Silent): 
               submit_s = "bsub -q " + queue + " -o " + fill_log_n + " " + fill_src_n
          else:
               submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fill_src_n

          print '\n[job #' + str(ijob) + '] :: ' + submit_s

          # actually submitting filling tasks
          submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
          output = submitJobs.communicate()
          print output

          # avoid overlapping submission
          time.sleep(1)

      checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
      datalines = (checkJobs.communicate()[0]).splitlines()

      print 'Waiting for filling jobs to be finished...'

      # Daemon cheking running jobs
      while len(datalines)>=2 :
          for entry in datalines:
              entry = entry.rstrip()
              entry = entry.split()[0]
              #print entry
              if(entry.find('JOBID')!=-1): continue
              i = int(entry)

          time.sleep(1)

          checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
          datalines = (checkJobs.communicate()[0]).splitlines()

      print 'Done with loop'
      print 'Now adding files...'

      # Computing Nunber of hadd
      inputlist_v = inputlistbase_v[:]
      NrelJob = float(len(inputlist_v)) / float(ijobmax)
      if( float(int(NrelJob) - NrelJob) < 0. ):
          NrelJob = int(NrelJob) + 1
      Nlist_flo = float(NrelJob/nHadd) + 1.
      Nlist = int(Nlist_flo)
      print Nlist
      # hadd to sum the epsilon histograms
      for nHadds in range(Nlist):
          Hadd_src_n = srcPath + "/hadd/HaddCfg_iter_" + str(iters) + "_job_" + str(nHadds) + ".sh"
          Hadd_log_n = logPath + "/HaddCfg_iter_" + str(iters) + "_job_" + str(nHadds) + ".log"
          Hsubmit_s = "bsub -q " + queue + " -o " + Hadd_log_n + " bash " + Hadd_src_n
          #Before each HADD we need ot check if the all the files in the list are present
          Grepcommand = "grep -i list " + Hadd_src_n + " | grep -v echo | awk '{print $4}'"
          myGrep = subprocess.Popen([Grepcommand], stdout=subprocess.PIPE, shell=True )
          FoutGrep = myGrep.communicate()
          FoutGrep_2 = str(FoutGrep)[3:]
          FoutGrep_2 = str(FoutGrep_2)[:-10]
          print 'Checking ' + str(FoutGrep_2)
          #Chech The size for each line
          f = open( str(FoutGrep_2) )
          lines = f.readlines()
          f.close()
          NumToRem = 0
          for filetoCheck in lines:
              if( NumToRem!=0 ):
                 Num = NumToRem - 1
                 f2 = open(str(FoutGrep_2) + str(Num))
                 lines = f2.readlines()
                 f2.close()
              filetoCheck2 = str(filetoCheck)[22:]
              CheckComm = 'cmsLs -l ' + str(filetoCheck2)
              myCheck =  subprocess.Popen([CheckComm], stdout=subprocess.PIPE, shell=True )
              Check_output = myCheck.communicate()
              #If file is not present, remove it from the list
              if "No such" in str(Check_output):
                 print 'HADD::MISSING: ' + str(filetoCheck2)
                 print 'removing from Hadd, in: ' + str(FoutGrep_2) + str(NumToRem)
                 f1 = open(str(FoutGrep_2) + str(NumToRem),"w")
                 NumToRem = NumToRem + 1
                 for line in lines:
                     if line!=str(filetoCheck):
                          f1.write(line)
                 f1.close()
              else:
                 Splitted =  str(Check_output).split( );
                 print "size: " + str(Splitted[1])
                 #If is corrupted (size too small), remove it from the list
                 if( int(Splitted[1])<10000 ):
                      print 'HADD::Bad size for: ' + str(filetoCheck2)
                      print 'removing from Hadd'
                      f1 = open(str(FoutGrep_2) + str(NumToRem),"w")
                      NumToRem = NumToRem + 1
                      lines1 = f1.readlines()
                      for line in lines:
                          if line!=str(filetoCheck):
                               f1.write(line)
                      f1.close()
          #moving the .list to the correct one
          if( NumToRem!=0 ):
              NumToRem = NumToRem - 1
              MoveComm = "cp " + str(FoutGrep_2) + str(NumToRem) + " " + str(FoutGrep_2)
              MoveC = subprocess.Popen([MoveComm], stdout=subprocess.PIPE, shell=True);
              mvOut = MoveC.communicate()
          #End of the check
          subJobs = subprocess.Popen([Hsubmit_s], stdout=subprocess.PIPE, shell=True);
          outJobs = subJobs.communicate()
          print outJobs
          time.sleep(1)

      checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
      datalines = (checkJobs.communicate()[0]).splitlines()

      print 'Waiting for all the hadd...'

      # Daemon cheking running jobs
      while len(datalines)>=2 :
          for entry in datalines:
              entry = entry.rstrip()
              entry = entry.split()[0]
              #print entry
              if(entry.find('JOBID')!=-1): continue
              i = int(entry)

          time.sleep(1)

          checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
          datalines = (checkJobs.communicate()[0]).splitlines()

      print 'Done with various hadd'

      print 'Now The Final One...'
      FHadd_src_n = srcPath + "/hadd/Final_HaddCfg_iter_" + str(iters) + ".sh"
      FHadd_log_n = logPath + "/Final_HaddCfg_iter_" + str(iters) + ".log"
      FHsubmit_s = "bsub -q " + queue + " -o " + FHadd_log_n + " bash " + FHadd_src_n
      FsubJobs = subprocess.Popen([FHsubmit_s], stdout=subprocess.PIPE, shell=True);
      FoutJobs = FsubJobs.communicate()
      print FoutJobs
      time.sleep(3)

      FcheckJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
      Fdatalines = (FcheckJobs.communicate()[0]).splitlines()
      print 'Waiting for the Final hadd...'
      # Daemon cheking running jobs
      while len(Fdatalines)>=2 :
          for entry in Fdatalines:
              entry = entry.rstrip()
              entry = entry.split()[0]
              #print entry
              if(entry.find('JOBID')!=-1): continue
              i = int(entry)

          time.sleep(1)

          FcheckJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
          Fdatalines = (FcheckJobs.communicate()[0]).splitlines()

      print 'Done with final hadd'

      print 'Done with staging the final epsilonPlots.root'

      # removing useles file
#      for nRm in range(Nlist):
#          remove_s = 'cmsRm ' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/'+ NameTag +'epsilonPlots_' + str(nRm) + '.root'
#          print '[Removed] :: ' + remove_s
#          removeFile = subprocess.Popen([remove_s],stdout=subprocess.PIPE, shell=True)
#          nFilesRemoved = 0
#          filesRemoved = (removeFile.communicate()[0]).splitlines()

   # N of Fit to send
   nEB = 61199/nFit
   if (61199%nFit != 0) :
       nEB = int(nEB) +1
   nEE = 14647/nFit
   if (14647%nFit != 0) :
       nEE = int(nEE) +1
   # For final hadd
   ListFinaHadd = list()
   # preparing submission of fit tasks (EB)
   print 'Submitting ' + str(nEB) + ' jobs to fit the Barrel'
   for inteb in range(nEB):
       #fit_log_n = logPath + "/fitEpsilonPlot_EB_" + str(inteb) + "_iter_" + str(iters) + ".log"
       fit_src_n = srcPath + "/Fit/submit_EB_" + str(inteb) + "_iter_"     + str(iters) + ".sh"
       #submit_s = "bsub -q " + queue + " -o " + fit_log_n + " source " + fit_src_n
       submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fit_src_n
       ListFinaHadd.append('root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName )
       if(onlyFinalHadd=='False'):
         print 'About to EB fit:'
         print 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName
         print submit_s
         # actually submitting fit tasks (EB)
         submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
         output = submitJobs.communicate()
         print output

   # preparing submission of fit tasks (EE)
   print 'Submitting ' + str(nEE) + ' jobs to fit the Endcap'
   for inte in range(nEE):
       #fit_log_n = logPath + "/fitEpsilonPlot_EE_" + str(inte) + "_iter_" + str(iters) + ".log"
       fit_src_n = srcPath + "/Fit/submit_EE_" + str(inte) + "_iter_"     + str(iters) + ".sh"
       # redirecting output on /dev/null is required to avoid LSF directory to be created automatically
       submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fit_src_n
       ListFinaHadd.append('root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName)
       if(onlyFinalHadd=='False'):
         print 'About to EE fit:'
         print 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName
         print submit_s
         # actually submitting fit tasks (EE)
         submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
         output = submitJobs.communicate()
         print output

   if(onlyFinalHadd=='False'):
      # checking number of running/pending jobs
      checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
      datalines = (checkJobs.communicate()[0]).splitlines()
   
      print 'Waiting for fit jobs to be finished...'
   
      #Daemon cheking running jobs
      while len(datalines)>=2 :
          for entry in datalines:
              entry = entry.rstrip()
              entry = entry.split()[0]
              #print entry
              if(entry.find('JOBID')!=-1): continue
              i = int(entry)
   
          time.sleep(1)
   
          checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
          datalines = (checkJobs.communicate()[0]).splitlines()
   
      print "Done with fitting! Now we have to merge all fits in one Calibmap.root"

   # Merge Final CalibMap
   from ROOT import *
   from PhysicsTools.PythonAnalysis import *
   gSystem.Load("libFWCoreFWLite.so")
   AutoLibraryLoader.enable()
   f = TFile('/tmp/' + NameTag + calibMapName, 'recreate')

   # Create a struct
   gROOT.ProcessLine(\
     "struct EBStruct{\
       Int_t rawId_;\
       Int_t hashedIndex_;\
       Int_t ieta_;\
       Int_t iphi_;\
       Int_t iSM_;\
       Int_t iMod_;\
       Int_t iTT_;\
       Int_t iTTeta_;\
       Int_t iTTphi_;\
       Int_t iter_;\
       Double_t coeff_;\
       Double_t Signal_;\
       Double_t Backgr_;\
       Double_t Chisqu_;\
       Double_t Ndof_;\
       Double_t fit_mean_;\
       Double_t fit_sigma_;\
       Double_t fit_Snorm_;\
       Double_t fit_b0_;\
       Double_t fit_b1_;\
       Double_t fit_b2_;\
       Double_t fit_b3_;\
       Double_t fit_Bnorm_;\
     };")
   gROOT.ProcessLine(\
     "struct EEStruct{\
       Int_t ix_;\
       Int_t iy_;\
       Int_t zside_;\
       Int_t sc_;\
       Int_t isc_;\
       Int_t ic_;\
       Int_t iquadrant_;\
       Int_t hashedIndex_;\
       Int_t iter_;\
       Double_t coeff_;\
       Double_t Signal_;\
       Double_t Backgr_;\
       Double_t Chisqu_;\
       Double_t Ndof_;\
       Double_t fit_mean_;\
       Double_t fit_sigma_;\
       Double_t fit_Snorm_;\
       Double_t fit_b0_;\
       Double_t fit_b1_;\
       Double_t fit_b2_;\
       Double_t fit_b3_;\
       Double_t fit_Bnorm_;\
     };")
   s = EBStruct()
   t = EEStruct()

   gROOT.ProcessLine(\
     "struct EB1Struct{\
       Int_t rawId;\
       Int_t hashedIndex;\
       Int_t ieta;\
       Int_t iphi;\
       Int_t iSM;\
       Int_t iMod;\
       Int_t iTT;\
       Int_t iTTeta;\
       Int_t iTTphi;\
       Int_t iter;\
       Double_t coeff;\
       Double_t Signal;\
       Double_t Backgr;\
       Double_t Chisqu;\
       Double_t Ndof;\
       Double_t fit_mean;\
       Double_t fit_sigma;\
       Double_t fit_Snorm;\
       Double_t fit_b0;\
       Double_t fit_b1;\
       Double_t fit_b2;\
       Double_t fit_b3;\
       Double_t fit_Bnorm;\
     };")
   gROOT.ProcessLine(\
     "struct EE1Struct{\
       Int_t ix;\
       Int_t iy;\
       Int_t zside;\
       Int_t sc;\
       Int_t isc;\
       Int_t ic;\
       Int_t iquadrant;\
       Int_t hashedIndex;\
       Int_t iter;\
       Double_t coeff;\
       Double_t Signal;\
       Double_t Backgr;\
       Double_t Chisqu;\
       Double_t Ndof;\
       Double_t fit_mean;\
       Double_t fit_sigma;\
       Double_t fit_Snorm;\
       Double_t fit_b0;\
       Double_t fit_b1;\
       Double_t fit_b2;\
       Double_t fit_b3;\
       Double_t fit_Bnorm;\
     };")
   s1 = EB1Struct()
   t1 = EE1Struct()

   calibMap_EB = TH2F("calibMap_EB", "EB calib coefficients: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
   calibMap_EEm = TH2F("calibMap_EEm", "EE- calib coefficients", 100,0.5,100.5,100,0.5,100.5)
   calibMap_EEp = TH2F("calibMap_EEp", "EE+ calib coefficients", 100,0.5,100.5,100,0.5,100.5)
   TreeEB = TTree("calibEB", "Tree of EB Inter-calibration constants")
   TreeEB.Branch('rawId_'      , AddressOf(s,'rawId_'),'rawId_/I')
   TreeEB.Branch('hashedIndex_', AddressOf(s,'hashedIndex_'),'hashedIndex_/I')
   TreeEB.Branch('ieta_'       , AddressOf(s,'ieta_'),'ieta_/I')
   TreeEB.Branch('iphi_'       , AddressOf(s,'iphi_'),'iphi_/I')
   TreeEB.Branch('iSM_'        , AddressOf(s,'iSM_'),'iSM_/I')
   TreeEB.Branch('iMod_'       , AddressOf(s,'iMod_'),'iMod_/I')
   TreeEB.Branch('iTT_'        , AddressOf(s,'iTT_'),'iTT_/I')
   TreeEB.Branch('iTTeta_'     , AddressOf(s,'iTTeta_'),'iTTeta_/I')
   TreeEB.Branch('iTTphi_'     , AddressOf(s,'iTTphi_'),'iTTphi_/I')
   TreeEB.Branch('iter_'       , AddressOf(s,'iter_'),'iter_/I')
   TreeEB.Branch('coeff_'      , AddressOf(s,'coeff_'),'coeff_/F')
   TreeEB.Branch('Signal_'     , AddressOf(s,'Signal_'),'Signal_/F')
   TreeEB.Branch('Backgr_'     , AddressOf(s,'Backgr_'),'Backgr_/F')
   TreeEB.Branch('Chisqu_'     , AddressOf(s,'Chisqu_'),'Chisqu_/F')
   TreeEB.Branch('Ndof_'       , AddressOf(s,'Ndof_'),'Ndof_/F')
   TreeEB.Branch('fit_mean_'   , AddressOf(s,'fit_mean_'),'fit_mean_/F')
   TreeEB.Branch('fit_sigma_'  , AddressOf(s,'fit_sigma_'),'fit_sigma_/F')
   TreeEB.Branch('fit_Snorm_'  , AddressOf(s,'fit_Snorm_'),'fit_Snorm_/F')
   TreeEB.Branch('fit_b0_'     , AddressOf(s,'fit_b0_'),'fit_b0_/F')
   TreeEB.Branch('fit_b1_'     , AddressOf(s,'fit_b1_'),'fit_b1_/F')
   TreeEB.Branch('fit_b2_'     , AddressOf(s,'fit_b2_'),'fit_b2_/F')
   TreeEB.Branch('fit_b3_'     , AddressOf(s,'fit_b3_'),'fit_b3_/F')
   TreeEB.Branch('fit_Bnorm_'  , AddressOf(s,'fit_Bnorm_'),'fit_Bnorm_/F')

   TreeEE = TTree("calibEE", "Tree of EE Inter-calibration constants")
   TreeEE.Branch('ix_'         , AddressOf(t,'ix_'),'ix_/I')
   TreeEE.Branch('iy_'         , AddressOf(t,'iy_'),'iy_/I')
   TreeEE.Branch('zside_'      , AddressOf(t,'zside_'),'zside_/I')
   TreeEE.Branch('sc_'         , AddressOf(t,'sc_'),'sc_/I')
   TreeEE.Branch('isc_'        , AddressOf(t,'isc_'),'isc_/I')
   TreeEE.Branch('ic_'         , AddressOf(t,'ic_'),'ic_/I')
   TreeEE.Branch('iquadrant_'  , AddressOf(t,'iquadrant_'),'iquadrant_/I')
   TreeEE.Branch('hashedIndex_', AddressOf(t,'hashedIndex_'),'hashedIndex_/I')
   TreeEE.Branch('iter_'       , AddressOf(t,'iter_'),'iter_/I')
   TreeEE.Branch('coeff_'      , AddressOf(t,'coeff_'),'coeff_/F')
   TreeEE.Branch('Signal_'     , AddressOf(t,'Signal_'),'Signal_/F')
   TreeEE.Branch('Backgr_'     , AddressOf(t,'Backgr_'),'Backgr_/F')
   TreeEE.Branch('Chisqu_'     , AddressOf(t,'Chisqu_'),'Chisqu_/F')
   TreeEE.Branch('Ndof_'       , AddressOf(t,'Ndof_'),'Ndof_/F')
   TreeEE.Branch('fit_mean_'   , AddressOf(t,'fit_mean_'),'fit_mean_/F')
   TreeEE.Branch('fit_sigma_'  , AddressOf(t,'fit_sigma_'),'fit_sigma_/F')
   TreeEE.Branch('fit_Snorm_'  , AddressOf(t,'fit_Snorm_'),'fit_Snorm_/F')
   TreeEE.Branch('fit_b0_'     , AddressOf(t,'fit_b0_'),'fit_b0_/F')
   TreeEE.Branch('fit_b1_'     , AddressOf(t,'fit_b1_'),'fit_b1_/F')
   TreeEE.Branch('fit_b2_'     , AddressOf(t,'fit_b2_'),'fit_b2_/F')
   TreeEE.Branch('fit_b3_'     , AddressOf(t,'fit_b3_'),'fit_b3_/F')
   TreeEE.Branch('fit_Bnorm_'  , AddressOf(t,'fit_Bnorm_'),'fit_Bnorm_/F')

   for thisfile_s in ListFinaHadd:
       thisfile_s = thisfile_s.rstrip()
       print thisfile_s
       thisfile_f = TFile.Open(thisfile_s)
       #Taking Interval and EB or EE
       h_Int = thisfile_f.Get("hint")
       init = h_Int.GetBinContent(1)
       finit = h_Int.GetBinContent(2)
       EEoEB = h_Int.GetBinContent(3)

       #TTree
       if EEoEB == 0:
          thisTree = thisfile_f.Get("calibEB")
          thisTree.SetBranchAddress( 'rawId',AddressOf(s1,'rawId'));
          thisTree.SetBranchAddress( 'hashedIndex',AddressOf(s1,'hashedIndex'));
          thisTree.SetBranchAddress( 'ieta',AddressOf(s1,'ieta'));
          thisTree.SetBranchAddress( 'iphi',AddressOf(s1,'iphi'));
          thisTree.SetBranchAddress( 'iSM',AddressOf(s1,'iSM'));
          thisTree.SetBranchAddress( 'iMod',AddressOf(s1,'iMod'));
          thisTree.SetBranchAddress( 'iTT',AddressOf(s1,'iTT'));
          thisTree.SetBranchAddress( 'iTTeta',AddressOf(s1,'iTTeta'));
          thisTree.SetBranchAddress( 'iTTphi',AddressOf(s1,'iTTphi'));
          thisTree.SetBranchAddress( 'iter',AddressOf(s1,'iter'));
          thisTree.SetBranchAddress( 'coeff',AddressOf(s1,'coeff'));
          thisTree.SetBranchAddress( 'Signal',AddressOf(s1,'Signal'));
          thisTree.SetBranchAddress( 'Backgr',AddressOf(s1,'Backgr'));
          thisTree.SetBranchAddress( 'Chisqu',AddressOf(s1,'Chisqu'));
          thisTree.SetBranchAddress( 'Ndof',AddressOf(s1,'Ndof'));
          thisTree.SetBranchAddress( 'fit_mean',AddressOf(s1,'fit_mean'));
          thisTree.SetBranchAddress( 'fit_sigma',AddressOf(s1,'fit_sigma'));
          thisTree.SetBranchAddress( 'fit_Snorm',AddressOf(s1,'fit_Snorm'));
          thisTree.SetBranchAddress( 'fit_b0',AddressOf(s1,'fit_b0'));
          thisTree.SetBranchAddress( 'fit_b1',AddressOf(s1,'fit_b1'));
          thisTree.SetBranchAddress( 'fit_b2',AddressOf(s1,'fit_b2'));
          thisTree.SetBranchAddress( 'fit_b3',AddressOf(s1,'fit_b3'));
          thisTree.SetBranchAddress( 'fit_Bnorm',AddressOf(s1,'fit_Bnorm'));
          for ntre in range(thisTree.GetEntries()):
              thisTree.GetEntry(ntre);
              if (ntre>=init and ntre<=finit):
                  s.rawId_ = s1.rawId
                  s.hashedIndex_ = s1.hashedIndex
                  s.ieta_ = s1.ieta
                  s.iphi_ = s1.iphi
                  s.iSM_ = s1.iSM
                  s.iMod_ = s1.iMod
                  s.iTT_ = s1.iTT
                  s.iTTeta_ = s1.iTTeta
                  s.iTTphi_ = s1.iTTphi
                  s.iter_ = s1.iter
                  s.coeff_ = s1.coeff
                  s.Signal_ = s1.Signal
                  s.Backgr_ = s1.Backgr
                  s.Chisqu_ = s1.Chisqu
                  s.Ndof_ = s1.Ndof
                  s.fit_mean_ = s1.fit_mean
                  s.fit_sigma_ = s1.fit_sigma
                  s.fit_Snorm_ = s1.fit_Snorm
                  s.fit_b0_ = s1.fit_b0
                  s.fit_b1_ = s1.fit_b1
                  s.fit_b2_ = s1.fit_b2
                  s.fit_b3_ = s1.fit_b3
                  s.fit_Bnorm_ = s1.fit_Bnorm
                  TreeEB.Fill()
       else:
          thisTree = thisfile_f.Get("calibEE")
          thisTree.SetBranchAddress( 'ix',AddressOf(t1,'ix'));
          thisTree.SetBranchAddress( 'iy',AddressOf(t1,'iy'));
          thisTree.SetBranchAddress( 'zside',AddressOf(t1,'zside'));
          thisTree.SetBranchAddress( 'sc',AddressOf(t1,'sc'));
          thisTree.SetBranchAddress( 'isc',AddressOf(t1,'isc'));
          thisTree.SetBranchAddress( 'ic',AddressOf(t1,'ic'));
          thisTree.SetBranchAddress( 'iquadrant',AddressOf(t1,'iquadrant'));
          thisTree.SetBranchAddress( 'hashedIndex',AddressOf(t1,'hashedIndex'));
          thisTree.SetBranchAddress( 'iter',AddressOf(t1,'iter'));
          thisTree.SetBranchAddress( 'coeff',AddressOf(t1,'coeff'));
          thisTree.SetBranchAddress( 'Signal',AddressOf(t1,'Signal'));
          thisTree.SetBranchAddress( 'Backgr',AddressOf(t1,'Backgr'));
          thisTree.SetBranchAddress( 'Chisqu',AddressOf(t1,'Chisqu'));
          thisTree.SetBranchAddress( 'Ndof',AddressOf(t1,'Ndof'));
          thisTree.SetBranchAddress( 'fit_mean',AddressOf(t1,'fit_mean'));
          thisTree.SetBranchAddress( 'fit_sigma',AddressOf(t1,'fit_sigma'));
          thisTree.SetBranchAddress( 'fit_Snorm',AddressOf(t1,'fit_Snorm'));
          thisTree.SetBranchAddress( 'fit_b0',AddressOf(t1,'fit_b0'));
          thisTree.SetBranchAddress( 'fit_b1',AddressOf(t1,'fit_b1'));
          thisTree.SetBranchAddress( 'fit_b2',AddressOf(t1,'fit_b2'));
          thisTree.SetBranchAddress( 'fit_b3',AddressOf(t1,'fit_b3'));
          thisTree.SetBranchAddress( 'fit_Bnorm',AddressOf(t1,'fit_Bnorm'));
          for ntre in range(thisTree.GetEntries()):
              thisTree.GetEntry(ntre);
              if (ntre>=init and ntre<=finit):
                  t.ix_ = t1.ix
                  t.iy_ = t1.iy
                  t.zside_ = t1.zside
                  t.sc_ = t1.sc
                  t.isc_ = t1.isc
                  t.ic = t1.ic
                  t.iquadrant_ = t1.iquadrant
                  t.hashedIndex_ = t1.hashedIndex
                  t.iter_ = t1.iter
                  t.coeff_ = t1.coeff
                  t.Signal_ = t1.Signal
                  t.Backgr_ = t1.Backgr
                  t.Chisqu_ = t1.Chisqu
                  t.Ndof_ = t1.Ndof
                  t.fit_mean_ = t1.fit_mean
                  t.fit_sigma_ = t1.fit_sigma
                  t.fit_Snorm_ = t1.fit_Snorm
                  t.fit_b0_ = t1.fit_b0
                  t.fit_b1_ = t1.fit_b1
                  t.fit_b2_ = t1.fit_b2
                  t.fit_b3_ = t1.fit_b3
                  t.fit_Bnorm_ = t1.fit_Bnorm
                  TreeEE.Fill()
       #TH2
       thisHistoEB = thisfile_f.Get("calibMap_EB")
       thisHistoEEm = thisfile_f.Get("calibMap_EEm")
       thisHistoEEp = thisfile_f.Get("calibMap_EEp")
       if EEoEB == 0:
          MaxEta = 85
          Init = int(init)
          Fin = int(finit+1)
          for nFitB in range(Init,Fin):
              if nFitB < 61200:
                 myRechit = EBDetId( EBDetId.detIdFromDenseIndex(nFitB) )
                 bin_x = myRechit.ieta()+MaxEta+1
                 bin_y = myRechit.iphi()
                 value = thisHistoEB.GetBinContent(bin_x,bin_y)
                 calibMap_EB.SetBinContent(bin_x,bin_y,value)
       else :
          Init1 = int(init)
          Fin1 = int(finit+1)
          for nFitE in range(Init1,Fin1):
              if nFitE < 14648:
                 myRechitE = EEDetId( EEDetId.detIdFromDenseIndex(nFitE) )
                 if myRechitE.zside() < 0 :
                    value = thisHistoEEm.GetBinContent(myRechitE.ix(),myRechitE.iy())
                    calibMap_EEm.SetBinContent(myRechitE.ix(),myRechitE.iy(),value)
                 if myRechitE.zside() > 0 :
                    value = thisHistoEEp.GetBinContent(myRechitE.ix(),myRechitE.iy())
                    calibMap_EEp.SetBinContent(myRechitE.ix(),myRechitE.iy(),value)

       thisfile_f.Close()
   f.cd()
   f.Write()
   f.Close()

   print 'Now staging calibMap.root on EOS'
   stage_s_fin = 'cmsStage /tmp/' + NameTag + calibMapName + ' ' + eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + NameTag + calibMapName
   print stage_s_fin
   stageCalibFile = subprocess.Popen([stage_s_fin], stdout=subprocess.PIPE, shell=True);
   print stageCalibFile.communicate()
   print 'Done with staging the final calibMap.root'

   # checking that calibMap.root is actually available on EOS
   print "Checking availabilty of " + NameTag + calibMapName
   checkFileAvailability_s = 'cmsLs ' + eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + NameTag + calibMapName
   print checkFileAvailability_s
   checkFileAvailability = subprocess.Popen([checkFileAvailability_s], stdout=subprocess.PIPE, shell=True);
   output = checkFileAvailability.communicate()[0]
   print output

   for iTrial in range(20):
       if('o such file' in output):
           print '[trial #' + str(iTrial) + '] ' + NameTag + calibMapName + ' is not available. Trying again in 30s...'
           time.sleep(30)
           checkFileAvailability = subprocess.Popen([checkFileAvailability_s], stdout=subprocess.PIPE, shell=True);
           output = checkFileAvailability.communicate()[0]
       else:
           break

   print "Done with iteration " + str(iters)
   onlyFIT='False'
   onlyFinalHadd='False'

print "---THE END---"
