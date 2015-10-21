#!/usr/bin/env python

import subprocess, time, sys, os
from methods import *

mode = str(sys.argv[1])
if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):  # Beacause in IIHE the pwd give a link to the area, and you don't want that
    pwd         = os.getenv('PWD')
    num         = 7 #For IIHE T2 you submit all in teh same queue, and the bjobs give an output of 6 for the daemon
else:
    pwd         = os.getcwd()
    num         = 2

if ( mode.find('CRAB')==-1 and mode.find('BATCH_RESU')==-1 ): # Batch system
   if len(sys.argv) != 3:
       print "usage thisPyton.py nITER queue"
       sys.exit(1)
elif ( mode.find('BATCH_RESU') != -1 ):                                   # Batch Resubmission
    if len(sys.argv) != 5:
       print "usage thisPyton.py BATCH_RESU nITER queue nJobs"
       sys.exit(1)
else:                                                                                 # CRAB
    if len(sys.argv) != 6:
       print "usage thisPyton.py CRAB currentITER queue (bsub -q " + queueForDaemon + " 'source " + pwd + "/" + dirname + "/CRAB_files/HaddSendafterCrab_XXX.sh')"
       print "or"
       print "Change CRAB with CRAB_RESU_FinalHadd in: (bsub -q " + queueForDaemon + " 'source " + pwd + "/" + dirname + "/CRAB_files/HaddSendafterCrab_XXX.sh')"
       print "or"
       print "Change CRAB with CRAB_RESU_FitOnly in: (bsub -q " + queueForDaemon + " 'source " + pwd + "/" + dirname + "/CRAB_files/HaddSendafterCrab_XXX.sh')"
       sys.exit(1)
#Selec what mode you are running
RunCRAB = True; RunBatch = True; RunResub = True;
if ( mode.find('CRAB') != -1 ):
     RunCRAB = True; RunBatch = False; RunResub = False;
elif ( mode.find('BATCH_RESU') != -1 ):
     RunCRAB = False; RunBatch = False; RunResub = True;
else:
     RunCRAB = False; RunBatch = True; RunResub = False;
ONLYHADD = False; ONLYFINHADD = False; ONLYFIT = False
if ( mode.find('ONLYHADD') != -1 ):
     ONLYHADD = True;
if ( mode.find('ONLYFINALHADD') != -1 ):
     ONLYFINHADD = True;
if ( mode.find('ONLYFIT') != -1 ):
     ONLYFIT = True;

Add_path = ''
Add_pathOLDIter = ''
ListPaths = []
if ( RunCRAB ):
    nIterations = 1
    njobs = 0
    ListPaths = sys.argv[4].replace('~',' ').split(' ')
    Add_path = ListPaths[0] #When you save the outputs you save them into the 1st path done by CRAB
    Add_pathOLDIter = sys.argv[5] #Path of the previous IC Map
    queue = sys.argv[3]
elif ( RunResub ):
    njobs = int(sys.argv[4])
    queue = sys.argv[3]
    nIterations = nIterations - int(sys.argv[2])
else:
    njobs = int(sys.argv[1])
    queue = sys.argv[2]

outputdir = pwd+'/'+dirname
logPath = outputdir + '/log'
srcPath  = outputdir + '/src'
cfgHaddPath  = outputdir + '/src/hadd'

# To compute the num of hadd
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = inputlist_f.readlines()

for iters in range(nIterations):
    if ( RunCRAB ):
        iters = int(sys.argv[2])
    if ( RunResub ):
        iters = iters + int(sys.argv[2])
    if ( not RunCRAB and not ONLYHADD and not ONLYFIT and not ONLYFINHADD ):
        print "\n*******  ITERATION " + str(iters) + "/" + str(nIterations-1) + "  *******"
        print "Submitting " + str(njobs) + " jobs"
        for ijob in range(njobs):
            #In case you want the stat. syst
            if ( mode.find('BATCH_RESU_SYST_1') != -1 ):
                 env_script_n = open(outputdir + "/cfgFile/Fill/fillEpsilonPlot_iter_" + str(iters) + "_job_" + str(ijob) + ".py", 'a')
                 SystParamLine = 'process.analyzerFillEpsilon.SystOrNot = cms.untracked.double(1)\n'
                 env_script_n.write(SystParamLine)
                 env_script_n.close()
            if ( mode.find('BATCH_RESU_SYST_2') != -1 ):
                 env_script_n = open(outputdir + "/cfgFile/Fill/fillEpsilonPlot_iter_" + str(iters) + "_job_" + str(ijob) + ".py", 'a')
                 SystParamLine = 'process.analyzerFillEpsilon.SystOrNot = cms.untracked.double(2)\n'
                 env_script_n.write(SystParamLine)
                 env_script_n.close()
            # preparing submission of filling tasks
            fill_log_n = logPath + "/fillEpsilonPlot_iter_" + str(iters) + "_job_" + str(ijob) + ".log"
            fill_src_n = srcPath + "/Fill/submit_iter_"     + str(iters) + "_job_" + str(ijob) + ".sh"
            submit_s=""
            if not(Silent):
                 submit_s = "bsub -q " + queue + " -o " + fill_log_n + " " + fill_src_n
            else:
                 submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fill_src_n

            print '\n[job #' + str(ijob) + '] :: ' + submit_s

            # actually submitting filling tasks
            submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
            output = submitJobs.communicate()
            print "Out: " + str(output)

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

            time.sleep(10)

            checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
            datalines = (checkJobs.communicate()[0]).splitlines()
            checkJobs2 = subprocess.Popen(['rm -rf ' + pwd + '/core.*'], stdout=subprocess.PIPE, shell=True);
            datalines2 = (checkJobs2.communicate()[0]).splitlines()
        print 'Done with the Fill part'
    #Crab start from HADD, but it need to rebuild the list of files. So he has this additional part
    if ( mode == 'CRAB' ):
        getGoodfile_str = ''
        if( storageSite=="T2_CH_CERN" ):
           for Extra_path in ListPaths:
               print 'LETS TRY: ' + Extra_path
               print 'Getting Good file: ' + "cmsLs " + eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + Extra_path + " | awk '{print $5}' | grep root | grep -v epsilonPlots | grep -v Barrel | grep -v Endcap | grep " + outputFile +"_"
               getGoodfile = subprocess.Popen(["cmsLs " + eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + Extra_path +  " | awk '{print $5}' | grep root | grep -v epsilonPlots | grep -v Barrel | grep -v Endcap | grep " + outputFile + "_" ], stdout=subprocess.PIPE, shell=True)
               getGoodfile_c = getGoodfile.communicate()
               getGoodfile_str += str(getGoodfile_c)
        if( isOtherT2 and storageSite=="T2_BE_IIHE" ):
           for Extra_path in ListPaths:
               print 'LETS TRY: ' + Extra_path
               print 'Getting Good file: ' + "ls /pnfs/iihe/cms/" + outLFN + "/iter_" + str(iters) + "/" + Extra_path + " | awk '{print $1}' | grep root | awk '{print\"dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + "/" + Extra_path + "/\"$0}'"
               getGoodfile = subprocess.Popen(["ls /pnfs/iihe/cms/" + outLFN + "/iter_" + str(iters) + "/" + Extra_path + " | awk '{print $1}' | grep root | awk '{print\"dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + "/" + Extra_path +"/\"$0}'"], stdout=subprocess.PIPE, shell=True)
               getGoodfile_c = getGoodfile.communicate()
               getGoodfile_str += str(getGoodfile_c)
        getGoodfile_str = getGoodfile_str.replace("\\n", " ")
        getGoodfile_str = getGoodfile_str.replace("'", "")
        getGoodfile_str = getGoodfile_str.replace("(", "")
        getGoodfile_str = getGoodfile_str.replace(")", "")
        getGoodfile_str = getGoodfile_str.replace("None", "")
        getGoodfile_str = getGoodfile_str.replace(",", "")
        getGoodfile_list = getGoodfile_str.split()
        NrelJob = float(len(getGoodfile_list))
        Nlist_flo = float(NrelJob/nHadd) + 1.
        if ( NrelJob%nHadd == 0 ): Nlist_flo -= 1
        Nlist = int(Nlist_flo)
        print "Number of Hadds in parallel (CRAB): " + str(Nlist)
        #Remove Old .list and create new ones
        rmOLDlist = subprocess.Popen(["rm -rf " + srcPath + "/hadd/hadd_iter_" + str(iters) + "_step_*.list" ], stdout=subprocess.PIPE, shell=True)
        rmOLDlist_c = rmOLDlist.communicate()
        rmOLDlist1 = subprocess.Popen(["rm -rf " + srcPath + "/hadd/hadd_iter_" + str(iters) + "_final.list" ], stdout=subprocess.PIPE, shell=True)
        rmOLDlist1_c = rmOLDlist1.communicate()
        haddSrc_n_s = list()
        haddSrc_f_s = list()
        haddSrc_final_n_s = srcPath + "/hadd/hadd_iter_" + str(iters) + "_final.list"
        haddSrc_final_f_s = open(  haddSrc_final_n_s, 'w')
        for num_list in range(Nlist):
            haddSrc_n_s.append( srcPath + "/hadd/hadd_iter_" + str(iters) + "_step_" + str(num_list)+ ".list")
            haddSrc_f_s.append( open(  haddSrc_n_s[num_list], 'w') )
            if(fastHadd):
               fileToAdd_final_n_s = eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
            else:
               fileToAdd_final_n_s = 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
            if( isOtherT2 and storageSite=="T2_BE_IIHE" ):
               fileToAdd_final_n_s = "dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + '/' + NameTag + 'epsilonPlots_' + str(num_list) + '.root\n'
            for nj in range(nHadd):
                nEff = num_list*nHadd+nj
                if( nEff < len(getGoodfile_list) ):
                    if(fastHadd):
                        fileToAdd_n_s = str(getGoodfile_list[nEff]) + '\n'
                    else:
                        fileToAdd_n_s = 'root://eoscms//eos/cms' + str(getGoodfile_list[nEff]) + '\n'
                    if( isOtherT2 and storageSite=="T2_BE_IIHE" ):
                       fileToAdd_n_s = str(getGoodfile_list[nEff]) + '\n'
                    haddSrc_f_s[num_list].write(fileToAdd_n_s)
            haddSrc_final_f_s.write(fileToAdd_final_n_s)
            haddSrc_f_s[num_list].close()
        haddSrc_final_f_s.close()
        #Remove Old .sh and create new ones
        rmOLDsh = subprocess.Popen(["rm -rf " + srcPath + "/hadd/HaddCfg_iter_" + str(iters) + "_job_*.sh" ], stdout=subprocess.PIPE, shell=True)
        rmOLDsh_c = rmOLDsh.communicate()
        dest = eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path
        if( isOtherT2 and storageSite=="T2_BE_IIHE" ):
           dest = "srm://maite.iihe.ac.be:8444/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + "/"
        for num_list in range(Nlist):   
            hadd_cfg_n = cfgHaddPath + "/HaddCfg_iter_" + str(iters) + "_job_" + str(num_list) + ".sh"
            hadd_cfg_f = open( hadd_cfg_n, 'w' )
            HaddOutput = NameTag + "epsilonPlots_" + str(num_list) + ".root"
            if(fastHadd):
                printParallelHaddFAST(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd, num_list )
            else:
                printParallelHadd(hadd_cfg_f, HaddOutput, haddSrc_n_s[num_list], dest, pwd )
            hadd_cfg_f.close()
            changePermission = subprocess.Popen(['chmod 777 ' + hadd_cfg_n], stdout=subprocess.PIPE, shell=True);
            debugout = changePermission.communicate()
        # print Final hadd
        Fhadd_cfg_n = cfgHaddPath + "/Final_HaddCfg_iter_" + str(iters) + ".sh"
        Fhadd_cfg_f = open( Fhadd_cfg_n, 'w' )
        if(fastHadd):
            printFinalHaddFAST(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd )
        else:
            printFinalHadd(Fhadd_cfg_f, haddSrc_final_n_s, dest, pwd )
        Fhadd_cfg_f.close()

    #HADD for batch and CRAB, if you do not want just the finalHADD or the FIT
    if ( mode != 'CRAB_RESU_FinalHadd' and mode != 'CRAB_RESU_FitOnly' and not ONLYFIT and not ONLYFINHADD ):
        print 'Now adding files...'
        if not( RunCRAB ):
           inputlist_v = inputlistbase_v[:]
           NrelJob = float(len(inputlist_v)) / float(ijobmax)
           if( float(int(NrelJob) - NrelJob) < 0. ):
               NrelJob = int(NrelJob) + 1
           Nlist_flo = float(NrelJob/nHadd) + 1.
           Nlist = int(Nlist_flo)
        print "Number of Hadd in parallel: " + str(Nlist)
        for nHadds in range(Nlist):
            Hadd_src_n = srcPath + "/hadd/HaddCfg_iter_" + str(iters) + "_job_" + str(nHadds) + ".sh"
            Hadd_log_n = logPath + "/HaddCfg_iter_" + str(iters) + "_job_" + str(nHadds) + ".log"
            Hsubmit_s = "bsub -q " + queue + " -o " + Hadd_log_n + " bash " + Hadd_src_n
            if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
               Hsubmit_s = "qsub -q localgrid@cream02 -o /dev/null -e /dev/null " +  Hadd_src_n
            #Before each HADD we need ot check if the all the files in the list are present
            #BUT we do that only if you are working on batch
            if not( RunCRAB ):
               if(fastHadd):
                  Grepcommand = "grep -i list " + Hadd_src_n + " | grep -v echo | grep -v bash | awk '{print $2}'"
               else:
                  Grepcommand = "grep -i list " + Hadd_src_n + " | grep -v echo | awk '{print $4}'"
               myGrep = subprocess.Popen([Grepcommand], stdout=subprocess.PIPE, shell=True )
               FoutGrep = myGrep.communicate()
               if(fastHadd):
                  FoutGrep_2 = str(FoutGrep)[2:]
               else:
                  FoutGrep_2 = str(FoutGrep)[3:]
               if(fastHadd):
                  FoutGrep_2 = str(FoutGrep_2)[:-11]
               else:
                  FoutGrep_2 = str(FoutGrep_2)[:-10]
               print 'Checking ' + str(FoutGrep_2)
               #Chech The size for each line
               f = open( str(FoutGrep_2) )
               lines = f.readlines()
               f.close()
               NumToRem = 0
               for filetoCheck in lines:
                   if(fastHadd):
                      filetoCheck = "root://eoscms//eos/cms" + filetoCheck
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
                           print 'removing from Hadd, in: ' + str(FoutGrep_2) + str(NumToRem)
                           f1 = open(str(FoutGrep_2) + str(NumToRem),"w+")
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
            #End of the check, sending the job
            print Hsubmit_s
            subJobs = subprocess.Popen([Hsubmit_s], stdout=subprocess.PIPE, shell=True);
            outJobs = subJobs.communicate()
            print outJobs
            time.sleep(5)
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
           checkJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
           datalines = (checkJobs.communicate()[0]).splitlines()
        else:
           checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
           datalines = (checkJobs.communicate()[0]).splitlines()

        print 'Waiting for all the hadd...'

        # Daemon cheking running jobs
        while len(datalines)>=num :
            if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
               time.sleep(5)
               checkJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
               datalines = (checkJobs.communicate()[0]).splitlines()
            else:
               for entry in datalines:
                   entry = entry.rstrip()
                   entry = entry.split()[0]
                   #print entry
                   if(entry.find('JOBID')!=-1): continue
                   i = int(entry)

               time.sleep(5)

               checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
               datalines = (checkJobs.communicate()[0]).splitlines()
        print 'Done with various hadd'

    if ( mode != 'CRAB_RESU_FitOnly' and not ONLYFIT ):
        print 'Now The Final One...'
        FHadd_src_n = srcPath + "/hadd/Final_HaddCfg_iter_" + str(iters) + ".sh"
        FHadd_log_n = logPath + "/Final_HaddCfg_iter_" + str(iters) + ".log"
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
             FHsubmit_s = "qsub -q localgrid@cream02 -o /dev/null -e /dev/null " + FHadd_src_n
        else:
             FHsubmit_s = "bsub -q " + queue + " -o " + FHadd_log_n + " bash " + FHadd_src_n
        FsubJobs = subprocess.Popen([FHsubmit_s], stdout=subprocess.PIPE, shell=True);
        FoutJobs = FsubJobs.communicate()
        print FoutJobs
        time.sleep(5)

        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
           FcheckJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
           Fdatalines = (FcheckJobs.communicate()[0]).splitlines()
        else:
           FcheckJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
           Fdatalines = (FcheckJobs.communicate()[0]).splitlines()
        print 'Waiting for the Final hadd...'
        # Daemon cheking running jobs
        while len(Fdatalines)>=num :
            if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
               time.sleep(5)
               FcheckJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
               Fdatalines = (FcheckJobs.communicate()[0]).splitlines()
            else:
               for entry in Fdatalines:
                   entry = entry.rstrip()
                   entry = entry.split()[0]
                   #print entry
                   if(entry.find('JOBID')!=-1): continue
                   i = int(entry)

               time.sleep(5)

               FcheckJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
               Fdatalines = (FcheckJobs.communicate()[0]).splitlines()
        print 'Done with final hadd'

    # N of Fit to send
    nEB = 61199/nFit
    if (61199%nFit != 0) :
        nEB = int(nEB) +1
    nEE = 14647/nFit
    if (14647%nFit != 0) :
        nEE = int(nEE) +1
    # For final hadd
    ListFinaHaddEB = list()
    ListFinaHaddEE = list()
    # preparing submission of fit tasks (EB)
    print 'Submitting ' + str(nEB) + ' jobs to fit the Barrel'
    for inteb in range(nEB):
        fit_src_n = srcPath + "/Fit/submit_EB_" + str(inteb) + "_iter_"     + str(iters) + ".sh"
        fit_cfg_n = outputdir + "/cfgFile/Fit/fitEpsilonPlot_EB_" + str(inteb) + "_iter_" + str(iters) + ".py"
        #if isCRAB change the path of the files
        if(isCRAB):
            #Correct path into the cfg
            if not( '#FILE_APPEPENDED' in open(fit_cfg_n).read()):
               with open(fit_cfg_n, 'a') as file:
                    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                        file.write("#FILE_APPEPENDED\n")
                        file.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + '/' + NameTag + "epsilonPlots.root')\n")
                        if(iters==0):
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters-1) + "/" + Add_path + "/" + NameTag + calibMapName + "')\n")
                        else:
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters-1) + "/" + Add_pathOLDIter + "/" + NameTag + calibMapName + "')\n")
                    else:
                        file.write("#FILE_APPEPENDED\n")
                        file.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + Add_path + "/" + NameTag + "epsilonPlots.root')\n")
                        if(iters==0):
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters-1) + "/" + Add_path + "/" + NameTag + calibMapName + "')\n")
                        else:
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters-1) + "/" + Add_pathOLDIter + "/" + NameTag + calibMapName + "')\n")
            #Correct path into the src
            logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EB_" + str(inteb) + "_iter_" + str(iters) + ".log"
            if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                 tmpsource = "$TMPDIR/" + NameTag + "Barrel_" + str(inteb) + "_" + calibMapName
                 destination = "srm://maite.iihe.ac.be:8443/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + "/" + NameTag + "Barrel_" + str(inteb)+ "_" + calibMapName
            else:
                 tmpsource = "/tmp/" + NameTag + "Barrel_" + str(inteb) + "_" + calibMapName
                 destination = eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + Add_path + "/" + NameTag + "Barrel_" + str(inteb)+ "_" + calibMapName
            if not( '#FILE_APPEPENDED' in open(fit_src_n).read()):
               with open(fit_src_n, 'a') as file2:
                    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                        file2.write("#FILE_APPEPENDED\n")
                        file2.write("echo 'srmcp file:///" + tmpsource + " " + destination + "' >> " + logpath  + "\n")
                        file2.write("srmcp file:///" + tmpsource + " " + destination + " >> " + logpath + " 2>&1 \n")
                    else:
                        file2.write("#FILE_APPEPENDED\n")
                        file2.write("echo 'cmsStage -f " + tmpsource + " " + destination + "' >> " + logpath  + "\n")
                        file2.write("cmsStage -f " + tmpsource + " " + destination + " >> " + logpath + " 2>&1 \n")
                    file2.write("echo 'rm -f " + tmpsource + "' >> " + logpath + " \n")
                    file2.write("rm -f " + tmpsource + " >> " + logpath + " 2>&1 \n")
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            submit_s = "qsub -q localgrid@cream02 -o /dev/null -e /dev/null " + fit_src_n
            ListFinaHaddEB.append("dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + '/' + NameTag + 'Barrel_'+str(inteb) + '_' + calibMapName)
        else:
            submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fit_src_n
            ListFinaHaddEB.append('root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName )
        print 'About to EB fit:'
        print 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Barrel_'+str(inteb)+'_' + calibMapName
        print submit_s
        # actually submitting fit tasks (EB)
        submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
        output = submitJobs.communicate()
        print output

    # preparing submission of fit tasks (EE)
    print 'Submitting ' + str(nEE) + ' jobs to fit the Endcap'
    for inte in range(nEE):        
        fit_src_n = srcPath + "/Fit/submit_EE_" + str(inte) + "_iter_"     + str(iters) + ".sh"
        fit_cfg_n = outputdir + "/cfgFile/Fit/fitEpsilonPlot_EE_" + str(inte) + "_iter_" + str(iters) + ".py"
        #if isCRAB change the path of the files
        if(isCRAB):
            #Correct path into the cfg
            if not( '#FILE_APPEPENDED' in open(fit_cfg_n).read()):
               with open(fit_cfg_n, 'a') as file:
                    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                        file.write("#FILE_APPEPENDED\n")
                        file.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + '/' + NameTag + "epsilonPlots.root')\n")
                        if(iters==0):
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters-1) + "/" + Add_path + "/" + NameTag + calibMapName + "')\n")
                        else:
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters-1) + "/" + Add_pathOLDIter + "/" + NameTag + calibMapName + "')\n")
                    else:
                        file.write("#FILE_APPEPENDED\n")
                        file.write("process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters) + "/" + Add_path + "/" + NameTag + "epsilonPlots.root')\n")
                        if(iters==0):
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters-1) + "/" + Add_path + "/" + NameTag + calibMapName + "')\n")
                        else:
                           file.write("process.fitEpsilon.calibMapPath = cms.untracked.string('root://eoscms//eos/cms" + eosPath + "/" + dirname + "/iter_" + str(iters-1) + "/" + Add_pathOLDIter + "/" + NameTag + calibMapName + "')\n")
            #Correct path into the src
            logpath = pwd + "/" + dirname + "/log/" + "fitEpsilonPlot_EE_" + str(inte) + "_iter_" + str(iters) + ".log"
            if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                tmpsource = "$TMPDIR/" + NameTag + "Endcap_" + str(inte) + "_" + calibMapName
                destination = "srm://maite.iihe.ac.be:8443/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + "/" + NameTag + "Endcap_" + str(inte)+ "_" + calibMapName
            else:
                tmpsource = "/tmp/" + NameTag + "Endcap_" + str(inte) + "_" + calibMapName
                destination = eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + Add_path + "/" + NameTag + "Endcap_" + str(inte)+ "_" + calibMapName
            if not( '#FILE_APPEPENDED' in open(fit_src_n).read()):
               with open(fit_src_n, 'a') as file2:
                    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
                        file2.write("#FILE_APPEPENDED\n")
                        file2.write("echo 'srmcp file:///" + tmpsource + " " + destination + "' >> " + logpath  + "\n")
                        file2.write("srmcp file:///" + tmpsource + " " + destination + " >> " + logpath + " 2>&1 \n")
                    else:
                        file2.write("#FILE_APPEPENDED\n")
                        file2.write("echo 'cmsStage -f " + tmpsource + " " + destination + "' >> " + logpath  + "\n")
                        file2.write("cmsStage -f " + tmpsource + " " + destination + " >> " + logpath + " 2>&1 \n")
                    file2.write("echo 'rm -f " + tmpsource + "' >> " + logpath + " \n")
                    file2.write("rm -f " + tmpsource + " >> " + logpath + " 2>&1 \n")
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            submit_s = "qsub -q localgrid@cream02 -o /dev/null -e /dev/null " + fit_src_n
            ListFinaHaddEE.append("dcap://maite.iihe.ac.be/pnfs/iihe/cms" + outLFN + "/iter_" + str(iters) + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName)
        else:
            submit_s = "bsub -q " + queue + " -o /dev/null -e /dev/null " + fit_src_n
            ListFinaHaddEE.append('root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName)
        print 'About to EE fit:'
        print 'root://eoscms//eos/cms' + eosPath + '/' + dirname + '/iter_' + str(iters) + '/' + Add_path + '/' + NameTag + 'Endcap_'+str(inte) + '_' + calibMapName
        print submit_s
        # actually submitting fit tasks (EE)
        submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
        output = submitJobs.communicate()
        print output

    # checking number of running/pending jobs
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        checkJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
        datalines = (checkJobs.communicate()[0]).splitlines()
    else:
        checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
        datalines = (checkJobs.communicate()[0]).splitlines()
    print 'Waiting for fit jobs to be finished...'

    #Daemon cheking running jobs
    while len(datalines)>=num :
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            time.sleep(5)
            checkJobs = subprocess.Popen(['qstat -u $USER localgrid@cream02'], stdout=subprocess.PIPE, shell=True);
            datalines = (checkJobs.communicate()[0]).splitlines()
        else:
            for entry in datalines:
                entry = entry.rstrip()
                entry = entry.split()[0]
                #print entry
                if(entry.find('JOBID')!=-1): continue
                i = int(entry)

            time.sleep(5)

            checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
            datalines = (checkJobs.communicate()[0]).splitlines()

    print "Done with fitting! Now we have to merge all fits in one Calibmap.root"

    # Merge Final CalibMap
    from ROOT import *
    from PhysicsTools.PythonAnalysis import *
    gSystem.Load("libFWCoreFWLite.so")
    AutoLibraryLoader.enable()
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        f = TFile('$TMPDIR/' + NameTag + calibMapName, 'recreate')
    else:
        f = TFile('/tmp/' + NameTag + calibMapName, 'recreate')
    #Run only on EB or EE if needed
    ListFinaHadd = list()
    if Barrel_or_Endcap=='ONLY_BARREL':
       ListFinaHadd = ListFinaHaddEB
    if Barrel_or_Endcap=='ONLY_ENDCAP':
       ListFinaHadd = ListFinaHaddEE
    if (Barrel_or_Endcap=='ALL_PLEASE'):
       ListFinaHadd = ListFinaHaddEB
       ListFinaHadd = ListFinaHadd + ListFinaHaddEE

    # Create a struct
    if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
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
           Double_t fit_mean_err_;\
           Double_t fit_sigma_;\
           Double_t fit_Snorm_;\
           Double_t fit_b0_;\
           Double_t fit_b1_;\
           Double_t fit_b2_;\
           Double_t fit_b3_;\
           Double_t fit_Bnorm_;\
         };")
       s = EBStruct()
    if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
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
           Double_t fit_mean_err_;\
           Double_t fit_sigma_;\
           Double_t fit_Snorm_;\
           Double_t fit_b0_;\
           Double_t fit_b1_;\
           Double_t fit_b2_;\
           Double_t fit_b3_;\
           Double_t fit_Bnorm_;\
         };")
       t = EEStruct()
    if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
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
           Double_t fit_mean_err;\
           Double_t fit_sigma;\
           Double_t fit_Snorm;\
           Double_t fit_b0;\
           Double_t fit_b1;\
           Double_t fit_b2;\
           Double_t fit_b3;\
           Double_t fit_Bnorm;\
         };")
       s1 = EB1Struct()
    if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
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
           Double_t fit_mean_err;\
           Double_t fit_sigma;\
           Double_t fit_Snorm;\
           Double_t fit_b0;\
           Double_t fit_b1;\
           Double_t fit_b2;\
           Double_t fit_b3;\
           Double_t fit_Bnorm;\
         };")
       t1 = EE1Struct()

    calibMap_EB = TH2F("calibMap_EB", "EB calib coefficients: #eta on x, #phi on y", 171,-85.5,85.5 , 360,0.5,360.5)
    calibMap_EEm = TH2F("calibMap_EEm", "EE- calib coefficients", 100,0.5,100.5,100,0.5,100.5)
    calibMap_EEp = TH2F("calibMap_EEp", "EE+ calib coefficients", 100,0.5,100.5,100,0.5,100.5)
    TreeEB = TTree("calibEB", "Tree of EB Inter-calibration constants")
    if(Barrel_or_Endcap=='ONLY_BARREL' or Barrel_or_Endcap=='ALL_PLEASE'):
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
       TreeEB.Branch('fit_mean_err_'   , AddressOf(s,'fit_mean_err_'),'fit_mean_err_/F')
       TreeEB.Branch('fit_sigma_'  , AddressOf(s,'fit_sigma_'),'fit_sigma_/F')
       TreeEB.Branch('fit_Snorm_'  , AddressOf(s,'fit_Snorm_'),'fit_Snorm_/F')
       TreeEB.Branch('fit_b0_'     , AddressOf(s,'fit_b0_'),'fit_b0_/F')
       TreeEB.Branch('fit_b1_'     , AddressOf(s,'fit_b1_'),'fit_b1_/F')
       TreeEB.Branch('fit_b2_'     , AddressOf(s,'fit_b2_'),'fit_b2_/F')
       TreeEB.Branch('fit_b3_'     , AddressOf(s,'fit_b3_'),'fit_b3_/F')
       TreeEB.Branch('fit_Bnorm_'  , AddressOf(s,'fit_Bnorm_'),'fit_Bnorm_/F')

    TreeEE = TTree("calibEE", "Tree of EE Inter-calibration constants")
    if(Barrel_or_Endcap=='ONLY_ENDCAP' or Barrel_or_Endcap=='ALL_PLEASE'):
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
       TreeEE.Branch('fit_mean_err_'   , AddressOf(t,'fit_mean_err_'),'fit_mean_err_/F')
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
           thisTree.SetBranchAddress( 'fit_mean_err',AddressOf(s1,'fit_mean_err'));
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
                   s.fit_mean_err_ = s1.fit_mean_err
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
           thisTree.SetBranchAddress( 'fit_mean_err',AddressOf(t1,'fit_mean_err'));
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
                   t.fit_mean_err_ = t1.fit_mean_err
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
    if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
        stage_s_fin = 'srmcp file:///$TMPDIR//' + NameTag + calibMapName + ' srm://maite.iihe.ac.be:8443/pnfs/iihe/cms' + outLFN + "/iter_" + str(iters) + "/" + NameTag + calibMapName
    else:
        stage_s_fin = 'cmsStage -f /tmp/' + NameTag + calibMapName + ' ' + eosPath + '/' + dirname + '/iter_' + str(iters) + "/" + Add_path + "/" + NameTag + calibMapName
    print stage_s_fin
    stageCalibFile = subprocess.Popen([stage_s_fin], stdout=subprocess.PIPE, shell=True);
    print stageCalibFile.communicate()
    if ( mode.find('CRAB') != -1 ):
        print 'Since we are in CRAB mode I will also copy ' + NameTag + calibMapName + 'on data/ to be accessible from CRAB' 
        if( isOtherT2 and storageSite=="T2_BE_IIHE" and isCRAB ):
            stage_s_fin2 = 'cp $TMPDIR//' + NameTag + calibMapName + " " + pwd + "/../FillEpsilonPlot/data/"
        else:
            stage_s_fin2 = 'cp /tmp/' + NameTag + calibMapName + " " + pwd + "/../FillEpsilonPlot/data/"
        print stage_s_fin2
        stageCalibFile2 = subprocess.Popen([stage_s_fin2], stdout=subprocess.PIPE, shell=True);
        print stageCalibFile2.communicate()
    print 'Done with staging the final ' + NameTag + calibMapName

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
    if( ONLYHADD or ONLYFINHADD or ONLYFIT):
       mode = "BATCH_RESU"
       ONLYHADD = False; ONLYFINHADD = False; ONLYFIT=False;

print "---THE END---"
