import math
import ROOT
from ROOT import *
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *
#import print_options

#print_options.set_float_precision(4)
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()

#import EcalDetId
#from DataFormats.EcalDetId import *
#
import sys,os

allRecHits=False
lumi=393
lumi=-1
lumi=466
eventNumber=-1
eventNumber=685233340
#lumi=376
#eventNumber=344404776
#eventNumber=552337102
#eventNumber=577412878

#eventNumber=347006644
#eventNumber=344911692

lumi=374
eventNumber=548668449
# kRecovered lumi=376 eventNumber=552337102
# kRecovered lumi=393 eventNumber=576932539
# kRecovered lumi=394 eventNumber=578483322
lumi=394
eventNumber=578490502
lumi=394
eventNumber=579700192
lumi=395
eventNumber=579843406
lumi=401
eventNumber=588810213
lumi=402
eventNumber=591275401
# kRecovered lumi=403 eventNumber=593293410
lumi=403
eventNumber=591888388
# kRecovered lumi=406 eventNumber=597401546
# kRecovered lumi=407 eventNumber=598290564
lumi=415
eventNumber=610541757
lumi=415
eventNumber=610541436
lumi=416
eventNumber=612542602
# kRecovered lumi=419 eventNumber=616682572
# kRecovered lumi=419 eventNumber=615876590
lumi=422
eventNumber=620689835

lumi=466
eventNumber=685900276
lumi=467
eventNumber=687572911
lumi=472
eventNumber=694966852
# kRecovered
lumi=466
eventNumber=685233340

eventNumber=-1
eventMin=-1
lumi=-1

# for now look for events in two files with a given lumi section
maxEvents=2
event_counter=0

file="alcaraw.root"



events = Events (file)
print file
handleElectrons = Handle('std::vector<reco::GsfElectron>')

handleRecHitsES = Handle('edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >')
handleRecHitsEB_ALCARECO = Handle('edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >')
handleRecHitsEE_ALCARECO = Handle('edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >')

handle_hltdigis_EB = Handle('EBDigiCollection')
handle_hltdigis_EE = Handle('EEDigiCollection')

handle_dummyDigis_EB = Handle('EBDigiCollection')
handle_dummyDigis_EE = Handle('EEDigiCollection')

handle_cleanedDigis_EB = Handle('EBDigiCollection')
handle_cleanedDigis_EE = Handle('EEDigiCollection')


EErecHitmap_ele1 = TH2F("EErecHitmap_ele1", "EErecHitmap_ele1",
                   100,0,100,
                   100,0,100)

EBrecHitmap_ele1 = TH2F("EBrecHitmap_ele1", "EBrecHitmap_ele1",
                   171,-85,85,
                   360,0,360)

EErecHitmap_ele2 = TH2F("EErecHitmap_ele2", "EErecHitmap_ele2",
                   100,0,100,
                   100,0,100)

EBrecHitmap_ele2 = TH2F("EBrecHitmap_ele2", "EBrecHitmap_ele2",
                   171,-85,85,
                   360,0,360)

#print file_format, file, electronTAG, processName, maxEvents

print "run\tlumi, event, energy, eSC, rawESC, e5x5, E_ES, etaEle, phiEle, etaSC, phiSC, clustersSize, nRecHits"
for event in events:

    if(maxEvents > 0 and event_counter > maxEvents):
        break
    #if(event.eventAuxiliary.run()== 145351895):
    if lumi > 0 and int(event.eventAuxiliary().luminosityBlock()) != lumi :
            continue

    if(eventNumber > 0 and event.eventAuxiliary().event()!= eventNumber ):
        continue

    event.getByLabel("hltAlCaEtaEBRechitsToDigis", "etaEBDigis", "HLT", handle_hltdigis_EB)
    event.getByLabel("hltAlCaEtaEERechitsToDigis", "etaEEDigis", "HLT", handle_hltdigis_EE)

    event.getByLabel("hltAlCaPi0EBRechitsToDigis", "pi0EBDigis", "HLT", handle_hltdigis_EB)
    event.getByLabel("hltAlCaPi0EERechitsToDigis", "pi0EEDigis", "HLT", handle_hltdigis_EE)


    event.getByLabel("dummyHits"   ,"dummyBarrelDigis"       ,   "analyzerFillEpsilon", handle_dummyDigis_EB)
    event.getByLabel("dummyHits"   ,"dummyEndcapDigis"       ,   "analyzerFillEpsilon", handle_dummyDigis_EE)

    event.getByLabel("digiCleaning","cleanedEBDigiCollection",   "analyzerFillEpsilon", handle_cleanedDigis_EB)
    event.getByLabel("digiCleaning","cleanedEEDigiCollection",   "analyzerFillEpsilon", handle_cleanedDigis_EE)


    if(handle_hltdigis_EB.isValid()):
        digis= handle_hltdigis_EB.product()
        for digi in digis:
            print "EBDigi:", digi.id()

    if(handle_hltdigis_EE.isValid()):
        digis= handle_hltdigis_EE.product()
        for digi in digis:
            print "EEDigi:", digi.id()

    print "End of Event"
    event_counter=event_counter+1
    continue

# edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "ecalRecHit"                "EcalRecHitsEB"   "analyzerFillEpsilon"   
# edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >    "ecalRecHit"                "EcalRecHitsEE"   "analyzerFillEpsilon"   
# edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit> >    "ecalMultiFitUncalibRecHit"   "EcalUncalibRecHitsEB"   "analyzerFillEpsilon"   
# edm::SortedCollection<EcalUncalibratedRecHit,edm::StrictWeakOrdering<EcalUncalibratedRecHit> >    "ecalMultiFitUncalibRecHit"   "EcalUncalibRecHitsEE"   "analyzerFillEpsilon"   
# edm::TriggerResults                   "TriggerResults"            ""                "analyzerFillEpsilon"   
# unsigned int                          "bunchSpacingProducer"      ""                "analyzerFillEpsilon"   


    event.getByLabel("ecalRecHit", "EcalRecHitsEB", "analyzerFillEpsilon", handleRecHitsEB_ALCARECO)
    event.getByLabel("ecalRecHit", "EcalRecHitsEE", "analyzerFillEpsilon", handleRecHitsEE_ALCARECO)

    print "#############################"    
    print "RECHIT:  rawId    energy"
    recHits_ALCARECO_EB = handleRecHitsEB_ALCARECO.product()
    recHits_ALCARECO_EE = handleRecHitsEE_ALCARECO.product()
        
    nRecHits_ALCARECO=0

    for recHit in recHits_ALCARECO_EB:
        nRecHits_ALCARECO=nRecHits_ALCARECO+1
        #               if(recHit.checkFlag(EcalRecHit.kTowerRecovered)):
        #           print recHit.id().rawId(), recHit.checkFlag(EcalRecHit.kTowerRecovered)
        print recHit.id().rawId(), recHit.energy()

    for recHit in recHits_ALCARECO_EE:
        nRecHits_ALCARECO=nRecHits_ALCARECO+1
        #               if(recHit.checkFlag(EcalRecHit.kTowerRecovered)):
        #           print recHit.id().rawId(), recHit.checkFlag(EcalRecHit.kTowerRecovered)
        print recHit.id().rawId(), recHit.energy()

print event_counter



