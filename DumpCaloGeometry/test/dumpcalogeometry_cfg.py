import FWCore.ParameterSet.Config as cms

#process = cms.Process("Demo")
#
#process.load("FWCore.MessageService.MessageLogger_cfi")
#
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#
#process.source = cms.Source("PoolSource",
#    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        'file:myfile.root'
#    )
#)
#
#process.demo = cms.EDAnalyzer('DumpCaloGeometry'
#)
#
#
#process.p = cms.Path(process.demo)


process = cms.Process("GeometryTest")
from Geometry.CaloEventSetup.CaloTopology_cfi import *
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MC_3XY_V18::All'

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )


process.dumper = cms.EDAnalyzer("DumpCaloGeometry")

process.p1 = cms.Path(process.dumper)
