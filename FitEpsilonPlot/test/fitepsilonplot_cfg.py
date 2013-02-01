import FWCore.ParameterSet.Config as cms

process = cms.Process("FitEpsilonPlot")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source =   cms.Source("EmptySource")

process.fitEpsilon = cms.EDAnalyzer('FitEpsilonPlot'
)

process.fitEpsilon.OutputFile = cms.untracked.string('calibMap.root')
process.fitEpsilon.CalibType = cms.untracked.string('xtal')
process.fitEpsilon.OutputDir = cms.untracked.string('/cmshome/mgrassi/calib/428-parallelCalib/src/CalibCode/calibtest')
process.fitEpsilon.CurrentIteration = cms.untracked.int32(1)
process.fitEpsilon.EpsilonPlotFileName = cms.untracked.string('epsilonPlots.root')


process.p = cms.Path(process.fitEpsilon)
